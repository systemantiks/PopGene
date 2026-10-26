//! Bounded Staden STHUFF: canonical literal Huffman with interlaced code sets.
//! Format/algorithms: io_lib deflate_interlaced.c and RFC1951 dynamic headers.
//! BSD source attribution: vendor-notices/staden-io-lib/deflate_interlaced-LICENSE.txt.
use super::{LIMIT, Result, bad};
use std::collections::BTreeMap;
#[path = "huffman_tables.rs"]
mod predefined;

struct Bits<'a> {
    data: &'a [u8],
    position: usize,
}
impl Bits<'_> {
    fn read(&mut self, count: usize) -> Result<usize> {
        if count > 16
            || self
                .position
                .checked_add(count)
                .is_none_or(|end| end > self.data.len() * 8)
        {
            return Err(bad("Truncated STHUFF bitstream"));
        }
        let mut value = 0;
        for i in 0..count {
            value |= usize::from((self.data[self.position / 8] >> (self.position % 8)) & 1) << i;
            self.position += 1;
        }
        Ok(value)
    }
}
#[derive(Clone, Default)]
struct Node {
    children: [Option<usize>; 2],
    symbol: Option<u16>,
}
#[derive(Clone)]
struct Tree {
    nodes: Vec<Node>,
}
impl Tree {
    fn new(lengths: &[u8]) -> Result<Self> {
        let mut counts = [0usize; 16];
        for &length in lengths {
            if length > 15 {
                return Err(bad("STHUFF code exceeds15 bits"));
            }
            counts[length as usize] += 1;
        }
        if counts[0] == lengths.len() {
            return Err(bad("Empty STHUFF code tree"));
        }
        counts[0] = 0;
        let mut space = 1isize;
        for count in &counts[1..] {
            space = space * 2 - *count as isize;
            if space < 0 {
                return Err(bad("Oversubscribed STHUFF code tree"));
            }
        }
        let mut next = [0usize; 16];
        let mut code = 0;
        for length in 1..16 {
            code = (code + counts[length - 1]) << 1;
            next[length] = code;
        }
        let mut result = Self {
            nodes: vec![Node::default()],
        };
        for (symbol, &length) in lengths.iter().enumerate() {
            if length == 0 {
                continue;
            }
            let code = next[length as usize];
            next[length as usize] += 1;
            let mut node = 0;
            for bit in (0..length).rev() {
                if result.nodes[node].symbol.is_some() {
                    return Err(bad("Conflicting STHUFF prefix codes"));
                }
                let direction = (code >> bit) & 1;
                node = match result.nodes[node].children[direction] {
                    Some(child) => child,
                    None => {
                        let child = result.nodes.len();
                        result.nodes.push(Node::default());
                        result.nodes[node].children[direction] = Some(child);
                        child
                    }
                };
            }
            if result.nodes[node].symbol.is_some()
                || result.nodes[node].children.iter().any(Option::is_some)
            {
                return Err(bad("Duplicate STHUFF code"));
            }
            result.nodes[node].symbol = Some(symbol as u16);
        }
        Ok(result)
    }
    fn symbol(&self, bits: &mut Bits<'_>) -> Result<usize> {
        let mut node = 0;
        for _ in 0..15 {
            node = self.nodes[node].children[bits.read(1)?]
                .ok_or_else(|| bad("Invalid STHUFF prefix"))?;
            if let Some(symbol) = self.nodes[node].symbol {
                return Ok(symbol as usize);
            }
        }
        Err(bad("Invalid STHUFF code length"))
    }
}
#[derive(Clone)]
struct Codes {
    trees: Vec<Tree>,
    start_bit: usize,
}
impl Codes {
    fn nodes(&self) -> usize {
        self.trees.iter().map(|tree| tree.nodes.capacity()).sum()
    }
}
fn tree(bits: &mut Bits<'_>) -> Result<Tree> {
    let literal_count = bits.read(5)? + 257;
    let distance_count = bits.read(5)? + 1;
    let count = bits.read(4)? + 4;
    let order = [
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
    ];
    let mut code_lengths = [0u8; 19];
    for &symbol in &order[..count] {
        code_lengths[symbol] = bits.read(3)? as u8;
    }
    let code_tree = Tree::new(&code_lengths)?;
    let total = literal_count + distance_count;
    let mut lengths = Vec::with_capacity(total);
    while lengths.len() < total {
        let (value, count) = match code_tree.symbol(bits)? {
            v @ 0..=15 => (v as u8, 1),
            16 => (
                *lengths
                    .last()
                    .ok_or_else(|| bad("STHUFF repeat before first length"))?,
                bits.read(2)? + 3,
            ),
            17 => (0, bits.read(3)? + 3),
            18 => (0, bits.read(7)? + 11),
            _ => return Err(bad("Invalid STHUFF length symbol")),
        };
        if lengths.len() + count > total {
            return Err(bad("STHUFF code-length repeat exceeds its alphabet"));
        }
        lengths.extend(std::iter::repeat_n(value, count));
    }
    if lengths[257..literal_count].iter().any(|v| *v != 0) {
        return Err(bad("STHUFF only supports byte literals and end-of-block"));
    }
    Tree::new(&lengths[..257])
}
fn codes(bits: &mut Bits<'_>) -> Result<(Codes, bool)> {
    let final_block = bits.read(1)? != 0;
    let count = match bits.read(2)? {
        2 => 1,
        3 => {
            let width = bits.read(4)? + 1;
            bits.read(width)? + 1
        }
        _ => return Err(bad("STHUFF requires dynamic or interlaced Huffman blocks")),
    };
    if count > 256 {
        return Err(bad("STHUFF exceeds256 interlaced code tables"));
    }
    let trees = (0..count).map(|_| tree(bits)).collect::<Result<Vec<_>>>()?;
    if !trees
        .iter()
        .any(|t| t.nodes.iter().any(|n| n.symbol == Some(256)))
    {
        return Err(bad("STHUFF code set has no end-of-block symbol"));
    }
    Ok((
        Codes {
            trees,
            start_bit: bits.position % 8,
        },
        final_block,
    ))
}

#[derive(Default)]
pub(crate) struct HuffmanTables {
    tables: BTreeMap<u8, Codes>,
    nodes: usize,
}
impl HuffmanTables {
    /// HUFF payload after ordinary compression has been removed: RAW,id,codes.
    pub(crate) fn insert(&mut self, data: &[u8]) -> Result<()> {
        if data.len() < 3 || data[0] != 0 || data[1] < 128 {
            return Err(bad("Invalid shared HUFF header"));
        }
        if self.tables.contains_key(&data[1]) {
            return Err(bad("Duplicate shared HUFF code-set identifier"));
        }
        let (codes, _) = codes(&mut Bits {
            data: &data[2..],
            position: 0,
        })?;
        let nodes = self.nodes + codes.nodes();
        if nodes > 1_000_000 {
            return Err(bad("Shared HUFF tables exceed the node budget"));
        }
        self.nodes = nodes;
        self.tables.insert(data[1], codes);
        Ok(())
    }
}
pub(super) fn decode(data: &[u8], shared: &HuffmanTables) -> Result<Vec<u8>> {
    if data.len() < 3 {
        return Err(bad("Truncated STHUFF header"));
    }
    let mut current = match data[1] {
        0 => None,
        id @ 1..=3 => Some(Codes {
            trees: vec![Tree::new(&predefined::PREDEFINED[id as usize - 1])?],
            start_bit: if id == 1 { 5 } else { 1 },
        }),
        id @ 128..=255 => Some(
            shared
                .tables
                .get(&id)
                .ok_or_else(|| bad("Missing shared HUFF code set"))?
                .clone(),
        ),
        _ => return Err(bad("Unknown predefined STHUFF code-set identifier")),
    };
    let mut bits = Bits {
        data: &data[2..],
        position: current.as_ref().map_or(0, |c| c.start_bit),
    };
    let mut output = Vec::new();
    for _ in 0..1024 {
        let (code, final_block) = match current.take() {
            Some(c) => (c, true),
            None => codes(&mut bits)?,
        };
        let mut column = 0;
        loop {
            let symbol = code.trees[column].symbol(&mut bits)?;
            if symbol == 256 {
                break;
            }
            if output.len() == LIMIT {
                return Err(bad("STHUFF expansion exceeds64 MiB"));
            }
            if output.len() == output.capacity() {
                let capacity = output.capacity().saturating_mul(2).clamp(1024, LIMIT);
                output.reserve_exact(capacity - output.len());
            }
            output.push(symbol as u8);
            column += 1;
            if column == code.trees.len() {
                column = 0;
            }
        }
        if final_block {
            // At most the unused tail of the final byte may remain.
            if bits.data.len() * 8 - bits.position >= 8 {
                return Err(bad("Trailing bytes after STHUFF final block"));
            }
            return Ok(output);
        }
    }
    Err(bad("STHUFF exceeds1024 blocks"))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_staden_huffman_vectors_and_shared_file() {
        fn hex(text: &str) -> Vec<u8> {
            text.as_bytes()
                .as_chunks::<2>()
                .0
                .iter()
                .map(|v| u8::from_str_radix(std::str::from_utf8(v).unwrap(), 16).unwrap())
                .collect()
        }
        let facts: serde_json::Value = serde_json::from_str(include_str!(
            "../../../../validation/fixtures/trace/extended/huffman-vectors.json"
        ))
        .unwrap();
        for vector in facts["vectors"].as_array().unwrap() {
            let mut tables = HuffmanTables::default();
            if let Some(table) = vector["table"].as_str() {
                tables.insert(&hex(table)).unwrap();
            }
            let actual = super::super::decode_with_tables(
                &hex(vector["input"].as_str().unwrap()),
                None,
                &tables,
            );
            assert!(actual.is_ok(), "{}: {:?}", vector["name"], actual.err());
            assert_eq!(
                actual.unwrap().0,
                hex(vector["expected"].as_str().unwrap()),
                "{}",
                vector["name"]
            );
        }
        let fixture = &facts["fixture"];
        let bytes =
            include_bytes!("../../../../validation/fixtures/trace/extended/shared-huffman.ztr");
        let doc = crate::import("shared-huffman.ztr", bytes).unwrap();
        assert_eq!(doc.sequence.as_ref(), fixture["sequence"].as_str().unwrap());
        let trace = doc.trace.as_ref().unwrap();
        assert_eq!(
            trace.channels,
            serde_json::from_value::<Vec<Vec<u16>>>(fixture["channels"].clone()).unwrap()
        );
        assert_eq!(
            trace.qualities,
            serde_json::from_value::<Vec<u8>>(fixture["qualities"].clone()).unwrap()
        );
    }
    #[test]
    fn rejects_missing_tables_and_malformed_codes() {
        assert!(decode(&[77, 128, 0], &HuffmanTables::default()).is_err());
        assert!(decode(&[77, 4, 0], &HuffmanTables::default()).is_err());
        assert!(decode(&[77, 0, 0], &HuffmanTables::default()).is_err());
        assert!(Tree::new(&[1, 1, 1]).is_err());
        assert!(Tree::new(&[0, 0]).is_err());
        assert!(decode(&[77, 0, 255, 255, 127], &HuffmanTables::default()).is_err());
        let facts: serde_json::Value = serde_json::from_str(include_str!(
            "../../../../validation/fixtures/trace/extended/huffman-vectors.json"
        ))
        .unwrap();
        let raw = facts["vectors"]
            .as_array()
            .unwrap()
            .iter()
            .find_map(|v| v["table"].as_str())
            .unwrap();
        let data = raw
            .as_bytes()
            .as_chunks::<2>()
            .0
            .iter()
            .map(|p| u8::from_str_radix(std::str::from_utf8(p).unwrap(), 16).unwrap())
            .collect::<Vec<_>>();
        let mut tables = HuffmanTables::default();
        tables.insert(&data).unwrap();
        assert!(tables.insert(&data).is_err());
    }
}

//! Explicit creation inputs; transformations are opt-in and invalid symbols are never dropped.
use super::{Document, Error, MAX_BASES, Molecule};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Mode {
    #[default]
    Sequence,
    Unknown,
    Placeholder,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Request {
    pub molecule: Molecule,
    pub mode: Mode,
    pub name: String,
    pub text: String,
    pub count: String,
    pub circular: bool,
    pub single_stranded: bool,
    pub reverse: bool,
    pub convert_tu: bool,
    pub ignore_numbers: bool,
}
impl Request {
    pub fn new(molecule: Molecule) -> Self {
        Self {
            molecule,
            mode: Mode::Sequence,
            name: "Untitled".into(),
            text: String::new(),
            count: "100".into(),
            circular: false,
            single_stranded: molecule == Molecule::Rna,
            reverse: false,
            convert_tu: false,
            ignore_numbers: true,
        }
    }
    pub fn create(&self) -> Result<Document, Error> {
        let name = self.name.trim();
        if name.is_empty() || name.chars().any(char::is_control) {
            return Err(Error::Operation(
                "Enter a sequence name without control characters.".into(),
            ));
        }
        let mut sequence = match self.mode {
            Mode::Sequence => parse_text(
                &self.text,
                self.molecule,
                self.convert_tu,
                self.ignore_numbers,
            )?,
            Mode::Unknown => {
                let count = self
                    .count
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .filter(|n| (1..=MAX_BASES).contains(n))
                    .ok_or_else(|| {
                        Error::Operation("Enter a length from 1 to 10,000,000.".into())
                    })?;
                if self.molecule == Molecule::Protein {
                    "X"
                } else {
                    "N"
                }
                .repeat(count)
            }
            Mode::Placeholder => String::new(),
        };
        if sequence.is_empty() && self.mode != Mode::Placeholder {
            return Err(Error::Operation(
                "Enter a sequence, or choose Empty placeholder.".into(),
            ));
        }
        let mut document = Document::new_typed(
            name,
            if sequence.is_empty() { "A" } else { &sequence },
            self.molecule,
        )?;
        if self.mode == Mode::Placeholder {
            document.sequence = std::mem::take(&mut sequence).into();
            document.placeholder = true;
        }
        document.circular = self.circular;
        document.single_stranded = self.single_stranded;
        document.validate()?;
        if self.reverse {
            document.revcomp()?;
        }
        Ok(document)
    }
}
/// One raw sequence or one FASTA record. No punctuation, gaps or invalid letters are discarded.
pub fn parse_text(
    text: &str,
    molecule: Molecule,
    convert_tu: bool,
    ignore_numbers: bool,
) -> Result<String, Error> {
    if text.len() > 32 * 1024 * 1024 {
        return Err(Error::Operation(
            "Pasted text exceeds 32 MiB. Open a sequence file instead.".into(),
        ));
    }
    let text = text.trim_start_matches('\u{feff}').trim();
    let body = if text.starts_with('>') {
        let (header, body) = text.split_once('\n').unwrap_or((text, ""));
        if header[1..].trim().is_empty() {
            return Err(Error::Operation("FASTA needs a record name.".into()));
        }
        if body.lines().any(|l| l.trim_start().starts_with('>')) {
            return Err(Error::Operation(
                "Multiple FASTA records found. Create one document per record.".into(),
            ));
        }
        body
    } else {
        text
    };
    let input: String = body
        .chars()
        .filter(|c| !(ignore_numbers && c.is_ascii_digit()))
        .map(|c| {
            if convert_tu {
                match (molecule, c) {
                    (Molecule::Rna, 'T' | 't') => 'U',
                    (Molecule::Dna, 'U' | 'u') => 'T',
                    _ => c,
                }
            } else {
                c
            }
        })
        .collect();
    molecule.normalize(&input)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn typed_creation_preserves_alphabets_and_topology() {
        let mut r = Request::new(Molecule::Rna);
        r.name = "Synthetic RNA".into();
        r.text = " 1 acguNry 7\n".into();
        r.circular = true;
        let mut d = r.create().unwrap();
        assert_eq!(&*d.sequence, "ACGUNRY");
        assert!(d.circular && d.single_stranded);
        let original = d.clone();
        d.revcomp().unwrap();
        assert_eq!(&*d.sequence, "RYNACGU");
        d.revcomp().unwrap();
        assert_eq!(d, original);
        r.text = "ATG".into();
        assert!(r.create().is_err());
        r.convert_tu = true;
        assert_eq!(&*r.create().unwrap().sequence, "AUG");
    }
    #[test]
    fn protein_has_no_dna_operations_and_edits_keep_its_alphabet() {
        let mut r = Request::new(Molecule::Protein);
        r.text = "MKWVTFISLL*BXZUO".into();
        let mut d = r.create().unwrap();
        let before = d.clone();
        assert!(d.revcomp().is_err());
        assert_eq!(d, before);
        d.replace(0, 1, "X").unwrap();
        assert!(d.find("TFI").unwrap().len() == 1);
        r.circular = true;
        assert!(r.create().is_err());
    }
    #[test]
    fn placeholders_and_unknown_sequences_are_explicit() {
        for molecule in [Molecule::Dna, Molecule::Rna, Molecule::Protein] {
            let mut r = Request::new(molecule);
            assert!(r.create().is_err());
            r.mode = Mode::Placeholder;
            let mut d = r.create().unwrap();
            assert!(d.placeholder && d.sequence.is_empty());
            d.replace(0, 0, "ACG").unwrap();
            assert!(!d.placeholder);
            r.mode = Mode::Unknown;
            r.count = "12".into();
            let d = r.create().unwrap();
            assert_eq!(d.sequence.len(), 12);
            assert!(d.sequence.chars().all(|c| c
                == if molecule == Molecule::Protein {
                    'X'
                } else {
                    'N'
                }));
            for n in ["0", "-1", "10000001", "oops"] {
                r.count = n.into();
                assert!(r.create().is_err());
            }
        }
    }
    #[test]
    fn pasted_fasta_is_one_record_and_invalid_text_is_rejected() {
        assert_eq!(
            parse_text(">example\nacgt n", Molecule::Dna, false, true).unwrap(),
            "ACGTN"
        );
        for input in [">a\nACGT\n>b\nAG", "ACGT!", "ACG-T", "LOCUS x"] {
            assert!(parse_text(input, Molecule::Dna, false, true).is_err());
        }
        assert!(parse_text("1 ACG", Molecule::Dna, false, false).is_err());
        assert!(parse_text("AU", Molecule::Dna, false, true).is_err());
        assert_eq!(parse_text("AU", Molecule::Dna, true, true).unwrap(), "AT");
    }
}

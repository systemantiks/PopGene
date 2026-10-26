use super::Node;
use popgene_core::metadata::Reference;

pub(super) fn read(notes: &Node) -> Vec<Reference> {
    notes
        .child("References")
        .into_iter()
        .flat_map(|refs| &refs.children)
        .filter(|node| node.name == "Reference")
        .map(|node| {
            let opt = |key: &str| node.attrs.get(key).cloned();
            Reference {
                title: opt("title").unwrap_or_default(),
                authors: opt("authors"),
                journal: opt("journal"),
                pubmed: opt("pubMedID"),
                extra: node
                    .attrs
                    .iter()
                    .filter(|(key, _)| {
                        !matches!(key.as_str(), "title" | "authors" | "journal" | "pubMedID")
                    })
                    .map(|(key, value)| (key.clone(), value.clone()))
                    .collect(),
                ..Default::default()
            }
        })
        .collect()
}

pub(super) fn write(refs: &[Reference]) -> Node {
    Node {
        name: "References".into(),
        children: refs
            .iter()
            .map(|r| {
                let mut node = Node {
                    name: "Reference".into(),
                    ..Default::default()
                };
                for (key, value) in &r.extra {
                    // These attributes were observed in original reference records.
                    // Other source formats' attributes are not guessed into native XML.
                    if matches!(
                        key.as_str(),
                        "type"
                            | "date"
                            | "issue"
                            | "journalName"
                            | "pages"
                            | "volume"
                            | "doi"
                            | "patentNumber"
                            | "patentApplicationNumber"
                            | "country"
                            | "assignees"
                            | "inventors"
                    ) {
                        node.attrs.insert(key.clone(), value.clone());
                    }
                }
                node.attrs
                    .entry("type".into())
                    .or_insert_with(|| "Journal Article".into());
                node.attrs.insert("title".into(), r.title.clone());
                for (key, value) in [
                    ("authors", &r.authors),
                    ("journal", &r.journal),
                    ("pubMedID", &r.pubmed),
                ] {
                    if let Some(value) = value {
                        node.attrs.insert(key.into(), value.clone());
                    }
                }
                node
            })
            .collect(),
        ..Default::default()
    }
}

pub(super) fn supplemental(refs: &[Reference]) -> Option<String> {
    let mut lines = vec![];
    for (index, r) in refs.iter().enumerate() {
        for (key, value) in [
            ("Consortium", &r.consortium),
            ("Reference remark", &r.remark),
            ("Original reference scope", &r.scope),
        ] {
            if let Some(value) = value {
                lines.push(format!("Reference {} — {key}: {value}", index + 1));
            }
        }
        for (key, value) in &r.extra {
            if !matches!(
                key.as_str(),
                "type"
                    | "date"
                    | "issue"
                    | "journalName"
                    | "pages"
                    | "volume"
                    | "doi"
                    | "patentNumber"
                    | "patentApplicationNumber"
                    | "country"
                    | "assignees"
                    | "inventors"
            ) {
                lines.push(format!("Reference {} — {key}: {value}", index + 1));
            }
        }
    }
    (!lines.is_empty()).then(|| lines.join("\n"))
}

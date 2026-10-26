//! Context digests must not depend on Cargo's serde_json feature unification.
use crate::Result;
use serde::Serialize;
use serde_json::Value;
use sha2::{Digest, Sha256};

pub(super) struct ContextHash {
    pub canonical: String,
    legacy_insertion: String,
}

impl ContextHash {
    pub fn matches(&self, expected: Option<&Value>) -> bool {
        expected
            .and_then(Value::as_str)
            .is_some_and(|expected| expected == self.canonical || expected == self.legacy_insertion)
    }

    #[cfg(test)]
    pub fn legacy_insertion(&self) -> &str {
        &self.legacy_insertion
    }
}

/// Callers provide typed structs in the exact historical field order. Direct
/// struct serialization reproduces the old preserve_order build's context on
/// every build; it is a read-only compatibility alternative, not a new writer.
/// Recursive sorting then defines the one canonical digest for all new writes.
pub(super) fn context<T: Serialize>(typed: &T) -> Result<ContextHash> {
    let legacy_insertion = format!("{:x}", Sha256::digest(serde_json::to_vec(typed)?));
    let mut value = serde_json::to_value(typed)?;
    value.sort_all_objects();
    let canonical = format!("{:x}", Sha256::digest(serde_json::to_vec(&value)?));
    Ok(ContextHash {
        canonical,
        legacy_insertion,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Serialize)]
    struct Nested<'a> {
        z: &'a str,
        a: Option<u8>,
    }
    #[derive(Serialize)]
    struct Context<'a> {
        z: Vec<Nested<'a>>,
        a: bool,
    }

    #[test]
    fn context_recursively_sorts_objects_but_preserves_arrays_and_scalars() {
        let input = Context {
            z: vec![
                Nested { z: "α\n", a: None },
                Nested {
                    z: "first",
                    a: Some(7),
                },
            ],
            a: false,
        };
        let hashes = context(&input).unwrap();
        let expected = br#"{"a":false,"z":[{"a":null,"z":"\u03b1\n"},{"a":7,"z":"first"}]}"#;
        // The literal above deliberately escapes alpha; serde_json emits its
        // UTF-8 bytes. Compare independently specified canonical bytes instead.
        let expected = std::str::from_utf8(expected)
            .unwrap()
            .replace("\\u03b1", "α");
        assert_eq!(
            hashes.canonical,
            format!("{:x}", Sha256::digest(expected.as_bytes()))
        );
        assert_ne!(hashes.canonical, hashes.legacy_insertion);
        assert!(hashes.matches(Some(&Value::from(hashes.canonical.as_str()))));
        assert!(hashes.matches(Some(&Value::from(hashes.legacy_insertion.as_str()))));
        assert!(!hashes.matches(Some(&Value::from("0".repeat(64)))));
        assert!(!hashes.matches(Some(&Value::Null)));
        assert!(!hashes.matches(None));
        let reversed = Context {
            z: input.z.into_iter().rev().collect(),
            a: input.a,
        };
        assert_ne!(hashes.canonical, context(&reversed).unwrap().canonical);
    }
}

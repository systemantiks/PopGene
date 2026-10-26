//! Read-only project directory inventory. Files are parsed lazily as tabs open.
use std::{
    fs,
    path::{Path, PathBuf},
};
#[derive(Clone, Debug)]
pub struct Entry {
    pub path: PathBuf,
    pub name: String,
    pub directory: bool,
    pub supported: bool,
    pub children: Vec<Entry>,
}
#[derive(Debug)]
pub struct Project {
    pub root: PathBuf,
    pub entries: Vec<Entry>,
    pub warnings: Vec<String>,
    pub files: usize,
}
pub fn supported(path: &Path) -> bool {
    super::SaveFormat::from_path(path).is_some()
        || path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("popgene"))
}
pub fn scan(root: &Path) -> std::io::Result<Project> {
    let root = fs::canonicalize(root)?;
    if !root.is_dir() {
        return Err(std::io::Error::other("Choose a project folder"));
    }
    let mut p = Project {
        root: root.clone(),
        entries: vec![],
        warnings: vec![],
        files: 0,
    };
    let mut visited = 0;
    p.entries = walk(&root, 0, &mut visited, &mut p)?;
    Ok(p)
}
fn walk(
    path: &Path,
    depth: usize,
    visited: &mut usize,
    p: &mut Project,
) -> std::io::Result<Vec<Entry>> {
    if depth > 32 {
        p.warnings
            .push(format!("Depth limit reached: {}", path.display()));
        return Ok(vec![]);
    }
    let mut result = vec![];
    for item in fs::read_dir(path)? {
        if *visited >= 5000 {
            p.warnings.push(
                "Project inventory capped at 5,000 entries; narrow the selected folder.".into(),
            );
            break;
        }
        let entry = match item {
            Ok(e) => e,
            Err(e) => {
                p.warnings.push(e.to_string());
                continue;
            }
        };
        let name = entry.file_name().to_string_lossy().into_owned();
        if name.starts_with('.') {
            continue;
        }
        *visited += 1;
        let path = entry.path();
        let kind = match entry.file_type() {
            Ok(k) => k,
            Err(e) => {
                p.warnings.push(format!("{name}: {e}"));
                continue;
            }
        };
        if kind.is_symlink() {
            p.warnings
                .push(format!("Link not followed: {}", path.display()));
            continue;
        }
        let directory = kind.is_dir();
        let children = if directory {
            match walk(&path, depth + 1, visited, p) {
                Ok(c) => c,
                Err(e) => {
                    p.warnings.push(format!("{}: {e}", path.display()));
                    vec![]
                }
            }
        } else {
            vec![]
        };
        if !directory {
            p.files += 1;
        }
        result.push(Entry {
            supported: kind.is_file() && supported(&path),
            path,
            name,
            directory,
            children,
        });
    }
    result.sort_by(|a, b| {
        b.directory
            .cmp(&a.directory)
            .then_with(|| a.name.to_lowercase().cmp(&b.name.to_lowercase()))
    });
    Ok(result)
}
impl Project {
    pub fn supported_files(&self) -> Vec<PathBuf> {
        fn collect(entries: &[Entry], out: &mut Vec<PathBuf>) {
            for e in entries {
                if e.supported {
                    out.push(e.path.clone());
                }
                collect(&e.children, out);
            }
        }
        let mut out = vec![];
        collect(&self.entries, &mut out);
        out
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nested_project_retains_unsupported_entries() {
        let t = tempfile::tempdir().unwrap();
        fs::create_dir(t.path().join("DNA Files")).unwrap();
        fs::write(t.path().join("DNA Files/one.dna"), b"").unwrap();
        fs::write(t.path().join("readme.pdf"), b"").unwrap();
        fs::write(t.path().join(".sglock"), b"").unwrap();
        let p = scan(t.path()).unwrap();
        assert_eq!(p.files, 2);
        assert_eq!(p.supported_files().len(), 1);
        assert_eq!(p.entries[0].name, "DNA Files");
    }
}

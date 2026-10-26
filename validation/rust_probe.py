"""Build validation probes with the same toolchain and artifacts as cargo-local."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def tool_environment():
    env = dict(os.environ)
    local = ROOT / ".tools/cargo/bin/cargo"
    if local.is_file() and os.access(local, os.X_OK):
        env["CARGO_HOME"] = str(ROOT / ".tools/cargo")
        env["RUSTUP_HOME"] = str(ROOT / ".tools/rustup")
        env["PATH"] = str(local.parent) + os.pathsep + env.get("PATH", "")
    return env


def build(packages, example=None):
    env = tool_environment()
    cargo = shutil.which("cargo", path=env.get("PATH"))
    if cargo is None:
        raise RuntimeError("Install Rust/Cargo or configure the optional .tools toolchain")
    command = [cargo, "build", "--locked", "--offline", "--message-format=json"]
    for package in packages:
        command += ["-p", package]
    if example:
        command += ["--example", example]
    result = subprocess.run(command, cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE)
    artifacts = [json.loads(line) for line in result.stdout.splitlines() if line.startswith("{")]
    for artifact in artifacts:
        if artifact.get("reason") == "compiler-message" and artifact["message"].get("rendered"):
            print(artifact["message"]["rendered"], file=sys.stderr, end="")
    result.check_returncode()
    return artifacts


def build_example(package, name):
    for artifact in build([package], name):
        if artifact.get("reason") == "compiler-artifact" and artifact["target"]["name"] == name:
            if artifact.get("executable"):
                return Path(artifact["executable"])
    raise RuntimeError(f"Cargo did not produce example {name}")


def compile_probe(source, binary, packages, optimize=False):
    artifacts = build(packages)
    env = tool_environment()
    rustc = env.get("RUSTC") or shutil.which("rustc", path=env.get("PATH"))
    if rustc is None:
        raise RuntimeError("The selected Cargo toolchain has no rustc compiler")
    command = [rustc, "--edition=2024", str(source), "-o", str(binary)]
    if optimize:
        command.append("-O")
    for package in packages:
        crate = package.replace("-", "_")
        library = next((Path(filename) for artifact in artifacts
                        if artifact.get("reason") == "compiler-artifact"
                        and artifact["target"]["name"] == crate
                        for filename in artifact["filenames"] if filename.endswith(".rlib")), None)
        if library is None:
            raise RuntimeError(f"Cargo did not produce the {crate} library")
        dependencies = library.parent if library.parent.name == "deps" else library.parent / "deps"
        command += ["--extern", f"{crate}={library}", "-L", f"dependency={dependencies}"]
    subprocess.run(command, cwd=ROOT, env=env, check=True)

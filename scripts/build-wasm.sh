#!/bin/sh
# Builds the npm package nec2-wasm in target/npm/nec2-wasm.
# It needs the Rust target wasm32-unknown-unknown and the program
# wasm-bindgen in the version of the wasm-bindgen crate in Cargo.lock:
#   rustup target add wasm32-unknown-unknown
#   cargo install wasm-bindgen-cli --version <that version> --locked
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
out=target/npm/nec2-wasm
version=$(sed -n '/^\[workspace.package\]/,/^\[/s/^version = "\(.*\)"/\1/p' Cargo.toml)

cargo build --release -p nec2-wasm --target wasm32-unknown-unknown
rm -rf "$out"
wasm-bindgen --target web --out-dir "$out" --out-name nec2 \
  target/wasm32-unknown-unknown/release/nec2_wasm.wasm
cp crates/nec2-wasm/README.md LICENSE "$out"/
cat > "$out/package.json" <<JSON
{
  "name": "nec2-wasm",
  "version": "$version",
  "description": "NEC-2 antenna modeling in WebAssembly: the method of moments for wire antennas",
  "keywords": ["nec2", "antenna", "electromagnetics", "method-of-moments", "wasm"],
  "license": "MIT",
  "repository": { "type": "git", "url": "git+https://github.com/CyberDev-space/nec2-rs.git" },
  "homepage": "https://github.com/CyberDev-space/nec2-rs/tree/main/crates/nec2-wasm",
  "type": "module",
  "main": "nec2.js",
  "types": "nec2.d.ts",
  "exports": {
    ".": { "types": "./nec2.d.ts", "default": "./nec2.js" },
    "./nec2_bg.wasm": "./nec2_bg.wasm"
  },
  "files": ["nec2.js", "nec2.d.ts", "nec2_bg.wasm", "nec2_bg.wasm.d.ts", "README.md", "LICENSE"],
  "sideEffects": false
}
JSON
echo "$out: nec2-wasm $version"
ls -l "$out"

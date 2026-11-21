#!/bin/sh
# Fetch nec2c 1.3.1 (public domain) from the Debian archive and build the two
# reference programs in reference/bin:
#   nec2c        the pristine sources
#   nec2c-fixed  the sources with scripts/patch-nec2c.py applied
# The sources are compiled directly; the build scripts of the archive are not run.
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
ref=$root/reference
archive=nec2c_1.3.1.orig.tar.bz2
sha256=271bd2ea6037f896912139a66e3abb75981c72315bab5221fcf384a1736d0337
sources="main.c calculations.c fields.c geometry.c ground.c input.c matrix.c misc.c network.c radiation.c shared.c somnec.c"

mkdir -p "$ref/dl" "$ref/bin"
if [ ! -f "$ref/dl/$archive" ]; then
  curl -fsSL -o "$ref/dl/$archive" "https://deb.debian.org/debian/pool/main/n/nec2c/$archive"
fi
echo "$sha256  $ref/dl/$archive" | sha256sum -c -

if [ ! -d "$ref/nec2c-1.3.1" ]; then
  tar -xjf "$ref/dl/$archive" -C "$ref"
fi
(cd "$ref/nec2c-1.3.1" && cc -O2 -w -DPACKAGE_STRING='"nec2c 1.3.1"' $sources -lm -o "$ref/bin/nec2c")

if [ ! -d "$ref/nec2c-1.3.1-fixed" ]; then
  cp -r "$ref/nec2c-1.3.1" "$ref/nec2c-1.3.1-fixed"
  python3 "$root/scripts/patch-nec2c.py" "$ref/nec2c-1.3.1-fixed"
fi
(cd "$ref/nec2c-1.3.1-fixed" && cc -O2 -w -DPACKAGE_STRING='"nec2c 1.3.1 with fixes"' $sources -lm -o "$ref/bin/nec2c-fixed")
echo "built $ref/bin/nec2c and $ref/bin/nec2c-fixed"

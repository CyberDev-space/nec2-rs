#!/bin/sh
# Regenerate the reference reports for every deck in tests/fixtures/decks:
#   tests/fixtures/nec2c        from the pristine nec2c 1.3.1
#   tests/fixtures/nec2c-fixed  from nec2c 1.3.1 with scripts/patch-nec2c.py applied
# Next to a report NAME.out there is NAME.plt if the deck asks for a plot file.
# Build the two binaries first with scripts/build-reference.sh.
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
# nec2c rejects file names longer than 75 characters, so it runs with relative paths.
cd "$root/tests/fixtures"
status=0
for variant in nec2c nec2c-fixed; do
  mkdir -p "$variant"
  for deck in decks/*.nec; do
    name=$(basename "$deck" .nec)
    if ! "$root/reference/bin/$variant" -i "$deck" -o "$variant/$name.out" >"$variant/$name.log" 2>&1; then
      # Decks named err_* are meant to stop with an error. The stock nec2c
      # crashes on the asymmetry check without networks and refuses a strip
      # of patches.
      case $variant:$name in
        *:err_*) ;;
        nec2c:s6_asymmetry) ;;
        nec2c:s8_patch_strip) ;;
        *) echo "FAILED ($variant): $name" >&2; status=1 ;;
      esac
    fi
    # The log is kept only when nec2c printed something to the terminal.
    [ -s "$variant/$name.log" ] || rm -f "$variant/$name.log"
    # A deck with a PL card: nec2c writes the plot file next to the deck.
    [ -f "$deck.plt" ] && mv "$deck.plt" "$variant/$name.plt"
  done
done
exit $status

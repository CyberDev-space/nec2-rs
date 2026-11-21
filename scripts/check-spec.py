#!/usr/bin/env python3
"""Mechanical check of a clean-room specification before the clean side reads it.

Usage: check-spec.py <cleanroom.md>

The script looks for text that reads like C or C++ source and reports only
line numbers and the kind of marker, never the text itself, so that running
it does not expose the content. Blocks that are NEC decks (lines that start
with a two-letter card mnemonic) and tables of numbers are accepted.

Exit status 0 means that nothing suspicious was found.
"""
import re
import sys

CARDS = {
    "CM", "CE", "GW", "GC", "GA", "GH", "GM", "GR", "GS", "GX", "GE", "GF", "SP", "SM", "SC",
    "FR", "LD", "GN", "GD", "EX", "NT", "TL", "XQ", "RP", "NX", "PT", "PQ", "KH", "NE", "NH",
    "EK", "CP", "PL", "EN", "WG", "XT",
}

MARKERS = [
    ("statement ending in ';'", re.compile(r";\s*$")),
    ("'->' operator", re.compile(r"->")),
    ("'::' scope operator", re.compile(r"\w::\w")),
    ("preprocessor directive", re.compile(r"^\s*#\s*(include|define|ifdef|ifndef|endif|pragma)\b")),
    ("brace", re.compile(r"[{}]")),
    ("C keyword with parenthesis", re.compile(r"\b(for|while|switch|sizeof)\s*\(.*\)\s*$")),
    ("increment or compound assignment", re.compile(r"(\+\+|--(?!-)|\+=|-=|\*=|/=)\s*\w|\w\s*(\+\+|\+=|-=|\*=|/=)")),
    ("pointer or address expression", re.compile(r"(\*\w+\s*=[^=])|(&\w+\[)")),
    ("C comment", re.compile(r"/\*|\*/|^\s*//")),
]


def is_deck_line(line):
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return True
    return stripped[:2].upper() in CARDS and (len(stripped) == 2 or stripped[2] in " \t,")


def is_numeric_line(line):
    tokens = line.replace("|", " ").replace(",", " ").split()
    numeric = sum(bool(re.fullmatch(r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?[j%]?", t)) for t in tokens)
    return bool(tokens) and numeric >= max(1, len(tokens) // 2)


def main(path):
    findings = []
    in_fence = False
    fence_start = 0
    fence_lines = []
    with open(path, encoding="utf-8") as spec:
        for number, line in enumerate(spec, 1):
            line = line.rstrip("\n")
            if line.lstrip().startswith("```"):
                if in_fence:
                    other = [n for n, text in fence_lines if not (is_deck_line(text) or is_numeric_line(text))]
                    if other:
                        findings.append((fence_start, f"fenced block with {len(other)} lines that are neither cards nor numbers"))
                else:
                    fence_start = number
                    fence_lines = []
                in_fence = not in_fence
                continue
            if in_fence:
                fence_lines.append((number, line))
                if is_deck_line(line) or is_numeric_line(line):
                    continue
            # Program and language names are not increments.
            checked = re.sub(r"(?i)\b(nec2|nec-2|c)\+\+", "", line)
            for name, pattern in MARKERS:
                if pattern.search(checked):
                    findings.append((number, name))
    for number, name in findings:
        print(f"line {number}: {name}")
    print(f"{len(findings)} findings")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

#!/usr/bin/env bash
# =============================================================================
# build_methods_pdf.sh — render docs/methods.md to docs/methods.pdf (pandoc + xelatex).
#
# The markdown keeps Unicode superscripts (10⁻⁴) so it reads well on GitHub; Times New Roman lacks
# those glyphs, so they are rewritten to pandoc super/subscript syntax for the PDF build only.
# Usage:  scripts/build_methods_pdf.sh [--font "Times New Roman"]
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../docs/methods.md"; OUT="$HERE/../docs/methods.pdf"
FONT="Times New Roman"
[[ "${1:-}" == "--font" ]] && FONT="$2"
command -v pandoc >/dev/null || { echo "pandoc not found" >&2; exit 1; }
tmpd="$(mktemp -d -t methods)"; tmp="$tmpd/methods.md"; trap 'rm -rf "$tmpd"' EXIT
# ⁻⁴ etc. -> ^−4^ (pandoc superscript, U+2212 minus); ᵥ -> ~v~ (pandoc subscript)
sed -e 's/⁻⁴/^−4^/g; s/⁻⁵/^−5^/g; s/⁻⁶/^−6^/g; s/⁻³/^−3^/g; s/⁻²/^−2^/g; s/ᵥ/~v~/g' "$SRC" > "$tmp"
pandoc "$tmp" -s --shift-heading-level-by=-1 --pdf-engine=xelatex \
    -V geometry:margin=1in -V mainfont="$FONT" -V monofont="Menlo" -V fontsize=11pt \
    -V colorlinks=true -V linkcolor=blue -V urlcolor=blue \
    -M author="hprv pipeline — github.com/jlanej/high_priority_rare_variant" \
    -M date="$(git -C "$HERE" log -1 --format=%cs 2>/dev/null || date +%F)" \
    -o "$OUT"
echo "wrote $OUT"

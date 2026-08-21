#!/usr/bin/env python3
"""Doxygen input filter : Markdown dollar-math -> Doxygen native math.

Converts `$...$` (inline) and `$$...$$` (display) LaTeX equations into
Doxygen's `\\f$...\\f$` and `\\f[...\\f]` commands, so that:

  1. formula bodies are protected from Doxygen's Markdown pass (a raw `$c_{mix}$`
     would have its underscore parsed as emphasis and reach MathJax broken);
  2. rendering goes through the stock Doxygen + MathJax machinery
     (USE_MATHJAX = YES in doc/Doxyfile), no client-side delimiter hacks.

This mirrors the convention of the user documentation (Hugo + MathJax), so the
same equation markup works in both doc trees. Applied to *.md files only, via
FILTER_PATTERNS in doc/Doxyfile — C/C++ comments keep using \\f$ natively
(and their `$MIN/$MAX/$DEFAULT` introspection tokens are never touched).

Rules, matching the usual Markdown-math conventions:
  - `$$...$$` may span several lines; `$...$` must stay on one line;
  - the opening `$` must be immediately followed by a non-space, the closing
    `$` immediately preceded by one (so "5$ and 10$" or "$5 vs $10" are left
    alone);
  - `\\$` is an escaped, literal dollar;
  - fenced code blocks (``` / ~~~) and inline code spans (`...`) are never
    converted.
"""

import re
import sys

# display math, possibly multiline ; escaped \$ do not delimit
_DISPLAY = re.compile(r"(?<!\\)\$\$(.+?)(?<!\\)\$\$", re.S)
# inline math, single line, snug delimiters ; body may contain escaped chars ;
# the closing $ must be preceded by a non-space that is not a backslash, so an
# escaped \$ can never close a formula (backtracking would otherwise allow it)
_INLINE = re.compile(r"(?<![\\$])\$(?=\S)((?:\\.|[^$\n])+?)(?<=\S)(?<!\\)\$(?!\$)")
# fenced code blocks, indented code blocks (4 spaces / tab, the old Markdown
# style), then inline code spans, all kept verbatim
_FENCED = re.compile(r"(```.*?```|~~~.*?~~~)", re.S)
_INDENTED = re.compile(r"((?:^(?:[ ]{4,}|\t).*\n?)+)", re.M)
_CODE_SPAN = re.compile(r"(`[^`\n]*`)")


def _convert_prose(text):
    text = _DISPLAY.sub(r"\\f[\1\\f]", text)
    return _INLINE.sub(r"\\f$\1\\f$", text)


def convert(text):
    out = []
    for i, block in enumerate(_FENCED.split(text)):
        if i % 2 == 1:  # fenced code : verbatim
            out.append(block)
            continue
        for j, chunk in enumerate(_INDENTED.split(block)):
            if j % 2 == 1:  # indented code : verbatim
                out.append(chunk)
                continue
            for k, span in enumerate(_CODE_SPAN.split(chunk)):
                out.append(span if k % 2 == 1 else _convert_prose(span))
    return "".join(out)


if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8", errors="replace") as f:
        sys.stdout.write(convert(f.read()))

#!/usr/bin/env python3
"""Check markdown for constructs that break pandoc -> pdflatex.

Pandoc's inline-math rules (a violation makes both dollars literal, which then
emits bare LaTeX commands outside math mode and fails to compile):
  1. an opening $ must not be followed by whitespace
  2. a closing $ must not be preceded by whitespace
  3. a closing $ must not be followed immediately by a digit
Plus: non-ASCII beyond em dash, en dash and the section sign.
"""
import re, sys

OK_NONASCII = {'—', '–', '§'}

def strip_code(text):
    """Blank out anything pandoc will not parse as markdown, preserving offsets.

    HTML comments included: pandoc passes them through as raw HTML and the LaTeX
    writer drops them, so their contents cannot break a render.
    """
    out = list(text)
    for m in re.finditer(r'<!--.*?-->', text, re.S):
        for i in range(m.start(), m.end()):
            if out[i] != '\n': out[i] = ' '

    for m in re.finditer(r'^```.*?^```', text, re.S | re.M):
        for i in range(m.start(), m.end()):
            if out[i] != '\n': out[i] = ' '
    for m in re.finditer(r'`[^`\n]*`', ''.join(out)):
        for i in range(m.start(), m.end()):
            if out[i] != '\n': out[i] = ' '
    return ''.join(out)

def check(path):
    raw = open(path, encoding='utf-8').read()
    body = strip_code(raw)
    problems = []

    for i, ch in enumerate(raw):
        if ord(ch) > 127 and ch not in OK_NONASCII:
            line = raw.count('\n', 0, i) + 1
            problems.append((line, f"non-ASCII U+{ord(ch):04X} {ch!r}"))

    # Pair dollars left to right, ignoring escaped ones.
    pos, opens = 0, []
    while pos < len(body):
        if body[pos] == '$' and (pos == 0 or body[pos-1] != '\\'):
            opens.append(pos)
        pos += 1
    for a, b in zip(opens[0::2], opens[1::2]):
        line = raw.count('\n', 0, a) + 1
        span = raw[a:b+1]
        if len(span) < 3 or '\n\n' in span:
            continue
        if body[a+1].isspace():
            problems.append((line, f"opening $ followed by space: {span[:40]!r}"))
        if body[b-1].isspace():
            problems.append((line, f"closing $ preceded by space: {span[:40]!r}"))
        if b + 1 < len(body) and body[b+1].isdigit():
            problems.append((line, f"closing $ followed by digit: {raw[a:b+6]!r}"))
    return problems

bad = 0
for path in sys.argv[1:]:
    ps = check(path)
    print(f"{path}: {'OK' if not ps else str(len(ps)) + ' problem(s)'}")
    for line, msg in sorted(set(ps)):
        print(f"  L{line}: {msg}")
    bad += len(ps)
sys.exit(1 if bad else 0)

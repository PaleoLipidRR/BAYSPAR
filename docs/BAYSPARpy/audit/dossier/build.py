"""Assemble dossier.html from template.html, the real source files, and data.json.

Run gen_data.py first. Code excerpts come from the MATLAB files at the repository root, from
BAYSPARpy/src/baysparpy, and from the installed brews/baysparpy (`import bayspar`).

Stage blocks in the template look like

    <!--stage
    id: std-thin
    title: Thin the draws
    note: <p>...</p>
    M: bayspar_tex.m:84-89 | ref | caption
    P: bayspar_tex.py:75 ; utils.py:55-59 | same | caption
    B: predict.py:178-185:hl=180 | diff | caption
    -->

and code ranges are pulled from the files with their real line numbers.
"""
import html, importlib.util, re
from pathlib import Path

from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import MatlabLexer, PythonLexer

ROOT = Path(__file__).resolve().parents[4]          # the BAYSPAR repository
HERE = Path(__file__).parent
SRC = {
    "M": (ROOT, "MATLAB"),
    "P": (ROOT / "BAYSPARpy/src/baysparpy", "BAYSPARpy"),
    "B": (Path(importlib.util.find_spec("bayspar").origin).parent, "brews/baysparpy"),
}
SHOWN_PATH = {"M": "{}", "P": "baysparpy/{}", "B": "bayspar/{}"}
STATUS = {
    "ref": "Reference", "same": "Identical math", "equiv": "Same result",
    "changed": "Changed on purpose", "diff": "Differs", "defect": "Defect",
    "missing": "Not implemented", "fixed": "Fixes the defect",
}
FMT = HtmlFormatter(nowrap=True)


def ranges(spec):
    out = []
    for part in spec.split("+"):
        a, _, b = part.partition("-")
        out.append((int(a), int(b or a)))
    return out


def code_block(key, spec):
    """One file excerpt: 'file:ranges[:hl=a-b+c]'."""
    bits = spec.strip().split(":")
    fname, rng = bits[0], bits[1]
    hl = set()
    for b in bits[2:]:
        if b.startswith("hl="):
            for a, z in ranges(b[3:]):
                hl.update(range(a, z + 1))
    base, _ = SRC[key]
    lines = (base / fname).read_text().splitlines()
    lexer = MatlabLexer(stripnl=False) if fname.endswith(".m") else PythonLexer(stripnl=False)
    chunks = []
    rs = ranges(rng)
    for a, z in rs:
        if not (1 <= a <= z <= len(lines)):
            raise ValueError(f"{key}:{fname} has {len(lines)} lines; asked for {a}-{z}")
    # dedent across every range of this excerpt so they line up
    sel = [lines[i - 1] for a, z in rs for i in range(a, z + 1)]
    indent = min((len(l) - len(l.lstrip()) for l in sel if l.strip()), default=0)
    for k, (a, z) in enumerate(rs):
        if k:
            chunks.append('<span class="gap" aria-label="lines omitted">⋯</span>')
        raw = "\n".join(l[indent:] for l in lines[a - 1:z]) + "\n"
        out = highlight(raw, lexer, FMT).rstrip("\n").split("\n")
        # pygments drops leading blank lines; pad back so numbering stays true
        out += [""] * ((z - a + 1) - len(out))
        for n, h in zip(range(a, z + 1), out):
            cls = "ln hl" if n in hl else "ln"
            chunks.append(f'<span class="{cls}" data-n="{n}">{h or " "}</span>')
    label = (f"L{rs[0][0]}" if rs[0][0] == rs[-1][1] and len(rs) == 1
             else ", ".join(f"L{a}–{z}" if a != z else f"L{a}" for a, z in rs))
    return (f'<figure class="code"><figcaption><span class="path">'
            f'{html.escape(SHOWN_PATH[key].format(fname))}</span><span class="lines">{label}</span>'
            f'</figcaption><pre><code>{"".join(chunks)}</code></pre></figure>')


def column(key, spec):
    parts = [p.strip() for p in spec.split("|")]
    codes, status, caption = parts[0], parts[1], parts[2] if len(parts) > 2 else ""
    _, name = SRC[key]
    body = ""
    if codes and codes != "none":
        body = "".join(code_block(key, c) for c in codes.split(";"))
    else:
        body = '<div class="nocode">No equivalent code</div>'
    cap = f'<p class="cap">{caption}</p>' if caption else ""
    return (f'<div class="col" data-impl="{key}"><div class="colhead"><span class="who">{name}</span>'
            f'<span class="st st-{status}">{STATUS[status]}</span></div>{body}{cap}</div>')


def stage(block):
    f = {}
    for line in block.strip().splitlines():
        k, _, v = line.partition(":")
        f[k.strip()] = v.strip()
    return (f'<article class="stage" id="{f["id"]}"><header class="stagehead">'
            f'<h3>{f["title"]}</h3></header>'
            f'<div class="stagenote">{f.get("note", "")}</div>'
            f'<div class="tri">{column("M", f["M"])}{column("P", f["P"])}{column("B", f["B"])}</div>'
            f'</article>')


tpl = (HERE / "template.html").read_text()
tpl = re.sub(r"<!--stage\n(.*?)-->", lambda m: stage(m.group(1)), tpl, flags=re.S)
tpl = re.sub(r"<!--code:([MPB]):(.*?)-->", lambda m: code_block(m.group(1), m.group(2)), tpl)
tpl = tpl.replace("/*DATA*/null", (HERE / "data.json").read_text())
(HERE / "dossier.html").write_text(tpl)
print("dossier.html", len(tpl) // 1024, "KB")

# CLAUDE.md

> **What is this file?** It is a briefing for AI coding assistants such as [Claude Code](https://docs.anthropic.com/en/docs/claude-code), which read it automatically when they work in this repository.

## This repository

A fork of Jess Tierney's BAYSPAR (MATLAB; Tierney & Tingley, 2014, 2015) that also holds BAYSPARpy, Ronnakrit Rattanasriampaipong's Python recode, in `BAYSPARpy/` (package `bayspar`, imported as `baysparpy`; not yet released as a package). The recode's design notes, porting record and audit material are in `docs/BAYSPARpy/`. Its tests read the calibration files in `ModelOutput/` and the saved reference output in `docs/BAYSPARpy/audit/`.

```bash
cd BAYSPARpy && uv sync && uv run python -m pytest     # 83 passed, 2 skipped
```

## Conventions

- Write every paragraph of prose as one line, with no line breaks inside it. Break lines only between paragraphs, list items, table rows, headings and code blocks. This applies to Markdown files, notebook markdown cells, and PR and issue descriptions; code comments, docstrings and commit messages keep their usual line-length conventions.
- Leave the upstream MATLAB files and `ReadMe.md` as Jess Tierney wrote them; changes for the recode go in `BAYSPARpy/` and `docs/BAYSPARpy/`.
- `docs/BAYSPARpy/PORTING.md` is checked by `tests/test_porting_doc.py`: keep its `### ID · STATUS · title` headings and register-table rows intact.
- After changing a ported function, rebuild and recommit `docs/BAYSPARpy/audit/dossier/dossier.html` (see that folder's README).

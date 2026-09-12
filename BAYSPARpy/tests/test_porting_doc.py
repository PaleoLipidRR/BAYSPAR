"""The record and the code must move together.

Asserts the correspondence promised in docs/BAYSPARpy/PORTING.md: every entry
has a test, every test has an entry, and every behavioural difference is in the
register a reviewer actually reads.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

TESTS = Path(__file__).parent
ENTRY = re.compile(r"^### ([A-Z]+-\d\d) · ([A-Z]+) · (.+)$", re.M)
REGISTER_ROW = re.compile(r"^\| \[`([A-Z]+-\d\d)`\]\(#[a-z]+-\d\d\) \| ([A-Z]+) \|", re.M)
TEST_DEF = re.compile(r"^def (test_[A-Z]+_\d\d)\(", re.M)
NEEDS_REGISTER = {"DEVIATION", "DEFECT"}


def _porting_md() -> Path:
    """Locate PORTING.md, wherever the package is checked out."""
    for parent in TESTS.resolve().parents:
        candidate = parent / "docs" / "BAYSPARpy" / "PORTING.md"
        if candidate.is_file():
            return candidate
    pytest.skip("PORTING.md not found; it lives in the BAYSPAR MATLAB repository")


@pytest.fixture(scope="module")
def doc() -> str:
    return _porting_md().read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def test_names() -> set[str]:
    return set(TEST_DEF.findall((TESTS / "test_port.py").read_text(encoding="utf-8")))


def test_every_entry_has_a_test(doc, test_names):
    """No entry may be documented without a test pinning it."""
    entries = {i for i, _, _ in ENTRY.findall(doc)}
    assert entries, "no entries parsed -- has the heading format changed?"
    missing = sorted(i for i in entries if f"test_{i.replace('-', '_')}" not in test_names)
    assert not missing, f"entries with no test in test_port.py: {missing}"


def test_every_test_has_an_entry(doc, test_names):
    """No test may claim an ID the record does not define."""
    entries = {i for i, _, _ in ENTRY.findall(doc)}
    orphans = sorted(n for n in test_names
                     if n.removeprefix("test_").replace("_", "-") not in entries)
    assert not orphans, f"tests naming an ID that PORTING.md does not define: {orphans}"


def test_deviations_are_registered(doc):
    """Every DEVIATION and DEFECT appears in the register at the top."""
    registered = {i for i, _ in REGISTER_ROW.findall(doc)}
    unregistered = sorted(i for i, status, _ in ENTRY.findall(doc)
                          if status in NEEDS_REGISTER and i not in registered)
    assert not unregistered, (
        f"behavioural differences missing from the register: {unregistered}")


def test_register_rows_are_real_entries(doc):
    """The register may not cite an ID that has no entry."""
    entries = {i for i, _, _ in ENTRY.findall(doc)}
    rows = REGISTER_ROW.findall(doc)
    assert rows, "no register rows parsed"
    ghosts = sorted(i for i, _ in rows if i not in entries)
    assert not ghosts, f"register rows with no entry: {ghosts}"


def test_register_statuses_match_their_entries(doc):
    """A row's status must match the entry it points at."""
    by_id = {i: status for i, status, _ in ENTRY.findall(doc)}
    mismatched = [(i, row_status, by_id[i]) for i, row_status in REGISTER_ROW.findall(doc)
                  if by_id.get(i) != row_status]
    assert not mismatched, f"register/entry status mismatch: {mismatched}"


def test_entry_ids_are_consecutive(doc):
    """IDs run 1..n per MATLAB file, so a gap means an entry was dropped."""
    from collections import defaultdict
    groups: dict[str, list[int]] = defaultdict(list)
    for ident, _, _ in ENTRY.findall(doc):
        prefix, num = ident.rsplit("-", 1)
        groups[prefix].append(int(num))
    for prefix, nums in groups.items():
        assert sorted(nums) == list(range(1, len(nums) + 1)), f"gap in {prefix} ids: {sorted(nums)}"

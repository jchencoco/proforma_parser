#!/usr/bin/env python3

"""

=============================================================================

Proteoform Compendium Processing Pipeline

=============================================================================
 
WHAT THIS SCRIPT DOES

---------------------

Reads an Excel or CSV file where each row describes a proteoform using a ProForma

string and a UniProt accession. For every row it:
 
  1. Parses the ProForma string to extract the amino acid sequence and all

     post-translational modifications (PTMs).

  2. Downloads the canonical protein sequence and gene metadata from

     UniProt (results are cached locally so each protein is only fetched once).

  3. Maps the proteoform sequence onto the full protein sequence and detects

     amino acid substitutions or initiator methionine removal by direct

     sequence comparison.

  4. Assigns a Proteoform Identification Level (1-5) following the Consortium

     for Top-Down Proteomics standard.
 
Fourteen new columns are added to your original data and the

result is saved as a new CSV file - your input file is never modified.
 
HOW TO RUN

----------

  python profoma_parser.py
 
That is all. Edit the "USER CONFIGURATION" section below to point at your

files and tell the script which columns contain the ProForma strings and

UniProt IDs.
 
REQUIREMENTS

------------

Install these Python packages before running (one-time setup):
 
  pip install pandas openpyxl pyteomics[proforma] requests requests_cache
 
  - pandas / openpyxl     : read CSV files and Excel input files

  - pyteomics[proforma]   : pure-Python ProForma 2.0 parser, plus the

                            `psims` dependency it uses to resolve PSI-MOD /

                            XL-MOD / GNO accessions

  - requests              : query the UniProt REST API

  - requests_cache        : cache UniProt responses on disk so you don't hit

                            rate limits when re-running the script
 
ABOUT pyteomics.proforma

--------------------------

pyteomics is a pure-Python proteomics toolkit; we use it for ProForma

parsing only: `pyteomics.proforma` parses ProForma 2.0 strings into a

structured object (sequence, terminal mods, positional mods, unlocalized

mods), and resolves modification names and formula deltas.

Mass calculation is handled locally by mash_mass.py, ported from

mash_reference/seqFunctions.cs. If an

accession genuinely can't be resolved by any

vocabulary, it is reported as a raw mass shift and flagged

identity-ambiguous (Level 2B), same as before.
 
=============================================================================

"""
 
# Makes type hints like `float | None` and `dict | None` safe to use even on

# Python versions older than 3.10 - see explanation further down where this

# actually matters.

from __future__ import annotations
 
# =============================================================================

# USER CONFIGURATION - edit this section

# =============================================================================
 
# -----------------------------------------------------------------------------
# 1) INPUT AND OUTPUT FILES
#
# Put your input file in this folder, or use a full path such as:
# r"C:\Users\jchen2422\Desktop\my_file.csv"
#
# Input can be CSV, TSV, XLS, XLSX, or XLSM.
# Output is written as both CSV and XLSX.
# -----------------------------------------------------------------------------

INPUT_FILE  = r"data/Protein ID Compendium published 9.csv"
OUTPUT_FILE = r"data/Protein ID Compendium published 9 output.csv"

# -----------------------------------------------------------------------------
# 2) INPUT COLUMN NAMES
#
# These defaults match the newer file format. Older columns ProForma and
# UniProt_ID are also accepted automatically.
# -----------------------------------------------------------------------------

PROFORMA_COLUMN = "proForma"

UNIPROT_COLUMN = "protein_accessions"

SOURCE_COLUMN = "Source"

ORGAN_COLUMN = "Organ"

EXPERIMENTAL_MONOISOTOPIC_MASS_COLUMN = "Experimental Monoisotopic Mass"

EXPERIMENTAL_AVERAGE_MASS_COLUMN = "Experimental Average Mass"

EXPERIMENTAL_MOST_ABUNDANT_MASS_COLUMN = "Experimental Most Abundant Mass"


# -----------------------------------------------------------------------------
# 3) OUTPUT COLUMN NAMES
#
# Usually you should not need to edit these.
# -----------------------------------------------------------------------------

PROFORMA_OUTPUT_COLUMN = "ProForma"

UNIPROT_OUTPUT_COLUMN = "UniProt ID"
 
# =============================================================================

# END OF USER CONFIGURATION

# =============================================================================
 
import re

import sys

from pathlib import Path

import mercury_mass
import mash_mass

SCRIPT_DIR = Path(__file__).resolve().parent

INPUT_FILE = Path(INPUT_FILE)
if not INPUT_FILE.is_absolute():
    INPUT_FILE = SCRIPT_DIR / INPUT_FILE

OUTPUT_FILE = Path(OUTPUT_FILE)
if not OUTPUT_FILE.is_absolute():
    OUTPUT_FILE = SCRIPT_DIR / OUTPUT_FILE


def _csv_output_path(path: Path) -> Path:

    return path if path.suffix.lower() == ".csv" else path.with_suffix(".csv")


OUTPUT_FILE = _csv_output_path(OUTPUT_FILE)
 
# -- third-party imports with friendly error messages -----------------------

try:

    import pandas as pd

except ImportError:

    sys.exit("[ERROR] pandas is not installed. Run: pip install pandas openpyxl")
 
try:

    from pyteomics import proforma as pf_module

except ImportError:

    sys.exit("[ERROR] pyteomics is not installed. Run: pip install pyteomics[proforma]")
 
# PSI-MOD (and XLMOD/GNO) resolution inside pyteomics.proforma depends

# entirely on the separate `psims` library - if it's missing, pyteomics

# doesn't fall back to anything, it raises an ImportError the first time a

# PSI-MOD tag is actually encountered, potentially hundreds of rows into a

# run. We check for it up front instead, since this pipeline's data relies

# on PSI-MOD tags. Unimod resolution alone does not strictly require psims.

try:

    import psims  # noqa: F401

except ImportError:

    sys.exit("[ERROR] The 'psims' package is required to resolve PSI-MOD "

              "([MOD:xxxxx]) modifications, which this dataset uses heavily.\n"

              "        Run: pip install psims")
 
try:

    import requests

except ImportError:

    sys.exit("[ERROR] requests is not installed. Run: pip install requests")
 
try:

    import requests_cache

    _SESSION = requests_cache.CachedSession(

        str(SCRIPT_DIR / "pipeline_cache"),  # pipeline_cache.sqlite

        expire_after=7 * 24 * 3600,  # cached responses expire after 7 days

    )

    print("[INFO] API responses will be cached in pipeline_cache.sqlite")

except ImportError:

    _SESSION = requests.Session()

    print("[WARN] requests_cache is not installed - API responses will not be "

          "cached. Install it with: pip install requests_cache")
 
# =============================================================================

# SECTION 1 - Modification database access
#
# pyteomics.proforma resolves ProForma tags into names and formula deltas.
# It is no longer used for mass arithmetic: monoiso and average masses are
# recalculated by mash_mass, using seqFunctions.cs as the source of truth.
#
# We only reach into pyteomics' PSI-MOD resolver for the pre-parse
# massless-term workaround below.

# =============================================================================
 
 
def _get_psimod_db():

    """The PSI-MOD database pyteomics.proforma itself uses."""

    return pf_module.PSIModModification.resolver.database


def _composition_to_formula(composition) -> dict[str, float] | None:

    """Convert a pyteomics/psims composition object to a plain formula dict.

    The composition is used only as the resolved chemical formula delta. Mass
    values are always recalculated by mash_mass.
    """

    if composition is None:

        return None

    if hasattr(composition, "items"):

        items = composition.items()

    else:

        try:

            items = dict(composition).items()

        except Exception:

            return None

    formula = {}

    for element, count in items:

        element = str(element)

        if element.startswith("isotope_string"):

            return None

        try:

            count = float(count)

        except (TypeError, ValueError):

            return None

        if count:

            formula[element] = formula.get(element, 0.0) + count

    return formula
 
 
_MOD_TAG_RE = re.compile(r'\[MOD:(\d+)\]')

_ACCESSION_TAG_SPACE_RE = re.compile(r'\[(MOD|UNIMOD|Unimod):\s*(\d+)\s*\]')


def _normalize_accession_tag_spacing(proforma_str: str) -> str:

    return _ACCESSION_TAG_SPACE_RE.sub(lambda m: f"[{m.group(1)}:{m.group(2)}]", proforma_str)
 
 
def _patch_massless_psimod(proforma_str: str):

    """

    Some PSI-MOD terms describe a *kind* of change rather than one exact

    chemistry (e.g. MOD:00677 'hydroxylated residue') and carry no mass

    value in the database at all. pyteomics needs every tag to resolve to

    a mass just to finish reading the string - so without this check, the

    WHOLE row fails to parse, even though pyteomics knows perfectly well

    what the modification is called.
 
    We look ahead for exactly this situation before handing the string to

    the real parser, and swap only the problem tag(s) for a harmless

    placeholder it can always read, remembering the real name so we can

    put it back afterward. Every other tag - including PSI-MOD tags that

    DO have a mass - is left completely untouched.
 
    Returns (patched_string, {sentinel_mass: real_name, ...})

    """

    placeholder_registry: dict = {}

    counter = [0]
 
    def _check(m: "re.Match") -> str:

        acc = m.group(1)

        try:

            term = _get_psimod_db()[f"MOD:{acc}"]

        except Exception:

            return m.group(0)  # unknown accession - let the real parser report that clearly

        if any(k in term for k in ("DiffMono", "MassMono")):

            return m.group(0)  # this one has a mass - nothing to patch

        counter[0] += 1

        sentinel = counter[0] * 1e-6  # a tiny, unique marker value, never a real observed mass

        placeholder_registry[sentinel] = getattr(term, "name", None) or f"MOD:{acc}"

        return f"[Obs:+{sentinel:.6f}]"
 
    patched = _MOD_TAG_RE.sub(_check, proforma_str)

    return patched, placeholder_registry
 
 
#

# pyteomics.proforma.ProForma.parse() reads a ProForma 2.0 string and gives

# us back a structured object:

#   pf.sequence   - list of (amino_acid, [tag, tag, ...]) pairs

#   pf.n_term / pf.c_term            - terminal modification tags

#   pf.unlocalized_modifications     - unlocalized modification groups

#   pf.intervals                     - bracketed/parenthesised sub-ranges,

#                                       including any tags attached to that

#                                       interval

#

# Each "tag" is an object describing one modification exactly as written

# (e.g. [Unimod:1], [MOD:00046], [Phospho], or a bare [+42.0106]). We resolve

# each tag into a small ResolvedMod(name, mono_mass, avg_mass, is_raw) record

# by asking the matching database for it. If nothing can resolve a tag, we

# fall back to reporting it as a raw mass shift, exactly as the previous

# version of this script did - this keeps the Level 2B logic unchanged.

# =============================================================================
 
class ResolvedMod:

    """A modification with its resolved display name and, where available,

    monoisotopic and average mass."""

    __slots__ = (
        "name", "mono_mass", "avg_mass", "formula", "numeric_shift",
        "is_raw", "isotope_approximate",
    )
 
    def __init__(
        self, name, mono_mass=None, avg_mass=None, formula=None,
        numeric_shift=None, is_raw=False, isotope_approximate=False,
    ):

        self.name = name

        self.mono_mass = mono_mass

        self.avg_mass = avg_mass

        self.formula = formula

        self.numeric_shift = numeric_shift

        self.is_raw = is_raw

        self.isotope_approximate = isotope_approximate
 
 
class _ParsedProteoform:

    """A small stand-in object exposing just the pieces the rest of this

    script needs, already resolved to ResolvedMod records."""
 
    def __init__(self, stripped_sequence, positions, n_term, c_term,

                 interval_mods, unlocalized_mods, has_unlocalized,
                 has_ambiguous_order):

        self.stripped_sequence = stripped_sequence

        self.positions = positions            # [(aa, [ResolvedMod, ...]), ...]

        self.n_term = n_term                  # [ResolvedMod, ...]

        self.c_term = c_term                  # [ResolvedMod, ...]

        self.interval_mods = interval_mods    # [(ResolvedMod, "posX-Y"), ...]

        self.unlocalized_mods = unlocalized_mods  # [ResolvedMod, ...]

        self.has_unlocalized = bool(has_unlocalized)

        self.has_ambiguous_order = has_ambiguous_order
 
 
def _resolve_tag(tag) -> ResolvedMod:

    """

    Resolve a single modification tag to a ResolvedMod.
 
    pyteomics.proforma resolves vocabulary tags into a display name and a

    chemical composition delta. We keep the name/formula, then calculate

    monoiso and average mass from that formula with mash_mass.
 
    MassModification (a bare [+34.5] or an explicit [Obs:+34.5] - both parse

    to this exact same class) has no formula. We add its numeric shift

    directly to mono/average mass and mark future isotope simulation as

    approximate. It is still Level 2B because the identity is ambiguous.

    """

    if isinstance(tag, pf_module.MassModification):

        return ResolvedMod(
            name=f"{tag.mass:+g} Da",
            mono_mass=tag.mass,
            avg_mass=tag.mass,
            numeric_shift=tag.mass,
            is_raw=True,
            isotope_approximate=True,
        )
 
    try:

        name = tag.name

        formula = _composition_to_formula(getattr(tag, "composition", None))

        if name and formula:

            mass = mash_mass.mass_from_formula(formula)

            return ResolvedMod(
                name=name,
                mono_mass=mass.mono,
                avg_mass=mass.average,
                formula=formula,
            )

    except Exception:

        pass  # not found in any vocabulary - fall through to raw/unresolved
 
    # Nothing matched: report the raw text the person wrote, flag Level 2B.

    return ResolvedMod(name=str(getattr(tag, "value", tag)), is_raw=True)


def _looks_like_mod_tag(value) -> bool:

    """Return True for pyteomics modification tag objects, False for counts

    and containers. This keeps unlocalized groups like ([MOD:...], 2) from

    accidentally resolving the count as if it were a modification."""

    if isinstance(value, pf_module.MassModification):

        return True

    return any(hasattr(value, attr) for attr in ("mass", "provider", "value"))


def _count_value(value) -> int:

    for attr in ("count", "multiplier"):

        raw = getattr(value, attr, None)

        if isinstance(raw, (int, float, str)):

            try:

                return int(raw)

            except (TypeError, ValueError):

                pass

    return 1


def _flatten_unlocalized_value(value):

    """Yield raw modification tag objects from pyteomics' unlocalized field.

    Different pyteomics releases have represented unlocalized modifications

    as direct tags, tuples, dict-like records, or small objects with tag/count

    attributes. This helper accepts those shapes and repeats tags when a count

    is present, while ignoring bookkeeping fields such as the count itself."""

    if value is None:

        return

    if _looks_like_mod_tag(value):

        yield value

        return

    if isinstance(value, dict):

        count = int(value.get("count") or value.get("multiplier") or 1)

        for key, item in value.items():

            if _looks_like_mod_tag(key) and isinstance(item, int):

                for _ in range(count * item):

                    yield key

        for key in ("tag", "mod", "modification", "tags", "mods", "modifications"):

            if key in value:

                for tag in _flatten_unlocalized_value(value[key]):

                    for _ in range(count):

                        yield tag

                return

        for item in value.values():

            yield from _flatten_unlocalized_value(item)

        return

    for attr in ("tag", "mod", "modification", "tags", "mods", "modifications"):

        if hasattr(value, attr):

            for tag in _flatten_unlocalized_value(getattr(value, attr)):

                for _ in range(_count_value(value)):

                    yield tag

            return

    if isinstance(value, (list, tuple)):

        local_count = _count_value(value)

        candidates = []

        for item in value:

            if isinstance(item, int):

                local_count *= item

            else:

                candidates.append(item)

        for item in candidates:

            for tag in _flatten_unlocalized_value(item):

                for _ in range(local_count):

                    yield tag


def _resolve_unlocalized_modifications(parsed, resolve_func) -> list[ResolvedMod]:

    raw_unlocalized = getattr(parsed, "unlocalized_modifications", None) or []

    return [resolve_func(tag) for tag in _flatten_unlocalized_value(raw_unlocalized)]


def _count_sequence_residues(text: str) -> int:

    count = 0

    in_tag = False

    for ch in text:

        if ch == "[":

            in_tag = True

        elif ch == "]":

            in_tag = False

        elif not in_tag and ch.isalpha():

            count += 1

    return count


def _parse_inline_tag_block(tag_block: str):

    """Parse one or more ProForma tags by attaching them to a dummy residue."""

    try:

        parsed = pf_module.ProForma.parse("A" + tag_block)

    except Exception:

        return []

    if not parsed.sequence:

        return []

    return parsed.sequence[0][1] or []


def _extract_tagged_interval_modifications(proforma_str: str, resolve_func):

    interval_mods = []

    for match in re.finditer(r"\(([A-Za-z]+)\)((?:\[[^\]]+\])+)", proforma_str):

        start = _count_sequence_residues(proforma_str[:match.start()]) + 1

        end = start + len(match.group(1)) - 1

        location = f"pos{start}-{end}"

        for tag in _parse_inline_tag_block(match.group(2)):

            interval_mods.append((resolve_func(tag), location))

    return interval_mods
 
 
def _parse_proforma(proforma_str: str) -> _ParsedProteoform:

    """Parse a ProForma string with pyteomics and resolve every modification

    tag to a ResolvedMod. Raises ValueError if the string itself is invalid

    ProForma (this is a syntax error, not a missing-database problem)."""

    normalized_str = _normalize_accession_tag_spacing(proforma_str)

    patched_str, placeholder_registry = _patch_massless_psimod(normalized_str)

    try:

        parsed = pf_module.ProForma.parse(patched_str)

    except Exception as e:

        raise ValueError(f"Could not parse ProForma string '{proforma_str}': {e}")
 
    def _resolve(tag) -> ResolvedMod:

        # A placeholder we inserted for a named-but-massless PSI-MOD term:

        # restore its real name. Since we DO know its identity, this is not

        # identity-ambiguous (not Level 2B) - just missing a mass number.

        if isinstance(tag, pf_module.MassModification) and tag.mass in placeholder_registry:

            return ResolvedMod(name=placeholder_registry[tag.mass], is_raw=False)

        return _resolve_tag(tag)
 
    positions = [(aa, [_resolve(t) for t in (tags or [])])

                 for aa, tags in parsed.sequence]

    parsed_interval_mods = []

    has_sequence_interval = False

    fallback_interval_mods = _extract_tagged_interval_modifications(normalized_str, _resolve)

    for interval in (getattr(parsed, "intervals", None) or []):

        # pyteomics interval coordinates are zero-based and end-exclusive.

        location = f"pos{interval.start + 1}-{interval.end}"

        tags = getattr(interval, "tags", None) or []

        if not tags:

            has_sequence_interval = True

        if fallback_interval_mods and tags:

            continue

        for tag in tags:

            parsed_interval_mods.append((_resolve(tag), location))

    interval_mods = []

    seen_interval_mods = set()

    for m, location in parsed_interval_mods + fallback_interval_mods:

        key = (m.name, location)

        if key in seen_interval_mods:

            continue

        interval_mods.append((m, location))

        seen_interval_mods.add(key)

    raw_unlocalized = getattr(parsed, "unlocalized_modifications", None)

    unlocalized_mods = _resolve_unlocalized_modifications(parsed, _resolve)
 
    return _ParsedProteoform(

        stripped_sequence="".join(aa for aa, _ in positions),

        positions=positions,

        n_term=[_resolve(t) for t in (parsed.n_term or [])],

        c_term=[_resolve(t) for t in (parsed.c_term or [])],

        interval_mods=interval_mods,

        unlocalized_mods=unlocalized_mods,

        has_unlocalized=bool(raw_unlocalized),

        has_ambiguous_order=has_sequence_interval,

    )
 
 
# =============================================================================

# SECTION 3 - UniProt REST API

# =============================================================================
 
_UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"

_uniprot_cache: dict = {}
 
def fetch_uniprot(uniprot_id: str) -> dict:

    if uniprot_id not in _uniprot_cache:

        resp = _SESSION.get(_UNIPROT_URL.format(uniprot_id), timeout=20)

        resp.raise_for_status()

        _uniprot_cache[uniprot_id] = resp.json()

    return _uniprot_cache[uniprot_id]


def _split_uniprot_accessions(uniprot_id: str) -> list[str]:

    """Split accession/protein-group fields while leaving normal IDs intact."""

    return [x.strip() for x in re.split(r"[;,|]", uniprot_id) if x.strip()]


def get_primary_gene_names(entry: dict) -> list[str]:

    names = []

    for gene in (entry.get("genes", []) or []):

        name = gene.get("geneName", {}).get("value")

        if name:

            names.append(name)

    return sorted(set(names))


def get_protein_name(entry: dict) -> str:

    protein_description = entry.get("proteinDescription", {}) or {}

    name_sources = [
        protein_description.get("recommendedName"),
        *((protein_description.get("submissionNames") or [])),
    ]

    for source in name_sources:

        if not source:

            continue

        full_name = source.get("fullName", {}) or {}

        value = full_name.get("value")

        if value:

            return value

    return entry.get("proteinName") or "N/A"


def has_gene_origin_ambiguity(uniprot_id: str, entry: dict) -> bool:

    # Multiple accessions in one input field means the protein/gene origin is
    # not uniquely assigned by the row.

    accessions = _split_uniprot_accessions(uniprot_id)

    if len(accessions) != 1:

        return True

    # Multiple primary gene entries within one UniProt record are not treated
    # as ambiguous here; only a missing/unclear gene field is considered 2D.

    return len(get_primary_gene_names(entry)) == 0
 
 
# =============================================================================

# SECTION 4 - Sequence alignment

# =============================================================================
 
_MAX_MISMATCH_FRACTION = 0.15   # up to 15 % of residues may differ

_MIN_MATCH_LENGTH      = 1      # always allow at least this many mismatches
 
 
def map_to_protein(query: str, full_seq: str) -> dict | None:

    q, f = len(query), len(full_seq)

    if q == 0 or f == 0 or q > f:

        return None
 
    max_mm = max(_MIN_MATCH_LENGTH, int(q * _MAX_MISMATCH_FRACTION))

    best   = None
 
    for start in range(f - q + 1):

        mm = [i for i in range(q) if full_seq[start + i] != query[i]]

        if len(mm) == 0:

            best = (0, start, mm)

            break

        if len(mm) <= max_mm and (best is None or len(mm) < best[0]):

            best = (len(mm), start, mm)
 
    if best is None:

        return None
 
    _, start, mm_idxs = best

    return {

        "start":     start + 1,

        "end":       start + q,

        "mismatches": [(start + i + 1, full_seq[start + i], query[i])

                       for i in mm_idxs],

    }
 
 
# =============================================================================

# SECTION 5 - Sequence alterations from direct sequence comparison

# =============================================================================
 
def build_sequence_alterations(mapping: dict | None, full_seq: str) -> list[str]:

    if mapping is None:

        return []
 
    alterations = []
 
    for prot_pos, uni_res, obs_res in mapping["mismatches"]:

        alterations.append(f"Substitution: {uni_res}{prot_pos}{obs_res}")

    if full_seq.startswith("M") and mapping["start"] == 2:

        alterations.append("Met removal")

    return alterations
 
 
# =============================================================================

# SECTION 6 - Proteoform Identification Level classification

# =============================================================================
 
def _has_sequence_ambiguity(proforma_str: str, pf) -> bool:

    if re.search(r"\([A-Za-z]{2,}\)(?!\[)", proforma_str):

        return True

    return bool(getattr(pf, "has_ambiguous_order", False))
 
 
def classify_level(mapping, ptm_localization_ambiguous, ptm_identity_ambiguous,

                   seq_ambiguous, gene_ambiguous=False):

    if mapping is None:

        return "Level 5"

    fail_A = bool(ptm_localization_ambiguous)

    fail_B = bool(ptm_identity_ambiguous)

    fail_C = bool(seq_ambiguous)

    fail_D = bool(gene_ambiguous)

    n = sum([fail_A, fail_B, fail_C, fail_D])

    if n == 0:

        return "Level 1"

    if n == 1:

        return "Level 2A" if fail_A else "Level 2B" if fail_B else "Level 2C" if fail_C else "Level 2D"

    if n == 2:

        return "Level 3"

    return "Level 4"
 
 
# =============================================================================

# SECTION 7 - Per-row processing

# =============================================================================
 
def _format_mod(m: "ResolvedMod", location: str) -> str:

    if location == "Unlocalized":

        return m.name

    return f"{m.name}@{location}"


def _mod_mass_inputs(mods: list[ResolvedMod]):

    formulas = []

    numeric_shifts = []

    for mod in mods:

        if mod.formula is not None:

            formulas.append(mod.formula)

        elif mod.numeric_shift is not None:

            numeric_shifts.append(mod.numeric_shift)

        else:

            return None

    return formulas, numeric_shifts


def _calculate_mash_proteoform_mass(sequence: str, mods: list[ResolvedMod]):

    inputs = _mod_mass_inputs(mods)

    if inputs is None:

        return None

    formulas, numeric_shifts = inputs

    try:

        return mash_mass.peptide_mass(sequence, formulas, numeric_shifts)

    except ValueError:

        return None


def _calculate_mercury_most_abundant_mass(sequence: str, mods: list[ResolvedMod]):

    inputs = _mod_mass_inputs(mods)

    if inputs is None:

        return None

    formulas, numeric_shifts = inputs

    try:

        formula = mash_mass.peptide_formula(sequence, formulas)

        return mercury_mass.calculate_most_abundant_mass(formula, numeric_shifts)

    except ValueError:

        return None
 
 
def process_row(proforma_str: str, uniprot_id: str) -> dict:

    # Guard against stray whitespace (including invisible non-breaking

    # spaces, e.g. "P13929\xa0" from a copy-paste) that would otherwise

    # break the UniProt API URL or the ProForma parser.

    proforma_str = proforma_str.strip()

    uniprot_id = uniprot_id.strip()

    accessions = _split_uniprot_accessions(uniprot_id)

    lookup_uniprot_id = accessions[0] if accessions else uniprot_id
 
    pf = _parse_proforma(proforma_str)

    clean_seq = pf.stripped_sequence
 
    entry    = fetch_uniprot(lookup_uniprot_id)

    full_seq = entry.get("sequence", {}).get("value", "")

    protein_name = get_protein_name(entry)

    col2 = full_seq
 
    mapping = map_to_protein(clean_seq, full_seq)

    col3    = f"{mapping['start']}-{mapping['end']}" if mapping else "Not Found"
 
    mod_parts    = []

    mod_entries  = []

    any_raw_mass = False
 
    for m in pf.n_term:

        location = "N-Term"

        mod_parts.append(_format_mod(m, location))

        mod_entries.append((m, location))

        any_raw_mass = any_raw_mass or m.is_raw
 
    for idx, (aa, mods) in enumerate(pf.positions, start=1):

        for m in mods:

            location = f"pos{idx}({aa})"

            mod_parts.append(_format_mod(m, location))

            mod_entries.append((m, location))

            any_raw_mass = any_raw_mass or m.is_raw
 
    for m, location in pf.interval_mods:

        mod_parts.append(_format_mod(m, location))

        mod_entries.append((m, location))

        any_raw_mass = any_raw_mass or m.is_raw

    for m in pf.c_term:

        location = "C-Term"

        mod_parts.append(_format_mod(m, location))

        mod_entries.append((m, location))

        any_raw_mass = any_raw_mass or m.is_raw

    for m in pf.unlocalized_mods:

        location = "Unlocalized"

        mod_parts.append(_format_mod(m, location))

        mod_entries.append((m, location))

        any_raw_mass = any_raw_mass or m.is_raw
 
    ptm_localization_ambiguous = pf.has_unlocalized or bool(pf.interval_mods)

    if pf.has_unlocalized and not pf.unlocalized_mods:

        mod_parts.append("Unlocalized modification")
 
    alteration_parts = mod_parts + build_sequence_alterations(mapping, full_seq)

    col4 = ", ".join(alteration_parts) if alteration_parts else "None"

    proteoform_mass = _calculate_mash_proteoform_mass(
        clean_seq,
        [m for m, _location in mod_entries],
    )

    col_proteoform_mono = f"{proteoform_mass.mono:.4f}" if proteoform_mass else "N/A"

    col_proteoform_avg = f"{proteoform_mass.average:.4f}" if proteoform_mass else "N/A"

    most_abundant_mass = _calculate_mercury_most_abundant_mass(
        clean_seq,
        [m for m, _location in mod_entries],
    )

    col_most_abundant = (
        f"{most_abundant_mass.most_abundant_mass:.4f}" if most_abundant_mass else "N/A"
    )

    seq_ambig    = _has_sequence_ambiguity(proforma_str, pf)

    # A repeated sequence within one UniProt entry is not gene-of-origin
    # ambiguity; Level 2D is reserved for ambiguous gene/protein origin.
    gene_ambig   = has_gene_origin_ambiguity(uniprot_id, entry)

    level        = classify_level(
        mapping,
        ptm_localization_ambiguous,
        any_raw_mass,
        seq_ambig,
        gene_ambig,
    )
 
    return {

        "Proteoform Sequence":          clean_seq,

        "Full Protein Sequence":        col2,

        "Protein Name":                 protein_name,

        "Sequence Range":                col3,

        "Modifications":                col4,

        "Experimental Monoisotopic Mass": "",

        "Experimental Average Mass": "",

        "Experimental Most Abundant Mass": "",

        "Theoretical Monoisotopic Mass": col_proteoform_mono,

        "Theoretical Average Mass": col_proteoform_avg,

        "Theoretical Most Abundant Mass": col_most_abundant,

        "ID Level":                     level,

        "Organ":                        "",

        "Source":                       "",

    }
 
 
# =============================================================================

# SECTION 8 - Inline sanity check

# =============================================================================
 
_MOCK_SEQ = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDL"
 
 
def _run_sanity_check() -> bool:

    print("[INFO] Running sanity check...")

    ok = True
 
    def _check(label, condition):

        nonlocal ok

        status = "PASS" if condition else "FAIL"

        print(f"  {label}: {status}")

        if not condition:

            ok = False
 
    pf_a = _parse_proforma("M[Unimod:1]VLSPADKTNVK")

    seq_a = pf_a.stripped_sequence

    map_a = map_to_protein(seq_a, _MOCK_SEQ)

    _check("Case A - sequence",   seq_a == "MVLSPADKTNVK")

    _check("Case A - range",      map_a is not None and map_a["start"] == 1 and map_a["end"] == 12)

    _check("Case A - no mismatch", map_a is not None and len(map_a["mismatches"]) == 0)
 
    pf_b = _parse_proforma("XYZINVALIDSEQUENCE")

    map_b = map_to_protein(pf_b.stripped_sequence, _MOCK_SEQ)

    _check("Case B - Not Found",  map_b is None)
 
    pf_c = _parse_proforma("MVLWPADKTNVK")

    seq_c = pf_c.stripped_sequence

    map_c = map_to_protein(seq_c, _MOCK_SEQ)

    alterations_c = build_sequence_alterations(map_c, _MOCK_SEQ)

    _check("Case C - range",       map_c is not None and map_c["start"] == 1 and map_c["end"] == 12)

    _check("Case C - S4W detected", any("S4W" in item for item in alterations_c))

    pf_d = _parse_proforma("VLSPADKTNVK")

    map_d = map_to_protein(pf_d.stripped_sequence, _MOCK_SEQ)

    alterations_d = build_sequence_alterations(map_d, _MOCK_SEQ)

    _check("Case D - methionine removal", "Met removal" in alterations_d)

    mash_seq = (
        "GREFGNLTRMRHVISYSLSPFEQRAYPHVFTKGIPNVLRRIRESFFRVVPQFVVFYLIYTWGTEEF"
        "ERSKRKNPAAYENDK"
    )

    mash_mass_check = mash_mass.peptide_mass(mash_seq)

    _check("Case E - MASH sequence mono mass",
           round(mash_mass_check.mono, 4) == 9769.0802)

    _check("Case E - MASH sequence average mass",
           round(mash_mass_check.average, 4) == 9775.1930)

    pf_f = _parse_proforma(
        "SEAEDASLLSFMQGYMKHATKTAKDALSSVQESQVAQQARGWVTDGFSSLKDYWSTVKDKFSEF"
        "WDLDPEVRPT[MOD:00528 ]SAVAA"
    )

    _check("Case F - MOD accession with trailing space",
           any(m.name == "Hex1HexNAc1NeuAc2 glycosylated residue"
               for _aa, mods in pf_f.positions for m in mods))

    return ok
 
 
def _read_input_table(path: Path) -> "pd.DataFrame":

    suffix = path.suffix.lower()

    if suffix == ".csv":

        return pd.read_csv(path)

    if suffix == ".tsv":

        return pd.read_csv(path, sep="\t")

    if suffix in {".xls", ".xlsx", ".xlsm"}:

        return pd.read_excel(path)

    raise ValueError(
        f"Unsupported input file type '{path.suffix}'. Use .csv, .tsv, .xls, .xlsx, or .xlsm."
    )


def _excel_safe_sequence_range(value):

    text = str(value)

    if re.fullmatch(r"\d+-\d+", text):

        return f'="{text}"'

    return value


def _write_excel_output(dataframe: "pd.DataFrame", path: Path) -> None:

    sheet_name = "Proteoforms"

    with pd.ExcelWriter(path, engine="openpyxl") as writer:

        dataframe.to_excel(writer, index=False, sheet_name=sheet_name)

        worksheet = writer.sheets[sheet_name]

        range_column = dataframe.columns.get_loc("Sequence Range") + 1

        for cells in worksheet.iter_cols(
            min_col=range_column,
            max_col=range_column,
            min_row=2,
        ):

            for cell in cells:

                cell.number_format = "@"


def _resolve_input_column(columns, preferred: str, aliases: list[str], label: str) -> str:

    candidates = [preferred] + aliases

    for candidate in candidates:

        if candidate in columns:

            return candidate

    lower_to_actual = {str(col).lower(): col for col in columns}

    for candidate in candidates:

        match = lower_to_actual.get(candidate.lower())

        if match is not None:

            return match

    raise ValueError(
        f"Column for {label} not found in {INPUT_FILE}.\n"
        f"        Tried: {candidates}\n"
        f"        Available columns: {list(columns)}"
    )


def _resolve_optional_input_column(columns, preferred: str, aliases: list[str]) -> str | None:

    try:

        return _resolve_input_column(columns, preferred, aliases, preferred)

    except ValueError:

        return None


def _copy_input_value(row, column: str | None) -> str:

    if column is None:

        return ""

    value = row[column]

    if pd.isna(value):

        return ""

    return str(value)


# =============================================================================

# SECTION 9 - Main entry point

# =============================================================================
 
def main():

    if not _run_sanity_check():

        print("\n[ERROR] Sanity Check FAILED: Pipeline logic, substitution tracking, "

              "or database mapping is broken! Aborting before file read.")

        sys.exit(1)

    print("[SUCCESS] Sanity Check passed.\n")
 
    print(f"[INFO] Reading input file: {INPUT_FILE}")

    try:

        df = _read_input_table(INPUT_FILE)

    except FileNotFoundError:

        sys.exit(f"[ERROR] Input file not found: {INPUT_FILE}\n"

                 f"        Edit the INPUT_FILE variable at the top of this script.")

    except ValueError as e:

        sys.exit(f"[ERROR] {e}")
 
    try:

        proforma_column = _resolve_input_column(
            df.columns,
            PROFORMA_COLUMN,
            ["ProForma", "proforma"],
            "ProForma string",
        )

        uniprot_column = _resolve_input_column(
            df.columns,
            UNIPROT_COLUMN,
            ["UniProt_ID", "uniprot id", "uniprot_id"],
            "UniProt accession",
        )

    except ValueError as e:

        sys.exit(f"[ERROR] {e}")

    source_column = _resolve_optional_input_column(
        df.columns,
        SOURCE_COLUMN,
        ["Source", "source"],
    )

    organ_column = _resolve_optional_input_column(
        df.columns,
        ORGAN_COLUMN,
        ["Organ", "organ", "Tissue", "tissue"],
    )

    experimental_mass_columns = {
        "Experimental Monoisotopic Mass": _resolve_optional_input_column(
            df.columns,
            EXPERIMENTAL_MONOISOTOPIC_MASS_COLUMN,
            ["experimental monoisotopic mass", "experimental_monoisotopic_mass"],
        ),
        "Experimental Average Mass": _resolve_optional_input_column(
            df.columns,
            EXPERIMENTAL_AVERAGE_MASS_COLUMN,
            ["experimental average mass", "experimental_average_mass"],
        ),
        "Experimental Most Abundant Mass": _resolve_optional_input_column(
            df.columns,
            EXPERIMENTAL_MOST_ABUNDANT_MASS_COLUMN,
            ["experimental most abundant mass", "experimental_most_abundant_mass"],
        ),
    }
 
    NEW_COLS = [

        "Protein Name", "Proteoform Sequence", "Full Protein Sequence",

        "Sequence Range",

        "Modifications",

        "Experimental Monoisotopic Mass",

        "Experimental Average Mass",

        "Experimental Most Abundant Mass",

        "Theoretical Monoisotopic Mass",

        "Theoretical Average Mass",

        "Theoretical Most Abundant Mass",
        
        "ID Level",

        "Organ",

        "Source",

    ]

    results = {c: [] for c in NEW_COLS}

    n_rows  = len(df)
 
    for idx, row in df.iterrows():

        proforma  = str(row[proforma_column])

        uniprot   = str(row[uniprot_column])

        row_label = f"Row {idx + 1}/{n_rows} (UniProt={uniprot})"

        try:

            res = process_row(proforma, uniprot)

        except Exception as e:

            print(f"[ERROR] {row_label}: {e}")

            res = {c: "ERROR" for c in NEW_COLS}

            res["ID Level"] = "ERROR"

        if source_column is not None and not pd.isna(row[source_column]):

            res["Source"] = str(row[source_column])

        else:

            res["Source"] = ""

        if organ_column is not None and not pd.isna(row[organ_column]):

            res["Organ"] = str(row[organ_column])

        else:

            res["Organ"] = ""

        for output_column, input_column in experimental_mass_columns.items():

            res[output_column] = _copy_input_value(row, input_column)

        for c in NEW_COLS:

            results[c].append(res[c])
 
        if (idx + 1) % 100 == 0 or (idx + 1) == n_rows:

            print(f"[INFO] Processed {idx + 1}/{n_rows} rows...")
 
    for c in NEW_COLS:

        df[c] = results[c]

    df[UNIPROT_OUTPUT_COLUMN] = df[uniprot_column]

    df[PROFORMA_OUTPUT_COLUMN] = df[proforma_column]

    requested_order = [

        UNIPROT_OUTPUT_COLUMN,

        "Protein Name",

        PROFORMA_OUTPUT_COLUMN,

        "Proteoform Sequence",

        "Full Protein Sequence",

        "Sequence Range",

        "Modifications",

        "Experimental Monoisotopic Mass",

        "Experimental Average Mass",

        "Experimental Most Abundant Mass",

        "Theoretical Monoisotopic Mass",

        "Theoretical Average Mass",

        "Theoretical Most Abundant Mass",

        "ID Level",

        "Organ",

    ]

    ordered_cols = [c for c in requested_order if c in df.columns]

    df = df[ordered_cols + ["Source"]]

    csv_df = df.copy()

    csv_df["Sequence Range"] = csv_df["Sequence Range"].map(_excel_safe_sequence_range)
 
    try:

        csv_df.to_csv(OUTPUT_FILE, index=False)

        saved_csv_file = OUTPUT_FILE

    except PermissionError:

        output_path = Path(OUTPUT_FILE)

        fallback_file = output_path.with_name(

            f"{output_path.stem}_updated{output_path.suffix}"

        )

        csv_df.to_csv(fallback_file, index=False)

        saved_csv_file = fallback_file

        print(f"\n[WARN] {OUTPUT_FILE} is open or locked. Saved to: {saved_csv_file}")

    excel_file = Path(saved_csv_file).with_suffix(".xlsx")

    try:

        _write_excel_output(df, excel_file)

        saved_excel_file = excel_file

    except PermissionError:

        fallback_excel_file = excel_file.with_name(f"{excel_file.stem}_updated.xlsx")

        _write_excel_output(df, fallback_excel_file)

        saved_excel_file = fallback_excel_file

        print(f"\n[WARN] {excel_file} is open or locked. Saved to: {saved_excel_file}")

    print(f"\n[SUCCESS] CSV results saved to:   {saved_csv_file}")

    print(f"[SUCCESS] Excel results saved to: {saved_excel_file}")

    print(f"          {n_rows} rows processed, {len(NEW_COLS)} columns added.")
 
 
if __name__ == "__main__":

    main()
 

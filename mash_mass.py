"""MASH sequence and formula mass calculations.

This module follows ``mash_reference/seqFunctions.cs`` for the peptide mass
path used by the MASH UI:

* residue mono/average masses are summed directly from ``AA_RES``
* peptide termini add ``H2O``
* residue formulas are still summed for downstream isotope work
* PTM formula deltas use the ``ELEMENTS`` masses from ``seqFunctions.cs``
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import re


@dataclass(frozen=True)
class MassResult:
    mono: float
    average: float
    isotope_approximate: bool = False


# seqFunctions.cs ELEMENTS.
MASH_ELEMENT_MASSES: dict[str, tuple[float, float]] = {
    "C": (12.0, 12.011),
    "H": (1.0078250, 1.00794),
    "N": (14.0030740, 14.00674),
    "O": (15.9949146, 15.9994),
    "F": (18.9984033, 18.9984),
    "P": (30.9737634, 30.97376),
    "S": (31.9720718, 32.066),
    "Cl": (34.9688527, 35.4527),
    "Br": (78.9183361, 79.904),
    "proton": (1.0078250 - 0.0005486, 1.00794 - 0.0005486),
    "el": (0.0005486, 0.0005486),
}


# seqFunctions.cs AA_RES: residual formula, residual mono mass, residual avg mass.
RESIDUES: dict[str, tuple[str, float, float]] = {
    "A": ("C3H5N1O1", 71.03711, 71.0788),
    "R": ("C6H12N4O1", 156.10111, 156.1876),
    "V": ("C5H9N1O1", 99.06841, 99.1326),
    "N": ("C4H6N2O2", 114.04293, 114.1039),
    "D": ("C4H5N1O3", 115.02694, 115.0886),
    "C": ("C3H5N1O1S1", 103.00919, 103.1448),
    "E": ("C5H7N1O3", 129.04259, 129.1155),
    "Q": ("C5H8N2O2", 128.05858, 128.1308),
    "G": ("C2H3N1O1", 57.02146, 57.0520),
    "H": ("C6H7N3O1", 137.05891, 137.1412),
    "I": ("C6H11N1O1", 113.08406, 113.1595),
    "L": ("C6H11N1O1", 113.08406, 113.1595),
    "K": ("C6H12N2O1", 128.09496, 128.1742),
    "M": ("C5H9N1O1S1", 131.04049, 131.1986),
    "F": ("C9H9N1O1", 147.06841, 147.1766),
    "P": ("C5H7N1O1", 97.05276, 97.1167),
    "S": ("C3H5N1O2", 87.03203, 87.0782),
    "T": ("C4H7N1O2", 101.04768, 101.1051),
    "W": ("C11H10N2O1", 186.07931, 186.2133),
    "Y": ("C9H9N1O2", 163.06333, 163.1760),
}

WATER_FORMULA = {"H": 2.0, "O": 1.0}
WATER_MASS = MassResult(
    mono=2 * 1.0078250 + 15.9949146,
    average=2 * 1.00794 + 15.9994,
)

_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]*|proton|el)\s*([0-9]*\.?[0-9]*)")


def add_formula(target: dict[str, float], source: Mapping[str, float], factor: float = 1.0) -> None:
    for element, count in source.items():
        target[element] = target.get(element, 0.0) + float(count) * factor


def parse_formula(formula: str) -> dict[str, float]:
    result: dict[str, float] = {}
    index = 0
    while index < len(formula):
        if formula[index] in " \t":
            index += 1
            continue
        match = _FORMULA_TOKEN_RE.match(formula, index)
        if not match:
            raise ValueError(f"Bad molecular formula at position {index + 1}: {formula!r}")
        element, raw_count = match.groups()
        if element not in MASH_ELEMENT_MASSES:
            raise ValueError(f"Unknown element in MASH seqFunctions formula: {element}")
        result[element] = result.get(element, 0.0) + (float(raw_count) if raw_count else 1.0)
        index = match.end()
    return result


RESIDUE_FORMULAS = {
    residue: parse_formula(formula)
    for residue, (formula, _mono, _average) in RESIDUES.items()
}


def mass_from_formula(formula: Mapping[str, float] | str) -> MassResult:
    composition = parse_formula(formula) if isinstance(formula, str) else formula
    mono = 0.0
    average = 0.0
    for element, count in composition.items():
        try:
            element_mono, element_average = MASH_ELEMENT_MASSES[element]
        except KeyError as exc:
            raise ValueError(f"Unknown element in MASH seqFunctions formula: {element}") from exc
        mono += element_mono * float(count)
        average += element_average * float(count)
    return MassResult(mono=mono, average=average)


def peptide_formula(sequence: str, modification_formulas: list[Mapping[str, float]] | None = None) -> dict[str, float]:
    formula: dict[str, float] = {}
    for residue in sequence:
        try:
            add_formula(formula, RESIDUE_FORMULAS[residue])
        except KeyError as exc:
            raise ValueError(f"Unsupported residue for MASH peptide formula: {residue}") from exc
    add_formula(formula, WATER_FORMULA)
    for mod_formula in modification_formulas or []:
        add_formula(formula, mod_formula)
    return {element: count for element, count in formula.items() if abs(count) > 1e-12}


def peptide_mass(
    sequence: str,
    modification_formulas: list[Mapping[str, float]] | None = None,
    numeric_shifts: list[float] | None = None,
) -> MassResult:
    mono = WATER_MASS.mono
    average = WATER_MASS.average
    for residue in sequence:
        try:
            _formula, residue_mono, residue_average = RESIDUES[residue]
        except KeyError as exc:
            raise ValueError(f"Unsupported residue for MASH peptide mass: {residue}") from exc
        mono += residue_mono
        average += residue_average

    for mod_formula in modification_formulas or []:
        mod_mass = mass_from_formula(mod_formula)
        mono += mod_mass.mono
        average += mod_mass.average

    isotope_approximate = False
    for shift in numeric_shifts or []:
        mono += float(shift)
        average += float(shift)
        isotope_approximate = True

    return MassResult(mono=mono, average=average, isotope_approximate=isotope_approximate)

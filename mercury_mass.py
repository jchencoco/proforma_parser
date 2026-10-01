"""Mercury-style most-abundant isotope mass calculation.

This ports the high-resolution path in
``mash_reference/MercuryIsotopeDistribution.h`` and reads isotope data from
``mash_reference/AtomicInformation.cpp``.  This module is only for isotope
simulation / most-abundant mass.  Sequence mono/average masses remain in
``mash_mass.py`` and are based on ``seqFunctions.cs``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
import math
import re

import numpy as np


PI = math.pi
ELECTRON_MASS = 0.00054858
DEFAULT_MERCURY_RESOLUTION = 100000.0


@dataclass(frozen=True)
class ElementalIsotopes:
    symbol: str
    average_mass: float
    mass_variance: float
    isotope_masses: tuple[float, ...]
    isotope_probs: tuple[float, ...]


@dataclass(frozen=True)
class MercuryResult:
    most_abundant_mass: float
    most_abundant_mz: float
    most_abundant_intensity: float
    mono_mass: float
    average_mass: float
    approximate: bool = False


_ATOMIC_CPP = Path(__file__).resolve().parent / "mash_reference" / "AtomicInformation.cpp"
_ATOMIC_RE = re.compile(
    r'strcpy\(isotopes\.marr_symbol, "([^"]+)"\).*?'
    r'isotopes\.mdbl_average_mass = ([0-9.]+);.*?'
    r'isotopes\.mdbl_mass_variance = ([0-9.]+) ;.*?'
    r'isotopes\.mint_num_isotopes = ([0-9]+);(?P<body>.*?)'
    r'mvect_elemental_isotopes\.push_back\(isotopes\)',
    re.DOTALL,
)
_MASS_RE = re.compile(r'isotopes\.marr_isotope_mass\[(\d+)\] = ([0-9.]+) ;')
_PROB_RE = re.compile(r'isotopes\.marr_isotope_prob\[(\d+)\] = ([0-9.]+) ;')


def _load_atomic_information() -> dict[str, ElementalIsotopes]:
    text = _ATOMIC_CPP.read_text()
    elements: dict[str, ElementalIsotopes] = {}
    for match in _ATOMIC_RE.finditer(text):
        symbol = match.group(1)
        if symbol in elements:
            continue
        average_mass = float(match.group(2))
        variance = float(match.group(3))
        num_isotopes = int(match.group(4))
        body = match.group("body")
        masses = {int(i): float(v) for i, v in _MASS_RE.findall(body)}
        probs = {int(i): float(v) for i, v in _PROB_RE.findall(body)}
        isotope_masses = tuple(masses[i] for i in range(num_isotopes))
        isotope_probs = tuple(probs[i] for i in range(num_isotopes))
        elements[symbol] = ElementalIsotopes(
            symbol=symbol,
            average_mass=average_mass,
            mass_variance=variance,
            isotope_masses=isotope_masses,
            isotope_probs=isotope_probs,
        )
    if not elements:
        raise RuntimeError(f"Could not load Mercury isotope data from {_ATOMIC_CPP}")
    return elements


ATOMIC_INFORMATION = _load_atomic_information()


def _mass_range_from_variance(variance: float, charge: int) -> int:
    if charge == 0:
        mass_range = int(math.sqrt(1 + variance) * 10)
    else:
        mass_range = int(math.sqrt(1 + variance) * 10.0 / float(charge))
    for value in (1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1):
        if value < mass_range:
            return value * 2
    return 1


def _parabolic_peak(x1: float, y1: float, x2: float, y2: float, x3: float, y3: float) -> float:
    denom = (y2 - y1) * (x3 - x2) - (y3 - y2) * (x2 - x1)
    if denom == 0:
        return x2
    return ((x1 + x2) - ((y2 - y1) * (x3 - x2) * (x1 - x3)) / denom) / 2.0


def _calculate_masses(formula: Mapping[str, float]) -> tuple[float, float]:
    mono = 0.0
    average = 0.0
    for symbol, count in formula.items():
        element = ATOMIC_INFORMATION[symbol]
        atomicity = int(count)
        if atomicity:
            mono += atomicity * element.isotope_masses[0]
            average += atomicity * element.average_mass
    return mono, average


def _calculate_variance(formula: Mapping[str, float]) -> float:
    variance = 0.0
    for symbol, count in formula.items():
        variance += ATOMIC_INFORMATION[symbol].mass_variance * int(count)
    return variance


def _calc_frequencies(formula: Mapping[str, float], charge: int, mass_range: int, num_points: int) -> np.ndarray:
    data = np.empty(num_points, dtype=np.complex128)
    half = num_points // 2
    items = [(ATOMIC_INFORMATION[symbol], int(count)) for symbol, count in formula.items() if int(count)]

    for idx in range(num_points):
        if idx < half:
            freq = idx / mass_range
        else:
            freq = (idx - num_points) / mass_range

        radius = 1.0
        theta = 0.0
        for element, atomicity in items:
            real = 0.0
            imag = 0.0
            for isotope_mass, isotope_prob in zip(element.isotope_masses, element.isotope_probs):
                wrap_freq = 0.0
                if len(element.isotope_masses) > 1:
                    wrap_freq = isotope_mass / charge - element.average_mass / charge
                    if wrap_freq < 0:
                        wrap_freq += mass_range
                x = 2 * PI * wrap_freq * freq
                real += isotope_prob * math.cos(x)
                imag += isotope_prob * math.sin(x)

            tempr = math.sqrt(real * real + imag * imag)
            radius *= tempr ** atomicity
            if real > 0:
                theta += atomicity * math.atan(imag / real)
            elif real < 0:
                theta += atomicity * (math.atan(imag / real) + PI)
            elif imag > 0:
                theta += atomicity * PI / 2
            else:
                theta -= atomicity * PI / 2

        data[idx] = complex(radius * math.cos(theta), radius * math.sin(theta))
    return data


def _apodize_gaussian(data: np.ndarray, subscript: float) -> None:
    num_points = len(data)
    if subscript == 0:
        return
    expdenom = (num_points / subscript) * (num_points / subscript)
    for idx in range(num_points):
        c_index = idx + 1
        if c_index <= num_points / 2:
            ap_val = math.exp(-(c_index - 1) * (c_index - 1) / expdenom)
        else:
            ap_val = math.exp(-(num_points - c_index - 1) * (num_points - c_index - 1) / expdenom)
        data[idx] *= ap_val


def calculate_most_abundant_mass(
    formula: Mapping[str, float],
    numeric_shifts: list[float] | None = None,
    charge: int = 1,
    resolution: float = DEFAULT_MERCURY_RESOLUTION,
    mercury_size: int = 8192,
    threshold: float = 0.0,
) -> MercuryResult:
    """Calculate Mercury most-abundant neutral mass for a formula."""
    if charge == 0:
        raise ValueError("Mercury charge 0 path is not valid for most-abundant mass output")
    clean_formula = {symbol: count for symbol, count in formula.items() if int(count)}
    unknown = sorted(set(clean_formula) - set(ATOMIC_INFORMATION))
    if unknown:
        raise ValueError(f"Unknown elements for Mercury isotope simulation: {', '.join(unknown)}")

    mono_mass, average_mass = _calculate_masses(clean_formula)
    variance = _calculate_variance(clean_formula)
    mass_range = _mass_range_from_variance(variance, charge)
    points_per_amu = mercury_size // mass_range
    num_points = mass_range * points_per_amu
    cc_mass = 1.00727638
    min_mz = average_mass / charge + (cc_mass - ELECTRON_MASS) - mass_range / 2.0

    ap_subscript = ((average_mass / (resolution * abs(charge))) * mercury_size * 2.0) / mass_range
    frequency_data = _calc_frequencies(clean_formula, charge, mass_range, num_points)
    _apodize_gaussian(frequency_data, ap_subscript)

    # Numerical Recipes Four1(..., isign=-1) uses the same sign convention as
    # numpy.fft.fft. Mercury normalizes later, so the missing scale factor is OK.
    transformed = np.fft.fft(frequency_data)
    raw_intensity = transformed.real
    max_intensity = float(np.max(raw_intensity))
    if max_intensity == 0:
        raise ValueError("Mercury isotope simulation produced zero intensity")
    intensity_by_index = raw_intensity * (100.0 / max_intensity)

    mz_values = []
    intensities = []
    half = num_points // 2
    for c_index in range(half + 1, num_points + 1):
        mz = (c_index - num_points - 1) / points_per_amu + average_mass / charge + (cc_mass - ELECTRON_MASS)
        mz_values.append(mz)
        intensities.append(float(intensity_by_index[c_index - 1]))
    for c_index in range(1, half + 1):
        mz = (c_index - 1) / points_per_amu + average_mass / charge + (cc_mass - ELECTRON_MASS)
        mz_values.append(mz)
        intensities.append(float(intensity_by_index[c_index - 1]))

    highest = -float("inf")
    peak_index = 0
    for idx, intensity in enumerate(intensities):
        if intensity > threshold and intensity > highest:
            highest = intensity
            peak_index = idx

    mz = mz_values[peak_index]
    y2 = intensities[peak_index]
    x2 = mz
    if peak_index > 0:
        x1 = mz_values[peak_index - 1]
        y1 = intensities[peak_index - 1]
    else:
        x1 = mz - 1.0 / (charge * points_per_amu)
        y1 = 0.0
    if peak_index < num_points - 1:
        x3 = mz_values[peak_index + 1]
        y3 = intensities[peak_index + 1]
    else:
        x3 = mz + 1.0 / (charge * points_per_amu)
        y3 = y2

    max_peak_mz = _parabolic_peak(x1, y1, x2, y2, x3, y3)
    most_intense_mw = max_peak_mz * charge - cc_mass * charge + ELECTRON_MASS * charge
    if abs(most_intense_mw - mono_mass) < 0.5 * 1.003 / charge:
        most_intense_mw = mono_mass
        max_peak_mz = (most_intense_mw + cc_mass * charge - ELECTRON_MASS * charge) / charge

    approximate = False
    for shift in numeric_shifts or []:
        most_intense_mw += float(shift)
        mono_mass += float(shift)
        average_mass += float(shift)
        approximate = True

    return MercuryResult(
        most_abundant_mass=most_intense_mw,
        most_abundant_mz=max_peak_mz,
        most_abundant_intensity=highest,
        mono_mass=mono_mass,
        average_mass=average_mass,
        approximate=approximate,
    )

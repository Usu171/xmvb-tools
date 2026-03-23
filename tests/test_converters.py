from pathlib import Path
import shutil

import numpy as np
import pytest

from xmvb_tools import gus2molden, molden2gus, no2molden, vb2molden
from xmvb_tools.molden import ReadMolden


DATA_DIR = Path(__file__).resolve().parent / 'data'
ATOL = 1e-8
RTOL = 1e-7


def normalize_coefficients(primitives: np.ndarray) -> np.ndarray:
    coefficients = primitives[:, 1]
    scale = np.max(np.abs(coefficients))
    if scale <= ATOL:
        return coefficients.copy()
    return coefficients / scale


def assert_molden_equal(actual_path: Path, expected_path: Path) -> None:
    actual = ReadMolden(str(actual_path))
    expected = ReadMolden(str(expected_path))

    assert actual['is_Cartesian'] == expected['is_Cartesian']
    assert actual['orbs'].shape == expected['orbs'].shape
    np.testing.assert_allclose(actual['orbs'], expected['orbs'], rtol=RTOL, atol=ATOL)

    assert actual['Atoms'].keys() == expected['Atoms'].keys()
    for atom_index in actual['Atoms']:
        actual_atom = actual['Atoms'][atom_index]
        expected_atom = expected['Atoms'][atom_index]
        assert actual_atom[0] == expected_atom[0]
        assert actual_atom[1] == pytest.approx(expected_atom[1], rel=RTOL, abs=ATOL)
        np.testing.assert_allclose(actual_atom[2], expected_atom[2], rtol=RTOL, atol=ATOL)

    assert actual['GTO'].keys() == expected['GTO'].keys()
    for atom_index in actual['GTO']:
        actual_shells = actual['GTO'][atom_index]
        expected_shells = expected['GTO'][atom_index]
        assert len(actual_shells) == len(expected_shells)

        for actual_shell, expected_shell in zip(actual_shells, expected_shells, strict=True):
            assert actual_shell[0] == expected_shell[0]
            assert len(actual_shell[1]) == len(expected_shell[1])
            actual_primitives = np.asarray(actual_shell[1], dtype=float)
            expected_primitives = np.asarray(expected_shell[1], dtype=float)
            np.testing.assert_allclose(
                actual_primitives[:, 0], expected_primitives[:, 0], rtol=RTOL, atol=ATOL
            )
            np.testing.assert_allclose(
                normalize_coefficients(actual_primitives),
                normalize_coefficients(expected_primitives),
                rtol=RTOL,
                atol=ATOL,
            )


@pytest.mark.parametrize('stem', ['C6H6_xmvb4', 'SF6_28'])
def test_vb2molden_matches_reference(tmp_path: Path, stem: str) -> None:
    shutil.copy(DATA_DIR / f'{stem}.orb', tmp_path / f'{stem}.orb')
    shutil.copy(DATA_DIR / f'{stem}.xmo', tmp_path / f'{stem}.xmo')

    prefix = tmp_path / stem
    vb2molden.WriteMolden(str(prefix), str(prefix), convert_basis=True)

    assert_molden_equal(tmp_path / f'{stem}.molden', DATA_DIR / f'{stem}_ori.molden')


@pytest.mark.parametrize('stem', ['C6H6_xmvb4', 'SF6_28'])
def test_gus2molden_matches_reference(tmp_path: Path, stem: str) -> None:
    shutil.copy(DATA_DIR / f'{stem}.xmo', tmp_path / f'{stem}.xmo')

    prefix = tmp_path / stem
    gus2molden.WriteMolden(str(prefix), convert_basis=True)

    assert_molden_equal(tmp_path / f'{stem}_vbgus.molden', DATA_DIR / f'{stem}_vbgus_ori.molden')


def test_molden2gus_matches_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    stem = 'SO4_lo'
    shutil.copy(DATA_DIR / f'{stem}.molden', tmp_path / f'{stem}.molden')
    inputs = (DATA_DIR / f'{stem}.txt').read_text(encoding='utf-8').splitlines()
    input_iter = iter(inputs)
    monkeypatch.setattr('builtins.input', lambda _prompt='': next(input_iter))

    molden2gus.main(str(tmp_path / stem))

    assert_molden_equal(tmp_path / f'{stem}_gus.molden', DATA_DIR / f'{stem}_gus_ori.molden')


def test_no2molden_reads_chunked_natural_orbitals(tmp_path: Path) -> None:
    sample = tmp_path / 'chunked.xmo'
    sample.write_text(
        '              ******   COMPUTED NATURAL ORBITALS   ******\n\n\n'
        '                          1          2          3          4          5\n'
        '                      2.000000   1.000000   0.500000   0.250000   0.125000\n'
        '    1  H  1  S        0.100000   0.200000   0.300000   0.400000   0.500000\n'
        '    2  H  1  PX       0.010000   0.020000   0.030000   0.040000   0.050000\n'
        '    3  H  1  PY       0.001000   0.002000   0.003000   0.004000   0.005000\n'
        '    4  H  1  PZ      -0.001000  -0.002000  -0.003000  -0.004000  -0.005000\n'
        '    5  H  2  S        0.110000   0.120000   0.130000   0.140000   0.150000\n'
        '    6  H  2  PX       0.210000   0.220000   0.230000   0.240000   0.250000\n'
        '\n\n'
        '                          6\n'
        '                      0.062500\n'
        '    1  H  1  S        0.600000\n'
        '    2  H  1  PX       0.060000\n'
        '    3  H  1  PY       0.006000\n'
        '    4  H  1  PZ      -0.006000\n'
        '    5  H  2  S        0.160000\n'
        '    6  H  2  PX       0.260000\n',
        encoding='utf-8',
    )

    occupations, eigenvectors = no2molden.ReadNaturalOrbitals(str(sample), 6)

    np.testing.assert_allclose(
        occupations,
        np.array([2.0, 1.0, 0.5, 0.25, 0.125, 0.0625]),
        rtol=RTOL,
        atol=ATOL,
    )
    np.testing.assert_allclose(
        eigenvectors[:, 0],
        np.array([0.1, 0.01, 0.001, -0.001, 0.11, 0.21]),
        rtol=RTOL,
        atol=ATOL,
    )
    np.testing.assert_allclose(
        eigenvectors[:, 5],
        np.array([0.6, 0.06, 0.006, -0.006, 0.16, 0.26]),
        rtol=RTOL,
        atol=ATOL,
    )


def test_no2molden_writes_molden_from_xmo_natural_orbitals(tmp_path: Path) -> None:
    stem = 'SF6_28'
    shutil.copy(DATA_DIR / f'{stem}.xmo', tmp_path / f'{stem}.xmo')

    prefix = tmp_path / stem
    n = no2molden.ReadNxmo(str(prefix.with_suffix('.xmo')))
    expected_occ, expected_orbs = no2molden.ReadNaturalOrbitals(str(prefix.with_suffix('.xmo')), n)
    _, expected_orbs = no2molden.ReadBasis(str(prefix.with_suffix('.xmo')), expected_orbs)

    no2molden.WriteMolden(str(prefix), convert_basis=True)

    actual = ReadMolden(str(tmp_path / f'{stem}_no.molden'))
    np.testing.assert_allclose(actual['occ'], expected_occ, rtol=RTOL, atol=ATOL)
    np.testing.assert_allclose(actual['orbs'], expected_orbs, rtol=RTOL, atol=ATOL)


# Copyright (C) 2024  Usu171

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse

import numpy as np

try:
    from .utils import ReadBasis, ReadGeoxmo, ReadNxmo, WriteSimpleMolden, float1
except ImportError:
    from utils import ReadBasis, ReadGeoxmo, ReadNxmo, WriteSimpleMolden, float1


NATURAL_ORBITALS_HEADER = 'COMPUTED NATURAL ORBITALS'


def is_orbital_index_line(parts):
    return bool(parts) and all(part.isdigit() for part in parts)


def ReadNaturalOrbitals(filename, n):
    occupations = np.zeros(n)
    eigenvectors = np.zeros((n, n))

    with open(filename, 'r') as file:
        lines = file.readlines()

    start = None
    for index, line in enumerate(lines):
        if NATURAL_ORBITALS_HEADER in line:
            start = index + 1
            break

    if start is None:
        raise ValueError(f'{NATURAL_ORBITALS_HEADER} section not found in {filename}')

    i = start
    while i < len(lines):
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        if not is_orbital_index_line(parts):
            i += 1
            continue

        columns = [int(part) - 1 for part in parts]
        if max(columns) >= n:
            raise ValueError(f'Natural orbital index exceeds basis size in {filename}')

        i += 1
        while i < len(lines) and not lines[i].split():
            i += 1
        if i >= len(lines):
            break

        occ_values = [float1(value) for value in lines[i].split()]
        if len(occ_values) != len(columns):
            raise ValueError(f'Occupation count mismatch in {filename}')
        for column, occupation in zip(columns, occ_values, strict=True):
            occupations[column] = occupation

        i += 1
        row_count = 0
        while i < len(lines):
            parts = lines[i].split()
            if not parts:
                i += 1
                break
            if is_orbital_index_line(parts):
                break
            if len(parts) < 4 + len(columns):
                raise ValueError(
                    f'Natural orbital row is too short in {filename}: {lines[i].rstrip()}'
                )

            row = int(parts[0]) - 1
            coeffs = [float1(value) for value in parts[4 : 4 + len(columns)]]
            for column, coefficient in zip(columns, coeffs, strict=True):
                eigenvectors[row, column] = coefficient
            row_count += 1
            i += 1

        if row_count != n:
            raise ValueError(
                f'Expected {n} natural orbital rows for block {columns[0] + 1}-{columns[-1] + 1}, got {row_count}'
            )

        if np.count_nonzero(occupations) == n:
            break

    return occupations, eigenvectors


def WriteMolden(file, convert_basis=True):
    n = ReadNxmo(f'{file}.xmo')
    occupations, eigenvectors = ReadNaturalOrbitals(f'{file}.xmo', n)
    gto_dict, eigenvectors = ReadBasis(f'{file}.xmo', eigenvectors)
    geo = ReadGeoxmo(f'{file}.xmo')
    WriteSimpleMolden(
        f'{file}_no.molden',
        geo,
        gto_dict,
        eigenvectors,
        occupations=occupations,
        convert_basis=convert_basis,
    )


def cli(argv=None):
    parser = argparse.ArgumentParser(
        description='Convert natural orbitals embedded in an XMVB .xmo file to a Molden file.'
    )
    parser.add_argument('xmo_prefix', help='Input filename without the .xmo suffix')
    parser.add_argument(
        '-b',
        action='store_true',
        help='Disable basis conversion.',
    )
    args = parser.parse_args(argv)
    WriteMolden(args.xmo_prefix, convert_basis=not args.b)


if __name__ == '__main__':
    cli()

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

import math
import re
import numpy as np

try:
    from .molden import WriteGTO
except ImportError:
    from molden import WriteGTO


basis_dict = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6, 'J': 7, 'K': 8}


def ReadNxmo(filename):
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(' NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS'):
                return int(line.split()[7])
    raise ValueError(f'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS not found in {filename}')


def ReadGeoxmo(filename):
    with open(filename, 'r') as file:
        geo = []
        geo_start = False
        count = 0
        for line in file:
            parts = line.split()
            if parts == ['CHARGE', 'X', 'Y', 'Z']:
                geo_start = True
                continue
            if geo_start:
                if parts:
                    count += 1
                    geo.append(
                        f'{parts[0]}{count:>6}'
                        + ''.join(f'{float(j):>15.8f}' for j in parts[1:5])
                        + '\n'
                    )
                else:
                    break
    return ''.join(geo)


def ReadOrb(filename, n):
    matrix = np.zeros((n, n))

    with open(filename, 'r') as file:
        next(file)
        col = None
        for line in file:
            if line.startswith('# ORBITAL'):
                parts = line.split()
                col = int(parts[2]) - 1
            elif col is not None:
                values = line.split()
                for i in range(0, len(values), 2):
                    value = float(values[i])
                    row = int(values[i + 1]) - 1
                    matrix[row, col] = value

    return matrix


def ReadOrbGus(filename, n):
    matrix = np.zeros((n, n))
    num_start = False
    item_start = False
    count1 = 0
    with open(filename, 'r') as file:
        counts = []
        for line in file:
            stripped = line.strip()
            if stripped == '---------------End of Guess--------------':
                break
            if stripped == '---------------Initial Guess---------------':
                num_start = True
                continue
            if not num_start:
                continue

            parts = line.split()
            if not parts:
                continue
            if '.' not in parts[0]:
                counts.append(parts)
                continue

            if not item_start:
                column_counts = [int(item) for row in counts for item in row]
                column_ends = np.cumsum(column_counts)
                item_start = True

            for i in range(0, len(parts), 2):
                value = float(parts[i])
                row = int(parts[i + 1]) - 1
                for col_index, end in enumerate(column_ends):
                    if count1 < end:
                        matrix[row, col_index] = value
                        break
                count1 += 1

    return matrix


def reorderD(matrix, i):
    temp = matrix[i + 1 : i + 6, :].copy()
    matrix[i + 1, :] = temp[2, :]  # XY 0  YY i+1
    matrix[i + 2, :] = temp[4, :]  # XZ 1  ZZ
    matrix[i + 3, :] = temp[0, :]  # YY 2  XY
    matrix[i + 4, :] = temp[1, :]  # YZ 3  XZ
    matrix[i + 5, :] = temp[3, :]  # ZZ 4  YZ
    return matrix


def reorderF(matrix, i):
    temp = matrix[i + 1 : i + 10, :].copy()
    matrix[i + 1, :] = temp[5, :]  # XXY 0  YYY i+1
    matrix[i + 2, :] = temp[8, :]  # XXZ 1  ZZZ
    # XYY 2  XYY
    matrix[i + 4, :] = temp[0, :]  # XYZ 3  XXY
    matrix[i + 5, :] = temp[1, :]  # XZZ 4  XXZ
    matrix[i + 6, :] = temp[4, :]  # YYY 5  XZZ
    matrix[i + 7, :] = temp[7, :]  # YYZ 6  YZZ
    matrix[i + 8, :] = temp[6, :]  # YZZ 7  YYZ
    matrix[i + 9, :] = temp[3, :]  # ZZZ 8  XYZ
    return matrix


def cartesian_ao_count(l_val):
    return (l_val + 1) * (l_val + 2) // 2


def ReadBasis(filename, matrix, atom_pattern=r'\b[A-Z][a-z]?\b'):
    atom_re = re.compile(atom_pattern)
    gto_dict = {}
    count1 = 0
    atom_index = 0
    basis_start = False

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if 'TOTAL NUMBER' in line:
                break

            if 'SHELL TYPE' in line:
                basis_start = True
                next(file)
                continue

            if not basis_start:
                continue

            if atom_re.match(line):
                atom_index += 1
                gto_dict[atom_index] = []
                continue

            if not line:
                continue

            parts = line.split()
            orbital_type = parts[1].upper()
            if orbital_type == 'L':
                count1 += 4
                s_primitives = []
                p_primitives = []
                while line and line.strip():
                    parts = line.split()
                    s_primitives.append([float(parts[3]), float(parts[4])])
                    p_primitives.append([float(parts[3]), float(parts[5])])
                    line = next(file).strip()

                gto_dict[atom_index].append([basis_dict['S'], s_primitives])
                gto_dict[atom_index].append([basis_dict['P'], p_primitives])
                continue

            l_val = basis_dict[orbital_type]
            primitives = []
            while line and line.strip():
                parts = line.split()
                primitives.append([float(parts[3]), float(parts[4])])
                line = next(file).strip()

            if orbital_type == 'D':
                reorderD(matrix, count1)
            elif orbital_type == 'F':
                reorderF(matrix, count1)
            count1 += cartesian_ao_count(l_val)
            gto_dict[atom_index].append([l_val, primitives])

    return gto_dict, matrix


def WriteSimpleMolden(
    filename,
    geo,
    gto_dict,
    matrix,
    occupations=None,
    title='qaq',
    convert_basis=False,
):
    num_orbitals = matrix.shape[1]
    if occupations is None:
        occupations = [2.0] * num_orbitals

    gto_to_write = convert_xmvb_basis(gto_dict) if convert_basis else gto_dict

    with open(filename, 'w') as mol:
        mol.write('[Molden Format]\n[Title]\n')
        mol.write(f'{title}\n')
        mol.write('[Atoms] AU\n')
        mol.write(geo)
        WriteGTO(mol, gto_to_write)
        mol.write('\n' * 3)
        mol.write('[6D10F]\n[MO]\n')
        for i in range(num_orbitals):
            mol.write(f'Sym=     1a\nEne= 1.0\nSpin= Alpha\nOccup= {occupations[i]:>15.10f}\n')
            for j in range(matrix.shape[0]):
                mol.write(f'{j + 1:>4}  {matrix[j, i]:>15.10f}\n')


def double_factorial(n):
    if n <= 0:
        return 1
    res = 1
    for i in range(n, 0, -2):
        res *= i
    return res


def convert_xmvb_basis(gto_dict):
    bse_data = {}

    for atom_idx, shells in gto_dict.items():
        new_shells = []
        for l_val, primitives in shells:
            dfact = double_factorial(2 * l_val - 1)

            k_values = []
            for alpha, coeff_xmvb in primitives:
                n_cart = (
                    ((2 * alpha / math.pi) ** 0.75)
                    * ((4 * alpha) ** (l_val / 2.0))
                    / math.sqrt(dfact)
                )

                k_i = coeff_xmvb / n_cart
                k_values.append(k_i)

            s_total = 0.0
            for i in range(len(primitives)):
                alpha_i = primitives[i][0]
                for j in range(len(primitives)):
                    alpha_j = primitives[j][0]

                    s_ij = ((2 * math.sqrt(alpha_i * alpha_j)) / (alpha_i + alpha_j)) ** (
                        l_val + 1.5
                    )
                    s_total += k_values[i] * k_values[j] * s_ij

            scale = math.sqrt(s_total)

            bse_primitives = []
            for i in range(len(primitives)):
                alpha = primitives[i][0]
                c_bse = k_values[i] / scale
                bse_primitives.append([alpha, c_bse])

            new_shells.append([l_val, bse_primitives])

        bse_data[atom_idx] = new_shells

    return bse_data


def float1(num):
    return float(num.replace('D', 'E'))

def ReadEig(filename, n):
    eigenvalues = []
    eigenvectors = np.zeros((n, n))
    row1, row2 = np.divmod(n, 5)
    col = -1
    with open(filename, 'r') as file:
        for line in file:
            if line:
                parts = line.split()
                eigenvalues.append(float1(parts[0]))
                col += 1
                for row in range(row1 + 1):
                    parts = next(file).split()
                    i = row * 5
                    if row2 == 0 and row == row1 - 1:
                        break
                    if row == row1:
                        eigenvectors[i : i + row2, col] = [float1(v) for v in parts]
                        break
                    eigenvectors[i : i + 5, col] = [float1(v) for v in parts]
            else:
                break

    return np.array(eigenvalues), eigenvectors


def SortEig(eigenvalues, eigenvectors):
    indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[indices]
    sorted_eigenvectors = eigenvectors[:, indices]
    return sorted_eigenvalues, sorted_eigenvectors


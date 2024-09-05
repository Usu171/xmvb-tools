import re
import sys

import numpy as np


def ReadNxmo(filename):
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(' NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS'):
                n = int(line.split()[7])
                break
    return n


def ReadGEOxmo(filename):
    with open(filename, 'r') as file:
        geo = []
        GEO_start = False
        count = 0
        for line in file:
            parts = line.split()
            if parts == ['CHARGE', 'X', 'Y', 'Z']:
                GEO_start = True
                continue
            elif GEO_start:
                if parts:
                    count += 1
                    geo.append(f'{parts[0]}{count:>6}' + ''.join(
                        f'{float(j):>15.8f}' for j in parts[1:5]) + '\n')
                else:
                    break
    return ''.join(geo)


def ReadOrb(filename, n):
    matrix = np.zeros((n, n))

    with open(filename, 'r') as file:

        next(file)
        col = None
        for line in file:

            if line.startswith("# ORBITAL"):
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
        c = []
        for line in file:
            if line.startswith(' --------------End of Guess--------------'):
                break
            if line.startswith(' --------------Initial Guess--------------'):
                num_start = True
                continue
            parts = line.split()
            if num_start:
                if '.' not in parts[0]:
                    c.append(line.split())
                    continue
                else:
                    if not item_start:
                        c1 = [int(item) for row in c for item in row]
                        column_ends = np.cumsum(c1)
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


def ChangeD(matrix, i):
    temp = matrix[i + 1:i + 6, :].copy()
    matrix[i + 1, :] = temp[2, :]  # XY 0  YY i+1
    matrix[i + 2, :] = temp[4, :]  # XZ 1  ZZ
    matrix[i + 3, :] = temp[0, :]  # YY 2  XY
    matrix[i + 4, :] = temp[1, :]  # YZ 3  XZ
    matrix[i + 5, :] = temp[3, :]  # ZZ 4  YZ
    return matrix


def ChangeF(matrix, i):
    temp = matrix[i + 1:i + 10, :].copy()
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


def ReadBasis(filename, matrix):
    atom_pattern = re.compile(r'\b[A-Z][A-Z]?\b')
    output = []
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

            if basis_start:
                if atom_pattern.match(line):
                    atom_index += 1
                    output.append(f'\n{atom_index}  0\n')
                    continue

                if line:
                    parts = line.split()
                    orbital_type = parts[1]
                    coefficient = []
                    count = 0
                    if orbital_type == 'L':
                        count1 += 4
                        coefficient_s = []
                        coefficient_p = []
                        while line and line.strip():
                            parts = line.split()
                            coefficient_s.append(
                                f'  {parts[3]:>14}    {parts[4]:>16}')
                            coefficient_p.append(
                                f'  {parts[3]:>14}    {parts[5]:>16}')
                            count += 1
                            line = next(file).strip()

                        output.append(f'S {count} 1.0\n' + ''.join(
                            f'{info}\n' for info in coefficient_s))
                        output.append(f'P {count} 1.0\n' + ''.join(
                            f'{info}\n' for info in coefficient_p))
                    else:
                        if orbital_type == 'S':
                            count1 += 1
                        elif orbital_type == 'P':
                            count1 += 3
                        elif orbital_type == 'D':
                            ChangeD(matrix, count1)
                            count1 += 6
                        elif orbital_type == 'F':
                            ChangeF(matrix, count1)
                            count1 += 10
                        while line and line.strip():
                            parts = line.split()
                            coefficient.append(
                                f'  {parts[3]:>14}    {parts[4]:>16}')
                            count += 1
                            line = next(file).strip()

                        output.append(f'{orbital_type} {count} 1.0\n' + ''.join(
                            f'{info}\n' for info in coefficient))

    return ''.join(output), matrix


def WriteMolden(file, file2):
    n = ReadNxmo(f'{file2}.xmo')

    matrix = ReadOrbGus(f'{file}.xdat', n)
    # matrix = ChangeBasis(matrix, 'INFO')

    basis, matrix = ReadBasis(f'{file2}.xmo', matrix)

    geo = ReadGEOxmo(f'{file2}.xmo')

    with open(f'{file}.molden', 'w') as mol:
        mol.write('''[Molden Format]
[Title]
qaq
[Atoms] AU
''')
        mol.write(geo)
        mol.write('[GTO]')
        mol.write(basis)
        mol.write('\n' * 3)
        mol.write('''[6D10F]\n[MO]\n''')
        for i in range(0, n):
            mol.write('''Sym=     1a
Ene= 1.0
Spin= Alpha
Occup= 2.000000
''')
            for j in range(0, n):
                mol.write(f'{j+1:>4}  {matrix[j,i]:>15.10f}\n')


if __name__ == "__main__":

    WriteMolden(sys.argv[1], sys.argv[2])

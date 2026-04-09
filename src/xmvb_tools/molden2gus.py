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
    from .molden import ReadMolden, WriteMolden
    from .reorder import reorder_shell
except ImportError:
    from molden import ReadMolden, WriteMolden
    from reorder import reorder_shell


def ExtractRows(matrix, atoms, column_index, position):
    result = np.zeros((matrix.shape[0], 1))

    atom_start = np.hstack((0, position[:-1]))

    for atom in atoms:
        atom -= 1
        start = atom_start[atom]
        end = position[atom]
        result[start:end, 0] = matrix[start:end, column_index]

    return result


def Write(filename, matrix):
    counts = []
    details = []

    for col in range(matrix.shape[1]):
        col_data = matrix[:, col]
        non_zeros = np.nonzero(col_data)[0]
        non_zero_values = col_data[non_zeros]

        counts.append(len(non_zero_values))

        data = [f'# ORBITAL{col + 1:>11}']
        for i in range(0, len(non_zero_values), 4):
            line = ''.join(
                f'{non_zero_values[i + j]:>13.10f} {non_zeros[i + j] + 1:>5}  '
                for j in range(4)
                if i + j < len(non_zero_values)
            )
            data.append(line)
        details.append(data)

    with open(filename, 'w') as file:
        file.write(''.join(f'   {i}' for i in counts) + '\n')

        for data in details:
            for line in data:
                file.write(line + '\n')


def Rotate(matrix, angle, col1, clo2):
    angle = np.deg2rad(angle)
    Rmatrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    cols = matrix[:, [col1, clo2]]
    matrix[:, [col1, clo2]] = np.dot(cols, Rmatrix)
    return matrix


def ParseInp(str):
    result = []
    parts = str.split(',')
    for i in parts:
        if '-' in i:
            start, end = map(int, i.split('-'))
            result.extend(range(start, end + 1))
        else:
            result.append(int(i))
    return result


def main(filename):
    mol = ReadMolden(f'{filename}.molden')
    natm = len(mol['Atoms'])
    c = mol['orbs']
    row_len, col_len = c.shape

    if mol['is_Cartesian']:
        basis_type = [
            (basis[0] + 1) * (basis[0] + 2) // 2 for atom in mol['GTO'].values() for basis in atom
        ]

        position = np.cumsum(
            [
                sum((basis[0] + 1) * (basis[0] + 2) // 2 for basis in atom)
                for atom in mol['GTO'].values()
            ]
        )
    else:
        basis_type = [basis[0] * 2 + 1 for atom in mol['GTO'].values() for basis in atom]

        position = np.cumsum(
            [sum(basis[0] * 2 + 1 for basis in atom) for atom in mol['GTO'].values()]
        )
    result = np.empty((row_len, 0))
    while True:
        inp = input("""Please input atom numbers and orbital numbers
for example: 2,4 6-8,9
Input \'a <cols>\' to add all atoms
Input \'r<angle> <col1>,<col2>\' to rotate two orbitals
Input \'m<num> <col>\' to multip orbital by <num>
Input \'q\' to write and exit\n""")
        try:
            if inp.lower() == 'q':
                mol['orbs'] = result
                WriteMolden(mol, f'{filename}_gus.molden')
                print(f'{filename}_gus.molden has been written')

                if mol['is_Cartesian'] and all(i <= 10 for i in basis_type):
                    count = 0
                    for i in basis_type:
                        if i == 6:
                            reorder_shell(result, count, 'D', inverse=True)
                        elif i == 10:
                            reorder_shell(result, count, 'F', inverse=True)
                        count += i

                    Write(f'{filename}.gus', result)
                    print(f'{filename}.gus has been written')

                break
            elif inp.lower().startswith('r'):
                angle, cols = inp.split()
                angle = float(angle[1:])
                cols = [int(x) - 1 for x in cols.split(',')]
                if len(cols) != 2 or max(cols) >= col_len or min(cols) < 0:
                    raise ValueError('Invalid number')
                c = Rotate(c, angle, cols[0], cols[1])
                print(
                    f'Orbital {cols[0] + 1} and Orbital {cols[1] + 1} have been rotated by {angle} degrees'
                )
            elif inp.lower().startswith('m'):
                num = float(inp.split()[0][1:])
                col = int(inp.split()[1]) - 1
                if col >= col_len or col < 0:
                    raise ValueError('Invalid number')
                c[:, col] *= num
                print(f'Orbital {col + 1} has been multiplied by {num}')
            elif inp.lower().startswith('a'):
                cols = [x - 1 for x in ParseInp(inp.split()[1])]
                if max(cols) >= col_len or min(cols) < 0:
                    raise ValueError('Invalid number')
                for col in cols:
                    result = np.hstack((result, c[:, col].reshape(-1, 1)))
            else:
                atoms, cols = inp.split()
                atoms = sorted(list(ParseInp(atoms)))
                cols = [x - 1 for x in ParseInp(cols)]
                if atoms[-1] > natm or max(cols) >= col_len or atoms[0] < 1 or min(cols) < 0:
                    raise ValueError('Invalid number')
                for col in cols:
                    result = np.hstack((result, ExtractRows(c, atoms, col, position)))
        except ValueError as e:
            print(e)


def cli(argv=None):
    parser = argparse.ArgumentParser(
        description='Read a Molden file and generate XMVB guess orbitals.'
    )
    parser.add_argument('filename', help='Input filename without the .molden suffix')
    args = parser.parse_args(argv)
    main(args.filename)


if __name__ == '__main__':
    cli()

import re
import sys
import numpy as np


def float1(num):
    return float(num.replace('D', 'E'))


def SortEig1(eigenvalues, eigenvectors):

    indices = np.argsort(eigenvalues)[::-1]
    
    sorted_eigenvalues = eigenvalues[indices]
    sorted_eigenvectors = eigenvectors[:, indices]
    
    return sorted_eigenvalues, sorted_eigenvectors


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


def ChangeD(matrix, i):
    temp = matrix[i + 1:i + 6, :].copy()
    matrix[i + 1, :] = temp[2, :]  # XY 0  YY i+1
    matrix[i + 2, :] = temp[4, :]  # XZ 1  ZZ
    matrix[i + 3, :] = temp[0, :]  # YY 2  XY
    matrix[i + 4, :] = temp[1, :]  # YZ 3  XZ
    matrix[i + 5, :] = temp[3, :]  # ZZ 4  YZ
    return matrix


def ChangeF(matrix, i):
    temp = matrix[i + 1:i + 10,].copy()
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
                    elif row == row1:
                        eigenvectors[i: i + row2, col] = [float1(v) for v in parts]
                        break
                    eigenvectors[i: i + 5, col] = [float1(v) for v in parts]

            else:
                break

    return np.array(eigenvalues), eigenvectors

# 示例用法


def WriteMolden(file):
    n = ReadNxmo(f'{file}.xmo')

    eigenvalues, eigenvectors = ReadEig(f'xmvb.no', n)

    eigenvalues, eigenvectors = SortEig1(eigenvalues, eigenvectors)

    basis, eigenvectors = ReadBasis(f'{file}.xmo', eigenvectors)
    

    print(np.sum(eigenvalues))
    geo = ReadGEOxmo(f'{file}.xmo')


    with open(f'{file}_no.molden', 'w') as mol:
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
            mol.write(f'''Sym=     1a
Ene= 1.0
Spin= Alpha
Occup= {eigenvalues[i]:>15.10f}
''')
            for j in range(0, n):
                mol.write(f'{j+1:>4}  {eigenvectors[j,i]:>15.10f}\n')


if __name__ == '__main__':

    WriteMolden(sys.argv[1])

import sys

import numpy as np
from pyscf.tools.molden import *


def _parse_mo(lines, envs):
    mol = envs['mol']
    if not mol._built:
        try:
            mol.build(0, 0)
        except RuntimeError:
            mol.build(0, 0, spin=1)

    irrep_labels = []
    mo_energy = []
    spins = []
    mo_occ = []
    mo_coeff_prim = []  # primary data, will be reworked for missing values
    coeff_idx = []
    mo_id = 0
    for line in lines[1:]:
        line = line.upper()
        if 'SYM' in line:
            irrep_labels.append(line.split('=')[1].strip())
        elif 'ENE' in line:
            mo_energy.append(float(_d2e(line).split('=')[1].strip()))
            mo_id = len(mo_energy) - 1
        elif 'SPIN' in line:
            spins.append(line.split('=')[1].strip())
        elif 'OCC' in line:
            mo_occ.append(float(_d2e(line.split('=')[1].strip())))
        else:
            ao_id, c = line.split()[:2]
            coeff_idx.append([int(ao_id) - 1, mo_id])
            mo_coeff_prim.append(float(c))

    coeff_idx = np.array(coeff_idx)
    number_of_aos, number_of_mos = coeff_idx.max(axis=0) + 1
    mo_coeff = np.zeros([number_of_aos, number_of_mos])
    mo_coeff[coeff_idx[:, 0], coeff_idx[:, 1]] = mo_coeff_prim

    mo_energy = np.array(mo_energy)
    mo_occ = np.array(mo_occ)

    return mol, mo_energy, mo_coeff, mo_occ, irrep_labels, spins


def orbital_coeff(mol, fout, mo_coeff, spin='Alpha', symm=None, ene=None,
                  occ=None, ignore_h=IGNORE_H):
    from pyscf.symm import label_orb_symm

    if ignore_h:
        mol, mo_coeff = remove_high_l(mol, mo_coeff)

    nmo = mo_coeff.shape[1]
    if symm is None:
        symm = ['A']*nmo
        if mol.symmetry:
            try:
                symm = label_orb_symm(mol, mol.irrep_name, mol.symm_orb,
                                      mo_coeff, tol=1e-5)
            except ValueError as e:
                logger.warn(mol, str(e))
    if ene is None or len(ene) != nmo:
        ene = np.arange(nmo)
    assert (spin == 'Alpha' or spin == 'Beta')
    if occ is None:
        occ = np.zeros(nmo)
        neleca, nelecb = mol.nelec
        if spin == 'Alpha':
            occ[:neleca] = 1
        else:
            occ[:nelecb] = 1

    if spin == 'Alpha':
        # Avoid duplicated [MO] session when dumping beta orbitals
        fout.write('[MO]\n')

    for imo in range(nmo):
        fout.write(' Sym= %s\n' % symm[imo])
        fout.write(' Ene= %15.10g\n' % ene[imo])
        fout.write(' Spin= %s\n' % spin)
        fout.write(' Occup= %10.5f\n' % occ[imo])
        for i in mo_coeff.shape[1]:
            fout.write(' %3d    %18.14g\n' % (i+1, mo_coeff[i]))


def ExtractRows(matrix, atoms, column_index, position):

    result = np.zeros((matrix.shape[0], 1))

    atom_start = np.hstack((0, position[:-1]))

    for atom in atoms:
        atom -= 1
        start = atom_start[atom]
        end = position[atom]
        result[start:end, 0] = matrix[start:end, column_index]

    return result


def ChangeDi(matrix, i):
    temp = matrix[i + 1:i + 6, :].copy()
    matrix[i + 1, :] = temp[2, :]  # YY 0  XY i+1
    matrix[i + 2, :] = temp[3, :]  # ZZ 1  XZ
    matrix[i + 3, :] = temp[0, :]  # XY 2  YY
    matrix[i + 4, :] = temp[4, :]  # XZ 3  YZ
    matrix[i + 5, :] = temp[1, :]  # YZ 4  ZZ
    return matrix


def ChangeFi(matrix, i):
    temp = matrix[i + 1:i + 10,].copy()
    matrix[i + 1, :] = temp[3, :]  # YYY 0  XXY i+1
    matrix[i + 2, :] = temp[4, :]  # ZZZ 1  XXZ
    # XYY 2  XYY
    matrix[i + 4, :] = temp[8, :]  # XXY 3  XYZ
    matrix[i + 5, :] = temp[5, :]  # XXZ 4  XZZ
    matrix[i + 6, :] = temp[0, :]  # XZZ 5  YYY
    matrix[i + 7, :] = temp[7, :]  # YZZ 6  YYZ
    matrix[i + 8, :] = temp[6, :]  # YYZ 7  YZZ
    matrix[i + 9, :] = temp[1, :]  # XYZ 8  ZZZ
    return matrix


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
            line = ''
            for j in range(4):
                if i + j < len(non_zero_values):
                    line += f'{non_zero_values[i + j]:>13.10f} {non_zeros[i + j] + 1:>5}  '
            data.append(line)
        details.append(data)

    with open(filename, 'w') as file:
        file.write(''.join(f'   {i}' for i in counts) + '\n')

        for data in details:
            for line in data:
                file.write(line + '\n')


def Rotate(matrix, angle, col1, clo2):
    angle = np.deg2rad(angle)
    Rmatrix = np.array([[np.cos(angle), -np.sin(angle)],
                        [np.sin(angle), np.cos(angle)]])
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
    a = read(f'{filename}.molden')
    natm = a[0].natm
    c = a[2]
    row_len, col_len = c.shape
    count = 0
    basis_type = [(i[0] + 1) * (i[0] + 2) // 2
                for key in a[0]._basis.values()
                for i in key]

    position = np.cumsum([
        sum((i[0] + 1) * (i[0] + 2) // 2
            for i in key)
        for key in a[0]._basis.values()
    ])
    result = np.empty((row_len, 0))
    while True:
        inp = input('''Please input atom numbers and orbital numbers
for example: 2,4 6-8,9
Input \'a <cols>\' to add all atoms
Input \'r<angle> <col1>,<col2>\' to rotate two orbitals
Input \'m<num> <col>\' to multip orbital by <num>
Input \'q\' to write and exit\n''')
        try:
            if inp.lower() == 'q':
                from_mo(a[0], f'{filename}_gus.molden', result)

                for i in basis_type:
                    if i == 6:
                        ChangeDi(result, count)
                    elif i == 10:
                        ChangeFi(result, count)
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
                    f'Orbital {cols[0]+1} and Orbital {cols[1]+1} have been rotated by {angle} degrees'
                )
            elif inp.lower().startswith('m'):
                num = float(inp.split()[0][1:])
                col = int(inp.split()[1]) - 1
                if col >= col_len or col < 0:
                    raise ValueError('Invalid number')
                c[:, col] *= num
                print(f'Orbital {col+1} has been multiplied by {num}')
            elif inp.lower().startswith('a'):
                cols = [x - 1 for x in ParseInp(inp.split()[1])]
                if max(cols) >= col_len or min(cols) < 0:
                    raise ValueError('Invalid number')
                for col in cols:
                    result = np.hstack(
                        (result, c[:, col].reshape(-1, 1)))
            else:
                atoms, cols = inp.split()
                atoms = sorted([x for x in ParseInp(atoms)])
                cols = [x - 1 for x in ParseInp(cols)]
                if atoms[-1] > natm or max(
                        cols) >= col_len or atoms[0] < 1 or min(cols) < 0:
                    raise ValueError('Invalid number')
                for col in cols:
                    result = np.hstack(
                        (result, ExtractRows(c, atoms, col, position)))
        except ValueError as e:
            print(e)


if __name__ == '__main__':

    main(sys.argv[1])

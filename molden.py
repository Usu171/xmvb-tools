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

import re
import numpy as np


basis_dict = {
    'S': 0,
    'P': 1,
    'D': 2,
    'F': 3,
    'G': 4,
    'H': 5,
    'I': 6,
    'J': 7,
    'K': 8
}

reversed_basis_dict = {v: k for k, v in basis_dict.items()}


def read_parts(file):
    current_part = []
    current_part_name = None
    in_part = False

    for line in file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        if match1 := re.match(r'^\[(.*)\]', line):
            if in_part:
                yield current_part_name, current_part
                current_part = []
            current_part_name = match1[1].upper()
            in_part = True
        elif in_part:
            current_part.append(line)

    if current_part_name is not None:
        yield current_part_name, current_part


def process_title(part, mol):
    if 'orca' in part[0]:
        mol['Title'] = 'orca'


def read_GTO_parts(part):
    current_part = []
    number = None
    in_part = False
    for line in part:
        if match1 := re.match(r'^(\d+)\s*0$', line):
            if in_part:
                yield number, current_part
                current_part = []
            number = int(match1[1])
            in_part = True
        elif in_part:
            current_part.append(line)
    if number is not None:
        yield number, current_part


def read_basis_parts(part):
    current_part = []
    basis_type = None
    in_part = False
    for line in part:
        parts = line.upper().split()
        if parts[0] in basis_dict:
            if in_part:
                yield basis_type, current_part
                current_part = []
            basis_type = basis_dict[parts[0]]
            in_part = True
        elif in_part:
            basis_list = [float(i) for i in parts]
            current_part.append(basis_list)

    yield basis_type, current_part


def read_MO_parts(part):
    orbs = []
    orb = []
    sym_list = []
    energy_list = []
    spin_list = []
    occup_list = []

    in_part = False
    for line in part:
        line = [i.strip() for i in line.upper().split('=')]
        if line[0][0].isalpha():
            if line[0].startswith('SYM'):
                sym_list.append(line[1])
            elif line[0].startswith('ENE'):
                energy_list.append(float(line[1]))
            elif line[0].startswith('SPIN'):
                spin_list.append(line[1])
            elif line[0].startswith('OCC'):
                occup_list.append(float(line[1]))
            if in_part:
                in_part = False
                orbs.append(orb)
                orb = []
        else:
            in_part = True
            orb.append(float(line[0].split()[1]))

    orbs.append(orb)
    return orbs, sym_list, energy_list, spin_list, occup_list


def process_Atoms(part, mol):
    atoms_dict = {}
    for line in part:
        element, number, charge, x, y, z = line.split()
        number = int(number)
        charge = int(charge)
        coordinate = np.array([float(x), float(y), float(z)])
        atoms_dict[number] = [element, charge, coordinate]
    mol['Atoms'] = atoms_dict


def process_GTO(part, mol):
    gto_dict = {
        number:
        [[basis_type, line] for basis_type, line in read_basis_parts(part1)]
        for number, part1 in read_GTO_parts(part)
    }
    mol['GTO'] = gto_dict


def process_MO(part, mol):

    orbs, sym_list, energy_list, spin_list, occup_list = read_MO_parts(part)
    mol['orbs'] = np.array(orbs).T
    mol['sym'] = sym_list
    mol['energy'] = energy_list
    mol['spin'] = spin_list
    mol['occ'] = occup_list


def ReadMolden(filename):
    mol = {'is_Cartesian': False, 'Title': 'qwq'}
    valid_part = {
        'ATOMS': process_Atoms,
        'GTO': process_GTO,
        'MO': process_MO,
        'TITLE': process_title
    }
    d5f7 = ('5D', '7F', '9G')
    d6f10 = ('6D', '10F', '15G')
    with open(filename, 'r') as file:
        for name, part in read_parts(file):
            if name in d5f7:
                mol['is_Cartesian'] = False
            elif name in d6f10:
                mol['is_Cartesian'] = True
            elif name in valid_part:
                valid_part[name](part, mol)
    d5_ao_number = sum(
        basis[0] * 2 + 1 for atom in mol['GTO'].values() for basis in atom
    )
    d6_ao_number = sum(
        (basis[0] + 1) * (basis[0] + 2) // 2
        for atom in mol['GTO'].values()
        for basis in atom
    )
    if mol['orbs'].shape[0] == d5_ao_number:
        mol['is_Cartesian'] = False
    elif mol['orbs'].shape[0] == d6_ao_number:
        mol['is_Cartesian'] = True
    return mol


def WriteHead(file, title):
    file.write('[Molden Format]\n[Title]\n')
    if title == 'orca':
        file.write('Molden file created by orca_2mkl for BaseName=1\n')
    else:
        file.write(f'{title}\n')


def WriteAtoms(file, Atoms):
    file.write('[Atoms] AU\n')
    for number, atom in Atoms.items():
        element, charge, coordinate = atom
        file.write(
            f'{element:>3} {number:>3} {charge:>3} {coordinate[0]:>17.10f} {coordinate[1]:>17.10f} {coordinate[2]:>17.10f}\n'
        )


def WriteGTO(file, GTO):
    file.write('[GTO]\n')
    for number, basis_list in GTO.items():
        file.write(f'{number} 0\n')
        for basis_type, basis_line in basis_list:
            file.write(
                f' {reversed_basis_dict[basis_type].lower()} {len(basis_line)} 1.0\n'
            )
            for basis in basis_line:
                file.write(f'  {basis[0]:>15.10f} {basis[1]:>15.10f}\n')
        file.write('\n')


def WriteBasisType(file, is_Cartesian):
    if is_Cartesian:
        file.write('[6D]\n[10F]\n[15G]\n')
    else:
        file.write('[5D]\n[7F]\n[9G]\n')


def WriteMO(file, orbs):
    file.write('[MO]\n')
    for i in range(orbs.shape[1]):
        file.write('''Sym=     1a
Ene= 1.0
Spin= Alpha
Occup= 2.000000
''')
        for j in range(orbs.shape[0]):
            file.write(f'{j+1:>4}  {orbs[j,i]:>15.12f}\n')


def WriteMolden(mol, filename):
    with open(filename, 'w') as file:
        WriteHead(file, mol['Title'])
        WriteAtoms(file, mol['Atoms'])
        WriteGTO(file, mol['GTO'])
        WriteBasisType(file, mol['is_Cartesian'])
        WriteMO(file, mol['orbs'])

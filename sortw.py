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

import sys


def sort_vb_weights(filename, type1='w'):
    wei_type = {
        'w': 'WEIGHTS OF STRUCTURES',
        'l': 'Lowdin Weights',
        'i': 'Inverse Weights',
        'r': 'Renormalized Weights',
        'c': 'COEFFICIENTS OF STRUCTURES',
        'lc': 'LOWDIN ORTHOGONALIZED COEFFICIENTS OF STRUCTURES',
    }
    wei_type1 = wei_type[type1]
    with open(f'{filename}.xmo', 'r') as file:
        str_data = []
        wei_start = False

        for line in file:
            if wei_type1 in line:
                wei_start = True
                next(file)
                continue

            if wei_start:
                if line.strip() == '':
                    break
                parts = line.split()
                if len(parts) > 1:
                    weight = float(parts[1])
                    num = int(parts[0])
                    info = ''.join(f'{i:>4}' for i in parts[3:])
                    str1 = (weight, num, info)
                    str_data.append(str1)

    str_data.sort(reverse=True, key=lambda x: x[0])

    one = 0
    for i in range(len(str_data)):
        one += str_data[i][0]
        print(
            f'{i+1:>5}{str_data[i][1]:>5}{str_data[i][0]:^20.8f}{str_data[i][2]}'
              )
    print(one)


if __name__ == '__main__':

    if len(sys.argv) > 2:
        sort_vb_weights(sys.argv[1], sys.argv[2])
    else:
        sort_vb_weights(sys.argv[1])

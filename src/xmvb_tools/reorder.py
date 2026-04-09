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

import numpy as np


REORDER_PERMUTATIONS = {
    # new index -> old index
    # D shell: XX, XY, XZ, YY, YZ, ZZ
    #        ->XX, YY, ZZ, XY, XZ, YZ
    'D': (0, 3, 5, 1, 2, 4),
    # new index -> old index
    # F shell: XXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ
    #       -> XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ
    'F': (0, 6, 9, 3, 1, 2, 5, 8, 7, 4),
}

# XX 0  XX
# XY 1  YY 
# XZ 2  ZZ
# YY 3  XY
# YZ 4  XZ
# ZZ 5  YZ

# XXX 0  XXX
# XXY 1  YYY
# XXZ 2  ZZZ
# XYY 3  XYY
# XYZ 4  XXY
# XZZ 5  XXZ
# YYY 6  XZZ
# YYZ 7  YZZ
# YZZ 8  YYZ
# ZZZ 9  XYZ


def apply_permutation(matrix, start, permutation):
    block = matrix[start : start + len(permutation), :].copy()
    matrix[start : start + len(permutation), :] = block[list(permutation), :]
    return matrix


def invert_permutation(permutation):
    inverse = [0] * len(permutation)
    for new_index, old_index in enumerate(permutation):
        inverse[old_index] = new_index
    return tuple(inverse)


def reorder_shell(matrix, start, shell_type, inverse=False):
    permutation = REORDER_PERMUTATIONS[shell_type]
    if inverse:
        permutation = invert_permutation(permutation)
    return apply_permutation(matrix, start, permutation)

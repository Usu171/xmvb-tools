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

try:
    from .utils import ReadBasis, ReadGeoxmo, ReadNxmo, ReadOrb, WriteSimpleMolden
except ImportError:
    from utils import ReadBasis, ReadGeoxmo, ReadNxmo, ReadOrb, WriteSimpleMolden


def WriteMolden(file, file2, convert_basis=True):
    n = ReadNxmo(f'{file2}.xmo')
    matrix = ReadOrb(f'{file}.orb', n)
    gto_dict, matrix = ReadBasis(f'{file2}.xmo', matrix)
    geo = ReadGeoxmo(f'{file2}.xmo')
    WriteSimpleMolden(f'{file}.molden', geo, gto_dict, matrix, convert_basis=convert_basis)


def cli(argv=None):
    parser = argparse.ArgumentParser(description='Convert XMVB VB orbitals to a Molden file.')
    parser.add_argument('orb_prefix', help='Input orbital filename without the .orb suffix')
    parser.add_argument('xmo_prefix', help='Input filename without the .xmo suffix')
    parser.add_argument(
        '-b',
        action='store_true',
        help='Disable basis conversion.',
    )
    args = parser.parse_args(argv)
    WriteMolden(args.orb_prefix, args.xmo_prefix, convert_basis=not args.b)


if __name__ == '__main__':
    cli()

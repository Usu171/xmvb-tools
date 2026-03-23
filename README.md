# xmvb-tools

`xmvb-tools` is a package for converting XMVB orbital data to and from Molden-related formats.

It currently provides five command-line tools:
- `vb2molden`
- `gus2molden`
- `molden2gus`
- `no2molden`
- `sortw`


## Installation

```bash
pip install xmvb-tools
```

## Requirements

- Python 3.10+
- `numpy`

## Scripts

### `vb2molden`

Convert XMVB VB orbitals from `.orb` , basis/geometry information from `.xmo` into a Molden file.

```bash
vb2molden <orb_prefix> <xmo_prefix>
```

Example:

```bash
vb2molden C6H6_xmvb4 C6H6_xmvb4
```

Output:
- `<orb_prefix>.molden`

Options:
- `-b`: disable basis conversion before writing the Molden GTO section


### `gus2molden`

Convert XMVB guess orbitals embedded in an `.xmo` file into a Molden file.

```bash
gus2molden <xmo_prefix>
```

Example:

```bash
gus2molden C6H6_xmvb4
```

Output:
- `<xmo_prefix>_vbgus.molden`

Options:
- `-b`: disable basis conversion before writing the Molden GTO section


### `molden2gus`

Read a Molden orbital file and interactively build an XMVB guess orbital file.

```bash
molden2gus <molden_prefix>
```

Example:

```bash
molden2gus SO4_lo
```

Output:
- `<molden_prefix>_gus.molden`
- `<molden_prefix>.gus`

Interactive commands:
- `<atoms> <orbitals>`: extract selected orbitals on selected atoms
- `a <orbitals>`: add whole orbitals for all atoms
- `r<angle> <orb1>,<orb2>`: rotate two orbitals by the given angle in degrees
- `m<number> <orbital>`: multiply one orbital by a scalar
- `q`: write output files and exit

Examples:
- `2,4 6-8,9`
- `a 1-5`
- `r45 1,2`
- `m-1 4`

You can also redirect commands from a text file:

```bash
molden2gus xx < xx.txt
```

### `no2molden`

Convert natural orbitals embedded directly in a modern XMVB `.xmo` file into a Molden file.

```bash
no2molden <xmo_prefix>
```

Example:

```bash
no2molden C6H6_xmvb4
```

Output:
- `<xmo_prefix>_no.molden`

Options:
- `-b`: disable basis conversion before writing the Molden GTO section


### `sortw`

Sort and print XMVB structure weights or coefficients from an `.xmo` file.

```bash
sortw <xmo_prefix> [weight_type]
```

Example:

```bash
sortw C6H6_xmvb4 w
```

Supported `weight_type` values:
- `w`: `WEIGHTS OF STRUCTURES`
- `l`: `Lowdin Weights`
- `i`: `Inverse Weights`
- `r`: `Renormalized Weights`
- `c`: `COEFFICIENTS OF STRUCTURES`

default: `w`


## License

GPL-3.0-or-later

# PseudoPotentialIO.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://azadoks.github.io/PseudoPotentialIO.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://azadoks.github.io/PseudoPotentialIO.jl/dev)
[![Build Status](https://github.com/azadoks/PseudoPotentialIO.jl/workflows/CI/badge.svg)](https://github.com/azadoks/PseudoPotentialIO.jl/actions)
[![Coverage](https://codecov.io/gh/azadoks/PseudoPotentialIO.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/azadoks/PseudoPotentialIO.jl)

PseudoPotentialIO aims to provide parsers for common pseudopotential file formats used by density functional theory codes and an interface for accesssing the quantities that they contain.

## Reading pseudopotentials

The following file formats are supported or have planned (or no planned) support.
If your favorite format does not appear in the table below, please file an [issue](https://github.com/azadoks/PseudoPotentialIO.jl/issues)!

| Format                   | Read        |
|--------------------------|-------------|
| UPF (old)                | ✅          |
| UPF (v2.0.1 with schema) | ✅          |
| PSP1                     | Not planned |
| PSP3/HGH                 | Partial     |
| PSP4/Teter               | Not planned |
| PSP5/"Phoney"            | Not planned |
| PSP6/FHI98PP             | Not planned |
| PSP7/ABINIT PAW          | Not planned |
| PSP8/ONCVPSP             | ✅          |
| PSP9/PSML                | Not planned |
| PSP17/ABINIT PAW XML     | Not planned |
| FPMD XML                 | Not planned |
| Vanderbilt USPP          | Not planned |

Pseudopotentials are read by calling `load_psp_file(path)`.
This returns one of the `[Format]File` structs which mirror very closely the contents of the file.
However, different file formats provide important physical quantities in slightly different forms.
Be sure to check the docstrings and comments here in PseudoPotentialIO.jl or the relevant documentation of the file format.

## Writing pseudopotentials
This package also supports writing pseudopotentials, currently only in the UPF v2.0.1 format.
This is done using `save_psp_file`:
```jl
save_psp_file("path/to/pseudo.upf", pseudo, "UPF v2.0.1")
```
Support for additional formats might be added in the future.

## Converting pseudopotentials
Conversion from PSP8 to UPF is implemented as well,
by calling the `UpfFile` constructor with a `Psp8File`.
Conversions between other formats might be added in the future.

Here is a full example of reading a PSP8 file, converting it to a UPF file and then writing it:
```jl
psp8 = load_psp_file("path/to/pseudo.psp8")
upf = UpfFile(psp8)
save_psp_file("path/to/pseudo.upf", upf, "UPF v2.0.1")
```

# Polymerization

[![C++23](https://img.shields.io/badge/C%2B%2B-23-00599C?logo=cplusplus&logoColor=white)](https://en.cppreference.com/w/cpp/23)
[![CMake](https://img.shields.io/badge/CMake-%E2%89%A53.20-064F8C?logo=cmake&logoColor=white)](https://cmake.org)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](LICENSE)

## Introduction

**A versatile simulation toolkit for polymerization and polymer-network formation.**

The way a polymer is made fixes its architecture, and with it the properties of the final
material. This toolkit follows that architecture as it forms, reaction by reaction, across a
wide range of polymerization chemistries and at several levels of description, from purely
statistical to lattice and coarse-grained models. In the spatially resolved models, molecules
move and react with their neighbors through the elementary steps of each mechanism, so effects
that simple theories leave out, such as cyclization, steric hindrance and diffusion control,
come out of the simulation instead of being assumed.

Our aim is one dependable code that researchers can pick up to simulate the polymerization in
their own work.

## Functions

| Model | Description | Module | Example | Status |
|:--|:--|:--|:--:|:--:|
| On-lattice: 3D bond fluctuation model | reactive Monte Carlo on a cubic lattice; free radical (FRP) and atom transfer radical (ATRP) polymerization | [**3dbfm**](src/3dbfm/README.md) | [example_1](examples/example_1/README.md) | ![available](https://img.shields.io/badge/available-brightgreen) |
| Random-walk packing | precursor chains (self-avoiding random walks on an FCC lattice) and cross-linkers packed into a box as the starting point of end-linking simulations in LAMMPS | [**randomwalk**](src/randomwalk/README.md) | [example_2](examples/example_2/README.md) | ![available](https://img.shields.io/badge/available-brightgreen) |
| On-lattice: dynamic lattice liquid model | reactive Monte Carlo with cooperative moves on a fully occupied lattice, with solvent; free radical (FRP) and atom transfer radical (ATRP) polymerization | [**dll**](src/dll/README.md) | [example_3](examples/example_3/README.md) | ![available](https://img.shields.io/badge/available-brightgreen) |
| Pure statistical model | coordinate-free stochastic simulation of the reaction kinetics | **statistical** | — | ![in preparation](https://img.shields.io/badge/in%20preparation-orange) |
| Coarse-grained model | reactive coarse-grained molecular dynamics as a fix for LAMMPS | **lammps** | — | ![in preparation](https://img.shields.io/badge/in%20preparation-orange) |

Each module keeps its code and README under `src/<module>/`, inputs under
`examples/example_<n>/` and tests under `tests/<module>/`. Code shared by the modules, such as the
random number generator, the configuration-file reader and the FRP/ATRP chemistry of the lattice
models, is in `src/common/`. All modules are built and tested together from the repository root.

## Build

Requires CMake ≥ 3.20 and a C++23 compiler, e.g. GCC 14+ or Clang 18+.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

Options: `-DPOLY_CXX_STANDARD=26` (build as C++26), `-DPOLY_NATIVE_ARCH=ON` (`-march=native`),
`-DPOLY_LTO=OFF` (no link-time optimization), `-DPOLY_BUILD_TESTS=OFF` (skip tests). Executables
are placed in `build/src/<module>/`.

Profile-guided optimization makes the simulations about 10 % faster. With Clang:

```bash
cmake -S . -B build-pgo -DPOLY_PGO=GENERATE && cmake --build build-pgo -j
LLVM_PROFILE_FILE=pgo/%p.profraw build-pgo/src/3dbfm/bfm_react examples/example_1/smoke_small.ini \
    stop_conversion=0.3 directory=pgo/run
llvm-profdata merge -o pgo/run.profdata pgo/*.profraw    # on macOS: xcrun llvm-profdata
cmake -S . -B build -DPOLY_PGO=USE -DPOLY_PGO_PROFILE=$PWD/pgo/run.profdata && cmake --build build -j
```

## How to cite

This code is a refined and restructured version of the code used in the following papers. Please
cite them when you use it:

- Z. Zhang, J. Krajniak, A. Z. Dookhith, Y. Tian, H. S. Sachar, N. Marioni, T. J. Duncan, J. Liu,
  G. E. Sanoja, V. Ganesan. Topology and Mechanical Properties of Polymer Networks Formed under
  Free Radical and Atom Transfer Radical Polymerizations. *Macromolecules* **2025**, 58,
  3168–3187. [10.1021/acs.macromol.4c03138](https://doi.org/10.1021/acs.macromol.4c03138)
- Z. Zhang, J. Krajniak, A. Z. Dookhith, H. Zhang, Y. Tian, H. S. Sachar, N. Marioni, T. J. Duncan,
  J. Liu, G. E. Sanoja, V. Ganesan. Fracture Behavior of Monodisperse and Polydisperse Polymer
  Networks: Effects of Topology and Strand-Length Dispersity. *Macromolecules* **2026**, 59,
  99–119. [10.1021/acs.macromol.5c01715](https://doi.org/10.1021/acs.macromol.5c01715)
- A. Z. Dookhith, Z. Zhang, V. Ganesan, G. E. Sanoja. Impact of Reversible Deactivation Radical
  Copolymerizations (RDRPs) on Gelation, Phase Separation, and Mechanical Properties of Polymer
  Networks. *Macromolecules* **2024**, 57, 8698–8711.
  [10.1021/acs.macromol.4c00905](https://doi.org/10.1021/acs.macromol.4c00905)
- A. Z. Dookhith, Z. Zhang, V. Ganesan, G. E. Sanoja. Designing soft and tough multiple-network
  elastomers: impact of reversible radical deactivation on filler network architecture and
  fracture toughness. *Soft Matter* **2025**, 21, 4029–4042.
  [10.1039/D5SM00045A](https://doi.org/10.1039/D5SM00045A)

## Acknowledgements

The original 3D bond fluctuation model code was obtained from the group of **Prof. Xuehao He,
School of Science, Tianjin University, China**, where it was developed for the following
studies:

1. L. Wang, X. He. Investigation of AB<sub>n</sub> (n = 2, 4) type hyperbranched polymerization
   with cyclization and steric factors: influences of monomer concentration, reactivity, and
   substitution effect. *J. Polym. Sci., Part A: Polym. Chem.* **2009**, 47, 523–533.
   [10.1002/pola.23169](https://doi.org/10.1002/pola.23169)
2. X. Yang, L. Wang, X. He. Kinetics of nonideal hyperbranched A<sub>2</sub> + B<sub>3</sub>
   polycondensation: simulation and comparison with experiments. *J. Polym. Sci., Part A: Polym.
   Chem.* **2010**, 48, 5072–5082. [10.1002/pola.24305](https://doi.org/10.1002/pola.24305)
3. X. He, J. Tang. Kinetics of self-condensing vinyl hyperbranched polymerization in
   three-dimensional space. *J. Polym. Sci., Part A: Polym. Chem.* **2008**, 46, 4486–4494.
   [10.1002/pola.22784](https://doi.org/10.1002/pola.22784)

The original random-walk packing code was obtained from **Mark A. Tschopp and Don K. Ward,
Center for Advanced Vehicular Systems, Mississippi State University, USA** (*Amorphous Polymer
Generator*, published on the center's
[Integrated Computational Materials Engineering website](https://www.cavs.msstate.edu/icme/code/lammps/tutorials/processing/amorphousGenerator.php)),
where it was used for the following studies:

4. D. Hossain, M. A. Tschopp, D. K. Ward, J. L. Bouvard, P. Wang, M. F. Horstemeyer. Molecular
   dynamics simulations of deformation mechanisms of amorphous polyethylene. *Polymer* **2010**,
   51, 6071–6083. [10.1016/j.polymer.2010.10.009](https://doi.org/10.1016/j.polymer.2010.10.009)
5. M. A. Tschopp, J. L. Bouvard, D. K. Ward, M. F. Horstemeyer. Atomic scale deformation
   mechanisms of amorphous polyethylene under tensile loading. In *Supplemental Proceedings:
   Materials Fabrication, Properties, Characterization, and Modeling* (TMS 2011); Wiley, **2011**;
   pp 789–794. [10.1002/9781118062142.ch95](https://doi.org/10.1002/9781118062142.ch95)

The dll module implements, from the published descriptions, the dynamic lattice liquid model of
**Prof. Tadeusz Pakula** and its application to polymerization by **Piotr Polanowski, Jeremiasz
K. Jeszka and co-workers, Lodz University of Technology, Poland**, with whom we used the model
in the following study (further references in the [module README](src/dll/README.md)):

6. J. Lyu, Y. Gao, Z. Zhang, U. Greiser, P. Polanowski, J. K. Jeszka, K. Matyjaszewski, H. Tai,
   W. Wang. Monte Carlo Simulations of Atom Transfer Radical (Homo)polymerization of Divinyl
   Monomers: Applicability of Flory–Stockmayer Theory. *Macromolecules* **2018**, 51, 6673–6681.
   [10.1021/acs.macromol.8b01630](https://doi.org/10.1021/acs.macromol.8b01630)

The reactive extensions and this version were developed in the group of Prof. Venkat Ganesan,
The University of Texas at Austin, with support from the Robert A. Welch Foundation (F-1599) and
the National Science Foundation (DMR-2225167).

## License

[PolyForm Noncommercial 1.0.0](LICENSE). Copyright (C) 2023–2026 Zidan Zhang.

The code is free for noncommercial use: research, teaching and personal study, and any use by
universities, public research institutions, government laboratories and other noncommercial
organizations, whatever their source of funding. Commercial use requires a separate license;
please contact Zidan Zhang (zhangzidan@gmail.com).

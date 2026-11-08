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
| Pure statistical model | coordinate-free kinetic Monte Carlo (Gillespie) of the reactions, with or without cyclization; free radical (FRP) and atom transfer radical (ATRP) polymerization | [**statistical**](src/statistical/README.md) | [example_4](examples/example_4/README.md) | ![available](https://img.shields.io/badge/available-brightgreen) |
| Coarse-grained model | reactive coarse-grained molecular dynamics as a fix for LAMMPS (`fix poly/react`): general reaction rules, starting with free radical (FRP) and atom transfer radical (ATRP) polymerization | [**lammps**](src/lammps/README.md) | [example_5](examples/example_5/README.md) | ![on-going](https://img.shields.io/badge/on--going-yellow) |

Each module keeps its code and README under `src/<module>/`, inputs under
`examples/example_<n>/` and tests under `tests/<module>/`; the LAMMPS fix, licensed separately, keeps
its tests in `src/lammps/`. Code shared by the modules, such as the
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

The structure of the LAMMPS fix borrows the ideas of **ChemLab**, the reactive coarse-grained
molecular dynamics code of **Jakub Krajniak, Zidan Zhang, Erik Nies and Giovanni Samaey, KU
Leuven, Belgium**
([cgchemlab/chemlab](https://github.com/cgchemlab/chemlab),
[10.5281/zenodo.817589](https://doi.org/10.5281/zenodo.817589)), developed for studies 8–10 below,
and of the earlier coarse-grained molecular dynamics – reactive Monte Carlo approach of study 7:

7. Z. Zhang, L. Wang, Z. Wang, X. He, Y. Chen, F. Müller-Plathe, M. C. Böhm. A coarse-grained
   molecular dynamics – reactive Monte Carlo approach to simulate hyperbranched polycondensation.
   *RSC Adv.* **2014**, 4, 56625–56636. [10.1039/C4RA10271A](https://doi.org/10.1039/C4RA10271A)
8. J. Krajniak, Z. Zhang, S. Pandiyan, E. Nies, G. Samaey. Reverse mapping method for complex
   polymer systems. *J. Comput. Chem.* **2018**, 39, 648–664.
   [10.1002/jcc.25129](https://doi.org/10.1002/jcc.25129)
9. J. Krajniak, Z. Zhang, S. Pandiyan, E. Nies, G. Samaey. Coarse-grained molecular dynamics
   simulations of polymerization with forward and backward reactions. *J. Comput. Chem.* **2018**,
   39, 1764–1778. [10.1002/jcc.25348](https://doi.org/10.1002/jcc.25348)
10. Z. Zhang, J. Krajniak, G. Samaey, E. Nies. A parallel multiscale simulation framework for
    complex polymerization: AB₂-type monomer hyperbranched polymerization as an example. *Adv.
    Theory Simul.* **2019**, 2, 1800102. [10.1002/adts.201800102](https://doi.org/10.1002/adts.201800102)

The reactive extensions and this version were developed in the group of Prof. Venkat Ganesan,
McKetta Department of Chemical Engineering, The University of Texas at Austin.

## License

[PolyForm Noncommercial 1.0.0](LICENSE). Copyright (C) 2023–2026 Zidan Zhang (zhangzidan@gmail.com) and Jakub
Krajniak (jkrajniak@gmail.com).

The code is free for noncommercial use: research, teaching and personal study, and any use by
universities, public research institutions, government laboratories and other noncommercial
organizations, whatever their source of funding. Commercial use requires a separate license;
please contact Zidan Zhang and Jakub Krajniak.

The LAMMPS fix in [`src/lammps`](src/lammps) is the exception: it is licensed under the GNU
General Public License, version 2, the license of LAMMPS ([src/lammps/LICENSE](src/lammps/LICENSE)),
and uses no code from the rest of the repository.

# Structured Analytic Mappings for Point Set Registration

This repository provides a C++ research prototype implementation of the paper:

**"Structured Analytic Mappings for Point Set Registration"**

The main registration algorithm implemented in this repository is referred to as **Analytic-ICP**. It embeds **Structured Analytic Mappings (SAM)** into an ICP-style hard-correspondence point set registration loop.

The current release focuses on the core algorithmic pipeline and selected reproducibility examples. It is intended as a research prototype accompanying the paper, not as a fully optimized general-purpose registration library.

---

## Overview

Structured Analytic Mappings (SAM) represent smooth point-set deformations by truncated multivariate Taylor mappings of vector-valued functions. Instead of using a point-indexed displacement field, SAM organizes deformation estimation in a finite-dimensional analytic function space whose coefficient dimension is controlled by the ambient dimension and analytic order.

In the Analytic-ICP implementation, the registration process follows an ICP-style hard-correspondence framework:

1. Estimate nearest-neighbor correspondences between the current moving point set and the fixed point set;
2. Fit a structured analytic mapping to the corresponding point pairs;
3. Update the moving point set by applying the estimated analytic map;
4. Repeat until convergence or until the maximum number of iterations is reached.

This framework provides a compact and interpretable deformation representation for smooth non-rigid point set registration.

---

## Project Structure

The current Visual Studio implementation follows a three-layer architecture:

- **`Analytic_ICP.dll`**  
  Core Analytic-ICP algorithm. It implements structured analytic mapping based 2D/3D registration and synthetic analytic deformation generation.

- **`SmoothAdjustment.dll`**  
  A lightweight file-level wrapper. It loads `Analytic_ICP.dll`, reads point sets from CSV files, invokes the 2D/3D registration interfaces, optionally generates analytic perturbations, and writes transformed point sets to CSV.

- **`Test.exe`**  
  A console frontend. It reads `experiment.ini`, initializes the DLLs, runs 2D or 3D registration, and writes output files.

A typical runtime directory should contain:

```text
Test.exe
SmoothAdjustment.dll
Analytic_ICP.dll
experiment.ini
```
---

## Build Environment

The current implementation is a Windows / Visual Studio research prototype.

Tested environment:

- **Operating system:** Windows
- **Compiler/IDE:** Microsoft Visual Studio
- **Build mode:** 64-bit Release build recommended
- **Dependencies:**
  - Eigen for matrix and vector computation
  - Boost for KD-tree based nearest-neighbor search
  - OpenCV for optional visualization and image output in the test/frontend modules

Linux, macOS, and CMake support are not yet provided in the current release.

---

## Implementation Notes

The ICP component is adapted from the public point-to-point ICP implementation by Andreas Geiger (2011). In this repository, the original code has been reorganized and converted into an Eigen-based implementation so that ICP-style nearest-neighbor registration and structured analytic mapping fitting can be evaluated within a consistent C++ numerical framework.

The nearest-neighbor component depends on Boost KD-tree functionality. Eigen is used for matrix and vector operations in the analytic fitting and registration modules. OpenCV is used mainly by the test and wrapper components for optional visualization and image output; the core analytic fitting logic itself is based on C++ and Eigen.

The implementation is single-threaded unless otherwise noted.

---

## Configuration File

The program uses `experiment.ini` to control the experiment mode, input paths, and output paths.

A typical configuration is:

```ini
[Param]
; 0: 2D registration, 1: 3D registration
RegistType=1

; 1: generate analytic perturbation before registration, 0: direct registration
AddPerturb=1

; input point sets
MovingPsPath=examples/3d/moving.csv
FixedPsPath=examples/3d/fixed.csv

; output point sets
PerturbedOutputPath=results/perturbed3d.csv
MovedOutputPath=results/moved3d.csv

; perturbation degrees used only when AddPerturb=1
PerturbDeg2D=8
PerturbDeg3D=2
```

---

## Parameter Description

| Parameter | Description |
|---|---|
| `RegistType` | `0` for 2D registration, `1` for 3D registration |
| `AddPerturb` | Whether to generate a synthetic analytic deformation before registration |
| `MovingPsPath` | Path to the moving/source point set |
| `FixedPsPath` | Path to the fixed/target point set |
| `PerturbedOutputPath` | Output path for the synthetically deformed point set |
| `MovedOutputPath` | Output path for the registered moving point set |
| `PerturbDeg2D` | Analytic degree for synthetic 2D deformation generation |
| `PerturbDeg3D` | Analytic degree for synthetic 3D deformation generation |

---

## Quick Start

1. Open the Visual Studio solution.

2. Build the following projects in **Release x64** mode:

   ```text
   Analytic_ICP.dll
   SmoothAdjustment.dll
   Test.exe
   ```

3. Place the following files in the same runtime directory:

   ```text
   Test.exe
   SmoothAdjustment.dll
   Analytic_ICP.dll
   experiment.ini
   ```

4. Edit `experiment.ini` to specify the fixed point set, moving point set, and output paths.

5. Run:

   ```bash
   Test.exe
   ```

After execution, the registered moving point set will be written to the path specified by `MovedOutputPath`.

If `AddPerturb=1`, the synthetically deformed point set will also be written to `PerturbedOutputPath`.

---

## Input Data Format

Input point sets are stored as CSV files without headers.

For 2D point sets:

```text
x,y
x,y
...
```

For 3D point sets:

```text
x,y,z
x,y,z
...
```

Each row corresponds to one point.

For pointwise error evaluation, the fixed and moving point sets should contain the same number of points and follow the same point ordering. For general registration without known pointwise correspondence, the output registered point set can still be visually inspected or evaluated using nearest-neighbor based metrics.

---

## Output Files

The program can produce the registered moving point set:

```text
results/moved3d.csv
```

or the path specified by `MovedOutputPath`.

When `AddPerturb=1`, the program can also produce the synthetically deformed point set:

```text
results/perturbed3d.csv
```

or the path specified by `PerturbedOutputPath`.

Optional visualization files may be generated by visualization-enabled modules when OpenCV output is enabled.

---

## Runtime Notes

Unless otherwise noted in the paper, experiments were run on a standard laptop with an Intel i5-6200U CPU and 8 GB RAM. The implementation is single-threaded C++ using Eigen, Boost, and OpenCV.

Wall-clock times reported in the paper should be compared within the same experimental subsection. Some experiments may have been executed on different hardware, which is explicitly indicated in the corresponding captions or subsections of the paper.

---

## Current Status

This repository is currently a research prototype accompanying the paper.

Current characteristics:

- C++ implementation
- Windows / Visual Studio project
- Single-threaded core implementation
- Eigen-based matrix and vector computation
- Boost-based KD-tree nearest-neighbor search
- OpenCV-enabled optional visualization frontend
- 2D and 3D Analytic-ICP registration support
- Synthetic analytic deformation generation for selected examples

The code is being organized for reproducibility. Additional examples, scripts, and documentation will be added progressively.

---

## Datasets

The paper uses public point-set and point-cloud datasets as well as generated deformation examples.

Due to dataset license restrictions, third-party datasets should be downloaded from their original sources when necessary. This repository provides selected example files and configuration templates for demonstrating the core algorithm.

Please follow the licenses and citation requirements of the original datasets.

---

## Relationship to Analytic-CPD

This repository implements the Structured Analytic Mapping framework and its hard-correspondence instantiation, Analytic-ICP.

A related follow-up project, Analytic-CPD, integrates Structured Analytic Mappings into the probabilistic posterior correspondence framework of Coherent Point Drift (CPD). In contrast to Analytic-ICP, which relies on nearest-neighbor hard correspondences, Analytic-CPD uses CPD-style soft posterior correspondences and formulates the M-step as posterior-weighted structured analytic mapping estimation.

Analytic-CPD repository:

```text
https://github.com/monge-ampere/Analytic-CPD
```

---

## Citation

If you use this code or the Structured Analytic Mapping deformation model, please cite the associated paper:

```bibtex
@article{feng2026structuredanalyticmappings,
  title={Structured Analytic Mappings for Point Set Registration},
  author={Feng, Wei and Wei, Tengda and Zheng, Haiyong},
  journal={arXiv preprint arXiv:2602.16753},
  year={2026},
  doi={10.48550/arXiv.2602.16753},
  note={Accepted for publication in SIAM Journal on Imaging Sciences}
}
```

If you use the Analytic-CPD extension, please also cite:

```bibtex
@article{feng2026structuredanalyticcpd,
  title={Structured Analytic Coherent Point Drift for Non-Rigid Point Set Registration},
  author={Feng, Wei and Zheng, Haiyong},
  journal={arXiv preprint arXiv:2605.00934},
  year={2026},
  doi={10.48550/arXiv.2605.00934}
}
```

---

## License

A formal open-source license will be added after the third-party dependency and code status is fully clarified.

Until then, please treat this repository as a research prototype made publicly available for academic reference. For redistribution, commercial use, or integration into other projects, please contact the authors and review all third-party code notices carefully.

The ICP component is adapted from the public point-to-point ICP implementation by Andreas Geiger (2011). Please review and respect the license terms of the original implementation and all third-party dependencies.

---

## Contact

For questions, comments, or collaboration related to Structured Analytic Mappings, Analytic-ICP, point set registration, or geometric deformation modeling, please contact the author through GitHub.

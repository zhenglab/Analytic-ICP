# Analytic-ICP

[Published paper](https://doi.org/10.1137/25M1752080) |
[arXiv preprint](https://arxiv.org/abs/2602.16753) |
[Analytic-CPD](https://github.com/monge-ampere/Analytic-CPD)

This repository provides the C++ research prototype accompanying:

> Wei Feng, Tengda Wei, and Haiyong Zheng, **Structured Analytic Mappings for Point Set Registration**, *SIAM Journal on Imaging Sciences*, 19 (2026), pp. 1655â€“1697.
> <https://doi.org/10.1137/25M1752080>

**Analytic-ICP** is the hard-correspondence point-set registration algorithm developed in the paper. It embeds **Structured Analytic Mappings (SAM)** into an ICP-style outer loop. SAM represents smooth deformations in a finite-dimensional analytic function space constructed from truncated multivariate Taylor expansions of vector-valued functions.

The repository contains the core 2D/3D implementation, a file-based DLL wrapper, a Windows console frontend, selected point sets, and one configured 3D example. It is a research prototype and source snapshot associated with the journal article, not a production registration library or a one-command reproduction package for every experiment in the paper.

## Method in Brief

Analytic-ICP alternates between nearest-neighbor correspondence estimation and structured analytic mapping estimation:

1. establish hard nearest-neighbor correspondences between the current moving set and the fixed set;
2. estimate a structured analytic map for the corresponding point pairs;
3. transform the moving set;
4. repeat until the internal stopping condition is reached.

The mapping is estimated through staged rigid, affine, and nonlinear analytic fitting. Its coefficient dimension is controlled by the ambient dimension and analytic order rather than by the number of input points.

The current hard-correspondence implementation is intended primarily for small and smooth non-rigid deformations. As discussed in the paper, large or nonlocal displacements can cause nearest-neighbor correspondence failures even when the analytic deformation model is sufficiently expressive. The related [Analytic-CPD](https://github.com/monge-ampere/Analytic-CPD) project combines SAM with soft posterior correspondences.

## Repository Layout

| Path | Purpose |
|---|---|
| `Analytic_ICP/` | Core 2D/3D Analytic-ICP fitting and synthetic analytic deformation functions |
| `SmoothAdjustment/` | File-level DLL wrapper around `Analytic_ICP.dll` |
| `Test/` | Console frontend, `experiment.ini`, one selected 3D example, and example outputs |
| `Data/` | Point sets retained from the research code; provenance review is still in progress |
| `SmoothAdjustment.sln` | Visual Studio solution containing all three projects |
| `CITATION.cff` | Machine-readable citation metadata for the published article |
| `API.md` | Exported DLL interfaces and call sequence |
| `REPRODUCIBILITY.md` | Build and selected-example reproducibility notes |
| `DATASETS.md` | Inventory and current provenance status of included data |
| `THIRD_PARTY_NOTICES.md` | Third-party source and licensing notices |

The runtime architecture is:

```text
Test.exe
  â””â”€ loads SmoothAdjustment.dll
       â””â”€ loads Analytic_ICP.dll
```

`Test.exe` resolves the wrapper functions dynamically with `GetProcAddress`; `SmoothAdjustment.dll` similarly resolves the core registration and deformation functions from `Analytic_ICP.dll`.

## Build Requirements

The checked-in project files describe the following legacy Windows environment:

| Component | Version or configuration encoded in the project |
|---|---|
| Operating system | Windows |
| Visual Studio solution | Visual Studio 14 / Visual Studio 2015 |
| Platform toolset | `v140` |
| Recommended configuration | `Release | x64` |
| Eigen | 3.3.8 path encoded |
| Boost | 1.60.0 path encoded |
| OpenCV | 3.4.16 libraries encoded as `opencv_world3416.lib` |

These values are read directly from the current `.sln` and `.vcxproj` files; a clean-machine build has not yet been revalidated for this publication snapshot.

> **Important:** the three `.vcxproj` files currently contain machine-specific absolute paths under `E:\ScientificComputingPackage\...`. Before building, replace those include and library paths with locations valid on your machine.

Although the analytic fitting code itself is based mainly on C++ and Eigen, the current projects include and/or link OpenCV. OpenCV is therefore a configured build dependency, not yet an optional dependency; release binaries may also require the matching OpenCV DLL when imports are retained. Boost is used by the KD-tree nearest-neighbor implementation.

Linux, macOS, CMake, and package-manager-based dependency setup are not provided in this snapshot.

## Quick Start

1. Clone the repository and open `SmoothAdjustment.sln`.

2. Update the Eigen, Boost, and OpenCV paths in the three `.vcxproj` files, or configure equivalent Visual Studio property sheets locally.

3. Select `Release | x64` and build:

   ```text
   Analytic_ICP
   SmoothAdjustment
   Test
   ```

4. Place the following files in one runtime directory:

   ```text
   Test.exe
   SmoothAdjustment.dll
   Analytic_ICP.dll
   experiment.ini
   ```

   If OpenCV is dynamically linked, ensure that `opencv_world3416.dll` is also in the runtime directory or available through `PATH`.

5. Use `Test/` as the working directory, or copy its `experiment.ini`, `examples/`, and `results/` layout to the runtime directory.

6. Run:

   ```text
   Test.exe
   ```

The selected example writes the registered moving point set to the path specified by `MovedOutputPath`. Additional build and verification details are provided in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Configuration

`Test.exe` reads `.\experiment.ini`. The checked-in configuration is:

```ini
[Param]
; 0: 2D registration, 1: 3D registration
RegistType=1

; 1: also generate a synthetic analytic deformation, 0: do not generate it
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

| Parameter | Current frontend behavior |
|---|---|
| `RegistType` | `0` selects the 2D wrapper and `1` selects the 3D wrapper |
| `AddPerturb` | Generates an additional synthetic deformation file when set to `1` |
| `MovingPsPath` | Moving/source point set passed to registration |
| `FixedPsPath` | Fixed/target point set passed to registration |
| `PerturbedOutputPath` | Destination of the synthetic analytic deformation |
| `MovedOutputPath` | Destination of the registered moving point set |
| `PerturbDeg2D` | Degree supplied to the 2D synthetic deformation function |
| `PerturbDeg3D` | Degree supplied to the 3D synthetic deformation function |

### Current `AddPerturb` Behavior

In the checked-in frontend, `AddPerturb=1` creates `PerturbedOutputPath`, but the registration call in the same run still reads `MovingPsPath`. The generated perturbation file is not automatically substituted as the registration input.

To register a generated perturbation without changing the code, first generate it, then run the program again with `MovingPsPath` set to that output file. This documents the current implementation; whether the single-run behavior should be changed will be decided separately from this documentation update.

The synthetic deformation functions seed their random-number generator from the current time, so perturbation output is expected to vary between runs.

## Point-Set File Format

Input files contain one point per row and no header. The current parser accepts comma-separated or whitespace-separated numeric coordinates.

2D examples:

```text
x,y
x,y
```

or

```text
x y
x y
```

3D examples:

```text
x,y,z
x,y,z
```

The output writer uses comma-separated coordinates. The registered output preserves the moving-set row order.

Known pointwise correspondence is not required by the nearest-neighbor registration loop. Equal point counts and matching row order are needed only when evaluating pointwise errors against known correspondences.

## Outputs

The checked-in example configuration uses:

```text
Test/results/moved3d.csv
Test/results/perturbed3d.csv
```

`moved3d.csv` is the registered version of `MovingPsPath`. When `AddPerturb=1`, `perturbed3d.csv` is also generated, subject to the behavior described above.

The console prints the scalar value returned by the registration function as `Final error`. This prototype does not yet emit a structured experiment log containing all parameters, timing, and evaluation metrics.

## Reproducibility Scope

The repository currently provides:

- the Windows C++ research implementation;
- 2D and 3D registration entry points;
- synthetic analytic deformation entry points;
- one configured 3D example with stored outputs;
- additional research point sets under `Data/`.

It does not yet provide:

- a portable dependency configuration;
- an automated clean-machine build;
- scripts reproducing every table and figure in the paper;
- automated numerical regression tests;
- completed provenance and redistribution records for every included point set.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the precise status and [DATASETS.md](DATASETS.md) before redistributing any data.

## Runtime and Performance Notes

The implementation is single-threaded unless otherwise noted. The paper reports experiments conducted on the hardware stated in the corresponding experimental subsections and captions. Wall-clock values should be compared within the same experimental setting rather than across different machines or separately described hardware.

The legacy prototype was also run on a laptop with an Intel i5-6200U CPU and 8 GB RAM. This machine description is retained as implementation history, not as a guarantee that all paper experiments used that hardware.

## Citation

If you use this code or the Structured Analytic Mapping deformation model, cite the published article:

```bibtex
@article{feng2026structuredanalyticmappings,
  author  = {Feng, Wei and Wei, Tengda and Zheng, Haiyong},
  title   = {Structured Analytic Mappings for Point Set Registration},
  journal = {SIAM Journal on Imaging Sciences},
  volume  = {19},
  number  = {3},
  pages   = {1655--1697},
  year    = {2026},
  doi     = {10.1137/25M1752080},
  url     = {https://doi.org/10.1137/25M1752080}
}
```

GitHub's citation interface reads the same journal metadata from [`CITATION.cff`](CITATION.cff). The [arXiv version](https://arxiv.org/abs/2602.16753) remains available as a preprint.

If you use the soft-correspondence extension, also cite its associated preprint:

```bibtex
@article{feng2026structuredanalyticcpd,
  author  = {Feng, Wei and Zheng, Haiyong},
  title   = {Structured Analytic Coherent Point Drift for Non-Rigid Point Set Registration},
  journal = {arXiv preprint arXiv:2605.00934},
  year    = {2026},
  doi     = {10.48550/arXiv.2605.00934},
  url     = {https://arxiv.org/abs/2605.00934}
}
```

## Licensing and Third-Party Code

No project-wide license has yet been declared for the original Analytic-ICP code.

The repository also contains adapted third-party ICP and KD-tree source with pre-existing license notices. In particular, the libicp-derived files state **GNU GPL version 2 or later**, while the KD-tree files refer to the **Academic Free License version 1.1** and additional provisions in a license file that is not currently present in this repository.

Before reuse or redistribution, read [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Public source availability should not be interpreted as permission to redistribute the entire repository under an arbitrary license.

## Contact

Questions, reproducibility reports, and research collaborations related to Structured Analytic Mappings, Analytic-ICP, point-set registration, or geometric deformation modeling are welcome through this repository's GitHub issue tracker.

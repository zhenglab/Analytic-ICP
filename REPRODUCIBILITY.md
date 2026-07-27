# Reproducibility Notes

This document describes what can currently be reproduced from the public Analytic-ICP repository and separates verified repository facts from items that still require confirmation.

The associated article is:

> Wei Feng, Tengda Wei, and Haiyong Zheng, *Structured Analytic Mappings for Point Set Registration*, SIAM Journal on Imaging Sciences, 19 (2026), pp. 1655–1697.
> <https://doi.org/10.1137/25M1752080>

## Current Reproducibility Level

The repository is a research-code snapshot with a selected runnable frontend. It is not yet a complete artifact for regenerating every figure and table in the paper.

Currently available:

- the core 2D and 3D Analytic-ICP source;
- a DLL wrapper with CSV input and output;
- an INI-driven Windows console frontend;
- one configured 3D example containing 4,706 fixed points and 4,706 moving points;
- stored registered and synthetically perturbed outputs;
- additional 2D and 3D research point sets.

Not yet available:

- a clean-machine build record;
- portable dependency discovery or CMake support;
- fixed-seed synthetic deformation generation;
- automated numerical regression tests;
- scripts mapping every paper figure and table to an executable command;
- complete data provenance and redistribution records.

## Environment Encoded in the Repository

The following information is taken from `SmoothAdjustment.sln` and the three `.vcxproj` files:

| Item | Encoded value |
|---|---|
| Visual Studio | Visual Studio 14 / Visual Studio 2015 |
| Platform toolset | `v140` |
| Primary configuration | `Release | x64` |
| Eigen include path | `eigen-3.3.8` |
| Boost include path | `boost_1_60_0` |
| OpenCV library | `opencv_world3416.lib` for Release |
| OpenCV debug library | `opencv_world3416d.lib` for Debug |

All three projects contain absolute paths rooted at:

```text
E:\ScientificComputingPackage\
```

Those paths must be replaced locally before another machine can build the solution. The documentation records this limitation; the first publication-documentation update does not alter the project files.

## Build Procedure

1. Install a compatible Visual Studio C++ toolchain.
2. Provide Eigen, Boost, and OpenCV.
3. Replace the absolute include and library paths in:

   ```text
   Analytic_ICP/Analytic_ICP.vcxproj
   SmoothAdjustment/SmoothAdjustment.vcxproj
   Test/Test.vcxproj
   ```

4. Open `SmoothAdjustment.sln`.
5. Select `Release | x64`.
6. Build the `Analytic_ICP`, `SmoothAdjustment`, and `Test` projects.
7. Collect `Analytic_ICP.dll`, `SmoothAdjustment.dll`, and `Test.exe` in one runtime directory.
8. Put `experiment.ini` in that directory and preserve the relative `examples/` and `results/` paths, or edit the paths in the INI file.
9. Make `opencv_world3416.dll` available beside the executables or through `PATH` when using the encoded dynamic OpenCV configuration.

This procedure reflects the source and project metadata. It still requires validation on a clean Windows machine.

## Selected 3D Example

The checked-in configuration is `Test/experiment.ini`:

```ini
[Param]
RegistType=1
AddPerturb=1
MovingPsPath=examples/3d/moving.csv
FixedPsPath=examples/3d/fixed.csv
PerturbedOutputPath=results/perturbed3d.csv
MovedOutputPath=results/moved3d.csv
PerturbDeg2D=8
PerturbDeg3D=2
```

Input inventory:

| File | Rows | Dimension |
|---|---:|---:|
| `Test/examples/3d/fixed.csv` | 4,706 | 3 |
| `Test/examples/3d/moving.csv` | 4,706 | 3 |

Run `Test.exe` with `Test/` as the working directory. A successful invocation should create or replace:

```text
Test/results/moved3d.csv
Test/results/perturbed3d.csv
```

Each output should contain 4,706 three-dimensional points.

### Important Execution Detail

With the current `Test.cpp`:

1. `AddPerturb=1` generates `PerturbedOutputPath` from `MovingPsPath`;
2. registration then still uses the original `MovingPsPath`;
3. `PerturbedOutputPath` is not automatically fed back into registration.

Consequently, the stored perturbation is an additional generated artifact, not the source automatically registered in the same run. To register it using the existing executable, perform a second run with `MovingPsPath` changed to the generated file.

### Randomness

The 2D and 3D synthetic deformation functions seed the C random-number generator with the current time. Therefore, `perturbed3d.csv` is not expected to be byte-for-byte identical across runs.

The current 2D perturbation implementation also attempts to write coefficients to the hard-coded path:

```text
E:\coef.csv
```

This is a known portability limitation and may prevent the 2D synthetic deformation example from running on machines without a writable `E:` drive. It is documented here but not changed in this documentation-only update.

## Interpreting Results

`Test.exe` prints the scalar returned by `p2d_regist` or `p3d_regist` as:

```text
Final error: ...
```

No machine-readable experiment log or asserted numerical tolerance is currently included. The stored output files provide a qualitative reference, but the repository does not yet define them as cross-platform golden files.

Timing values in the paper should be interpreted using the hardware and comparison scope stated in the corresponding experimental subsection. The legacy i5-6200U/8 GB laptop description in the README is implementation history and should not be assigned to paper experiments that explicitly state different hardware.

## Remaining Steps for a Strong Reproducibility Release

1. validate `Release | x64` on a clean Windows machine;
2. replace personal dependency paths with a portable configuration;
3. record exact compiler and dependency versions used for the release;
4. decide and test the intended `AddPerturb` workflow;
5. replace the hard-coded 2D coefficient-output path;
6. introduce an explicit random seed;
7. add numerical smoke tests and expected tolerances;
8. map selected paper figures and tables to datasets and configurations;
9. complete dataset and third-party licensing records.

No item in this list changes the scientific claims of the paper; it concerns packaging, verification, and reproducible execution.

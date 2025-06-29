# Smooth Adjustment

This repository provides the official C++ implementation of our paper:

**"Structured Analytic Mappings for Point Set Registration"**

It consists of three main components:

- **Analytic-ICP DLL** – Implements the core nth-order analytic ICP algorithm  
- **SmoothAdjustment DLL** – Provides post-registration tools such as surface smoothing, visualization, and measurements  
- **Test Application** – An executable frontend for running registration experiments and visualizing results

---

## Overview

The `SmoothAdjustment` module is the primary framework in this project. It builds upon the `Analytic-ICP` library to support advanced operations on registered surfaces, such as:

- Analytic warping of surfaces
- Measurement and evaluation tools
- Visualization (2D/3D)

While currently relying on `Analytic-ICP`, the design is modular and can interface with other registration backends.

---

## Prerequisites

To build and run this project, the following are required:

- **Operating System**: Windows  
- **Compiler/IDE**: Visual Studio 2013  
- **Dependencies**:
  - [Eigen](https://eigen.tuxfamily.org/)
  - [Boost](https://www.boost.org/) (for KD-tree support)
  - [OpenCV](https://opencv.org/) (for image display in the Test app)

> ⚠️ Note: This project currently supports only Visual Studio 2013. Later versions may require project file adjustments.

---

## Implementation Notes

- The ICP module is based on Prof. **Andreas Geiger**’s point-to-point ICP implementation.
- We have extended the original code by:
  - Integrating **Eigen** for all matrix and vector operations
  - Replacing the KD-tree implementation with **Boost KD-tree**
- These improvements enable the handling of large-scale point clouds in both 2D and 3D.
- The `SmoothAdjustment` library uses OpenCV for interactive image-based diagnostics and rendering.
- The `Test` application serves as a lightweight frontend to configure and run experiments.

---

## Datasets

The following datasets are included in the package:

- `SHREC’19` benchmark dataset
- `3DScanRep` real-world scans
- Synthetic examples (2D & 3D)

All datasets used in the experiments are bundled with the release. No external downloads are required.

---

## Quick Start

1. Download or clone this repository.
2. Navigate to the `Release/` folder.
3. Place the required OpenCV DLLs (e.g., `opencv_world300.dll`) into the `Release/` directory.
4. Edit the `experiment.ini` configuration file. Specify:
    - Path to input point sets
    - Registration parameters
    - Output path
5. Run `Test.exe`.

> A working example is provided. After setup, simply double-click `Test.exe` to launch the default experiment.

---

## File Structure

```plaintext
SmoothAdjustment/
├── AnalyticICP/           # Core registration module
├── SmoothAdjustment/      # Analytic operations and deformation
├── Test/                  # Executable interface
├── data/                  # Sample datasets (2D & 3D)
├── experiment.ini         # Configuration file
├── Release/               # Compiled executables and runtime DLLs
└── README.md              # This file

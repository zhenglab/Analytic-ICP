# Third-Party Source and Licensing Notices

This document summarizes license notices found in the checked-in source. It is a factual inventory, not legal advice and not a project-wide license grant.

## Andreas Geiger `libicp`

The following files contain copyright and GNU General Public License notices attributed to Andreas Geiger and the Institute of Measurement and Control Systems, Karlsruhe Institute of Technology:

```text
Analytic_ICP/algorithm/ICP/icp.h
Analytic_ICP/algorithm/ICP/icp.cpp
Analytic_ICP/algorithm/ICP/icpPointToPoint.h
Analytic_ICP/algorithm/ICP/icpPointToPoint.cpp
Analytic_ICP/algorithm/ICP/icpPointToPlane.h
Analytic_ICP/algorithm/ICP/icpPointToPlane.cpp
```

Their headers state that `libicp` may be redistributed and modified under the **GNU General Public License, version 2 or any later version**.

Original project page:

<https://www.cvlibs.net/software/libicp/>

The Analytic-ICP repository reorganizes and modifies this component, including integration with Eigen and the structured analytic fitting pipeline. Those modifications do not remove the original notices or their obligations.

The repository currently does not include a complete copy of the GNU GPL text. A compliant release should restore the applicable license text and retain the original notices.

## Matthew B. Kennel KD-tree

The following files contain a copyright notice attributed to Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004):

```text
Analytic_ICP/algorithm/ICP/kdtree.h
Analytic_ICP/algorithm/ICP/kdtree.cpp
```

Their headers state:

```text
Licensed under the Academic Free License version 1.1 found in file LICENSE
with additional provisions in that same file.
```

The referenced `LICENSE` file and its additional provisions are not present in the current repository. The original license material must be recovered and reviewed before the licensing status of this component can be documented completely.

## External Dependencies

The build also depends on external installations of:

- Eigen;
- Boost;
- OpenCV;
- Microsoft Visual C++ runtime components.

These dependencies are not vendored in the repository. Users and redistributors must comply with the terms of the versions they install and distribute.

## Original Analytic-ICP Code

No project-wide license has yet been declared for the original Structured Analytic Mapping, wrapper, and frontend code. Consequently, the repository should not currently be described as freely redistributable under a particular open-source license.

Before declaring a project-wide license:

1. recover all missing third-party license texts and additional provisions;
2. identify the exact upstream source/version or commit for each adapted component;
3. determine whether the intended project-wide license is compatible with the GPL-covered component;
4. alternatively, replace or isolate components whose terms conflict with the intended distribution model;
5. complete the dataset provenance review in `DATASETS.md`;
6. add a root `LICENSE` file and preserve all third-party notices.

Until those steps are complete, public availability of the source should not be interpreted as permission to redistribute the entire combined work under an arbitrary license.

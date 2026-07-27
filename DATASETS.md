# Included Point-Set Inventory

This file records what is physically present in the repository. It does **not** assert that redistribution rights have been verified for every point set.

The source URL, original publication, creator, license, and any transformations applied to each third-party dataset still need to be confirmed before a formal archival release. Until that review is complete, do not assume that files under `Data/` may be redistributed independently of this research snapshot.

## `Data/`

| File | Points | Dimension | Delimiter |
|---|---:|---:|---|
| `8.csv` | 248 | 2 | comma |
| `B.csv` | 248 | 2 | comma |
| `body1.csv` | 167,741 | 3 | comma |
| `body4.csv` | 200,735 | 3 | comma |
| `body_1.csv` | 4,706 | 3 | comma |
| `body_26.csv` | 10,050 | 3 | comma |
| `body_44.csv` | 6,890 | 3 | comma |
| `bunny.csv` | 14,290 | 3 | comma |
| `cowhead.csv` | 2,036 | 3 | comma |
| `curve.csv` | 78 | 2 | whitespace |
| `face.csv` | 392 | 3 | comma |
| `fish with outliers.csv` | 127 | 2 | whitespace |
| `fish.csv` | 91 | 2 | comma |
| `trash can.csv` | 359 | 2 | comma |

The counts above were derived from the checked-in files and are not dataset-origin claims.

## Configured 3D Example

| File | Points | Dimension | Note |
|---|---:|---:|---|
| `Test/examples/3d/fixed.csv` | 4,706 | 3 | Fixed/target input in `experiment.ini` |
| `Test/examples/3d/moving.csv` | 4,706 | 3 | Moving/source input in `experiment.ini`; byte-identical to `Data/body_1.csv` in the current snapshot |
| `Test/results/moved3d.csv` | 4,706 | 3 | Stored registration output |
| `Test/results/perturbed3d.csv` | 4,706 | 3 | Stored random analytic-deformation output |

The generated perturbation is time-seeded and can differ across runs.

## Information Needed to Complete Provenance

For each non-original dataset, record:

1. canonical dataset or project name;
2. original download URL;
3. creator or maintaining institution;
4. paper or technical report to cite;
5. dataset license or redistribution terms;
6. exact original filename or object identifier;
7. preprocessing, subsampling, normalization, deformation, or format conversion performed for this repository;
8. the paper figure or table in which it was used.

Files created entirely by the authors should instead be marked as original synthetic data and accompanied by the generation procedure and intended license.

## Recommended Archival Policy

After provenance is established:

- retain redistributable data with its license and citation;
- replace restricted data with download instructions and checksums;
- clearly separate original inputs, generated deformations, and algorithm outputs;
- avoid using a `.csv` extension for whitespace-delimited files unless the format is intentionally documented;
- add machine-readable manifests containing point count, dimension, checksum, source, and license.

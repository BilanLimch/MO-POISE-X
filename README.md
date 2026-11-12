# MO-POISE

**Assay-aware multi-objective evolutionary feature selection for compact multi-omics panels**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.13%20%7C%203.14-blue)
![pymoo](https://img.shields.io/badge/pymoo-0.6.2-orange)

MO-POISE (Multi-Objective Plug-and-play Optimizer for Imbalance-aware Subtyping of
Endometrial cancer) searches for small multi-omics panels that reconstruct the four
TCGA molecular subtypes of endometrial carcinoma (*POLE*, MSI, copy-number-low,
copy-number-high). Every panel is judged by three objectives at once:

| Objective | Minimized quantity |
|---|---|
| f1, error | 1 - macro-averaged F1 of logistic regression in inner 3-fold cross-validation |
| f2, assay cost | number of omics assays the panel needs (or any cost vector, see `--assay-cost`) |
| f3, size | number of selected features / N |

Any multi-objective evolutionary algorithm (MOEA) can be plugged in. The pipeline
with MOEA X is called MO-POISE-X, for example MO-POISE-NSGA-II. This repository
contains the code, the per-run results, and the scripts that reproduce every table and
analysis of the accompanying article.

---

## Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Quick start](#quick-start)
4. [Data](#data)
5. [Reproducing the benchmark](#reproducing-the-benchmark)
6. [Reproducing the additional analyses](#reproducing-the-additional-analyses)
7. [Command-line options](#command-line-options)
8. [Algorithms](#algorithms)
9. [Result files](#result-files)
10. [Repository structure](#repository-structure)
11. [Citation](#citation)
12. [License](#license)

---

## Overview

```
 TCGA-UCEC, 369 patients, 6 data blocks
          |
          v
 5 outer folds (stratified, one fixed partition)            <- leakage-controlled nested CV
          |
          +-- outer training fold
          |      ANOVA prefilter: top 500 features per omics block + 5 clinical (N = 2,505)
          |      MOEA search over binary panels, objectives f1, f2, f3
          |        shared operators: block-sparse initialization,
          |                          uniform crossover, block-aware mutation
          |      training-side choice: one-SE rule, per assay budget, per size cap
          |      filter reference fronts: ANOVA, mRMR, L1 over all 63 block subsets
          |
          +-- outer test fold
                 every panel of the final front scored once
                 test-side hypervolume, generalization gap, paired tests vs filter fronts
```

Main features

* **Assay-aware objective.** The number of omics assays is optimized directly, next to
  accuracy and panel size.
* **Plug-and-play.** Thirteen MOEAs share one problem definition, one evaluation, and
  one budget; new algorithms only need a constructor in `mopoise/search.py` or a module
  registered in `mopoise/moeas/`.
* **Block-sparse initialization.** Each initial panel samples its omics assays first and
  its features second, so the initial population covers every number of assays.
* **No selection bias.** The prefilter, the scaling, the search, and the choice of the
  reported panel use the outer training fold only; outer-test scores never choose a panel.
* **Filter reference fronts.** Every run is compared with univariate (ANOVA),
  redundancy-aware (mRMR), and L1-penalized panels built on the same training fold.

## Installation

Download `MOPOISE.zip` from https://github.com/BilanLimch/MO-POISE-X and unpack it.

```bash
unzip MOPOISE.zip && cd MOPOISE
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export PYTHONPATH=$PWD        # the scripts in scripts/ import the package from the root
```

The results were produced with Python 3.14, NumPy 2.5, SciPy 1.18, pandas 3.0,
scikit-learn 1.9, pymoo 0.6.2, and numba 0.68 (`requirements.txt` gives the accepted
ranges).

## Quick start

A run on synthetic data with the shape of the cohort needs no download and finishes in
seconds:

```bash
UCEC_RESULTS_DIR=/tmp/mopoise_smoke python -m mopoise.nested_cv \
    --algorithm NSGA2 --run 0 --synthetic --max-evals 600
```

The run prints the size of the final front, the inner and test-side hypervolume, and the
panel chosen by the one-standard-error rule, and writes
`/tmp/mopoise_smoke/NSGA2/fold0_seed2026.json`.

## Data

All data are public. Download the following TCGA-UCEC files from the GDC hub of UCSC
Xena (https://xenabrowser.net) into one directory, referred to as `$UCEC_RAW_DIR`:

| File | Content |
|---|---|
| `TCGA-UCEC.star_fpkm.tsv.gz` | gene expression, log2(FPKM + 1) |
| `TCGA-UCEC.mirna.tsv.gz` | miRNA expression, log2(RPM + 1) |
| `TCGA-UCEC.methylation450.tsv.gz` | DNA methylation, 450K beta values |
| `TCGA-UCEC.gene-level_ascat3.tsv.gz` | gene-level copy number |
| `TCGA-UCEC.somaticmutation_wxs.tsv.gz` | somatic mutations |
| `TCGA-UCEC.clinical.tsv.gz` | clinical data |
| `HM450.hg38.manifest.gencode.v36.probeMap` | probe annotation (additional analyses only) |

Add `TCGA-UCEC.molecular_subtype.tsv`, a tab-separated file with a header line, the
patient barcode in the first column, and the molecular subtype of the PanCancer Atlas
endometrial study in cBioPortal (https://www.cbioportal.org) in the second column.

Then build the cohort and the label-free data cache:

```bash
export UCEC_RAW_DIR=/path/to/raw
python scripts/build_cohort.py            # writes $UCEC_RAW_DIR/final_cohort.tsv
python -m mopoise.data --build            # writes the cache to ./data
```

## Reproducing the benchmark

```bash
JOBS=8 OMP_NUM_THREADS=2 ./run_all.sh            # filter fronts + 13 algorithms x 20 runs
JOBS=4 OMP_NUM_THREADS=4 ./run_all.sh NSGA2      # a subset of algorithms
python -m mopoise.aggregate                      # tables and statistics -> results/summary/
```

Run i uses outer fold i // 4 and seed 2026 + 98 i, so runs with the same index are
paired across algorithms. A unit whose result file already exists is skipped; pass
`--force` to `mopoise.nested_cv` to recompute it.

**Settings of the main runs.** Population 50 (45 for TELSO and IM-C-MOEA/D), 5,000
evaluations, 500 candidates per omics block, initial per-block density log-uniform in
[0.005, 0.2], block-sparse initialization, uniform crossover with probability 0.9,
bit-flip mutation with rate 1/N followed by the removal of one active block with
probability 0.1, equal assay costs, logistic regression with the scikit-learn defaults,
three inner folds. Algorithm-specific parameters are set in `mopoise/search.py` and
`mopoise/moeas/`.

## Reproducing the additional analyses

Unless stated otherwise, the analyses below use SPEA2, SMSEMOA, and OMNI, the three
algorithms with the highest mean test-side hypervolume, and runs 0 to 19. Each command
is shown for one run or fold; loop over `--run 0..19` or `--fold 0..4`.

| Analysis | Command | Output |
|---|---|---|
| Initial density up to 1.0 | `python -m mopoise.nested_cv --algorithm SPEA2 --run 0 --density-max 1.0 --tag _dens1` | `results/<ALG>_dens1/` |
| 20,000 evaluations | `python -m mopoise.nested_cv --algorithm SPEA2 --run 0 --max-evals 20000 --tag _20k` | `results/<ALG>_20k/` |
| Flat and uniform initialization | `python -m mopoise.nested_cv --algorithm SPEA2 --run 0 --init flat --tag _flat` (or `--init uniform --tag _unif`) | `results/<ALG>_flat/`, `results/<ALG>_unif/` |
| Flat initialization, other seven competitive algorithms | `python -m mopoise.nested_cv --algorithm NSGA2 --run 0 --init flat --out results_rev --tag _flat` | `results_rev/<ALG>_flat/` |
| 250 or 1,000 candidates per block | `python -m mopoise.nested_cv --algorithm SPEA2 --run 0 --k-per-block 250 --out results_rev --tag _k250` | `results_rev/<ALG>_k250/`, `results_rev/<ALG>_k1000/` |
| Filter fronts for these pools | `python -m mopoise.baselines --fold 0 --method anova --k-per-block 250 --out results_rev/k250` (methods `anova`, `mrmr`, `l1`) | `results_rev/k250/`, `results_rev/k1000/` |
| Assay-specific costs | `python -m mopoise.nested_cv --algorithm SPEA2 --run 0 --assay-cost mut=1.9,rna=1.1,meth=0.9,mirna=0.8,cnv=0.3 --out results_rev --tag _wcost` | `results_rev/<ALG>_wcost/` |
| Readouts of the clinical surrogate in the candidate pool | `python scripts/build_aggregates.py --out data/ucec_aggregates.npz`, then `python -m mopoise.nested_cv --algorithm SPEA2 --run 0 --aggregates data/ucec_aggregates.npz --out results_rev --tag _agg` and `python -m mopoise.baselines --fold 0 --method anova --aggregates data/ucec_aggregates.npz --out results_rev/agg` | `results_rev/<ALG>_agg/`, `results_rev/agg/` |

The following scripts read the stored panels and run no search:

| Analysis | Command | Output |
|---|---|---|
| MLH1 probe list | `python scripts/extract_mlh1_probes.py $UCEC_RAW_DIR/HM450.hg38.manifest.gencode.v36.probeMap` | `data/mlh1_probes.tsv` |
| Random forest and SVM refits, per-class precision and recall, Cohen's kappa, molecular surrogate of the ProMisE classification | `python scripts/revision_posthoc.py --fold 0` | `results_rev/posthoc/` |
| Tables of the post hoc analyses | `python scripts/revision_tables.py` | printed |
| Fold-wise results and selection stability | `python scripts/revision_stability.py` | `results_rev/stability/` |
| Re-scoring of the final fronts on repeated inner partitions | `python scripts/revision_rescore.py --fold 0`, then `python scripts/revision_rescore_summary.py` | `results_rev/rescore/` |
| Summaries | `python scripts/revision_sensitivity.py`, `python scripts/revision_ablation.py`, `python scripts/revision_aggregates.py` | `results_rev/sensitivity/`, `results_rev/ablation/`, printed |

`revision_posthoc.py` and `build_aggregates.py` read the somatic-mutation file from
`$UCEC_RAW_DIR`.

## Command-line options

`python -m mopoise.nested_cv` runs one (algorithm, fold, seed) unit.

| Option | Default | Meaning |
|---|---|---|
| `--algorithm` | required | one of the names in [Algorithms](#algorithms) |
| `--run` | | run index 0..19; sets the outer fold and the seed |
| `--fold`, `--seed` | | set both explicitly instead of `--run` |
| `--max-evals` | 5000 | evaluation budget |
| `--pop-size` | 50 | population size |
| `--k-per-block` | 500 | candidates kept per omics block by the prefilter |
| `--classifier` | LR | LR, SVM, RF, or KNN |
| `--init` | block | initialization: `block` (block-sparse), `flat` (sparse without the block step), `uniform` (probability 0.5) |
| `--density-max` | 0.2 | upper end of the initial per-block density |
| `--assay-cost` | equal | block=cost pairs, e.g. `mut=1.9,rna=1.1` |
| `--aggregates` | | npz of label-free aggregate features, always kept in the candidate pool |
| `--tag` | | suffix of the result directory |
| `--out` | `results` | result root (also `$UCEC_RESULTS_DIR`) |
| `--force` | off | recompute an existing unit |
| `--synthetic` | off | use synthetic data |

`python -m mopoise.baselines --fold K --method {anova,mrmr,l1}` builds one filter front
and accepts `--k-per-block`, `--aggregates`, and `--out` with the same meaning.

## Algorithms

| Name | Algorithm | Implementation | Operators |
|---|---|---|---|
| `SPEA2` | SPEA2 | pymoo | shared |
| `NSGA2` | NSGA-II | pymoo | shared |
| `SMSEMOA` | SMS-EMOA | pymoo | shared |
| `OMNI` | Omni-optimizer | pymoo | shared |
| `NSGA3` | NSGA-III | pymoo | shared |
| `UNSGA3` | U-NSGA-III | pymoo | shared |
| `RVEA` | RVEA | pymoo | shared |
| `AGEMOEA` | AGE-MOEA | pymoo | shared |
| `AGEMOEA2` | AGE-MOEA-II | pymoo | shared |
| `NRVMOEA` | NRV-MOEA | `mopoise/moeas/nrvmoea.py` | shared |
| `MOBCA` | MOBCA | `mopoise/moeas/mobca.py` | shared initialization, own real-coded moves |
| `IMCMOEAD` | IM-C-MOEA/D | `mopoise/moeas/imcmoead.py` | shared initialization, own real-coded moves |
| `TELSO` | TELSO | `mopoise/moeas/telso.py` | own binary operators |

## Result files

Each run is stored as `results/<ALG>/fold<k>_seed<s>.json`:

| Field | Content |
|---|---|
| `protocol` | all settings of the run |
| `front` | every panel of the final front, the set of panels non-dominated in the inner objectives (for runs of RVEA and the Omni-optimizer stored with dominated members of the final population, `mopoise.aggregate.load_run` removes these members and recomputes the test-side hypervolume, as all analysis scripts do): features, blocks, size, assay cost, inner error and its standard error, and outer-test macro-F1, macro-AUC, accuracy, and per-class F1 (runs made with the current version also store Cohen's kappa) |
| `selected` | index of the one-standard-error panel and of the panel per assay budget |
| `hypervolume` | inner and test-side hypervolume |
| `hv_curve` | inner hypervolume against the number of evaluations |
| `search_seconds`, `total_seconds`, `distinct_panels_scored` | run time and number of distinct panels evaluated |

Filter fronts are stored per fold in `results/baselines/`, `results/baselines_mrmr/`,
and `results/baselines_l1/`; `results/summary/` holds the tables and statistics of
`mopoise.aggregate`. The directories of the additional analyses are listed in the tables
above. The result files contain no patient identifiers.

## Repository structure

```
mopoise/
  data.py           cohort loader, label-free data cache, aggregate features
  prefilter.py      in-fold ANOVA top-k per block
  objectives.py     panel evaluator: inner CV, three objectives, cache
  search.py         problem, shared operators, pymoo algorithms, hypervolume
  moeas/            NRV-MOEA, MOBCA, IM-C-MOEA/D, TELSO
  selection.py      non-dominated filter, one-SE rule, per-budget rule
  nested_cv.py      one (algorithm, fold, seed) unit
  baselines.py      ANOVA, mRMR, and L1 filter fronts
  aggregate.py      tables and statistics
  synthetic.py      synthetic data for smoke runs
scripts/            data preparation and the additional analyses (see above)
results/            per-run results of the benchmark
results_rev/        per-run results of the additional analyses
run_all.sh          all filter fronts and all algorithm runs
```

## Citation

If you use this code, please cite the accompanying article:

> MO-POISE-X: Multi-Objective Evolutionary Search with Assay-Aware Block-Sparse
> Initialization for Compact Multi-Omics Panels of Endometrial Cancer Molecular Subtypes.

The full reference will be added here upon publication.

## Acknowledgements

The results are based on data generated by the TCGA Research Network
(https://www.cancer.gov/tcga), obtained from the GDC hub of UCSC Xena and from
cBioPortal.

## License

MIT, see [LICENSE](LICENSE).

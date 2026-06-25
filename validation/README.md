# Validation Report

**Date:** 22-06-2026
**Tool:** NDV Genotyper v1.4.1
**Author:** Antoine Cazin

---

## Overview

This report documents the accuracy evaluation of the NDV F Gene Genotyper against a labeled holdout dataset. The goal was to produce a leakage-free estimate of genotyping performance on sequences the classifier had never seen.

---

## Holdout Dataset

### Construction

A labeled dataset of NDV F gene sequences with known genotypes was assembled from publicly available sources (NCBI - Genbank). Prior to evaluation, the holdout was audited for **data leakage**: any sequence byte-identical to a sequence already present in the reference database was removed.

This audit revealed that a first version of the test set contained **84.6% leaked sequences** - sequences that were directly present in the reference DB, making any accuracy measured on them an overestimate. All matching sequences were removed from the test set.

Additionally, internal duplicates within the test set itself were removed.

### Final composition

| Total sequences | 151 |
| Leakage (identical to any reference) | 0 |
| Internal duplicates | 0 |
| Classes represented | Class I + Class II |

---

## Method

**Similarity metric:** Pairwise (identity over aligned positions)

Accuracy was measured at three levels of the NDV nested nomenclature:

| Level | Definition |
|---|---|
| **Exact** | Predicted sub-genotype matches label exactly (e.g. `VII.1.1 == VII.1.1`) |
| **Top genotype** | Correct at the top level (e.g. any `VII.*` counted as `VII`) |
| **Wrong** | fully wrong |

---

## Results

| Metric | Value |
|---|---|
| Exact (sub-genotype) accuracy | **97.4%** |
| Top genotype accuracy | 100.0% |
| Macro accuracy | **92.9%** |
| Sequences evaluated | 151 |

(Macro: Mean of per-genotype recall, rare genotypes weighted equally to common ones)

### Misclassifications

The 4 misclassifications were **"Genotype"-level** errors: the top-level genotype was correct but the sub-lineage differed. No "Wrong"-class errors were observed.

| Header | Real Genotype | Predicted Genotype | 
| --- | --- | --- |
| KF767468.1 | V | V.2 |
| KF767469.1 | V | V.2 |
| OR253530.1 | I.1.2.1 | I.2 | 
| OR253549.1 | I.2 | I.1.2.1 |

---

## Interpretation

A 97.4% exact sub-genotype accuracy on a **leakage-free** holdout is a strong result for a nearest-neighbour classifier. The main source of remaining error is sub-lineage ambiguity at the boundaries of closely related clades, which is a known challenge in NDV nomenclature and not a limitation of the method per se.

The pairwise similarity metric outperformed Hamming distance on this dataset, likely because it handles gap-rich regions more robustly after MAFFT alignment. 

---

## Limitations

- The reference database composition directly determines performance; sequences from underrepresented genotypes will have fewer neighbours and may be mis-assigned.
- The holdout covers sequences available as of June 2026; performance on future divergent strains may differ.
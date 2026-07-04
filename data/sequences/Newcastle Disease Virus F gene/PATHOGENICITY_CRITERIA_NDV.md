# NDV - Pathogenicity Criteria

Based on: **Dimitrov et al. (2019)** and **Wang et al. (2017)**

---

## Overview

Pathogenicity is determined by analyzing the **F protein cleavage site** - a sequence of 6 amino acids (positions 113-118) in the fusion protein. This tool searches for this motif within a ±15 nucleotide window around position 333 bp to account for insertions/deletions.

---

## Classification

### Likely Virulent
A known virulent motif (VFcs) is found in the cleavage region.

All virulent motifs share:
- Multiple basic residues (R/K) at positions P5-P1
- **F (Phenylalanine) at position P1'**

| Category | Motif | Genotypes |
|---|---|---|
| VFcs-1 | RRQKRF | II, V-VIII, XII-XIV, XVI-XVIII |
| VFcs-2 | RRQRRF | I, III-IV, VII, IX, XIII, XVII |
| VFcs-3 | KRQKRF | V-VII |
| VFcs-4 | GRQKRF | VI |
| VFcs-5 | RRKKRF | VI-VII |
| VFcs-6 | RRRKRF | VI-VII, XIII-XIV, XVII-XVIII |
| VFcs-7 | KRKKRF | VI |
| VFcs-8 | RRRRRF | XI |

### Likely Low-virulence
A known avirulent motif (AFcs) is found in the cleavage region.

All avirulent motifs share:
- Mono or multi-basic residues at positions P5-P1
- **L (Leucine) at position P1'** (except AFcs-5)

| Category | Motif | Genotypes |
|---|---|---|
| AFcs-1 | GRQGRL | Class I, I-II |
| AFcs-2 | GKQGRL | I-II, X |
| AFcs-3 | RRQGRL | I |
| AFcs-4 | ERQERL | Class I |
| AFcs-5 | RRQGRF | I (avirulent despite F) |
| AFcs-6 | ERQGRL | Class I, I |
| AFcs-7 | RKQGRL | I |
| AFcs-8 | EKQGRL | I, X |
| AFcs-9 | EQQERL | Class I |
| AFcs-10 | RRQRRL | I |

### Undetermined
No known motif was found in the expected region. The sequence may be incomplete, contain too many indels, or represent an unusual strain requiring further biological characterization.

---

## Important Notes

- This tool uses **motif-based detection** - results depend entirely on whether the cleavage site matches a known pattern from Wang et al. (2017)
- Results are reported as **"Likely"** Virulent or Low-virulence - molecular analysis alone is not sufficient for definitive pathogenicity classification
- The **ICPI (Intracerebral Pathogenicity Index)** remains the gold standard for official virulence determination (OIE/WOAH)
- Sequences shorter than ~320 bp may not contain the cleavage site region and will return Undetermined

---

## References

**Wang Y. et al. (2017)** - Characterization of Newcastle disease viruses isolated from chicken and goose flocks between 2013 and 2015 in some regions of China. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0183923

**Dimitrov K.M. et al. (2019)** - Updated unified phylogenetic classification system and revised nomenclature for Newcastle disease virus. *Infect. Genet. Evol.* https://doi.org/10.1016/j.meegid.2019.103917
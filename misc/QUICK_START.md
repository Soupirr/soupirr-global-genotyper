# NDV Sequence Analyzer - Quick Start Guide

---

## ⚠️ IMPORTANT - Adding New Sequences to the Reference Database

Before adding any new sequences to the reference database, **every sequence header must strictly follow this format:**

```
>NUMBER_GENOTYPE_INFORMATION_YEAR
```

**Example:**
```
>142_VII.1.1_MF100730.1_chicken_China_GD_390_2016
```

| Field | Description | Example |
|---|---|---|
| `NUMBER` | Unique sequence ID | `142` |
| `GENOTYPE` | NDV genotype | `VII.1.1` |
| `INFORMATION` | Accession, host, location | `MF100730.1_chicken_China_GD_` |
| `YEAR` | year | `2016` |

**Fields must be separated by underscores `_`.** The genotype must always be the second field. Sequences that do not follow this format will not be correctly identified and may cause errors in futur analysis. Year should also be in the last position.

A full list of the Gentoypes name can be found inside the "genotypes.txt" file.
Also note that adding too much sequences may make analysis slower over time.

To add new sequences:
1. Make sure your sequences follow the above format
2. Place your `.fas` file inside the `data/sequences/` folder
3. Restart the app

---

## What This App Does

The NDV Sequence Analyzer is a bioinformatics tool designed to identify the genotype and predict the pathogenicity of Newcastle Disease Virus (NDV) F gene sequences. It compares your input sequence against a curated reference database of 2944 sequences spanning 53 sub-genotypes.

---

## How to Use

### 1. Sequence Analysis Tab

This is the main tab of the app. It allows you to identify the genotype and pathogenicity of one or multiple NDV F gene sequences.

**Step 1 - Choose your similarity method**

In the left configuration panel, select one of two methods:

- **Hamming (fast)** - compares sequences position by position. Very fast but requires sequences of similar length. Results below 95% similarity may not be reliable due to insertions/deletions.
- **Pairwise Alignment (accurate)** - performs a global pairwise alignment before comparing. Handles insertions and deletions correctly. Recommended for sequences that may differ in length from the reference. Takes approximately 30 seconds per sequence.

**Step 2 - Choose the number of top matches**

Use the slider to select how many top matching genotypes you want to see in the results (1 to 10, default 5).

**Step 3 - Input your sequence**

You have two options:
- **Paste FASTA** - paste one or multiple sequences directly in FASTA format
- **Upload File** - upload a `.fasta`, `.fas`, `.fa`, or `.txt` file

Your input must be in FASTA format:
```
>sequence_name
ATGGGCTCCAGACCT...
```

Multiple sequences are supported - each sequence will be analyzed individually and results displayed in separate tabs.

**Step 4 - Analyze**

Click **Analyze Sequences** to start the analysis. A progress bar will show the status for each sequence.

---

### 2. Results

Results are displayed in individual tabs for each sequence analyzed.

**Sequence Information** - header and length of your input sequence.

**Genotype Identification** - shows:
- Detected class (Class I or Class II)
- Best matching genotype
- Average similarity percentage
- Number of reference sequences for that genotype
- Best individual match score and sequence name

A warning will appear if similarity is below 95% when using the Hamming method.

**Pathogenicity Analysis** - analyzes the F protein cleavage site region (positions ~333–357 bp) using a ±15 nucleotide window to account for insertions/deletions. Motifs are classified according to Dimitrov et al. 2019:

| Result | Meaning |
|---|---|
| **Virulent** | A known virulent cleavage motif (VFcs-1 to VFcs-8) was found |
| **Low-virulence** | A known avirulent cleavage motif (AFcs-1 to AFcs-10) was found |
| **Undetermined** | No known motif was found - sequence may be incomplete or unusual |

> ⚠️ Pathogenicity prediction is based solely on the F protein cleavage site motif. This is a strong indicator but not the only criterion for virulence. Results should always be interpreted in the context of additional biological data.

**Export** - at the bottom of the page, a summary table shows all analyzed sequences. You can download it as a CSV file.

---

### 3. Phylogenetic Tree Tab

This tab automatically builds a small phylogenetic tree showing your query sequence in context of its (by default) 20 closest neighbours from the reference database.

The tree is built using:
- **MAFFT** for multiple sequence alignment
- **FastTree** for tree construction (GTR model, SH-like bootstrap support)

Your query sequence is highlighted in **green**. Bootstrap values ≥ 50% are shown on internal nodes.

The tree is only generated after running an analysis in the Sequence Analysis tab.

---

### 4. Map Tab

Shows the geographic distribution of NDV genotypes from the reference database. You can filter by genotype to explore where specific strains have been reported.

---

## Input Requirements

- Format: FASTA (`.fasta`, `.fas`, `.fa`, `.txt`)
- Sequence type: NDV F gene nucleotide sequence
- Recommended length: ~1662 bp (full F gene)
- Shorter or longer sequences are accepted but may give less reliable results
- Sequences with insertions/deletions are handled better by the Pairwise method

---

## References

- Dimitrov et al. (2019) - Updated unified phylogenetic classification system and revised
nomenclature for Newcastle disease virus, https://doi.org/10.1016/j.meegid.2019.103917

- Wang et al. (2017) - Comprehensive analysis of amino acid sequence diversity at the F protein cleavage site of Newcastle disease virus in fusogenic activity, https://doi.org/10.1371/journal.pone.0183923.

- FastTree: Price et al. (2009/2010), PLOS ONE

- MAFFT: Katoh et al., multiple versions

---

## Newcastle Disease Virus Dataset Citation

Dimitrov, K.M., Abolnik, C., Afonso, C.L., Albina, E., Bahl, J., Berg, M., Briand, F.X., Brown, I.H., Choi, K.S., Chvala, I., Diel, D.G., Durr, P.A., Ferreira, H.L., Fusaro, A., Gil, P., Goujgoulova, G.V., Grund, C., Hicks, J.T., Joannis, T.M., Kim Torchetti, M., Kolosov, S., Lambrecht, B., Lewis, N.S., Liu, H., Liu, H., McCullough, S., Miller, P.J., Monne, I., Muller, C.P., Munir, M., Reischak, D., Sabra, M., Samal, S.K., Servan de Almeida, R., Shittu, I., Snoeck, C.J., Suarez, D.L., Van Borm, S., Wang, Z., Wong, F.Y.K., 2019. Updated unified phylogenetic classification system and revised nomenclature for Newcastle disease virus. *Infect. Genet. Evol.*, 103917.
https://doi.org/10.1016/j.meegid.2019.103917

## Datasets Used in This Tool

### Class I Dataset
- **Filename:** NDV_F_class_I_619_May_09_2022.fas
- **Contains:** 619 sequences from Class I genotypes
- **Source:** NDV Consortium GitHub
- **URL:** https://github.com/NDVconsortium/NDV_Sequence_Datasets

### Class II Dataset
- **Filename:** NDV_F_class_II_2157_May_09_2022.fas
- **Contains:** 2157 sequences from Class II genotypes (I-XXI)
- **Source:** NDV Consortium GitHub
- **URL:** https://github.com/NDVconsortium/NDV_Sequence_Datasets


---

**Last Updated:** June 17, 2026

**Source:** https://github.com/Soupirr/NDV-genotyper

*NDV Sequence Analyzer - developed as part of a bioinformatics internship project*

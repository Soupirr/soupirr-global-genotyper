# Sequence Analyzer Toolbox - Quick Start Guide

---

## Getting Started — Selecting an Entry

On the left sidebar, select an **entry** from the dropdown. Each entry is a reference dataset for a specific pathogen. If no entry exists yet, you need to add one first.

### Adding a New Entry

Use the **"Add new references datasets"** expander in the sidebar:

1. Enter a name for the entry (e.g. `Avian Influenza H5N1`)
2. Upload one or more reference FASTA files, directly as exported from **NCBI Virus** — the app automatically converts NCBI headers to its internal format, drops duplicate sequences and sequences missing a genotype, and shows a migration report (kept / dropped counts) after import
3. Optionally configure pathogenicity analysis:
   - **bp before the virulence motif area** — 0-indexed nucleotide position where the cleavage/virulence site starts. The app scans a fixed ±29-codon window around this position to tolerate indels.
   - **Upload a motif file** — CSV with columns `motif,label,type` (e.g. `RRQKRF,VFcs-1,virulent`). The `type` column can contain any category name you define — there is no limit on the number of categories.
4. Click **Add to the registry**

The entry folder will be created automatically under `data/sequences/`.

### Optional Files Per Entry

You can enrich an entry by placing additional files in its folder:

| File | Effect |
|---|---|
| `{entry_name}_motifs.csv` | Enables pathogenicity analysis |
| `{entry_name}_pathogenicity.csv` | Pre-computed pathogenicity stats (auto-generated via the Stats tab) |
| Any `.md` file | Adds a documentation tab in the Help — Statistics section |

---

## Sequence Header Format

Every reference sequence header must strictly follow this format, fields separated by pipes `|`:

```
>VIRUS|ACCESSION|GENOTYPE|HOST|COUNTRY|REGION|YEAR
```

**Example:**
```
>NDV|MH169357.1|VII.1.1|Chicken|Iran|?|2020
```

| Field | Description | Example |
|---|---|---|
| `VIRUS` | Virus/species name | `NDV` |
| `ACCESSION` | GenBank accession (version suffix stripped) | `MH169357.1` |
| `GENOTYPE` | Genotype, clade, or serotype identifier | `VII.1.1` |
| `HOST` | Host species | `Chicken` |
| `COUNTRY` | Country of isolation | `Iran` |
| `REGION` | State/province/region, or `?` if unknown or same as country | `?` |
| `YEAR` | Isolation year | `2020` |

If you upload sequences straight from an NCBI Virus export, you don't need to format headers by hand — the app's built-in migration step (see above) converts NCBI's own export format into this one automatically. If you're building a reference FASTA by hand, the genotype must always be the 3rd field and the year the last field; sequences that don't follow this format are skipped rather than mis-parsed.

---

## What This App Does

The Sequence Analyzer Toolbox is a bioinformatics tool for identifying genotypes and predicting pathogenicity from nucleotide sequences. It compares input sequences against a curated reference database and, when configured, analyzes a defined cleavage or virulence site region.

The tool supports any pathogen for which you can provide:
- A reference FASTA database with correctly formatted headers
- Optionally, a cleavage site position and a motif classification file

---

## How to Use

### 1. Sequence Analysis Tab

**Step 1 - Choose your similarity method**

- **Hamming (fast)** — compares sequences position by position. Very fast but requires sequences of similar length. Results below 95% similarity may not be reliable when insertions/deletions are present.
- **Pairwise Alignment (accurate)** — performs a global alignment before comparing. Handles insertions and deletions correctly. Takes approximately 30 seconds per sequence.

**Step 2 - Choose the number of top matches**

Use the slider to select how many top matching genotypes to display (1 to 10, default 5).

**Step 3 - Input your sequence**

- **Paste FASTA** — paste one or multiple sequences in FASTA format
- **Upload File** — upload a `.fasta`, `.fas`, `.fa`, or `.txt` file

Multiple sequences are supported — each will be analyzed individually and results displayed in separate tabs.

**Step 4 - Analyze**

Click **Analyze Sequences**. A progress bar will show the status for each sequence.

---

### 2. Results

**Sequence Information** — header and length of your input sequence.

**Genotype Identification** — shows:
- Best matching genotype and average similarity score
- Number of reference sequences for that genotype
- Best individual match score

A warning appears if similarity is below 95% when using the Hamming method.

**Pathogenicity Analysis** — analyzes the configured cleavage/virulence site region using a ±30 nucleotide window to account for insertions/deletions. Results are based on the motif file associated with the entry. The analysis is run across three reading frames (Main, +1, -1).

| Result | Meaning |
|---|---|
| **Configured type** (e.g. virulent) | A motif matching this category was found in the cleavage region |
| **Undetermined** | No known motif was found — sequence may be incomplete or unusual |

> ⚠️ Pathogenicity prediction is based on cleavage site motif matching only. Results should always be interpreted alongside additional biological data.

**Export** — a summary table at the bottom shows all analyzed sequences and can be downloaded as CSV.

---

### 3. Phylogenetic Tree Tab

Builds a phylogenetic tree placing your query sequence among its closest neighbours from the reference database. Configure before building:

- **Tree type** — *Per-query* (one tree per analyzed sequence) or *Combined* (all queries in a single tree)
- **Tree Mode** — *Cladogram* (uniform branch lengths) or *Phylogram* (real evolutionary distances)
- **Tree Method** — *FastTree* (fast, approximate ML) or *IQ-TREE2* (full ML with ModelFinder + ultrafast bootstrap — slower, publication-quality)

MAFFT is used for alignment in both cases. Your query sequence is highlighted in colour. Bootstrap support values are shown on internal nodes. Each tree can be downloaded in Newick format. The tree is only generated after running an analysis in the Sequence Analysis tab.

---

### 4. Help — Statistics Tab

**Statistics** — database overview: sequences per genotype, temporal distribution, host distribution, geographic distribution, and coverage health per genotype.

**Map** — geographic distribution of reference sequences, nested as a sub-tab here (shares the same "Load Statistics" click, so it loads at the same time). Filter by genotype to explore where specific strains have been reported. Location is parsed from the sequence header.

**Genotypes Pathogenicity** — shows pathogenicity statistics across all reference sequences. If no pre-computed data exists yet, click **Generate Pathogenicity Data** to automatically analyze all reference sequences and save the results. This only requires a motif file to be configured for the entry.

**Pathogenicity Documentation** *(only shown if a `.md` file is present in the entry folder)* — entry-specific documentation such as motif criteria, references, or dataset description.

---

### 5. Precision Validation Tab *(hidden — visit the app with `?dev=1` at the end of the URL)*

Sanity-checks the genotyper's own accuracy against the currently selected entry's reference dataset: it repeatedly holds out a random subset of reference sequences, removes them from the matching pool, and tests them against what's left. Because a single holdout draw can be lucky or unlucky, it reports the **mean accuracy across all runs** (± standard deviation) rather than a one-off score, along with a pooled confusion matrix and the most frequent genotype confusions. Configure the number of runs, the holdout size per run, and the similarity method before running.

---

## Tools Used

- **MAFFT**: Katoh et al. — multiple sequence alignment
- **FastTree**: Price et al. (2009/2010), PLOS ONE — approximate maximum-likelihood trees
- **IQ-TREE2**: Minh et al. (2020), MBE — full maximum-likelihood trees with model selection and ultrafast bootstrap

---

**Last Updated:** July 3, 2026

**Source:** https://github.com/Soupirr/NDV-genotyper

*Sequence Analyzer Toolbox — developed as part of a bioinformatics internship project*

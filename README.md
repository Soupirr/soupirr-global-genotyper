<div align="center">
  
![Python](https://img.shields.io/badge/python-3.9+-3776AB?logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/streamlit-1.57.0-FF4B4B?logo=streamlit&logoColor=white)
![Windows](https://img.shields.io/badge/Windows-0078D4?logo=windows&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-FCC624?logo=linux&logoColor=black)
![GitHub release](https://img.shields.io/github/v/release/Soupirr/soupirr-global-genotyper)
![License](https://img.shields.io/badge/license-MIT-green)
</div>

# Soupirr's Genotyper

A Streamlit web/local application for identifying pathogen genotypes and predicting pathogenicity from nucleotide sequences.

The tool is **pathogen-agnostic**: each "entry" is a self-contained reference dataset (a FASTA database plus optional pathogenicity/motif configuration) for a given pathogen and gene. Bundled entries include Newcastle Disease Virus (F gene), Avian Influenza (HA gene), and Bluetongue/Epizootic Hemorrhagic Disease Virus (VP2 gene), but adding a new pathogen only requires a correctly formatted reference FASTA and no code changes.

For more information on how to use the app and how to add new entries, see [**QUICK_START.md**](https://github.com/Soupirr/NDV-genotyper/blob/main/misc/QUICK_START.md).

---

## Download

Pre-built Windows executable available on the [Releases page](../../releases/latest).

Download `Genotyper_Windows.zip`, unzip, and double-click `Genotyper.exe`.

---

## Run from source (Linux / WSL)

### Requirements

Python 3.9+ is required. Install system dependencies:

```bash
sudo apt install python3 python3-pip mafft fasttree iqtree
```

Install Python dependencies:

```bash
pip install -r requirements.txt
```

Then run:

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

---


## References

Tools used for alignment and tree construction:

- **MAFFT** - Katoh K, Standley DM. *MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability.* Molecular Biology and Evolution, 30(4):772–780, 2013. https://doi.org/10.1093/molbev/mst010

- **FastTree** - Price MN, Dehal PS, Arkin AP. *FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments.* PLOS ONE, 5(3):e9490, 2010. https://doi.org/10.1371/journal.pone.0009490

- **IQ-TREE2** - Minh BQ, Schmidt HA, Chernomor O, et al. *IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era.* Molecular Biology and Evolution, 37(5):1530–1534, 2020. https://doi.org/10.1093/molbev/msaa015

- Dimitrov et al. (2019) - Updated unified phylogenetic classification system and revised
nomenclature for Newcastle disease virus, https://doi.org/10.1016/j.meegid.2019.103917 

- Wang et al. (2017) - Comprehensive analysis of amino acid sequence diversity at the F protein cleavage site of Newcastle disease virus in fusogenic activity, https://doi.org/10.1371/journal.pone.0183923.

### Newcastle Disease Virus Dataset Citation

Dimitrov, K.M., Abolnik, C., Afonso, C.L., Albina, E., Bahl, J., Berg, M., Briand, F.X., Brown, I.H., Choi, K.S., Chvala, I., Diel, D.G., Durr, P.A., Ferreira, H.L., Fusaro, A., Gil, P., Goujgoulova, G.V., Grund, C., Hicks, J.T., Joannis, T.M., Kim Torchetti, M., Kolosov, S., Lambrecht, B., Lewis, N.S., Liu, H., Liu, H., McCullough, S., Miller, P.J., Monne, I., Muller, C.P., Munir, M., Reischak, D., Sabra, M., Samal, S.K., Servan de Almeida, R., Shittu, I., Snoeck, C.J., Suarez, D.L., Van Borm, S., Wang, Z., Wong, F.Y.K., 2019. Updated unified phylogenetic classification system and revised nomenclature for Newcastle disease virus. *Infect. Genet. Evol.*, 103917.
https://doi.org/10.1016/j.meegid.2019.103917

---

## Gallery

<br><br>

<img width="700" alt="figtree_example" src="https://github.com/user-attachments/assets/e4b8cc8f-7cb3-46f6-829b-4beeada857a3" />

<br><br>

<img width="900" alt="Capture d&#39;écran 2026-06-25 101409" src="https://github.com/user-attachments/assets/1f65e62d-9504-4cb6-b958-0e0da67cc254" />

<br><br>

<img width="900" alt="Capture d&#39;écran 2026-06-25 101518" src="https://github.com/user-attachments/assets/2fca2686-8dfc-4afd-804c-959a5308477f" />

<br><br>

<img width="900" alt="Capture d&#39;écran 2026-06-19 150757" src="https://github.com/user-attachments/assets/b7a9c792-5714-4ee4-bb68-3660d6255980" />

<br><br>

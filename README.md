<div align="center">

<img width="100" alt="icon" src="https://github.com/user-attachments/assets/d48fb844-fed7-4c78-ad1d-e1d721a50f13" />

<br><br>
  
![Python](https://img.shields.io/badge/python-3.9+-3776AB?logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/streamlit-1.57.0-FF4B4B?logo=streamlit&logoColor=white)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-lightgrey)
![GitHub release](https://img.shields.io/github/v/release/Soupirr/NDV-genotyper)
![License](https://img.shields.io/badge/license-MIT-green)
</div>

# NDV Genotyper

A Streamlit web application for Newcastle Disease Virus F gene genotyping.

---

## Download

Pre-built Windows executable available on the [Releases page](../../releases/latest).

Download `NDVGenotyper_Windows.zip`, unzip, and double-click `NDVGenotyper.exe`.

---

## Run from source (Linux / WSL)

### Requirements

Python 3.9+ is required. Install system dependencies:

```bash
sudo apt install python3 python3-pip mafft fasttree
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

## Build the Windows executable

The Windows executable requires the MAFFT and FastTree Windows binaries placed in a `tools/` folder (not tracked by git).

1. Download [MAFFT for Windows](https://mafft.cbrc.jp/alignment/software/windows.html) and extract it to `tools/mafft-win/`
2. Download [FastTree](http://www.microbesonline.org/fasttree/) (`FastTree.exe`) and place it in `tools/`
3. Install PyInstaller: `pip install pyinstaller`
4. Run:

```bash
pyinstaller NDVGenotyper.spec
```

Output is in `dist/NDVGenotyper/`. Zip that folder to distribute.

---

For more information on how to use the app and how to add sequences to the database, see **QUICK_START.md**.

---

## References

- **MAFFT** - Katoh K, Standley DM. *MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability.* Molecular Biology and Evolution, 30(4):772–780, 2013. https://doi.org/10.1093/molbev/mst010

- **FastTree** - Price MN, Dehal PS, Arkin AP. *FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments.* PLOS ONE, 5(3):e9490, 2010. https://doi.org/10.1371/journal.pone.0009490

- **Dimitrov et al. (2019)** - Updated unified phylogenetic classification system and revised
nomenclature for Newcastle disease virus, https://doi.org/10.1016/j.meegid.2019.103917

- **Wang et al. (2017)** - Comprehensive analysis of amino acid sequence diversity at the F protein cleavage site of Newcastle disease virus in fusogenic activity, https://doi.org/10.1371/journal.pone.0183923.

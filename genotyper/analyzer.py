"""Core analysis module"""

from typing import Dict, List, Tuple
import subprocess
from Bio import Phylo
import io
import os
import platform
import shutil
import zlib
from genotyper.config import PALETTE

# ============================================================================
# ============================================================================

# Fonctions qui permet d'importer les package FastTree, Mafft et IQ-TREE2 sur Linux peut import ou ils sont


def get_mafft_cmd():
    system_mafft = shutil.which("mafft")
    if system_mafft:
        return system_mafft
    else:  # Windows packaged
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        return os.path.join(base_dir, "tools", "mafft-win", "mafft.bat")


def get_fasttree_cmd():
    system_fasttree = shutil.which("FastTree") or shutil.which("fasttree")
    if system_fasttree:
        return system_fasttree
    else:  # Windows packaged
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        return os.path.join(base_dir, "tools", "FastTree.exe")


def get_iqtree_cmd():
    system_iqtree = shutil.which("iqtree2") or shutil.which("iqtree")
    if system_iqtree:
        return system_iqtree
    else:  # Windows packaged
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        return os.path.join(
            base_dir, "tools", "iqtree-2.4.0-Windows", "bin", "iqtree2.exe"
        )


# ============================================================================
# ============================================================================


# parse_file lit un fichier .fas sur le disque et le transforme en dictionnaire {header: séquence}.
# parse_text même chose mais à partir d'un texte collé. Filtre les caractères illégaux en plus.
class FASTAParser:
    @staticmethod
    def parse_file(filepath: str) -> Dict[str, str]:
        sequences = {}  # Dictionnaire qui va stocker toutes les séquences : {header: séquence}
        current_header = None
        current_seq = []  # Liste qui accumule les lignes de la séquence en cours

        try:
            with open(filepath, "r") as f:  # Ouvre le fichier en mode lecture
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith(">"):
                        if current_header:
                            sequences[current_header] = "".join(current_seq)

                        # New header
                        current_header = line[
                            1:
                        ]  # On stocke le nouveau header sans le caractère >
                        current_seq = []
                    else:
                        # Add to current sequence
                        current_seq.append(line.upper())

                # Save last sequence
                if current_header:
                    sequences[current_header] = "".join(current_seq)

        except Exception as e:
            print(f"Error parsing FASTA file: {e}")
            return {}

        return sequences

    @staticmethod
    def parse_text(fasta_text: str) -> Dict[str, str]:
        sequences = {}
        lines = fasta_text.strip().split("\n")

        current_header = None
        current_seq = []

        for line in lines:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_seq)

                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(
                    "".join(c for c in line.upper() if c in "ATCGNRYSWKMBDHV")
                )

        if current_header:
            sequences[current_header] = "".join(current_seq)

        return sequences


# Class contenant les deux fonction d'analyse (Hamming et Pairwise)
class SequenceSimilarity:
    # Hamming :
    # Compare deux séquences position par position et compte les correspondances.
    # Très rapide mais si les séquences n'ont pas la même longueur, elle tronque simplement à la plus courte
    # ce qui peut fausser les résultats en cas d'insertions/délétions.

    @staticmethod
    def hamming_distance(seq1: str, seq2: str) -> float:
        if len(seq1) != len(seq2):
            # Align to shorter length
            min_len = min(len(seq1), len(seq2))
            seq1 = seq1[:min_len]
            seq2 = seq2[:min_len]

        if len(seq1) == 0:
            return 0.0

        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        similarity = (matches / len(seq1)) * 100

        return round(similarity, 2)

    # Pairwise :
    # Fait un alignement global avant de comparer, ce qui gère correctement les insertions/délétions.
    # Plus lent (~30 secondes) mais biologiquement plus correct.
    # Le score est calculé avec match=1, mismatch=0, gap=-1, puis normalisé sur la longueur maximale des deux séquences.

    @staticmethod
    def pairwise_similarity(seq1: str, seq2: str) -> float:
        from Bio.Align import PairwiseAligner

        if len(seq1) == 0 and len(seq2) == 0:
            return 100.0
        if len(seq1) == 0 or len(seq2) == 0:
            return 0.0

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = 0
        aligner.open_gap_score = -1
        aligner.extend_gap_score = 0

        score = aligner.score(seq1, seq2)
        max_score = max(len(seq1), len(seq2))
        similarity = (score / max_score) * 100

        return round(similarity, 2)


# Class d'analyse du site de cleavage pour l'analyse de la virulence
class CleavageSiteAnalyzer:
    # Cleavage site positions (approximate in nucleotides)
    # Position ~334-357 bp in F gene (codons 112-119)
    CLEAVAGE_START = 333  # 0-indexed
    CLEAVAGE_LENGTH = 24  # 8 codons = 24 nucleotides
    WINDOW = 30  # pour prendre en compte les mutations d'addition

    CODON_TABLE = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGC": "C",
        "TGT": "C",
        "TGA": "*",
        "TGG": "W",
    }
    # https://gist.github.com/juanfal/09d7fb53bd367742127e17284b9c47bf

    # Traduit une séquence ADN en protéine codon par codon (3 nucléotides à la fois).
    # Si le codon contient un N (nucléotide ambigu) ou est incomplet, il retourne X.
    @staticmethod
    def translate_to_protein(dna_sequence: str) -> str:
        protein = []
        for i in range(0, len(dna_sequence) - 2, 3):
            codon = dna_sequence[i : i + 3].upper()
            # Handle ambiguous nucleotides
            if "N" in codon or len(codon) < 3:
                protein.append("X")
            else:
                aa = CleavageSiteAnalyzer.CODON_TABLE.get(codon, "X")
                protein.append(aa)
        return "".join(protein)

    @staticmethod
    def analyze(
        sequence: str,
        cleavage_start: int,
        motifs_by_type: dict = None,
    ):
        motifs_by_type = motifs_by_type or {}

        WINDOW = 29 * 3  # ± tolerance for indels

        result = {
            "cleavage_region_found": False,
            "cleavage_nucleotides": None,
            "cleavage_protein": None,
            "pathogenicity": "Undetermined",
            "motif_type": "No known motif found in cleavage region",
            "motif_category": None,
        }

        result_plus_one = result.copy()
        result_minus_one = result.copy()

        # Vérifie que la séquence est assez longue
        if len(sequence) < cleavage_start - WINDOW:
            return result, result_plus_one, result_minus_one

        # Extrait la région autour du site de clivage avec la fenêtre de tolérance
        region_start = max(0, cleavage_start - WINDOW)
        region_end = min(len(sequence), cleavage_start + WINDOW)
        cleavage_region_nuc = sequence[region_start:region_end]
        # Extraction des cadres de lecture +1 -1
        region_start_plus_one = max(0, cleavage_start + 1 - WINDOW)
        region_end_plus_one = min(len(sequence), cleavage_start + 1 + WINDOW)
        cleavage_region_nuc_plus_one = sequence[
            region_start_plus_one:region_end_plus_one
        ]
        region_start_minus_one = max(0, cleavage_start - 1 - WINDOW)
        region_end_minus_one = min(len(sequence), cleavage_start - 1 + WINDOW)
        cleavage_region_nuc_minus_one = sequence[
            region_start_minus_one:region_end_minus_one
        ]

        # Traduit la région principale en protéine
        cleavage_region_prot = CleavageSiteAnalyzer.translate_to_protein(
            cleavage_region_nuc
        )
        cleavage_region_prot_plus_one = CleavageSiteAnalyzer.translate_to_protein(
            cleavage_region_nuc_plus_one
        )
        cleavage_region_prot_minus_one = CleavageSiteAnalyzer.translate_to_protein(
            cleavage_region_nuc_minus_one
        )

        result["cleavage_region_found"] = True
        result["cleavage_nucleotides"] = cleavage_region_nuc
        result["cleavage_protein"] = cleavage_region_prot

        result_plus_one["cleavage_region_found"] = True
        result_plus_one["cleavage_nucleotides"] = cleavage_region_nuc_plus_one
        result_plus_one["cleavage_protein"] = cleavage_region_prot_plus_one

        result_minus_one["cleavage_region_found"] = True
        result_minus_one["cleavage_nucleotides"] = cleavage_region_nuc_minus_one
        result_minus_one["cleavage_protein"] = cleavage_region_prot_minus_one

        for type_name, motifs in motifs_by_type.items():
            for frames, result_dict in [
                (cleavage_region_prot, result),
                (cleavage_region_prot_plus_one, result_plus_one),
                (cleavage_region_prot_minus_one, result_minus_one),
            ]:
                for motif, label in motifs.items():
                    if motif in frames:
                        result_dict["pathogenicity"] = type_name
                        result_dict["motif_type"] = motif
                        result_dict["motif_category"] = label
                        return result, result_plus_one, result_minus_one

        return result, result_plus_one, result_minus_one


# Class qui comparer une séquence inconnue contre toute la base de données
# et retourner les génotypes les plus similaires avec leurs statistiques.
class GenotypeIdentifier:
    def __init__(self, references: Dict[str, str]):
        # Initialise la classe avec le dictionnaire de séquences de référence
        self.references = references

    def identify(self, input_sequence: str, method="hamming", top_n=3) -> List[Tuple]:
        results = []

        # Choisit la méthode de comparaison (Hamming ou Pairwise)
        if method == "hamming":
            similarity_func = SequenceSimilarity.hamming_distance
        else:
            similarity_func = SequenceSimilarity.pairwise_similarity

        # Regroupe les scores par génotype et calcule la moyenne
        genotype_scores = {}
        genotype_best_match = {}

        for header, ref_sequence in self.references.items():
            similarity = similarity_func(input_sequence, ref_sequence)

            # Header : ID_GENOTYPE_...
            parts = header.split("|")
            genotype = parts[2] if len(parts) >= 3 else "Unknown"

            if genotype not in genotype_scores:
                genotype_scores[genotype] = []
                genotype_best_match[genotype] = (0, header)

            genotype_scores[genotype].append(similarity)

            # Garde en mémoire le meilleur match individuel pour chaque génotype
            if similarity > genotype_best_match[genotype][0]:
                genotype_best_match[genotype] = (similarity, header)

        # Calculate averages and prepare results
        for genotype, scores in genotype_scores.items():
            avg_score = sum(scores) / len(scores)
            best_score, best_header = genotype_best_match[genotype]

            results.append(
                (
                    genotype,
                    round(avg_score, 2),
                    len(scores),
                    best_header,
                    round(best_score, 2),
                )
            )

        # Trie les génotypes par score moyen décroissant
        results.sort(key=lambda x: x[1], reverse=True)

        return results[:top_n]


# ============================================================================
# ============================================================================


# fonction qui orchestre l'ensemble du pipeline :
# parse la séquence d'entrée,
# identifie le génotype via GenotypeIdentifier,
# analyse le site de clivage via CleavageSiteAnalyzer,
# et retourne tous les résultats dans un seul dictionnaire.
def analyze_sequence(
    input_fasta: str,
    reference_sequences: Dict[str, str],
    top_matches: int = 3,
    similarity_method: str = "hamming",
    pathogenicity_config: dict = None,  # utilisation de Hamming par default
) -> Dict:

    # Parse input
    parsed_input = FASTAParser.parse_text(input_fasta)

    if not parsed_input:
        return {"error": "Could not parse input FASTA"}

    # obtenir la première séquence
    input_header = list(parsed_input.keys())[0]
    input_sequence = parsed_input[input_header]

    # Identifier le génotype
    identifier = GenotypeIdentifier(reference_sequences)
    genotype_matches = identifier.identify(
        input_sequence, method=similarity_method, top_n=top_matches
    )

    cfg = pathogenicity_config or {}
    pathogenicity_configured = "cleavage_start" in cfg

    if pathogenicity_configured:
        cleavage_main, cleavage_plus_one, cleavage_minus_one = (
            CleavageSiteAnalyzer.analyze(
                input_sequence,
                cleavage_start=cfg["cleavage_start"],
                motifs_by_type=cfg.get("motifs_by_type"),
            )
        )
    else:
        _empty = {
            "cleavage_region_found": False,
            "cleavage_nucleotides": None,
            "cleavage_protein": None,
            "pathogenicity": "Not configured",
            "motif_type": None,
            "motif_category": None,
        }
        cleavage_main = cleavage_plus_one = cleavage_minus_one = _empty.copy()

    return {
        "input_header": input_header,
        "sequence_length": len(input_sequence),
        "genotype_matches": genotype_matches,
        "pathogenicity_configured": pathogenicity_configured,
        "cleavage_main": cleavage_main,
        "cleavage_plus_one": cleavage_plus_one,
        "cleavage_minus_one": cleavage_minus_one,
        "error": None,
    }


# fonction qui convertit le tuple brut retourné par GenotypeIdentifier.identify() en dictionnaire
def unpack_top_match(top_match):
    return {
        "genotype": top_match[0],
        "avg_similarity": top_match[1],
        "sample_count": top_match[2],
        "best_header": top_match[3],
        "best_score": top_match[4],
    }


def tree_to_newick(tree) -> str:
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue()


# fonction de construction de l'arbre avec FastTree
def build_tree_fasttree(aln_file):
    win_flags = (
        {"creationflags": subprocess.CREATE_NO_WINDOW}
        if platform.system() == "Windows"
        else {}
    )
    result = subprocess.run(
        [get_fasttree_cmd(), "-nt", "-gtr", aln_file],
        capture_output=True,
        text=True,
        **win_flags,
    )
    newick_str = result.stdout
    if not newick_str or not newick_str.strip():
        return None
    tree = Phylo.read(io.StringIO(newick_str), "newick")
    tree.root_at_midpoint()  # ultra important ça me pose des problèmes depuis 1 mois
    return tree


# fonction de construction de l'arbre avec IQ-TREE2 (ModelFinder + ultrafast bootstrap)
# même principe que pour fasttree
def build_tree_iqtree2(aln_file):
    win_flags = (
        {"creationflags": subprocess.CREATE_NO_WINDOW}
        if platform.system() == "Windows"
        else {}
    )
    subprocess.run(
        [
            get_iqtree_cmd(),
            "-s",
            aln_file,
            "-m",
            "MFP",
            "-bb",
            "1000",
            "-nt",
            "AUTO",
            "-redo",  # nécessaire parce que le pipeline "Per-query" réutilise le même chemin de fichier temporaire (tmp_aligned.fasta) à chaque itération de la boucle - sans -redo, IQ-TREE2 refuse de tourner une 2e fois sur un fichier de sortie déjà existant
        ],
        capture_output=True,
        text=True,
        **win_flags,
    )
    treefile = aln_file + ".treefile"
    if not os.path.exists(treefile):
        return None
    tree = Phylo.read(treefile, "newick")
    tree.root_at_midpoint()  # pareil que pour fasttree
    # IQ-TREE2 exprime le bootstrap ultrafast en 0-100, FastTree en 0-1 (voir tree_to_plotly) - on normalise ici
    for clade in tree.find_clades():
        if clade.confidence is not None:
            clade.confidence = clade.confidence / 100
    return tree


# alignement des séquences avec Mafft
def align_sequences_mafft(input_fasta_path, output_fasta_path):
    win_flags = (
        {"creationflags": subprocess.CREATE_NO_WINDOW}
        if platform.system() == "Windows"
        else {}
    )
    with open(output_fasta_path, "w") as out:
        subprocess.run(
            [get_mafft_cmd(), "--localpair", input_fasta_path],
            stdout=out,
            stderr=subprocess.PIPE,
            **win_flags,
        )


# trouve les voisins les plus proche de notre input (par default: 20)
def find_closest_neighbours(query_sequence, reference_sequences, n=20):
    similarities = []

    for header, seq in reference_sequences.items():
        score = SequenceSimilarity.hamming_distance(query_sequence, seq)
        similarities.append((header, seq, score))

    # Sort by similarity descending, take top N
    similarities.sort(key=lambda x: x[2], reverse=True)
    return similarities[:n]


# fonction qui permet de nettoyé notre input
def clean_sequence(seq):
    valid = set("ATCGNatcgnRYSWKMBDHVryswkmbdhv")
    return "".join(c for c in seq if c in valid)


# écrit les fichier fasta temporaires qui seront utilisé par FastTree
def write_temp_fasta(query_header, query_sequence, neighbours, output_path):
    with open(output_path, "w") as f:
        f.write(f">QUERY_{query_header}\n{clean_sequence(query_sequence)}\n")
        for header, seq, score in neighbours:
            f.write(f">{header}\n{clean_sequence(seq)}\n")


# ============================================================================
# ============================================================================


# Color mapping (generated)
def get_color(name):
    if not name:
        return "#888888"
    parts = name.split("|")
    if len(parts) < 3:
        return "#888888"
    genotype = parts[2]
    if genotype in ("?", "UNKNOWN", "") or genotype.startswith("UNCL"):
        return "#888888"
    idx = zlib.crc32(genotype.encode()) % len(PALETTE)
    return PALETTE[idx]

"""
Newcastle Disease Virus F Gene Genotyper
Core analysis module
"""

from typing import Dict, List, Tuple
import subprocess
from Bio import Phylo
import io
import os

# ============================================================================
# ============================================================================

# Fonctions qui permet d'importer les package FastTree et Mafft soit à partir de Windows soit à partir de WSL


def get_mafft_cmd():
    if os.path.exists("/usr/bin/mafft"):  # WSL/Linux
        return "mafft"
    else:  # Windows
        base_dir = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(base_dir, "tools", "mafft-win", "mafft.bat")


def get_fasttree_cmd():
    if os.path.exists("/usr/bin/FastTree"):  # WSL/Linux
        return "FastTree"
    else:  # Windows
        base_dir = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(base_dir, "tools", "FastTree.exe")


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

    # Known virulent motifs (VFcs) from Yanhong Wang et al. 2017
    VIRULENT_MOTIFS = {
        "RRQKRF": "VFcs-1",
        "RRQRRF": "VFcs-2",
        "KRQKRF": "VFcs-3",
        "GRQKRF": "VFcs-4",
        "RRKKRF": "VFcs-5",
        "RRRKRF": "VFcs-6",
        "KRKKRF": "VFcs-7",
        "RRRRRF": "VFcs-8",
    }

    # Known avirulent motifs (AFcs) from Yanhong Wang et al. 2017
    AVIRULENT_MOTIFS = {
        "GRQGRL": "AFcs-1",
        "GKQGRL": "AFcs-2",
        "RRQGRL": "AFcs-3",
        "ERQERL": "AFcs-4",
        "RRQGRF": "AFcs-5",  # avirulent même si F
        "ERQGRL": "AFcs-6",
        "RKQGRL": "AFcs-7",
        "EKQGRL": "AFcs-8",
        "EQQERL": "AFcs-9",
        "RRQRRL": "AFcs-10",
    }

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
    def analyze(sequence: str) -> Dict:

        result = {
            "cleavage_region_found": False,
            "cleavage_nucleotides": None,
            "cleavage_protein": None,
            "pathogenicity": None,
            "confidence": "Low",
            "motif_type": None,
            "motif_category": None,
        }

        CLEAVAGE_START = 333
        CLEAVAGE_LENGTH = 24
        WINDOW = 30  # ± tolerance for indels

        # Vérifie que la séquence est assez longue
        if len(sequence) < CLEAVAGE_START - WINDOW + CLEAVAGE_LENGTH:
            return result

        # Extrait la région autour du site de clivage avec la fenêtre de tolérance
        region_start = max(0, CLEAVAGE_START - WINDOW)
        region_end = min(len(sequence), CLEAVAGE_START + CLEAVAGE_LENGTH + WINDOW)
        cleavage_region_nuc = sequence[region_start:region_end]

        # Traduit cette région en protéine
        cleavage_region_prot = CleavageSiteAnalyzer.translate_to_protein(
            cleavage_region_nuc
        )

        result["cleavage_region_found"] = True
        result["cleavage_nucleotides"] = cleavage_region_nuc
        result["cleavage_protein"] = cleavage_region_prot

        # Cherche d'abord les motifs virulents,
        for motif, category in CleavageSiteAnalyzer.VIRULENT_MOTIFS.items():
            if motif in cleavage_region_prot:
                result["pathogenicity"] = "Likely Virulent"
                result["motif_type"] = motif
                result["motif_category"] = category
                result["confidence"] = "High"
                return result

        # puis les avirulents
        for motif, category in CleavageSiteAnalyzer.AVIRULENT_MOTIFS.items():
            if motif in cleavage_region_prot:
                result["pathogenicity"] = "Likely Low-virulence"
                result["motif_type"] = motif
                result["motif_category"] = category
                result["confidence"] = "High"
                return result

        # Si rien n'est trouvé:
        result["pathogenicity"] = "Undetermined"
        result["motif_type"] = "No known motif found in cleavage region"
        result["motif_category"] = None
        result["confidence"] = "Low"
        return result
        # Retourne le résultat avec le motif trouvé, sa catégorie et la pathogénicité prédite


# Class qui comparer une séquence inconnue contre toute la base de données
# et retourner les génotypes les plus similaires avec leurs statistiques.
class GenotypeIdentifier:
    def __init__(self, references: Dict[str, str]):
        # Initialise la classe avec le dictionnaire de séquences de référence
        self.references = references
        self.genotype_map = self._build_genotype_map()

    def _build_genotype_map(self) -> Dict[str, List[str]]:
        genotype_map = {}

        for header in self.references.keys():
            # Parcourt tous les headers de la base de données et les regroupe par génotype.
            # Exemple: 36_I.1.1_I_a_AF217084_chicken_Australia_Queensland_V_4_1966
            parts = header.split("_")
            if len(parts) >= 2:
                genotype = parts[1]  # e.g., "I.1.1", "II", "VII.1.1"

                if genotype not in genotype_map:
                    genotype_map[genotype] = []
                genotype_map[genotype].append(header)

        return genotype_map

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
            parts = header.split("_")
            genotype = parts[1] if len(parts) >= 2 else "Unknown"

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
def analyze_newcastle_sequence(
    input_fasta: str,
    reference_sequences: Dict[str, str],
    top_matches: int = 3,
    similarity_method: str = "hamming",  # utilisation de Hamming par default
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

    # lancer l'analyse du site de cleavage
    cleavage_analysis = CleavageSiteAnalyzer.analyze(input_sequence)

    return {
        "input_header": input_header,
        "sequence_length": len(input_sequence),
        "genotype_matches": genotype_matches,
        "cleavage_analysis": cleavage_analysis,
        "error": None,
    }


# fonction qui permet d'obtenir la Class d'un génotype en fonction de son nom
def get_class(genotype):
    if genotype.startswith("1"):
        return "Class I"
    else:
        return "Class II"


# fonction qui convertit le tuple brut retourné par GenotypeIdentifier.identify() en dictionnaire
def unpack_top_match(top_match):
    return {
        "genotype": top_match[0],
        "avg_similarity": top_match[1],
        "sample_count": top_match[2],
        "best_header": top_match[3],
        "best_score": top_match[4],
    }


# fonction de construction de l'arbre avec FastTree
def build_tree_fasttree(aln_file):
    result = subprocess.run(
        [get_fasttree_cmd(), "-nt", "-gtr", aln_file], capture_output=True, text=True,
        creationflags=subprocess.CREATE_NO_WINDOW,
    )
    newick_str = result.stdout
    tree = Phylo.read(io.StringIO(newick_str), "newick")
    return tree


# alignement des séquences avec Mafft
def align_sequences_mafft(input_fasta_path, output_fasta_path):
    with open(output_fasta_path, "w") as out:
        subprocess.run(
            [get_mafft_cmd(), "--auto", input_fasta_path],
            stdout=out,
            stderr=subprocess.PIPE,
            creationflags=subprocess.CREATE_NO_WINDOW,
        )


# trouve les 20 voisins les plus proche de notre input
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
CLASS_I = {"1.1", "1.1.1", "1.1.2", "1.2", "I.1.1", "I.1.2.1", "I.1.2.2", "I.2"}

CLASS_II_V = {
    "II": "#1a9850",
    "III": "#52b788",
    "IV": "#74c69d",
    "V": "#95d5b2",
    "V.1": "#b7e4c7",
    "V.2": "#d8f3dc",
    "V.3": "#40916c",
}

CLASS_VI = {
    "VI.1": "#e85d04",
    "VI.2.1": "#f48c06",
    "VI.2.1.1.1": "#faa307",
    "VI.2.1.1.2.1": "#ffba08",
    "VI.2.1.1.2.2": "#ffd60a",
    "VI.2.1.1.2.2.": "#ffdd00",
    "VI.2.1.2": "#f46036",
    "VI.2.2.1": "#e2711d",
    "VI.2.2.2": "#cc5803",
}

CLASS_VII = {
    "VII": "#d00000",
    "VII.1": "#e85252",
    "VII.1.1": "#ff6b6b",
    "VII.1.2": "#ff8fa3",
    "VII.2": "#a4133c",
}

CLASS_VIII_XIV = {
    "VIII": "#7b2d8b",
    "IX": "#9b4dca",
    "X": "#b57bee",
    "X.1": "#c77dff",
    "X.2": "#d4a5ff",
    "XI": "#6a0572",
    "XII": "#8338ec",
    "XII.1": "#9d4edd",
    "XII.2": "#c77dff",
    "XIII.1": "#5a189a",
    "XIII.1.1": "#6d23b6",
    "XIII.1.2": "#7b2d8b",
    "XIII.2": "#480ca8",
    "XIII.2.1": "#3a0ca3",
    "XIII.2.2": "#3c096c",
    "XIV": "#e0aaff",
    "XIV.1": "#c8b6ff",
    "XIV.2": "#b8c0ff",
}

CLASS_XV_PLUS = {
    "XVI": "#f72585",
    "XVII": "#ff4d6d",
    "XVIII.1": "#ff758f",
    "XVIII.2": "#ff85a1",
    "XIX": "#ffb3c1",
    "XX": "#fb6f92",
    "XXI": "#ff0a54",
    "XXI.1.1": "#ff477e",
    "XXI.1.2": "#ff5c8a",
    "XXI.2": "#c9184a",
}


# Fonction qui assigne chaque class à une couleur différente
def get_color(name):
    if not name:
        return "#888888"
    parts = name.split("_")
    if len(parts) < 2:
        return "#888888"
    genotype = parts[1]

    if genotype in CLASS_I:
        return "#4A90D9"
    elif genotype.startswith("UNCL"):
        return "#888888"
    elif genotype in CLASS_II_V:
        return CLASS_II_V[genotype]
    elif genotype in CLASS_VI:
        return CLASS_VI[genotype]
    elif genotype in CLASS_VII:
        return CLASS_VII[genotype]
    elif genotype in CLASS_VIII_XIV:
        return CLASS_VIII_XIV[genotype]
    elif genotype in CLASS_XV_PLUS:
        return CLASS_XV_PLUS[genotype]
    else:
        return "#888888"


# ============================================================================

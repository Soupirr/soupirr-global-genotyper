"""NCBI Virus -> Toolbox FASTA header migration."""

import pandas as pd
from genotyper.analyzer import FASTAParser
from genotyper.config import HOSTS_FOLDER, LOCATION_FOLDER
import os


def load_migration_data():
    df_hosts = pd.read_csv(os.path.join(HOSTS_FOLDER, "host_normalize.csv"))
    host_normalize = dict(
        zip(df_hosts["raw_host"].str.lower(), df_hosts["normalized_host"])
    )

    known_regions = {}
    df_us = pd.read_csv(os.path.join(LOCATION_FOLDER, "US_capitals.csv"))
    for name in df_us["state"]:
        known_regions[name.lower()] = name

    df_cn = pd.read_csv(os.path.join(LOCATION_FOLDER, "china_provinces.csv"))
    for name in df_cn["province"]:
        known_regions[name.lower()] = name

    df_ru = pd.read_csv(os.path.join(LOCATION_FOLDER, "russia_regions.csv"))
    for name in df_ru["region_ru"]:
        canonical = name.replace("_", " ")
        known_regions[name.lower()] = canonical
        known_regions[canonical.lower()] = canonical

    df_world = pd.read_csv(os.path.join(LOCATION_FOLDER, "countries.csv"))
    known_countries = {}
    for name in df_world["country"]:
        canonical = name.replace("_", " ")
        known_countries[name.lower()] = canonical
        known_countries[canonical.lower()] = canonical
    # overrides manuels
    known_countries["usa"] = "United States"
    known_countries["korea"] = "South Korea"
    known_countries["viet nam"] = "Vietnam"
    known_countries["czech republic"] = "Czechia"

    return host_normalize, known_regions, known_countries


def extract_year(date_str):
    return date_str.strip().split("-")[0]


def find_region(geo_location, known_regions):
    tokens = geo_location.replace(":", " ").replace(",", " ").split()
    for token in tokens:
        if token.lower() in known_regions:
            return known_regions[token.lower()].replace(" ", "_")
    return "?"


def convert_header(header, host_normalize, known_regions, known_countries):
    parts = [p.strip() for p in header.split("|")]
    if len(parts) < 7:
        return None

    virus = parts[0].replace(" ", "_")
    accession = parts[1].split(".")[0]
    genotype = parts[2].replace(" ", "_")
    host_raw = parts[3].lower()
    host = host_normalize.get(host_raw, parts[3]).replace(" ", "_")
    country_raw = parts[4].strip().lower()
    country = known_countries.get(country_raw, parts[4]).replace(" ", "_")
    region = find_region(parts[5], known_regions)
    year = extract_year(parts[6])

    return f"{virus}|{accession}|{genotype}|{host}|{country}|{region}|{year}"


# migre un texte FASTA brut NCBI vers le format toolbox :
# conversion des headers, suppression des séquences sans génotype,
# dédup par contenu de séquence identique (accession différente, séquence identique)
def migrate_fasta_text(raw_text: str):
    host_normalize, known_regions, known_countries = load_migration_data()
    parsed = FASTAParser.parse_text(raw_text)

    stats = {
        "input": len(parsed),
        "converted": 0,
        "no_genotype": 0,
        "duplicates": 0,
        "malformed": 0,
    }
    seen_sequences = set()
    lines = []

    for header, sequence in parsed.items():
        new_header = convert_header(
            header, host_normalize, known_regions, known_countries
        )
        if new_header is None:
            stats["malformed"] += 1
            continue

        genotype = new_header.split("|")[2]
        if not genotype:
            stats["no_genotype"] += 1
            continue

        if sequence in seen_sequences:
            stats["duplicates"] += 1
            continue
        seen_sequences.add(sequence)

        lines.append(f">{new_header}\n{sequence}\n")
        stats["converted"] += 1

    return "".join(lines), stats

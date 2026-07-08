# QClair_Natural_Products_Database (QC NPDB)
Qclairvoyance Natural Product Database (QClair NPDB): A Curated and Comprehensive Natural Product Repository for Researchers

---

# Natural Product Database Construction, Curation, and Enrichment Pipeline

This repository contains a collection of Python scripts developed for the construction, curation, validation, and enrichment of the **QClair Natural Product Database (QC-NPDB)**. The pipeline integrates multiple publicly available natural product repositories into a unified, standardized, and quality-controlled database.

The current release of the database contains:

* **1,084,361 natural product compounds**
* **Comprehensive molecular, structural, taxonomic, and classification annotations**

## Project Objectives

The primary objectives of this project are to:

* Integrate compounds from multiple public natural product databases
* Remove duplicate and inconsistent records
* Standardize molecular representations
* Fill missing molecular and taxonomic information
* Generate additional molecular descriptors using RDKit
* Validate molecular structures and identifiers
* Assign unique QClair identifiers (QIDs) for every compound

---

# 1. `missing_value_01.py`

## Overview

This script performs the initial preprocessing and exploratory analysis of the integrated dataset.

### Key Functions

* Loads the raw integrated dataset
* Performs exploratory data analysis (EDA)
* Identifies missing values across all columns
* Validates SMILES structures using RDKit
* Generates missing molecular formulas
* Generates missing SELFIES representations
* Attempts Murcko Framework generation
* Extracts missing IUPAC names using multiple sources
* Saves the cleaned intermediate dataset

---

# 2. Scientific Name Cleaning and Database Integration

The `scientific_name` field contained inconsistent, incomplete, and heterogeneous values originating from multiple source databases.

A dedicated cleaning and standardization pipeline was implemented to improve taxonomic consistency.

## Integrated Data Sources

* COCONUT
* LOTUS
* CMAUP
* NPASS
* NPAtlas
* StreptomeDB
* Seaweed Metabolite Database

These curated datasets were subsequently used for metadata enrichment throughout the pipeline.

---

# 3. `missing_value_02.py`

## Overview

This module performs automated enrichment using the cleaned source databases.

### Major Tasks

* Standardizes column names across repositories
* Cleans and harmonizes source datasets
* Matches compounds using molecular identifiers
* Fills missing information for:

  * Source database identifiers
  * IUPAC names
  * Species
  * Genus
  * Scientific names
  * Murcko frameworks
  * SELFIES

---

# 4. Web Scraping Modules

Several scripts were developed to retrieve missing information from external databases when it was unavailable in the integrated datasets.

| Script                     | Purpose                                                                    |
| -------------------------- | -------------------------------------------------------------------------- |
| `beautiful_soup.py`        | Retrieves molecular and taxonomic information using PubChem and ChemSpider |
| `iupac_name_extraction.py` | Extracts missing IUPAC names                                               |
| `scrape_iupac.py`          | Retrieves standardized IUPAC nomenclature                                  |
| `scrape_genus_species.py`  | Retrieves taxonomy information                                             |
| `web_scrapping.py`         | Collects additional metadata from public NP resources                      |
| `web_scrapping_1.py`       | Backup scraping pipeline for missing annotations                           |

---

# 5. `missing_value_05.py`

## Overview

This script enriches the database by generating additional molecular descriptors and annotations.

### Newly Generated Features

* Canonical SMILES
* Isomeric SMILES
* Common Name
* Chemical Class
* Direct Parent Classification
* NP Classifier Pathway
* NP Classifier Superclass
* NP Classifier Class
* NP Classifier Glycoside Annotation

### Processing Steps

* Generates standardized molecular representations using RDKit
* Computes missing descriptors
* Integrates newly generated annotations into the master database

---

# 6. `missing_value_06.py`

## Overview

This module performs comprehensive validation and correction of molecular structures.

### Validation Pipeline

* Validates SMILES
* Validates Canonical SMILES
* Validates Isomeric SMILES
* Detects invalid or corrupted molecular structures

### Recovery Strategy

Invalid structures are corrected using a two-stage approach:

1. PubChem lookup using the compound InChIKey
2. RDKit conversion from Standard InChI to SMILES (fallback)

This process ensures that missing or invalid molecular representations are recovered wherever possible.

---

# 7. Database Standardization

Additional quality-control procedures include:

* Duplicate detection and removal
* Identifier validation
* Molecular formula verification
* Canonicalization of molecular structures
* Taxonomic data standardization
* Molecular classification integration
* Consistency checks across all source databases

---

# Final Database Summary

After integration, cleaning, enrichment, validation, and descriptor generation, the latest release contains:

* **Total Compounds:** **1,084,361**
* **Unique QClair IDs (QIDs):** **1,084,361**
* **Integrated from seven public natural product repositories**
* **Validated molecular structures and standardized annotations**

---

# QID Assignment System

Each compound is assigned a unique **QClair Identifier (QID)**.

The repository prefixes are:

| Prefix | Source Database                                   |
| ------ | ------------------------------------------------- |
| QCAA   | COCONUT                                           |
| QCAB   | LOTUS                                             |
| QCAC   | CMAUP                                             |
| QCAD   | NPASS                                             |
| QCAE   | NPAtlas                                           |
| QCAF   | StreptomeDB                                       |
| QCAG   | Seaweed Metabolite Database                       |
| QCAH   | Other integrated sources / newly assigned records |

This QID system ensures global uniqueness, traceability, and consistent identification across the integrated database.

---

# Tools and Libraries

* Python 3
* Pandas
* NumPy
* RDKit
* Requests
* BeautifulSoup
* PubChemPy
* ChemSpiPy
* CSV processing and data integration utilities

---

# Final Output

The final processed database includes:

* Validated molecular structures
* Standardized SMILES representations
* Canonical and Isomeric SMILES
* Standard InChI and InChIKey
* Molecular formula
* Taxonomic information
* Molecular classifications
* NP Classifier annotations
* Murcko frameworks
* SELFIES representations
* Additional curated metadata

---

# QC-NPDB Website

https://npdb.qclairvoyance.com/

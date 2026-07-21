# QClair Natural Product Database

A Python-based cheminformatics pipeline for constructing a curated natural product database through automated molecular preprocessing, structural validation, descriptor computation, and database generation.

Database URL: https://npdb.qclairvoyance.com/

---

## Overview

Large-scale molecular datasets frequently contain duplicate structures, inconsistent molecular representations, missing annotations, and heterogeneous file formats, limiting their applicability in computational chemistry and drug discovery.

This repository provides a complete and reproducible workflow for processing molecular datasets into a chemically consistent, analysis-ready natural product database. The pipeline automates molecular standardization, structural validation, descriptor generation, missing value reconstruction, and database preparation for downstream cheminformatics applications.

The project is designed to support researchers, computational chemists, bioinformaticians, and developers working in natural product research and AI-assisted drug discovery.

---

## Key Features

- Automated molecular data preprocessing
- Multi-format molecular data parsing
- Structural validation using RDKit
- Molecular standardization
- Duplicate detection and removal
- Missing value reconstruction
- Molecular descriptor computation
- Physicochemical property calculation
- Chemical classification integration
- Molecular scaffold generation
- SMARTS generation
- SELFIES generation
- Database-ready output generation
- PostgreSQL and RDKit compatible workflow
- Modular Python implementation

---

## Workflow

```text
Raw Molecular Data
        │
        ▼
Data Collection
        │
        ▼
Format Standardization
        │
        ▼
Structural Validation
        │
        ▼
Molecular Standardization
        │
        ▼
Duplicate Detection
        │
        ▼
Missing Value Reconstruction
        │
        ▼
Descriptor Computation
        │
        ▼
Quality Validation
        │
        ▼
Curated Molecular Database
```

---

## Molecular Features Generated

The pipeline supports the generation and validation of a comprehensive set of molecular representations, descriptors, and annotations, including:

### Molecular Representations

- Standard InChI
- InChIKey
- Canonical SMILES
- Isomeric SMILES
- SMARTS
- SELFIES
- Murcko Framework

### Physicochemical Descriptors

- Molecular Formula
- Molecular Weight
- Exact Molecular Weight
- Heavy Atom Count
- Total Atom Count
- Bond Count
- Rotatable Bond Count
- Aromatic Ring Count
- Hydrogen Bond Donor (HBD)
- Hydrogen Bond Acceptor (HBA)
- LogP
- Topological Polar Surface Area (TPSA)
- Formal Charge
- van der Waals Surface Area
- van der Waals Volume

### Chemical Annotation

- Chemical Classification
- Parent Classification
- Biosynthetic Pathway
- Superclass
- Chemical Class
- Glycoside Classification

### Drug Discovery Descriptors

- QED Drug-Likeness
- Natural Product Likeness (NP-Likeness)

---

## Technologies Used

- Python 3.x
- RDKit
- PostgreSQL
- Open Babel
- Pandas
- NumPy
- SciPy
- PubChemPy
- Scikit-learn

---

## Applications

The generated database can be applied to:

- Cheminformatics
- Drug Discovery
- Computational Chemistry
- Bioinformatics
- Virtual Screening
- QSAR Modeling
- Molecular Similarity Search
- Chemical Space Analysis
- Scaffold Analysis
- Machine Learning
- Deep Learning
- AI-assisted Molecular Discovery

---

## Reproducibility

The workflow has been developed with reproducibility in mind. Each processing stage is modular, enabling independent execution, validation, and extension for different molecular datasets.

Researchers can adapt individual modules for their own database construction workflows or integrate additional cheminformatics analyses.

---


## License

This project is licensed under QCLAIR NPDB PROPRIETARY DATABASE

---

## Acknowledgements

This repository was developed as part of the QClair Natural Product Database project to support reproducible cheminformatics research and natural product based computational drug discovery.

---

## Contact

For questions, suggestions, or collaborations, please open an Issue in this repository.

Qclairvoyance Quantum Labs Private Limited, 
191, Hi-Tension Road, Sainikpuri, Secunderabad, Hyderabad, Telangana, India - 500094.

Email: info@qclairvoyance.in
Website: https://www.qclairvoyance.com

---

© 2026 Qclairvoyance Quantum Labs. All rights reserved.

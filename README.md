# RHDV Database

> An open, FAIR-compliant primer resource for **Rabbit Haemorrhagic Disease Virus (RHDV)** detection, integrating literature-derived primers (**AROLit**) and in silico designed primers (**iSOP**) to support PCR-based detection and environmental surveillance.

---

## Table of Contents

- [Background](#background)
- [Repository Structure](#repository-structure)
- [Datasets](#datasets)
- [Workflows](#workflows)
- [Validation](#validation)
- [How to Use](#how-to-use)
- [Citation](#citation)
- [Authors](#authors)
- [License](#license)

---

## Background

Rabbit haemorrhagic disease virus (RHDV) impacts the European rabbit (*Oryctolagus cuniculus*) and was first discovered in China in 1984. Since then, it has spread globally, causing significant rabbit population declines in both wild and domestic contexts.

To address the growing volume of RHDV genomic data, this project developed two complementary primer resources:

- **AROLit** — an automated computational method to retrieve and curate primers from the scientific literature.
- **iSOP** — a Python-based workflow to design primers *in silico* for RHDV detection.

Both workflows support optimal primer-pair selection for PCR applications and are transferable to other viral or bacterial genomes.

---

## Repository Structure

```
rhdv-primers-identification-db/
├── data/
│   ├── AROLit.xlsx           # Literature-derived primer dataset
│   └── iSOP/                 # In silico designed primer dataset (unzipped)
├── scripts/
│   └── AROLit_and_iSOP.py    # Main Python script for AROLit and iSOP workflows
├── workflows/
│   └── Workflow RHDV Orange.ows  # Orange Data Mining workflow file
├── LICENSE
└── README.md
```

---

## Datasets

The two core datasets are archived and publicly available on Zenodo:

| Dataset | Description | Link |
|---|---|---|
| **AROLit** | Literature-derived primers for RHDV detection | [zenodo.org/records/15112821](https://zenodo.org/records/15112821) |
| **iSOP** | In silico designed primers for RHDV detection | [zenodo.org/records/15113269](https://zenodo.org/records/15113269) |

- `data/AROLit.xlsx` — structured spreadsheet of primers curated from the scientific literature.
- `data/iSOP/` — primer sequences and associated metadata generated via the iSOP in silico design workflow.

---

## Workflows

### `scripts/AROLit_and_iSOP.py`

The main Python script implementing both workflows:

- **AROLit**: automates the retrieval and curation of RHDV primers from published literature into a structured format.
- **iSOP**: designs primers *in silico* from RHDV genomic sequences, applying computational primer design criteria optimised for PCR specificity and sensitivity.

Both workflows are freely accessible, follow **FAIR** principles, and are designed to be adaptable to other viral or bacterial genomes.

### `workflows/Workflow RHDV Orange.ows`

An [Orange Data Mining](https://orangedatamining.com/) workflow file for visual, no-code exploration and analysis of the primer datasets. Open this file directly in Orange Data Mining to interactively explore and compare AROLit and iSOP primer candidates.

---

## Validation

Five top primers from **AROLit** and six from **iSOP** were experimentally validated in the laboratory to evaluate their specificity against several RHDV strains. The best-performing primers demonstrated promising specificity, supporting the combined use of literature-derived and computationally designed primer resources for accurate pathogen surveillance.

---

## How to Use

### Requirements

- Python 3.x
- [Orange Data Mining](https://orangedatamining.com/) (for `.ows` workflow)
- Required Python packages (add your `requirements.txt` or list dependencies here)

### Run the Python workflow

```bash
# Clone the repository
git clone https://github.com/joaomiguelsov/rhdv-primers-identification-db.git
cd rhdv-primers-identification-db

# Run the main script
python scripts/AROLit_and_iSOP.py
```

### Explore with Orange Data Mining

1. Open [Orange Data Mining](https://orangedatamining.com/).
2. Go to **File → Open** and select `workflows/Workflow RHDV Orange.ows`.
3. Explore and compare the AROLit and iSOP primer datasets visually.

---

## Citation

If you use this repository or the associated datasets in your research, please cite:

**Preprint:**

> Carneiro, F., Cardeano, M., Lopes, A.M., Abrantes, J., & Carneiro, J. (2025). *Rabbit Haemorrhagic Disease Virus Database: Enhancing Environmental Surveillance with In Silico and Literature-Derived Primers for PCR Applications.* DOI: [10.21203/rs.3.rs-6364213/v1](https://doi.org/10.21203/rs.3.rs-6364213/v1)

**Datasets:**

> *AROLit – Literature-Derived Primers for RHDV Detection.* Zenodo. https://zenodo.org/records/15112821

> *iSOP – In Silico Designed Primers for RHDV Detection.* Zenodo. https://zenodo.org/records/15113269

---

## Authors

- F. Carneiro
- M. Cardeano
- A.M. Lopes
- J. Abrantes
- J. Carneiro

---

## License

This repository is licensed under the terms described in the [LICENSE](./LICENSE) file.

# RHDV Database

> An open, FAIR-compliant primer resource for **Rabbit Haemorrhagic Disease Virus (RHDV)** detection, integrating literature-derived primers (**AROLit**) and in silico designed primers (**iSOP**) to support PCR-based detection and environmental surveillance.

---

## Table of Contents

- [Background](#background)
- [Datasets](#datasets)
- [Workflows](#workflows)
- [Validation](#validation)
- [Repository Structure](#repository-structure)
- [How to Use](#how-to-use)
- [Citation](#citation)
- [Authors](#authors)
- [License](#license)

---

## Background

Rabbit haemorrhagic disease virus (RHDV) impacts the European rabbit (*Oryctolagus cuniculus*) and was first discovered in China in 1984. Since then, it has spread globally, leading to significant rabbit population declines in both wild and domestic populations.

To manage the growing volume of RHDV genomic data, this project developed two complementary primer resources:

- **AROLit** — an automated computational method to retrieve and curate primers from the scientific literature.
- **iSOP** — a Python-based workflow to design primers *in silico* for RHDV detection.

Both workflows are designed to support optimal primer-pair selection for PCR applications and are applicable beyond RHDV to any viral or bacterial genome.

---

## Datasets

The two core datasets are archived and publicly available on Zenodo:

| Dataset | Description | DOI / Link |
|---|---|---|
| **AROLit** | Literature-derived primers for RHDV detection | [zenodo.org/records/15112821](https://zenodo.org/records/15112821) |
| **iSOP** | In silico designed primers for RHDV detection | [zenodo.org/records/15113269](https://zenodo.org/records/15113269) |

---

## Workflows

### AROLit — Automated Retrieval of Literature Primers

AROLit is an automated pipeline that systematically retrieves primers from published scientific literature related to RHDV detection. It curates and organises primer sequences into a structured, reusable database suitable for direct PCR applications.

### iSOP — In Silico Primer Design

iSOP is a Python-based workflow for designing primers *in silico* against RHDV genomic sequences. It applies computational primer design principles to generate candidate primer pairs optimised for PCR sensitivity and specificity.

Both workflows:

- are freely accessible and adhere to **FAIR** (Findable, Accessible, Interoperable, Reusable) principles;
- produce structured outputs ready for downstream PCR validation;
- are transferable to other viral or bacterial genomes.

---

## Validation

Five top primers from **AROLit** and six from **iSOP** were experimentally validated in the laboratory to evaluate their specificity against several RHDV strains. The best-performing primers demonstrated promising specificity, supporting the combined use of literature-derived and computationally designed primer resources for accurate pathogen surveillance.

---

## Repository Structure

```
.
├── data/
│   ├── arolit/          # AROLit literature-derived primer data
│   └── isop/            # iSOP in silico designed primer data
├── scripts/             # Processing and analysis scripts
├── notebooks/           # Jupyter notebooks for exploration and validation
├── results/             # Processed outputs, tables, and figures
├── docs/                # Additional documentation
└── README.md
```

---

## How to Use

1. **Download the datasets** from the Zenodo records linked above.
2. **Place raw files** in `data/arolit/` and `data/isop/` respectively.
3. **Use the scripts** in `scripts/` or `notebooks/` to parse, filter, rank, or compare primer candidates.
4. **Save processed outputs** to `results/` for downstream PCR assay development.
5. **Cite** both datasets and the associated preprint in any derived work (see [Citation](#citation) below).

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

Please add the appropriate license for this repository. If datasets and code are released under different terms, document each separately in the relevant subdirectory.

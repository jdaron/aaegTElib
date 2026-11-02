# aaegTElib

A transposable element (TE) consensus library for the *Aedes aegypti* AaegL5 reference genome, together with the data and code of the paper *XXX* (Daron et al., 2026).

## Contents

`aaegTElib` is a library of TE consensus sequences identified in the *Aedes aegypti* AaegL5 reference genome. It consists of four files:

| File | Description |
|------|-------------|
| `aaegTElib.fna` | FASTA file of the TE consensus sequences. Each header contains the sequence ID, followed by a description giving the TE classification in three levels: class, order and superfamily. |
| `aaegTElib.classification.txt` | Tab-separated file listing the classification (class, order, superfamily) of each consensus sequence. |
| `aaegTElib.completness.txt` | Tab-separated file indicating whether each consensus sequence is `complete` or `fragmented`. |
| `aaegTElib.LTRcoordinate.txt` | Coordinates of the long terminal repeats (LTRs) within each LTR-containing consensus sequence. |

## Repository structure

The rest of the repository contains the data and code associated with the paper *XXX* (Daron et al., 2026):

| Path | Description |
|------|-------------|
| `figures/` | Input files and R scripts used to generate the figures of the paper. |
| `ressource/` | Resources shared across several analyses of the paper. |
| `scripts/` | Python scripts and bash command lines used to reproduce the analyses of the paper. |
| `TEfam.fna.gz` | Compressed FASTA file of the TEfam TE library, which contains previously published consensus sequences for *Ae. aegypti*. |

## License

- **Code** (`scripts/` and the R scripts in `figures/`): released under the [MIT License](https://opensource.org/licenses/MIT).
- **Data and TE library** (`aaegTElib.*`, `TEfam.fna.gz`, input files in `figures/` and `ressource/`): released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

Please add the corresponding `LICENSE` file (code) and `LICENSE-data` file (data) to the root of the repository.

## Citation

If you use the `aaegTElib` library, the code or the data in this repository, please cite:

> Daron J., et al. (2026). XXX. *JOURNAL*, VOLUME, PAGES. DOI: XXX

```bibtex
@article{daron2026,
  author  = {Daron, Josquin and others},
  title   = {XXX},
  journal = {JOURNAL},
  year    = {2026},
  doi     = {XXX}
}
```

## Contact

For questions or issues, please open a GitHub issue or contact Josquin Daron.


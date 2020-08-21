# eQTL mapping in SPIROMICS bronchial samples for genes related to COVID-19 response: candidate genes for COVID-19

This repository documents the code used to select genes from the Blanco-Melo et al. 2020 study to our COVID-19 candidate genes list.

### __COVID-19 candidate genes__

Candidate genes from different sources:

* Hoffmann et al. 2020, Cell, _SARS-CoV-2 Cell Entry Depends on ACE2 and TMPRSS2 and Is Blocked by a Clinically Proven Protease Inhibitor_, [link](https://doi.org/10.1016/j.cell.2020.02.052)
    + Identified ACE2 as the receptor to be exploited by the SARS-CoV-2 for cellular entry, and proteases TMPRSS2 and cathepsin B/L both to be used by SARS-CoV-2 for S protein priming, whilst only TMPRSS2 is essential for viral entry and viral spread
* Blanco-Melo et al. 2020, Cell, _Imbalanced host response to SARS-CoV-2 drives development of COVID-19_, [link](https://doi.org/10.1016/j.cell.2020.04.026)
    + Explored the transcriptional response to SARS-CoV-2 _in vitro_, _ex vivo_, and _in vivo_ models
* Gordon et al. 2020, Nature, _A SARS-CoV-2 protein interaction map reveals targets for drug repurposing_, [link](https://doi.org/10.1038/s41586-020-2286-9)
    + Identified 332 high-confidence SARS-CoV-2-human protein-protein interactions (PPIs)
* Gassen et al. 2020, bioRxiv, _Analysis of SARS-CoV-2-controlled autophagy reveals spermidine, MK-2206, and niclosamide as putative antiviral therapeutics_, [link](https://doi.org/10.1101/2020.04.15.997254)
    + Showed the role of SARS-CoV-2 infection in restricting AMPK/mTORC1 activation and autophagy
* Wang, K. et al. 2020, bioRxiv, _SARS-CoV-2 invades host cells via a novel route: CD147-spike protein_, [link](https://doi.org/10.1101/2020.03.14.988345)
    + Reported a mediating role of CD147 (also known as BSG) in SARS-CoV-2 viral invasion
* [COVID-19 Cell Atlas](https://www.covid19cellatlas.org)
    + Highlights 17 genes including cathepsins and other viral receptors or receptor associated enzymes

Overview of the COVID-19-related genes

```bash
candidate_genes_final.Rmd
```

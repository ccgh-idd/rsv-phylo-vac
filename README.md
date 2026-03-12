# RSV Vaccination and Evolutionary Dynamics

Does RSV vaccination change the evolutionary rate of the virus? This project uses publicly available RSV-A and RSV-B sequences to test whether the introduction of RSV vaccines (from 2023 onwards) is associated with a detectable shift in viral evolutionary dynamics.

## Approach

Sequences were downloaded from NCBI and filtered to countries with sufficient sampling. For each country, a maximum-likelihood phylogenetic tree was built from amino acid sequences aligned to the fusion (F) protein, and root-to-tip (RTT) regression was used to estimate the evolutionary rate before and after vaccine introduction. A meta-analysis across countries was used to assess whether vaccination was associated with a consistent change in rate.

## Data

- **Sequences**: RSV-A and RSV-B F-protein sequences from NCBI (downloaded March 2026)
- **Metadata**: Collection date, country, subtype
- **Vaccination dates**: National programme start years compiled from regulatory approvals (FDA, EMA, Health Canada, PMDA, etc.) for 21 countries

## Structure

```
R/                  analysis scripts
data/               vaccination dates and BEAST input files
raw_data/           raw sequence downloads
summaries/          metadata summaries
results/            phylogenetic trees and plots (generated)
beast_xml/          BEAST2 XML configs and post-processing
mtbd_rsv/           C++ MTBD model code
```

## Requirements

R (≥ 4.2), with packages: `ape`, `phangorn`, `ggtree`, `tidyverse`, `metafor`, `here`

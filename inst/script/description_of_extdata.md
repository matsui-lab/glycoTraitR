---
editor_options: 
  markdown: 
    wrap: 72
---

## Description of extdata (toy example)

The example files included in the extdata directory were generated from
a subset of the intact N-glycoproteomics dataset reported in
([https://www.mcponline.org/article/S1535-9476(22)00241-9/fulltext](https://www.mcponline.org/article/S1535-9476(22)00241-9/fulltext){.uri})
The complete dataset is publicly accessible under ProteomeXchange
identifier PXD032219, which contains LC–MS/MS measurements of enriched
N-glycoproteins obtained from postmortem human frontal cortex tissues
representing three groups:

-   individuals without cognitive impairment and without pathological AD

-   asymptomatic AD subjects

-   symptomatic AD subjects

For the examples, we used only the HILIC-enriched glycopeptide runs from
PXD032219. The raw LC–MS/MS files were re-analyzed using pGlyco3 and
GlycoDecipher under the default recommended settings for intact
N-glycopeptide identification, and the parameter settings are as
follows.

#### Parameter Setting

GlycoDecipher (v1.0.5) searches were performed using its built-in
N-glycan database containing 10,936 glycan entries.pGlyco3 (v3.1)
searches used the pGlyco3-N-Human glycan library (2,922 human
N-glycans). Both search engines used the same UniProtKB human reference
proteome (UP000005640, downloaded on March 26, 2025) as the protein
sequence database.

#### Complete files

The resulting glycopeptide spectral matches (GPSMs) from the full search
were archived in Zenodo (<https://zenodo.org/records/17759790>).

The original search engine output GPSM file(s) are too big to fit in the
R package. To construct the toy dataset provided in extdata, a small
subset of GPSM files was selected from this Zenodo archive, representing
typical outputs from both search engines. These files were extracted
without modification and formatted in the minimal structure required for
demonstrating the functionality of the glycoTraitR workflow, including:

-   a toy pGlyco3 GPSM matrix (pGlyco3_gpsm_toyexample.txt).
-   a toy GlycoDecipher GPSM folder (decipher_toyexample/).

This curated subset is intended for demonstration purposes. It preserves
the structure and essential characteristics of the original full dataset
while substantially reducing file size and complexity.

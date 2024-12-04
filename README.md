# Translating the Transcriptome: A Connectomics Approach for Gene-Network Mapping and Clinical Application

## Introduction

This repository contains the code and the preprocessed dataset for the preprint **Neudorfer et al. 2024**: *"Translating the Transcriptome: A Connectomics Approach for Gene-Network Mapping and Clinical Application.”* 

This repository can be used to generate gene-network maps, disease-network maps, and replicate the voxel-wise linear modeling analyses that relate gene expression data to functional connectivity patterns. Additionally, the code supports the calculation of correlations to compare disease and lesion network maps and associate clinical outcomes with DBS connectivity profiles.

---

## Abstract of the Paper

Genetic variation profoundly influences the structural and functional architecture of the human brain, yet the mechanisms by which gene expression patterns shape distributed brain networks remain poorly understood. Here, we introduce **gene network mapping**, a novel framework that integrates spatially resolved transcriptomic data with the normative functional connectome of the human brain. This approach enables the identification of functional brain networks associated with specific gene expression patterns, advancing our understanding of how genetic factors contribute to both physiological brain function and neurological disorders.

Leveraging data from the **Allen Human Brain Atlas (N=20,000 genes)** and **resting-state functional MRI** from a large cohort of healthy individuals (N=1,000), we generated gene-network maps that reveal the distributed connectivity patterns linked to individual gene expressions. By aggregating gene-network maps, we constructed disease-network maps that capture the cumulative impact of multiple genes associated with movement disorders such as Parkinson’s disease and dystonia on global brain networks. These maps were validated through comparisons with Lesion Network Mapping and deep brain stimulation (DBS) connectivity profiles, demonstrating a robust association between gene-driven network perturbations and clinical outcomes.

Our findings show that gene network mapping not only recapitulates known pathophysiological networks but also predicts therapeutic outcomes in patients undergoing DBS. This approach offers a powerful tool for dissecting the molecular underpinnings of brain network dysfunction and opens new avenues for precision medicine in neurology and psychiatry. By linking gene expression to brain connectivity, this work lays the foundation for the development of network-based diagnostics and targeted therapeutic interventions.

---

## Installation

### Prerequisites

- **Recommended RAM size**: 32GB or more
- **MATLAB version**: R2022b or later
- The following MATLAB toolboxes:
  - MATLAB Image Processing Toolbox
  - MATLAB Signal Processing Toolbox
  - MATLAB Statistics and Machine Learning Toolbox
  - MATLAB Curve Fitting Toolbox (optional)
  - MATLAB Parallel Computing Toolbox (optional)
- The **SPM12 toolbox**
- **abagen**: a toolbox for the Allen Brain Atlas genetics data
- **1000 GSP healthy subjects Connectome**, accessible via the repository of the Harvard Dataverse: [https://doi.org/10.7910/DVN/KKTJQC](https://doi.org/10.7910/DVN/KKTJQC)
- **Surf Ice**: [https://www.nitrc.org/projects/surfice/](https://www.nitrc.org/projects/surfice/)

### Normal Installation

- **Lead-DBS** can be downloaded from our website in fully functional form.
- **abagen** requires Python 3.6+. Clone the abagen repository from GitHub.

### Development Installation

Alternatively, especially if you wish to modify and contribute to Lead-DBS:

1. Ensure you meet the prerequisites.
2. Clone the Lead-DBS repository from GitHub.
3. Download the necessary data and unzip it into the cloned Git repository.
4. Clone the abagen repository from GitHub.

---

## Getting Started

You can run Lead-DBS by typing `lead dbs` into the MATLAB prompt, which will open up the main GUI.

---

## Minimal Dataset

The minimal dataset **Translating_The_Transcriptome** required to run the repository can be downloaded from Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14279062.svg)](https://doi.org/10.5281/zenodo.14279062)
.

- **Lead Connectome Mapper** was used to generate functional connectivity fingerprint maps. A detailed description and user guide for the Lead Connectome Mapper can be obtained from the [Lead-DBS User Guide](https://netstim.gitbook.io/leaddbs).
  
- To execute `gene_network_mapping.py`:
  - Install the abagen toolbox.
  - Copy the content in the folder `Translating_The_Transcriptome/Gene_Network_Mapping/`:
    - `Atlas_parcellation.nii.gz`: Atlas parcellations used by abagen to assign sampling locations.
    - `Brain_mask.nii.gz`: Brain mask for masking voxels outside brain tissue during gene-network map generation.
    - `Atlas_parcellation_annotation_table.csv`: Atlas annotations required by abagen to assign labels to atlas parcellations.
    - `Translating_The_Transcriptome/Gene_Network_Mapping/Gene_lists/`: List of genes associated with different diseases/symptoms/syndromes as reported in the manuscript.

- **Lead-DBS Network Mapping Explorer** was used to perform spatial correlations, cross-validations, and outcome predictions. Refer to *Neudorfer et al. 2024* for a detailed description. A step-by-step tutorial can be found in the [Lead-DBS User Guide](https://netstim.gitbook.io/leaddbs).

- **Note**: Patient imaging data cannot be publicly shared to comply with current data protection regulations. These are available from the principal investigators of the collecting sites upon reasonable request within the framework of a data-sharing agreement.

For inquiries or further information, please contact:
- **C.N.**: [cneudorfer@mgh.harvard.edu](mailto:cneudorfer@mgh.harvard.edu)
- **A.H.**: [ahorn1@bwh.harvard.edu](mailto:ahorn1@bwh.harvard.edu)

---

## Compatibility

This repository has been tested successfully with:
- **Python**: v3.10.10
- **MATLAB**: versions 2022b and 2023a
- **Operating System**: MacOS

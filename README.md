# Morphological Profiling Dataset of EU-OPENSCREEN Bioactive Compounds Over Multiple Imaging Sites and Cell Lines

Code authors: Carsten Beese and Christopher Schmied

Please cite:

Christopher Wolff, Martin Neuenschwander, Carsten Joern Beese, Divya Sitani, Maria C. Ramos, Alzbeta Srovnalova, Maria Jose Varela, Pavel Polishchuk, Katholiki E. Skopelitou, Ctibor Skuta, Bahne Stechmann, Jose Brea, Mads Hartvig Clausen, Petr Dzubak, Rosario Fernandez-Godino, Olga Genilloud, Marian Hajduch, Maria Isabel Loza, Martin Lehmann, Jens Peter von Kries, Han Sun, Christopher Schmied; Morphological Profiling Dataset of EU-OPENSCREEN Bioactive Compounds Over Multiple Imaging Sites and Cell Lines; bioRxiv 2024.08.27.609964; doi: https://doi.org/10.1101/2024.08.27.609964

# Resources:

Aggregated and processed profiles are hosted on a Zenodo repository: [https://doi.org/10.5281/zenodo.13309566](https://doi.org/10.5281/zenodo.13309566)

The raw image data is hosted on the AWS Cell Painting Gallery under cpg0036-EU-OS-bioactives: [https://cellpainting-gallery.s3.amazonaws.com/index.html#cpg0036-EU-OS-bioactives/](https://cellpainting-gallery.s3.amazonaws.com/index.html#cpg0036-EU-OS-bioactives/)

Information about the compounds: [https://www.probes-drugs.org/compounds/standardized#compoundset=353@AND](https://www.probes-drugs.org/compounds/standardized#compoundset=353@AND)

# Description:

Code shows how to load, process and analyse aggregated Cell Painting data. Each folder contains the analysis for a dataset from one of the four sources (FMP, IMTM, MEDINA, USC). The notebooks were used to create the figure panels of the associated publication.

# Notebooks:

* 1_Collect: Collects aggregated Cell Painting data into a single dataframe.
* 2_Normalization: Processing of profiles with normalization.
* 3_Feature-Selection: Feature selection, QC analysis and computation of consensus profiles.

# QC analysis:

* Number of toxic compounds.
* Number of low active compounds.
* Percent replication.

# Further analysis:

* Comparison between U2OS and HepG2 cell line from the FMP (FMP/4_Comparison_Cell_Lines.ipynb).
* UMAP for visualization of FMP U2OS and HepG2 datasets (UMAP_Viz.ipynb).
* Analysis of Batch effects using UMAPs (Batch_QCViz.ipynb).
* Overall cell numbers and cell numbers per control compound per dataset (CellNumber.ipynb).
* Characterization of Bioactive compounds (Characterize_Bioactive.ipynb).

# Tutorial data access and analysis:

## Access & analysis of Profiles

The profiles are hosted on Zenodo: [https://doi.org/10.5281/zenodo.13309565](https://doi.org/10.5281/zenodo.13309565). For the performed analysis please have a look at our article [https://doi.org/10.1101/2024.08.27.609964](https://doi.org/10.1101/2024.08.27.609964). In brief we extracted the profiles using a Cell Profiler based pipeline. This yields single cell profiles that were then aggregated using a median per well. 

### Aggregated profiles

You can get access to the per well aggregated profiles: [Aggregated_Profiles.zip](https://zenodo.org/records/13309566/files/Aggregated_Profiles.zip?download=1)

Read in aggregated profiles (one example FMP U2OS)
Processing: Normalization

### Normalized profiles

Access to normalized profiles Profile_Analysis_Results.zip

Perform Feature selection

### Profile aggregation and analysis

Access to normalized and reduced profiles
Perform basic analysis > Replication, Induction

## Image data

The image data is hosted on the Amazon Web Services (AWS) Cell Painting gallery ([Weisbart et al. 2024](https://doi.org/10.1038/s41592-024-02399-z)): [AWS Cell Painting Gallery](https://github.com/broadinstitute/cellpainting-gallery)

The dataset name is: cpg0036-EU-OS-bioactives

The dataset can be viewed and navigated here: [https://cellpainting-gallery.s3.amazonaws.com/index.html#cpg0036-EU-OS-bioactives/](https://cellpainting-gallery.s3.amazonaws.com/index.html#cpg0036-EU-OS-bioactives/)

The download process from the Cell Painting Gallery is documented here: [https://broadinstitute.github.io/cellpainting-gallery/download_instructions.html](https://broadinstitute.github.io/cellpainting-gallery/download_instructions.html)

The image data can be downloaded using the Amazon Web Services Command Line Interface ([AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html)). You will first need to install these tools: [https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)

Listing the dataset:

`DATASET=cpg0036-EU-OS-bioactives aws s3 ls s3://cellpainting-gallery/${DATASET}/ --no-sign-request`

### Download of entire dataset.

The entire image dataset is 3.5 TB in size.

Download the data will then be: `aws s3 cp --recursive "CPG_LOCATION" "LOCAL_DESTINATION"`

Thus would look like this: 

`DATASET=cpg0036-EU-OS-bioactives aws s3 cp --recursive s3://cellpainting-gallery/${DATASET}/ FOLDER/TO/LOCAL/ --no-sign-request`

If you just want to test the command before an actual download then use:

`DATASET=cpg0036-EU-OS-bioactives aws s3 cp --recursive s3://cellpainting-gallery/${DATASET}/ FOLDER/TO/LOCAL/ --no-sign-request --dryrun`

For the actual download just remove the flag --dryrun.

### Notes for image data download

Use the --dryrun flag before executing the download to test if the commands and as well as the source and destination locations are correct.

Please review the Cell Painting data structure guide to understand the structure of the data  provided: [https://broadinstitute.github.io/cellpainting-gallery/data_structure.html](https://broadinstitute.github.io/cellpainting-gallery/data_structure.html)

You do not need an AWS account for download of the files. If you get and error with the AWS CLI command add --no-sign-request to the end of the command.

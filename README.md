# Morphological Profiling Dataset of EU-OPENSCREEN Bioactive Compounds Over Multiple Imaging Sites and Cell Lines

Code authors: Carsten Beese and Christopher Schmied

Please cite: 

Christopher Wolff, Martin Neuenschwander, Carsten Joern Beese, Divya Sitani, Maria C. Ramos, Alzbeta Srovnalova, Maria Jose Varela, Pavel Polishchuk, Katholiki E. Skopelitou, Ctibor Skuta, Bahne Stechmann, Jose Brea, Mads Hartvig Clausen, Petr Dzubak, Rosario Fernandez-Godino, Olga Genilloud, Marian Hajduch, Maria Isabel Loza, Martin Lehmann, Jens Peter von Kries, Han Sun, Christopher Schmied; Morphological Profiling Dataset of EU-OPENSCREEN Bioactive Compounds Over Multiple Imaging Sites and Cell Lines; bioRxiv 2024.08.27.609964; doi: https://doi.org/10.1101/2024.08.27.609964

# Resources:

Aggregated and processed profiles are hosted on a Zenodo repository: [https://doi.org/10.5281/zenodo.13309566](https://doi.org/10.5281/zenodo.13309566)

The raw image data is hosted on the AWS Cell Painting Gallery under cpg0036-EU-OS-bioactives: [https://cellpainting-gallery.s3.amazonaws.com/index.html#cpg0036-EU-OS-bioactives/](https://cellpainting-gallery.s3.amazonaws.com/index.html#cpg0036-EU-OS-bioactives/)

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

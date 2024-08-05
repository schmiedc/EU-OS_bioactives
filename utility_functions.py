### functions for processing and analysis

### imports
import pandas as pd # type: ignore
import pycytominer # type: ignore
import numpy as np
import random
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore

##### get feature vector
def get_feature_vector(df):
    ######
    ### get metadata columns and define feature columns after deleting "0" Features

    Feature_Ident = ["Nuc_", "Cells_", "Cyto_"]
    feat = []
    
    for col in df.columns.tolist():
        if any([col.startswith(x) for x in Feature_Ident]):
                feat.append(col)
    return feat

#### get metadata vector
def get_metadata_vector(df):
    ######
    ### get metadata columns 

    metadata_Ident = ["Metadata"]
    meta = []
    
    for col in df.columns.tolist():
        if any([col.startswith(x) for x in metadata_Ident]):
                meta.append(col)
    return meta

###
##### Feature reduction
###

def feature_reduction(Data, variance_freq_cut=0.1, variance_unique_cut=0.1, outlier_cutoff=500, corr_threshold = 0.8, print_stats = True):
    """
    Reduces the number of featues by removing NaN columns, low variance, outlier and correlating features
    :param Data: DataFrame of well-averaged Profiles
    :param variance_freq_cut: frequencey cut off for the varience filter, defaults to 0.1
    :param variance_unique_cut: unique cut argument for the varience filter, defaults to 0.1
    :param outlier_cutoff: the threshold for the maximum or minimum value of a feature, defaults to 500
    :param corr_threshold: threshold for pearson correlation to determine correlating features, defaults to 0.8
    :param print_stats: boolean, if True prints statistics for each step
    :return: DataFarame with reduced features
    """

    ### 1) Remove NaN columns
    Data_DropNA = Data.drop(columns = Data.columns[(Data.isnull()).any()])

    ### 2) Remove Columns with low varience and thus low information

    Features = get_feature_vector(Data_DropNA)
        
    Data_VarianceThresh = pycytominer.feature_select(
            profiles = Data_DropNA, 
            features = Features, 
            samples='all', 
            operation='variance_threshold', 
            freq_cut=variance_freq_cut, # 2nd most common feature val / most common [default: 0.1]
            unique_cut=variance_unique_cut # float of ratio (num unique features / num samples) [default: 0.1] 
    )
    
    ### 3) Remove Columns with Outliers
    Features = get_feature_vector(Data_VarianceThresh)
     
    Data_dropOutliers = pycytominer.feature_select(
            profiles = Data_VarianceThresh, 
            features = Features, 
            samples='all', 
            operation='drop_outliers', 
            outlier_cutoff=outlier_cutoff # [default: 15] the threshold at which the maximum or minimum value of a feature
    )

    ### 4) Remove correlating features
    Features = get_feature_vector(Data_dropOutliers)
    Meta_Features = list(set(Data_dropOutliers.columns) - set(Features))
    
    Data_Reduced = pycytominer.feature_select(
            profiles = Data_dropOutliers, 
            features = Features, 
            samples='all', 
            operation='correlation_threshold',
            corr_threshold = corr_threshold, 
            corr_method='pearson'
    )
    
    if print_stats == True:
        Features_df = Feature_Vis(get_feature_vector(Data), "Original Features")
        Features_df = Features_df.merge(Feature_Vis(get_feature_vector(Data_VarianceThresh), "Variance Threshold"), on = "Category")
        Features_df["% Variance"] = round((100/Features_df["Original Features"])*Features_df["Variance Threshold"],1)
        Features_df = Features_df.merge(Feature_Vis(get_feature_vector(Data_dropOutliers), "Outlier Threshold"), on = "Category")
        Features_df["% Outlier"] = round((100/Features_df["Original Features"])*Features_df["Outlier Threshold"],1)
        Features_df = Features_df.merge(Feature_Vis(get_feature_vector(Data_Reduced), "Correlation Threshold"), on = "Category")
        Features_df["% Correlation"] = round((100/Features_df["Original Features"])*Features_df["Correlation Threshold"],1)
        print(Features_df.to_markdown(index=False))
        print()
        
        Features_df = Feature_Vis_Compartment(get_feature_vector(Data), "Original Features")
        Features_df = Features_df.merge(Feature_Vis_Compartment(get_feature_vector(Data_VarianceThresh), "Variance Threshold"), on = "Category")
        Features_df["% Variance"] = round((100/Features_df["Original Features"])*Features_df["Variance Threshold"],1)
        Features_df = Features_df.merge(Feature_Vis_Compartment(get_feature_vector(Data_dropOutliers), "Outlier Threshold"), on = "Category")
        Features_df["% Outlier"] = round((100/Features_df["Original Features"])*Features_df["Outlier Threshold"],1)
        Features_df = Features_df.merge(Feature_Vis_Compartment(get_feature_vector(Data_Reduced), "Correlation Threshold"), on = "Category")
        Features_df["% Correlation"] = round((100/Features_df["Original Features"])*Features_df["Correlation Threshold"],1)
        print(Features_df.to_markdown(index=False))
         
       
    return Data_Reduced

        
###
### Feature Reduction
### Which features are lost

def Feature_Vis(Feature_list, Ident):
    import re

    ## only intensity features
    r = re.compile(".*Intensity.*")
    Intensity_Features = list(filter(r.match, Feature_list))

    ## only Correlation features
    r = re.compile("Nuc_Correlation|Cells_Correlation|Cyto_Correlation.*")
    Correlation_Features = list(filter(r.match, Feature_list))

    ## only AreaShape features
    r = re.compile(".*AreaShape.*")
    AreaShape_Features = list(filter(r.match, Feature_list))

    ## only Granularity features
    r = re.compile(".*Granularity.*")
    Granularity_Features = list(filter(r.match, Feature_list))

    ## only Neighbors features
    r = re.compile(".*Neighbors.*")
    Neighbors_Features = list(filter(r.match, Feature_list))

    ## only RadialDistribution features
    r = re.compile(".*RadialDistribution.*")
    RadialDistribution_Features = list(filter(r.match, Feature_list))

    ## only Texture features
    r = re.compile(".*Texture.*")
    Texture_Features = list(filter(r.match, Feature_list))

    ## mito skeleton
    r = re.compile(".*ObjectSkeleton.*")
    Skeleton_Features = list(filter(r.match, Feature_list))

    Found_features = []
    Found_features.extend(Intensity_Features)
    Found_features.extend(Correlation_Features)
    Found_features.extend(AreaShape_Features)
    Found_features.extend(Granularity_Features)
    Found_features.extend(Neighbors_Features)
    Found_features.extend(RadialDistribution_Features)
    Found_features.extend(Texture_Features)
    Found_features.extend(Skeleton_Features)

    
    df_features = pd.DataFrame()
    df_features["Category"] = ["Total Features",
                               "Intensity",
                               "Correlation",
                               "AreaShape",
                               "Granularity",
                               "Neighbors",
                               "RadialDistribution",
                               "Texture",
                               'MitoSkeleton'] 
    df_features[Ident] = [len(Feature_list),
                          len(Intensity_Features),
                          len(Correlation_Features),
                          len(AreaShape_Features),
                          len(Granularity_Features),
                          len(Neighbors_Features),
                          len(RadialDistribution_Features),
                          len(Texture_Features),
                          len(Skeleton_Features)]
   
    return df_features


    ### Feature visualization by compartment
def Feature_Vis_Compartment(Feature_list, Ident):
    import re

    ## only nuclear features
    r = re.compile(".*Nuc.*")
    Nucleus_Features = list(filter(r.match, Feature_list))

    ## only Cell features
    r = re.compile(".*Cell.*")
    Cell_Features = list(filter(r.match, Feature_list))

    ## only Cytoplasm features
    r = re.compile(".*Cyto.*")
    Cyto_Features = list(filter(r.match, Feature_list))

    
    Found_features = []
    Found_features.extend(Nucleus_Features)
    Found_features.extend(Cell_Features)
    Found_features.extend(Cyto_Features)
    
    df_features = pd.DataFrame()
    df_features["Category"] = ["Total Features",
                               "Nucleus",
                               "Cell",
                               "Cytoplasm"] 
    df_features[Ident] = [len(Feature_list),
                          len(Nucleus_Features),
                          len(Cell_Features),
                          len(Cyto_Features)]
   
    return df_features 

###
##### removing toxic conditions
###
### we will remove all wells where the cell count is < median - 2.5 std of the population

def remove_tox(Data, key_col = ["Metadata_EOS", "Metadata_Plate", "Metadata_Concentration"], SD_Threshold = 2.5,  plot_distribution = True):
    """
    removes toxic conditions from aggregated CellProfiler Profiles
    :param Data: aggregated CellProfiler Profiles
    :param key_col: list of column name used as key to identify treatments
    :param plot_distribution: boolean, if True histograms of cell count distribution will be plotted
    :return: new DataFrame without toxic conditions, set of toxic conditions
    """
    Features = key_col.copy()
    Features.append("Metadata_Object_Count")

    Median_CellCount = pycytominer.consensus(
        profiles = Data.loc[:,Features], # A file or pandas DataFrame of profile data
        replicate_columns = key_col, # Metadata columns indicating which replicates to collapse, defaults to [“Metadata_Plate”, “Metadata_Well”]
        operation = "median", # (str) – The method used to form consensus profiles, defaults to “median”
        features = ["Metadata_Object_Count"], # (str, list) – The features to collapse, defaults to “infer”
    )

    tox_threshold = Median_CellCount["Metadata_Object_Count"].median() -(SD_Threshold*Median_CellCount["Metadata_Object_Count"].std())
    
    if tox_threshold < 50:
        tox_threshold = 50
    
    Tox_cond = Median_CellCount.loc[Median_CellCount["Metadata_Object_Count"] < tox_threshold, key_col]

    i1 = Data.set_index(key_col).index
    i2 = Tox_cond.set_index(key_col).index
    Data_tox = Data[~i1.isin(i2)]

    Data_tox = Data_tox.reset_index(drop = True)
    
    print("Toxic conditions removed with threshold", round(tox_threshold, 2))
    print("Old shape", Data.shape)
    print("New shape", Data_tox.shape)
    
    if plot_distribution == True:
        fig, axs = plt.subplots(ncols=2)

        sns.histplot(
            ax=axs[0],
            data = Data, 
            x = "Metadata_Object_Count",
            #hue= "Metadata_Plate"
        )
        axs[0].set_xlim([0, None])
        
        sns.histplot(
            ax=axs[1],
            data = Data_tox, 
            x = "Metadata_Object_Count",
            #hue= "Metadata_Plate"
        )
        axs[1].set_xlim([0, None])
        
        fig.suptitle('Cell Count')
        axs[0].set_title('Before')
        axs[1].set_title('After')
        
    
    return Data_tox, Tox_cond


###
##### Reproducibility
###

def corr_between_replicates_CP(df, group_by_feature):
    """
        Correlation between replicates
        Parameters:
        -----------
        df: pd.DataFrame
        group_by_feature: Feature name to group the data frame by
        Returns:
        --------
        list-like of correlation values
     """
    replicate_corr = []
    names = []
    replicate_grouped = df.groupby(group_by_feature)
    for name, group in replicate_grouped:
        group_features = group.loc[:, get_feature_vector(group)]
        corr = np.corrcoef(group_features)
        if len(group_features) == 1:  # If there is only one replicate on a plate
            replicate_corr.append(np.nan)
            names.append(name)
        else:
            np.fill_diagonal(corr, np.nan)
            replicate_corr.append(np.nanmedian(corr))  # median replicate correlation
            names.append(name)
    return replicate_corr, names
        

def corr_between_non_replicates_CP(df, n_samples, n_replicates, metadata_compound_name):
    """
        Null distribution between random "replicates".
        Parameters:
        ------------
        df: pandas.DataFrame
        n_samples: int
        n_replicates: int
        metadata_compound_name: Compound name feature
        Returns:
        --------
        list-like of correlation values, with a  length of `n_samples`
    """
    df.reset_index(drop=True, inplace=True)
    null_corr = []
    while len(null_corr) < n_samples:
        compounds = random.choices([_ for _ in range(len(df))], k = n_replicates)
        sample = df.loc[compounds].copy()
        if len(sample[metadata_compound_name].unique()) == n_replicates:
            sample_features = sample.loc[:, get_feature_vector(sample)]
            corr = np.corrcoef(sample_features)
            np.fill_diagonal(corr, np.nan)
            null_corr.append(np.nanmedian(corr))  # median non-replicate correlation
    return null_corr

def percent_score(null_dist, corr_dist, how):
    """
    Calculates the Percent strong or percent recall scores
    :param null_dist: Null distribution
    :param corr_dist: Correlation distribution
    :param how: "left", "right" or "both" for using the 5th percentile, 95th percentile or both thresholds
    :return: proportion of correlation distribution beyond the threshold
    """
    if how == 'right':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        return np.mean(above_threshold.astype(float))*100, perc_95
    if how == 'left':
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return np.mean(below_threshold.astype(float))*100, perc_5
    if how == 'both':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return (np.mean(above_threshold.astype(float)) + np.mean(below_threshold.astype(float)))*100, perc_95, perc_5


##### Roproducibility score

def remove_non_reproducible(Data, n_samples = 100, n_replicates = 4, ID_col = "Metadata_Gene_ID", cntrls = ["DMSO", "Nocodazole", "Tetrandrine"], description = "Data"):

    corr_replicating_df = pd.DataFrame()

    replicating_corr, names = list(corr_between_replicates_CP(Data, ID_col))
    null_replicating = list(corr_between_non_replicates_CP(Data, n_samples=n_samples, n_replicates=n_replicates, metadata_compound_name = ID_col))

    prop_95_replicating, value_95_replicating = percent_score(null_replicating, replicating_corr, how='right')
    
    ### this only works well for big data sets with bigger 3 repetitions. Otherwise the null distribution is very wide! Thus we set a threshold for these cases
    
    if value_95_replicating > 0.6:
        value_95_replicating = np.float64(0.6)
        above_threshold = replicating_corr > value_95_replicating
        prop_95_replicating = np.mean(above_threshold.astype(float))*100


    corr_replicating_df = corr_replicating_df.append({'Description': description,
                                                       # 'Modality':f'{modality}',
                                                        #'Cell':f'{cell}',
                                                        #'time':f'{time}',
                                                        'Replicating':replicating_corr,
                                                        'Null_Replicating':null_replicating,
                                                        'Percent_Replicating':'%.1f'%prop_95_replicating,
                                                        'Value_95':value_95_replicating}, ignore_index=True)


    print(corr_replicating_df[['Description','Percent_Replicating']].to_markdown(index=False))
    
    ### remove non replicating conditions
    df = pd.DataFrame(list(zip(names, replicating_corr)),
               columns =['Name', 'replicating_corr'])
    
    replicating = list(df.loc[df["replicating_corr"] > value_95_replicating, "Name"])
    replicating.extend(cntrls)
    replicating = list(set(replicating))

    Data_replicating = Data.loc[Data[ID_col].isin(replicating)].copy()
    
    print("Nonreplicating conditions removed with threshold", round(value_95_replicating, 2))
    if len(set(names) - set(replicating)) == 0:
        print("No conditions below threshold")
    #else:
        #print(set(names) - set(replicating))
    print("Old shape", Data.shape)
    print("New shape", Data_replicating.shape)
    
    return Data_replicating, corr_replicating_df


# Induction
# Calculated directly from the fingerprints by a feature-by-feature comparison of Z scores.
# From Christoforow et al., 2019; Foley et al., 2020; Laraia et al., 2020; Schneidewind et al., 2020; Zimmermann et al., 2019.  
# A significant change was defined as a deviation of more than three times the MAD from the median of the DMSO controls. 
# The induction value is then determined for every compound as the fraction of significantly changed features (as a percentage). 
# An induction of 5% or higher was considered a valid indication that the morphological change produced by the compound is meaningful. 

def remove_low_active(df: pd.DataFrame, 
                      key_col = ["Metadata_EOS", "Metadata_Plate", "Metadata_Concentration", "Metadata_Partner"], 
                      feature_activity_threshold=3.0, 
                      induction_threshold=5):
    """
    removes compounds via induction threshold
    :param df: consensus, feature selected CellProfiler Profiles
    :param key_col = ["Metadata_EOS", "Metadata_Plate", "Metadata_Concentration"]
    :param feature_activity_threshold: z-score where feature is considered active
    :param induction_threshold: % of active features where compound passes threshold
    :return: new DataFrame with active compounds, new DataFrame with non active compounds
    """
    # removes key columns
    feature_df = df.drop(columns=key_col)
 
    # percent of features equal or higher than activity threshold
    induction = (feature_df >= feature_activity_threshold).sum(axis=1) / len(feature_df.columns) * 100

    # treatments with induction equal or higher than induction threshold
    Data_active = df[(induction >= induction_threshold)]
    
    # treatments with induction lower than induction threshold
    Data_non_active = df[(induction < induction_threshold)]
    
    return Data_active, Data_non_active

# %Pairing
# Equivalent to % Matching 
# The signal distribution, which is the median pairwise correlation between each cell pair
# The null distribution, which is the median pairwise correlation of Compound-CRISPRs or Compound-ORF that target different genes, is computed for 1000 combinations of Compound-CRISPRs or Compound-ORFs.
# Percent Matching is computed as the percentage of the signal distribution that is the greater than the 95th percentile of null distribution
# The signal and noise distributions and the Percent Matching values are plotted and the table of Percent Matching is printed.
# Copied from: https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted/blob/main/benchmark/old_notebooks/2.percent_matching_across_modalities.ipynb
# Functions from: https://github.com/jump-cellpainting/2021_Chandrasekaran_submitted/blob/main/benchmark/old_notebooks/utils.py

def get_featuredata(df):
    """return dataframe of just featuredata columns"""
    return df[get_featurecols(df)]

def get_featurecols(df):
    """returna  list of featuredata columns"""
    return [c for c in df.columns if not c.startswith("Metadata")]


def correlation_between_modalities(modality_1_df, 
                                   modality_2_df, 
                                   modality_1, 
                                   modality_2, 
                                   metadata_common, 
                                   metadata_perturbation):
    """
    Compute the correlation between two different modalities.
    :param modality_1_df: Profiles of the first modality
    :param modality_2_df: Profiles of the second modality
    :param modality_1: "Compound", "ORF" or "CRISPR"
    :param modality_2: "Compound", "ORF" or "CRISPR"
    :param metadata_common: feature that identifies perturbation pairs
    :param metadata_perturbation: perturbation name feature
    :return: list-like of correlation values
    """
    list_common_perturbation_groups = list(np.intersect1d(list(modality_1_df[metadata_common]), 
                                                          list(modality_2_df[metadata_common])))

    merged_df = pd.concat([modality_1_df, modality_2_df], ignore_index=False, join='inner')

    modality_1_df = merged_df.query('Metadata_Cell_type==@modality_1')
    modality_2_df = merged_df.query('Metadata_Cell_type==@modality_2')

    corr_modalities = []

    for group in list_common_perturbation_groups:
        modality_1_perturbation_df = modality_1_df.loc[modality_1_df[metadata_common] == group]
        modality_2_perturbation_df = modality_2_df.loc[modality_2_df[metadata_common] == group]

        for sample_1 in modality_1_perturbation_df[metadata_perturbation].unique():
            for sample_2 in modality_2_perturbation_df[metadata_perturbation].unique():
                modality_1_perturbation_sample_df = modality_1_perturbation_df.loc[modality_1_perturbation_df[metadata_perturbation] == sample_1]
                modality_2_perturbation_sample_df = modality_2_perturbation_df.loc[modality_2_perturbation_df[metadata_perturbation] == sample_2]

                modality_1_perturbation_profiles = get_featuredata(modality_1_perturbation_sample_df)
                modality_2_perturbation_profiles = get_featuredata(modality_2_perturbation_sample_df)

                corr = np.corrcoef(modality_1_perturbation_profiles, modality_2_perturbation_profiles)
                corr = corr[0:len(modality_1_perturbation_profiles), len(modality_1_perturbation_profiles):]
                corr_modalities.append(np.nanmedian(corr))  # median replicate correlation

    return corr_modalities

def null_correlation_between_modalities(modality_1_df, modality_2_df, modality_1, modality_2, metadata_common, metadata_perturbation, n_samples):
    """
    Compute the correlation between two different modalities.
    :param modality_1_df: Profiles of the first modality
    :param modality_2_df: Profiles of the second modality
    :param modality_1: "Compound", "ORF" or "CRISPR"
    :param modality_2: "Compound", "ORF" or "CRISPR"
    :param metadata_common: feature that identifies perturbation pairs
    :param metadata_perturbation: perturbation name feature
    :param n_samples: int
    :return:
    """
    list_common_perturbation_groups = list(np.intersect1d(list(modality_1_df[metadata_common]), list(modality_2_df[metadata_common])))

    merged_df = pd.concat([modality_1_df, modality_2_df], ignore_index=False, join='inner')

    modality_1_df = merged_df.query('Metadata_Cell_type==@modality_1')
    modality_2_df = merged_df.query('Metadata_Cell_type==@modality_2')

    null_modalities = []

    while len(null_modalities) < n_samples:
        perturbations = random.choices(list_common_perturbation_groups, k=2)
        modality_1_perturbation_df = modality_1_df.loc[modality_1_df[metadata_common] == perturbations[0]]
        modality_2_perturbation_df = modality_2_df.loc[modality_2_df[metadata_common] == perturbations[1]]

        for sample_1 in modality_1_perturbation_df[metadata_perturbation].unique():
            for sample_2 in modality_2_perturbation_df[metadata_perturbation].unique():
                modality_1_perturbation_sample_df = modality_1_perturbation_df.loc[modality_1_perturbation_df[metadata_perturbation] == sample_1]
                modality_2_perturbation_sample_df = modality_2_perturbation_df.loc[modality_2_perturbation_df[metadata_perturbation] == sample_2]

                modality_1_perturbation_profiles = get_featuredata(modality_1_perturbation_sample_df)
                modality_2_perturbation_profiles = get_featuredata(modality_2_perturbation_sample_df)

                corr = np.corrcoef(modality_1_perturbation_profiles, modality_2_perturbation_profiles)
                corr = corr[0:len(modality_1_perturbation_profiles), len(modality_1_perturbation_profiles):]
                null_modalities.append(np.nanmedian(corr))  # median replicate correlation

    return null_modalities


def percent_score(null_dist, corr_dist, how):
    """
    Calculates the Percent strong or percent recall scores
    :param null_dist: Null distribution
    :param corr_dist: Correlation distribution
    :param how: "left", "right" or "both" for using the 5th percentile, 95th percentile or both thresholds
    :return: proportion of correlation distribution beyond the threshold
    """
    if how == 'right':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        return np.mean(above_threshold.astype(float))*100, perc_95
    if how == 'left':
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return np.mean(below_threshold.astype(float))*100, perc_5
    if how == 'both':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return (np.mean(above_threshold.astype(float)) + np.mean(below_threshold.astype(float)))*100, perc_95, perc_5
    
def distribution_plot(df, metric):
    """
    Generates the correlation distribution plots
    Parameters:
    -----------
    df: pandas.DataFrame
        dataframe containing the data points of replicate and null correlation distributions, description, Percent score and threshold values.
    output_file: str
        name of the output file. The file will be output to the figures/ folder.
    metric: str
        Percent Pairing
    Returns:
    -------
    None
    """

    if metric == 'Percent Replicating':
        metric_col = 'Percent_Replicating'
        null = 'Null_Replicating'
        null_label = 'non-replicates'
        signal = 'Replicating'
        signal_label = 'replicates'
        x_label = 'Replicate correlation'
    elif metric == 'Percent Matching':
        metric_col = 'Percent_Matching'
        null = 'Null_Matching'
        null_label = 'non-matching perturbations'
        signal = 'Matching'
        signal_label = 'matching perturbations'
        x_label = 'Correlation between Compounds'

    n_experiments = len(df)

    plt.rcParams['figure.facecolor'] = 'white'  # Enabling this makes the figure axes and labels visible in PyCharm Dracula theme
    plt.figure(figsize=[12, n_experiments * 6])

    for i in range(n_experiments):
        plt.subplot(n_experiments, 1, i + 1)
        plt.hist(df.loc[i, f'{null}'], label=f'{null_label}', density=True, bins=20, alpha=0.5)
        plt.hist(df.loc[i, f'{signal}'], label=f'{signal_label}', density=True, bins=20, alpha=0.5)
        plt.axvline(df.loc[i, 'Value_95'], label='95% threshold')
        plt.legend(fontsize=20)
        plt.title(
            f"{df.loc[i, 'Description']}\n" +
            f"{metric} = {df.loc[i, f'{metric_col}']}",
            fontsize=25
        )
        plt.ylabel("density", fontsize=25)
        plt.xlabel(f"{x_label}", fontsize=25)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        sns.despine()
    plt.tight_layout()


###
##### Visualizations
###
##### UMAP

def UMAP_proj(Data, dim = 2):
    """
    UMAP projection of profiles
    :param Data: DataFrame of Profiles
    :param dim: Dimensions of the UMAP dimensionality reduction
    :return: DataFarame UMAP axis and original Metadata
    """
    from sklearn.preprocessing import StandardScaler # type: ignore
    from umap import UMAP # type: ignore
    
    Data = Data.reset_index(drop = True)
    Features = get_feature_vector(Data)
    Meta_Features = list(set(Data.columns) - set(Features))
    
    # Normalization
    x = Data.loc[:, Features]
    x = StandardScaler().fit_transform(x) # normalizing the features. Actually the data already is normalized by pyCytoMiner!!!!!
    x.shape

    #Let's convert the normalized features into a tabular format with the help of DataFrame.
    feat_cols = ['feature'+str(i) for i in range(x.shape[1])]
    normalised_Features = pd.DataFrame(x,columns=feat_cols)

    umap = UMAP(n_components=dim, init='random', random_state=0)

    proj = umap.fit_transform(normalised_Features)

    ## Next, let's create a DataFrame that will have the principal component values
    if dim == 2:
        proj_df = pd.DataFrame(data = proj,
                               columns = ['Axis 1', 'Axis 2'])
        
    if dim == 3:
        proj_df = pd.DataFrame(data = proj,
                               columns = ['Axis 1', 'Axis 2', 'Axis 3'])
        
        
    UMAP_Data = pd.merge(Data.loc[:, Meta_Features], proj_df, left_index = True, right_index = True)
    print("UMAP projection performed")
    
    return UMAP_Data

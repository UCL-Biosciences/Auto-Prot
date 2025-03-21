### Read in and pre-process data ###

## Note this should be two files: protein abundance and metadata

import pandas as pd
import numpy as np
import logging
import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import sklearn
import sklearn.ensemble
from sklearn.experimental import enable_iterative_imputer
from sklearn.preprocessing import StandardScaler

# specify location of errors to standard output
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')

# to write to log file, add e.g. filename='data_processing.log'. filemode = 'a' would append instead of overwriting

def load_data(file_path):
    """
    Load data from a CSV or TSV file into a pandas DataFrame.

    Parameters:
    - file_path (str): Path to the file.

    Returns:
    - pd.DataFrame: Loaded data.
    """
    try:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        elif file_path.endswith('.tsv'):
            df = pd.read_csv(file_path, sep='\t')
        else:
            raise ValueError("Unsupported file format. Please use CSV or TSV.")
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def normalise_column_names(df, file_path=None, metadata = None):
    """
    Standardise column names by making them lowercase and replacing spaces with underscores.
    Note, need to do this in metadata cols too - see preprocess_data

    Parameters:
    - df (pd.DataFrame): The dataframe to modify.

    Returns:
    - pd.DataFrame: Dataframe with normalized column names.
    """
    df.columns = [col.lower().replace(' ', '_') for col in df.columns]
    df.columns = df.columns.astype(str)
    ### For phosphoproteomic data, there are abundances for phosphorylated proteins
    ### Each protein can be present multiple times - once per phosphorylation state
    ### For the analysis to proceed, we need a unique ID for the protein-phosphorylation state combination
    ## we append the phosporylation state (in column PTM.ModificationTitle and PTM.SiteAA) to the pg.genes column
    if ('ptm.modificationtitle' in df) and ('ptm.siteaa' in df) and ('ptm.sitelocation' in df) and ('pg.genes' in df):
        print("Gene names and phosphorylation state present. Combining to make unique gene names")
        df["pg.genes"] = df["pg.genes"] + '__' + df["ptm.modificationtitle"] + '_' + df["ptm.siteaa"] + '_' + df["ptm.sitelocation"].astype(str)
    # If 'proteindata' is in the file path, set a column containing 'genes' as the index
    if 'proteindata' in file_path:
        genes_columns = [col for col in df.columns if 'genes' in col.lower()]
        if genes_columns:  # If any column contains 'genes'
            df = df.set_index(genes_columns[0])
    return df

def clean_data(df, file_path = None, metadata = None, outPath = None):
    """
    Perform basic cleaning on the data: first filter protein data to include only protein abundance cols 
    and make sure they are all present.
    Non-numeric protein abundances are converted to NaNs, then duplicate rows or rows with NAs are removed.

    Parameters:
    - df (pd.DataFrame): The dataframe to clean.
    - metadata (pd.DataFrame): Optional metadata dataframe containing 'protein_abundance_name'.

    Returns:
    - pd.DataFrame: Cleaned dataframe.
    
    Raises:
    - ValueError: If metadata is not provided or does not contain 'protein_abundance_name'.

    """
    if 'metadata' in file_path:
        df['protein_abundance_name'] = df['protein_abundance_name'].str.lower().str.replace(' ', '_')
        df['sample_rep'] = (df['sample_id'] + '_' + df['replicate'].astype(str) ) # for a unique ID for each sample_replicate
        
        ### for plotting, assign a colour to each unique treatment group
        colours = sns.color_palette("colorblind", len(df['treatment'].unique())).as_hex() # list of colour-blind-friendly colours as long as number of unique treatment levels
        # colours = px.colors.qualitative.Alphabet[:len(df['treatment'].unique())] 
        col_map = dict(zip(df['treatment'].unique(), colours)) # zip (combine) the treatments to the colours, store as a dictionary
        df[ 'colours' ] = df['treatment'].map(col_map) # create a new var, mapping treatment in metadata to the colours in the dictionary
    if 'proteindata' in file_path:
        if metadata is None:
            raise ValueError("Error: Metadata is required but not provided.")   
        # Check if the required column is in metadata
        elif 'protein_abundance_name' not in metadata.columns:
            raise ValueError("Error: 'protein_abundance_name' column is missing in the metadata.")
        else:
            # Filter columns based on metadata['protein_abundance_name']
            # extract column names to keep
            valid_columns = metadata['protein_abundance_name'].tolist()
            # filter df to keep prot abundance columns
            df = df.loc[:, df.columns.isin(valid_columns)]
            # protein columns should have only numeric data
            # convert non-numeric values to NaN and print warning message
            if not df.equals(df.select_dtypes(include=[np.number])):
                print("Warning: DataFrame contains non-numeric values! Converting to NaN: these proteins will be imputed!")
                df = df.apply(pd.to_numeric, errors="coerce")
            ### protein columns can have long names - better to have just sample name
            # Create a mapping of old column names to new column names
            rename_mapping = dict(zip(metadata['protein_abundance_name'], metadata['sample_rep'] ) ) 
            # Rename columns in data_in based on the mapping
            df = df.rename(columns=rename_mapping)
            nrow_original = len(df.index)
#####
#####
            #### impute and normalise from AlphaStats functions #####
            ### log2 all vals
            df_log2 = np.log2(df)
            # Note It has been shown that normalizing the data first and then imputing the data performs better, than the other way around (Karpievitch et al. 2012). 
            # This preprocessing order is also acquired in AlphaStats (unless preprocessing is done in several steps).
            #### normalise
            minmax = sklearn.preprocessing.MinMaxScaler()
            scaler = sklearn.preprocessing.PowerTransformer()
            minmaxed_array = minmax.fit_transform(df_log2.values.transpose())
            normalized_array = scaler.fit_transform(minmaxed_array).transpose()
            df_norm = pd.DataFrame(normalized_array, index=df_log2.index, columns=df_log2.columns )
            #### impute ####
            ## using random forest approach, which the AlphStats paper says performs best
            imp = sklearn.ensemble.HistGradientBoostingRegressor(max_depth=10, max_iter=100, random_state=0)
            imp = sklearn.impute.IterativeImputer(random_state=0, estimator=imp)
            imputation_array = imp.fit_transform(df_norm.values)
            df_imp = pd.DataFrame(
                        imputation_array, index=df.index, columns=df.columns
                    )
            ### view distributions after different steps
            # List of dataframes and titles for the subplots
            dfs = [df, df_log2, df_norm, df_imp]
            titles = ['Raw Intensities', 'Log₂ Transformed', 'Normalised', 'Imputed']
            # Get unique treatment groups from metadata
            treatment_groups = metadata['treatment'].unique()
            # Set up the figure with a 2x2 grid
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            axes = axes.flatten()
            for ax, data, title in zip(axes, dfs, titles):
                # Convert the DataFrame to long format: one row per measurement with sample id and intensity.
                long_df = data.melt(var_name='sample_rep', value_name='intensity')
                # Merge with metadata to obtain treatment information for each sample.
                # Make sure metadata has 'sample_id' and 'treatment' columns.
                long_df = long_df.merge(metadata[['sample_rep', 'treatment']], on='sample_rep', how='left')
                # Plot boxplots: each sample's distribution is shown on the x-axis.
                # The boxes are colored by treatment.
                sns.boxplot(x='sample_rep', y='intensity', data=long_df, hue='treatment', ax=ax)
                ax.set_title(title)
                ax.set_xlabel("Sample ID")
                ax.set_ylabel("Intensity")
                # Rotate x tick labels for better readability
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
                # Remove redundant legends in subplots (optional)
                if ax != axes[0]:
                    ax.get_legend().remove()

            # Add a single legend for the entire figure
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(handles, labels, title='Treatment', loc='upper right')
            plt.tight_layout()
            plot_path = os.path.join(outPath, 'plots', 'boxplots_preProcessing_all_samples_plot.png')
            plt.savefig(plot_path, dpi=300)
            plt.close()
#####
#####
            # some proteins do not produce any associated genes. these values are left blank in the index
            # we replace the NaNs with Unknown-Gene-X, where X is a unique number for each unknown gene.\
            # Convert index to a Series to manipulate NaNs
            index_series = df.index.to_series()
            # Find NaN values in index
            nan_mask = index_series.isna()
            # Replace NaNs with "Unknown-gene-N"
            index_series[nan_mask] = [f"Unknown-gene-{i+1}" for i in range(nan_mask.sum())]
            # Set updated index
            df.index = index_series
            # drop rows with duplicated index values <<<< This is a weird quirk of the phosphoproteomics data. Need to find a way of uniquely identifying rows. Added to github issue
            df = df[df.index.duplicated(keep=False) == False]
            #### We also want to remove low abundance proteins. Explained in the docs
            #### But in short, zero abundances are hard to interpret - is the protein really absent or did you just not see it by chance?
            # 1. Remove proteins where abundance is 0 in all samples for any treatment
            for treatment in metadata['treatment'].unique():
                treatment_samples = metadata.loc[metadata['treatment'] == treatment, 'sample_rep'].tolist()  # Get sample names
                mask_all_zero = (df[treatment_samples] == 0).all(axis=1)
                df = df.loc[~mask_all_zero]
            # 2. Remove proteins where abundance is 0 in ≥50% of all samples
            mask_half_zero = (df == 0).sum(axis=1) >= (df.shape[1] / 2)
            df = df.loc[~mask_half_zero]
    df = df.drop_duplicates()
    df = df.dropna(axis=0) # where there are missing values, we will remove the protein (rows), not the sample (cols)   
    if 'proteindata' in file_path:
        return df, nrow_original
    else:
        return df

###########################################################################
#### This is the main function for processing data, combining the above ###
###########################################################################

def preprocess_data(file_path,
                    metadata=None,
                    json_out = None,
                    outPath = None):
    """
    Full preprocessing pipeline: load, clean, and normalize column names.

    Parameters:
    - file_path (str): Path to the data file.

    Returns:
    - pd.DataFrame: Preprocessed dataframe.
    """
    df = load_data(file_path)
    if df is not None:
        df = normalise_column_names(df, file_path=file_path, metadata=metadata)       
    #     # Clean data with or without metadata
    #     if metadata is None:
    #         df = clean_data(df, file_path = file_path)
    #     else:
    #         df, row_nan_count = clean_data(df, file_path = file_path, metadata=metadata) 
    
    if 'metadata' in file_path:
        df = clean_data(df, file_path = file_path)
#### df['protein_abundance_name'] = df['protein_abundance_name'].str.lower().str.replace(' ', '_')
        
        # generate some info to go into html report
        # Compute treatment-wise sample counts
        treatment_counts = df['treatment'].value_counts().to_dict()

        # Format as "group 1: n1, group 2: n2"
        treatment_summary = ", ".join([f"{treatment}: {n}" for treatment, n in treatment_counts.items()])

        metadata_values = {
            "NUM_SAMPLES": df['sample_rep'].nunique(),
            "NUM_TREATMENTS": df['treatment'].nunique(),
            "TREATMENTS": treatment_summary
        }

        # Save to a JSON file
        with open(json_out, "w") as f:    
            json.dump(metadata_values, f)

        ### run function to validate metadata
        validate_metadata(df)

        return df # if just metadata

    
    if 'proteindata' in file_path:
        df, nrow_original = clean_data(df, file_path = file_path, metadata=metadata, outPath = outPath) 
        
        #### protein summary
        NUM_PROTS = len(df.index)
        NUM_PROTS_REMOVED = nrow_original - NUM_PROTS
        
        # Compute mean value per row (per protein across all samples)
        mean_abundance = df.iloc[:, 1:].mean(axis=1)  # Exclude protein column

        # Compute min, max, and mean of these mean abundances
        abundance_stats = {
            "NUM_PROTS_OG": nrow_original,
            "NUM_PROTS_REMOVED": NUM_PROTS_REMOVED,
            "NUM_PROTS": f"{NUM_PROTS:,.0f}",
            "MIN_AVERAGE_ABUNDANCE": f"{mean_abundance.min():,.0f}",
            "MAX_AVERAGE_ABUNDANCE": f"{mean_abundance.max():,.0f}",
            "MEDIAN_AVERAGE_ABUNDANCE": f"{mean_abundance.median():,.0f}"
        }

        # read data from json file
        with open(json_out, "r") as f:
            existing_data = json.load(f)

        # Append new data
        existing_data.update(abundance_stats)

        # Write back to JSON file
        with open(json_out, "w") as f:
            json.dump(existing_data, f, indent=4)

        #### Standardising
        # PCAs find ways of maximising the distance between samples.
        # That means if we have a variable that has a much higher (or lower) range/scale than others
        # the PCA will be biased by that. We can scale the data so that they contribute to the PCA based
        # on the variation among samples within each variable rather than biasing towards those that are particularly variable
        # Standardising the features
        # standardscaler standardises: standardised score = ( sample - mean ) / SD
        # column-wise, i.e. by protein, not sample.
        df_T = df.T
        df_T.columns = df_T.columns.astype(str)
        df_scaled = StandardScaler().fit_transform(df_T)  # scale the data
        df_standardised = pd.DataFrame(df_scaled, index = df_T.index, columns = df_T.columns) # add additional data
        pd.DataFrame(df).to_csv(os.path.join(outPath, 'data/protAbundance.csv'), index=False)
        pd.DataFrame(df_standardised).to_csv(os.path.join(outPath, 'data/protAbundance_standardised.csv'), index=False)
        
        ### run function to validate protein abundance data
        validate_proteindata(data=df,
                             data_standardised = df_standardised, 
                             metadata = metadata)
        
        return df, df_standardised


### make outdir
def make_outdir(out_path,
                make_subdirs = True):
    """
    Make output dir, including checks

    Parameters:
    - out_path (str): Path to where project output should be stored.

    Returns:

    """
    if not os.path.exists(out_path):

        try:
            os.makedirs(out_path)              
            print(f"Directory '{out_path}' created successfully.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{out_path}'.")
        except Exception as e:
            print(f"An error occurred: {e}")   

    if (make_subdirs == True):
        out_subDirs = ['data', 'plots']
        for subDir in out_subDirs:
            path = os.path.join(out_path, subDir)
            if not os.path.exists(path):
                os.mkdir(path)
 


### Validation metadata function ###
 # checks column names, sample id duplicates, replicate is numeric, sample data columns are unique
def validate_metadata(metadata):
    required_columns = ['sample_id', 'treatment', 'replicate', 'protein_abundance_name']

    # Check if all required columns are present
    missing_columns = [col for col in required_columns if col not in metadata.columns]
    if missing_columns:
        raise ValueError(f"Error: Missing required columns: {missing_columns}")
    
    # Check for missing values (NaNs) in the entire DataFrame
    if metadata.isna().any().any():
        raise ValueError("Error: The metadata contains missing (NaN) values!")

    # Check if the combination of 'sample_id' and 'replicate' is unique
    if metadata[['sample_id', 'replicate']].duplicated().any():
        raise ValueError("Error: Each (sample_id, replicate) pair must be unique.")

    # Check if 'replicate' is numeric
    if not pd.api.types.is_numeric_dtype(metadata['replicate']):
        raise ValueError("Error: 'replicate' must be numeric.")

    # Check for duplicate protein_abundance_name values
    if metadata['protein_abundance_name'].duplicated().any():
        raise ValueError("Error: 'protein_abundance_name' must not contain duplicate values.")
    
### Validation protein abundance data ###
def validate_proteindata(data,
                         data_standardised,
                         metadata):
    ##### checks for raw protein abundance data #####
    # Check if data is empty
    if data.empty:
        raise ValueError("Error: The protein abundance data is empty.")
    # Check if row indices (proteins) are unique
    if not data.index.is_unique:
        raise ValueError("Error: Protein identifiers (row indices) must be unique.")
    # Check if column names (samples) are unique
    if not data.columns.is_unique:
        raise ValueError("Error: Sample identifiers (columns) must be unique.")
    # Validate protein abundance columns from metadata
    # names of abundance columns:
    abundance_columns = metadata['sample_rep']  # Columns to check. Note this depends on whether cols have been renamed
    # Ensure one-to-one mapping: all sample_id values in columns, and all columns are sample_ids
    if set(abundance_columns) != set(data.columns):
        raise ValueError("Error: Mismatch between metadata samples and protein abundance columns. Ensure a strict one-to-one mapping.")
    # If all protein_abundance_name columns are present,
    # Check numeric values only in specified abundance columns
    non_numeric_columns = [
        col for col in abundance_columns
        if not pd.api.types.is_numeric_dtype(data[col])
    ]
    if non_numeric_columns:
        raise ValueError(f"Error: The following columns contain non-numeric values: {non_numeric_columns}")
    # Check for missing values (NaNs) in the entire DataFrame
    if data.isna().any().any():
        raise ValueError("Error: The protein abundance df contains missing (NaN) values!")
    ##### checks for standardised protein abundance data #####
    # Check for missing values (NaNs) in the entire DataFrame
    if data_standardised.isna().any().any():
        raise ValueError("Error: The standardised protein abundance df contains missing (NaN) values!")
    ### check mean ~0 and variance ~1
    # Calculate column-wise mean and variance
    mean_vals = data_standardised.mean()
    var_vals = data_standardised.var(ddof=0) # ddof = 0 divides deviance by N-1 for sample var. Default divides by N and gives population var.
    # Define tolerance for deviation from expected values
    mean_tolerance = 1e-6  # Mean should be very close to 0
    var_tolerance = 0.1    # Variance should be around 1, allowing small deviations
    # Check mean is close to 0
    if not np.all(np.abs(mean_vals) < mean_tolerance):
        raise ValueError("Error: Standardised data mean is not sufficiently close to 0.")
    # Check variance is close to 1
    if not np.all(np.abs(var_vals - 1) < var_tolerance):
        raise ValueError("Error: Standardised data variance is not sufficiently close to 1.")

#!/bin/bash -l
# Request one hour of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=1:0:0
# Request 16 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=16G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set the name of the job.
#$ -N auto-prot

# Set the working directory to somewhere in your scratch space.  
#  This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
#$ -wd /home/<your_UCL_id>/Scratch/workspace

# Your work should be done here:
DIR=/path/to/working/directory

cd $DIR

########################
##### Load modules #####
########################

#### load conda
module load python/miniconda3
source $UCL_CONDA_PATH/etc/profile.d/conda.sh 

#######################
##### Set Options #####
#######################

### here we set parameters that determine behaviour later in the script - make sure these are correct!

### If you have not previously cloned the repository, set this to YES
clone_repo="NO"

### set branch to use
# stable - the main branch of the package and most recent, stable version
# staging - the staging branch of the package, which may contain new features and bug fixes that are not yet in the stable branch
# others are available, but not recommended for general use. See github for more details.
branch_name="dev-package-myriad-script"

### If you have not previously created the conda environments
### Set this to YES
create_environments="NO"

### if you want to generate random data, set this to YES
### note, this will overwrite input/data/proteindata.csv and input/data/metadata.csv
generate_data="YES"

###############################
##### Prepare Environment #####
###############################

#### Clone the repo
if [ $clone_repo = "YES" ]; then
    echo "Cloning Auto-Prot repository..."
    git clone https://github.com/jdgilbert245/Auto-Prot
    ### change to Auto-Prot directory
    cd Auto-Prot
fi

### if you want to use a branch of the package that is not the main branch ("stable"), set this to the branch name
git checkout $branch_name

### create environments (optional, see above)
if [ $create_environments = "YES" ]; then
    echo "Creating conda environments..."
    conda env create -f configs/auto-prot-env-myriad.yml
    conda env create -f configs/auto-prot-env-limma-myriad.yml
else
    echo "Skipping environment creation, assuming environments already exist."
fi

### activate environment
echo "Activating conda environment..."
conda activate auto-proteomics

##########################
##### Preparing Data #####
##########################

#### generate random data (optional, see above)
if generate_data == "YES"; then
    print( "Generating random data...")
    python ${DIR}/autoprot/utils/gen_random_data.py
fi

### if you have your own data, save the protein data and metadata in the input data folder
### the file names need to match those in the config file (configs/auto-prot-config.yaml)
print( "running Auto-Prot...")
python main.py


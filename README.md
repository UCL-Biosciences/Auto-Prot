# Proteomics Report
This repository contains the Auto-Proteomics analysis pipeline for generating standardised reports from protein intensity datasets. Welcome.

The pipeline is still in testing. Staging is the best branch to use but there might still be some things to fix. Would only suggest trying this if you know python and/or proteomics reasonably well :') clone the repo using the instructions below, then run `git checkout staging` to change to the staging branch.

This tool produces common outputs from MS proteomic and phosphoproteomic data. We want to make this possible for people with limited coding experience by requiring only some preparation steps and one line of code to run the pipeline. Parameters are determined in a "config" file so you don't need to look through the code to change things.

By producing common exploratory outputs automatically, we hope you can spend more time thinking about your data and results. Note, this pipeline produces results automatically without thinking about the idiosyncracies of your beautiful, unique dataset. **This should never be used without a full evaluation and validation of the results.** Best practice would be to run a similar analysis with your favourite software or collaborator to confirm our results are sensible.

If you think the results look reasonable and you want to make use of them, please review the tool carefully. We hope the docs give you an idea of what we are doing and why. We have tried to annotate the code so that you can also understand how the pipeline works. If you would like further info or to suggest an improvement to the pipeline, please get in touch. Contact details below.

## Features
- Processes mass spec-generated protein abundance data.
- Generates an HTML report with common proteomics outputs.
- Customisable templates.

## 🚀 Quick Start
### Branches
`stable` branch is the most tested and should work reliably. If not, please let me know! It would be best to quick start using that branch. `staging` has new features that haven't been fully tested. Other branches are for developing new features.

### Prerequisites
- a [conda](https://github.com/conda-forge/miniforge) distribution.
- [git for mac](https://git-scm.com/install/mac) or [Git bash for windows](https://git-scm.com/install/windows).
- Somewhere to run command line commands (e.g. [VS Code](https://code.visualstudio.com/)).

This tool has been tested on Windows and macOS.

### 1. Clone the repository
`git clone https://github.com/UCL-Biosciences/Auto-Prot`

### 2. Navigate into the project
`cd Auto-Prot`

**Optional** switch to staging branch: `git checkout staging`.

This is important - it allows you to make changes to the config file withour your setup being sent to the main repo: `git update-index --skip-worktree configs/auto-prot-config.yaml`

### 3. Install dependencies for your operating system
Create the environment and download dependencies using the information in configuration files: `conda env create -f configs/auto-prot-env-YOUR-OS-HERE.yml`. Which version you use will depend on whether you are using a Windows machine or a Mac. You will need to create environments for the general pipeline and the R functions:

- For general processing:  
  `conda env create -f configs/auto-prot-env-windowsOS.yml`
- For R-based differential expression:  
  `conda env create -f configs/auto-prot-env-limma-windowsOS.yml`

Activate the environment with: `conda activate auto-proteomics`. If you keep getting a message about running `git init`, try running `source activate`. You might have to do that at the start of every session.

*Optional* After creating the conda environment, run: `pre-commit install`. This sets up automatic code formatting and linting before each commit.

### 4. Generate outputs
#### Input
If you don't add your files to the `data` folder, you can generate an example dataset by running `python autoprot/utils/gen_random_data.py`. This will give a quick picture of how the input data should be formatted. Make sure your data are in `data/proteindata.csv` and `data/metadata.csv` (requirements given in /docs/Workflow.md). 

#### Generate output
To generate the output, you can run `python main.py`.

#### Outputs
Running `python main.py` will generate:
- Normalised protein abundance tables
- Summary statistics
- PCA and clustering plots
- Differential expression tables
- An html report summarising the results. An example generated from simulated data is available in `../docs/`

These are saved in `/output`. An example of the report generated is in [`report/example-report-out.html`](https://github.com/UCL-Biosciences/Auto-Prot/blob/staging/report/example-report-out.html) - if it doesn't show properly in the repo, try downloading and opening locally.

## 📂 Project Structure
```
/docs                 # Documentation files  
/report            # report template
/autoprot             # Processing & analysis scripts  
/main.py           # Python script for generating outputs  
/input/data                # input data. By default: `human_genes.txt`. Optional: `run autoprot/utils/gen_random_data.py` to generate additional proteindata.csv and metadata.csv in input/data
/output              # plots and tables generated by `main.py`. report saved as output/report-out.html
```

## 📝 Analysis Description
For more details on the analysis and calculations, see `docs/Workflow.md` and `docs/stats_details.md`/

## Contact and feedback
We want this tool to be useful and accessible for as many people as possible - comments and feedback are very welcome. Please send to [james.d.gilbert@ucl.ac.uk](mailto:james.d.gilbert@ucl.ac.uk). 

For help, please first open an issue and send an email if required.

## ⚙️ Configuration
There is a configuration (config) file that allows the user to control lots of parts of the analysis without editing any code. For instructions, see the workflow document. Changes to the config file are ignored (config file in .gitignore) so the default config doesn't get overwritten and continues to work when downloaded.

## 🤝 Contributing

Contributions are welcome!

If you’d like to contribute a bug fix, feature, or improvement:
1. Fork this repository and create a new branch from `staging`
2. Make your changes
3. Ensure code passes `pre-commit` checks (run `pre-commit run --all-files`)
4. Open a pull request to the `staging` branch with a short description of your changes

Please note:
- All code should be commented and use informative function/variable names

#### Project Structure and Packaging
This project uses a standard Python package layout and is installable via `pip`. The core code lives in the `autoprot/` directory and is defined as a package using `pyproject.toml`. Note: The package is not published on PyPI. Users must clone the repository and install it locally. It can also be cloned and used via `python main.py` - pip installing is not required.

#### Building for distribution
If you update the package structure or plan to share the repo: `python -m build` will generate `.tar.gz` and `.whl` files in the `dist/` directory. These are needed for `pip install .`

#### Testing
The project uses pytest for testing. Test files are located in the tests/ directory and follow standard pytest conventions.

To run all tests: `pytest`. If you’re developing new features or modifying existing functionality, please include corresponding tests to ensure the code remains correct and maintainable. Tests should cover both typical and edge cases where applicable.

To see the test coverage, run `pytest --cov=autoprot --cov-report=term-missing` from the main folder. This will show coverage for the autoprot package and highlight any lines not exercised by tests.

## Further Analysis
The outputs generated by this tool are exploratory only. We recommend a thorough examination of the data using the code in this repository or with other available tools. E.g.:
- [MS-DAP](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00513) R package
- [proDA](https://github.com/const-ae/proDA?tab=readme-ov-file) seems a good option if you have a lot of missing values or want to explicitly include missingness in your calculations.
- [AlphaPepStats](https://github.com/MannLabs/alphapept) Python package
- [Perseus](https://maxquant.net/perseus/) platform has a GUI with lots of functionality and doesn't require any coding.
  


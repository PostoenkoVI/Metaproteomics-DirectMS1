# MetaDirectMS1
`MetaDirectMS1` is a python-based full-workflow pipeline for metaproteomics analysis. It is designed to identify the most represented organisms in a sample based on data from an ultrafast proteomic (LC-MS) experiment and quantify its proteins. This is achieved by a three-stage search based on a `.fasta` database of potentially present organisms in the sample. For example, for microbiome researches full Uniprot database of bacteria is used.

## Installation

    pip install git+https://github.com/PostoenkoVI/Metaproteomics-DirectMS1

## Usage
Workflow consists of 5 major stages:

#### 1) Initial database reduction
Here input `.fasta` database is parsed to reduce it approximately 10-times, while keeping the most well-described organisms in it. It is the most time- and resource- consuming stage so it is not recommended to re-run it without necessity. Results of this stage are stored in `outdir/db_parcnig_folder`. It contains two subfolders with individual organisms `.fasta` files and could use a lot of disk space.

#### 2) Input files processing
Alongside with `.fasta` to search into, algorithm needs LC-MS data to analyze. There are three formats of accepted input files: `.raw` files, `.mzML` files, and `features.tsv` tab-separated tables with [Biosaur2](https://github.com/markmipt/biosaur2) output. These formats are then consistently converted into Biosaur2 features to use in further stages. Conversion also could be consuming depending on the complexity and sizes of experimental data. Generation of the `.mzML` files is using [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) with `--msLevel=1` setting by default to reduce operating time and used disk space, though any `.mzML` file is acceptable. Results of this stage are stored in `outdir/mzml`, `outdir/features`.

#### 3) Blind search
On this stage algorithm searches for the top N (15 by default) most represented organisms in each file of the input data. Then combines their proteins into the file-specific `.fasta` databases and runs [ms1searchpy](https://github.com/markmipt/ms1searchpy) - ultrafast proteomics search engine for LC-MS1 spectra - to identify as much proteins from those top N organisms as it can. All intermediate results for this stage are stored in `outdir/blind_search`. Combined tables and barplots with numbers of proteins found in each file for each organism or taxonomy group are stored in `outdir/results`. 

#### 4) Precise search
On this stage, based on the previous searches, algorithm unites individual `.fasta` files into one for another search on each file to be able to compare results between files and use it in quantitation. The criteria for the organism to be included in united database is simple - at least in one file at least 2% from all found (in that file on previous stage) proteins should be relate to the organism or taxonomy group. The results tables from previous stage show which organisms are written in united database. Additionally, it is possible for user to provide custom `.fasta` database not to unite databases automatically. Again, all intermediate results for this stage are stored in `outdir/precise_search`. Combined tables and barplots with numbers of proteins found in each file for each organism or taxonomy group are stored in `outdir/results`. 

#### 5) Quantitation
Quantitation stage conducts quantitative analysis on existing identification results from Precise search stage using [DirectMS1quantmulti](https://github.com/markmipt/ms1searchpy?tab=readme-ov-file#multi-condition-protein-profiling-using-directms1quantmulti). It needs sample file with details for all project files in format of a tab-separated table. More details on quantitation are described [there](https://github.com/markmipt/ms1searchpy?tab=readme-ov-file#multi-condition-protein-profiling-using-directms1quantmulti).

### Examples
#### Basic usage
    MetaDirectMS1 -full_uniprot_fasta /path/to/uniprot_fasta_file.fasta -mzml_folder /path/to/folder/with/mzml_files -ms1searchpy /path/to/ms1searchpy -biosaur2 /path/to/biosaur2
Though it is highly recommended to specify `-outdir`, by default `MetaDirectMS1` creates folder `current_date_MetaDirectMS1` in the input files folder.
#### Config file
    MetaDirectMS1 -cfg /path/to/config_file

OR

    MetaDirectMS1 -cfg /path/to/config_file -cfg_category my_settings
Each category of config file is independent set of parameters and only one is read from file. By default it uses DEFAULT category. Example of config file presented in the project folder on github.
#### Modes
`MetaDirectMS1` has `-mode` option that allows to run algorithm starting with every stage, rewriting results or continuing previous analysis. However, generally, to rewrite starting with `3) Blind search` it needs intermediate results from both previous stages, to unite protein database in `4) Precise search` individual `.fasta` files are needed and so on.

So if in any case analysis stopped before sucessfull ending, `-mode 0` option with exact same inputs and output folder should help to continue workflow from last completed stage. By default `-mode 1` is using to rewrite any old results in output folder.

### Full cmd options description
<...>
## Contacts
Valeriy Postoenko - v.i.postoenko@gmail.com

Elizaveta Kazakova - kazakovaem@gmail.com

Mark Ivanov - markmipt@gmail.com

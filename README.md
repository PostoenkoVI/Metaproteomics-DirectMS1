# MetaDirectMS1
`MetaDirectMS1` is a python-based full-workflow pipeline for metaproteomics analysis. It is designed to identify the most represented organisms in a sample based on data from an ultrafast proteomic (LC-MS) experiment and quantify its proteins. This is achieved by a three-stage search based on a `.fasta` database of potentially present organisms in the sample. For example, for microbiome researches full Uniprot database of bacteria is used.

## Installation
It is recommended to additionally install retention time prediction tool DeepLC version 1.1.2.2 (unofficial fork with small changes). Newer version has some issues right now. It is highly important to install DeepLC before ms1searchpy for outdated packages compatibility! Thus, the recommended way:

Create and activate virtual environment with python version 3.10.11 (for example using `pyenv` package, for detailed guide see [link](https://akrabat.com/creating-virtual-environments-with-pyenv/)) :

    pyenv install 3.10.11
    pyenv virtualenv 3.10.11 metadirectms1_env
    pyenv activate metadirectms1_env

Next, install DeepLC in clean environment to avoid version collisions:

    pip install https://github.com/markmipt/DeepLC/archive/refs/heads/alternative_best_model.zip

After that, install ms1searchpy search engine:

    pip install git+https://github.com/PostoenkoVI/ms1searchpy

And then, install MetaDirectMS1 freely:

    pip install git+https://github.com/PostoenkoVI/Metaproteomics-DirectMS1

It is also recommended to update the NCBI taxonomy database before using MetaDirectMS1. To do so, run the following commands in the same Python environment and wait for a successful download.

    from ete3 import NCBITaxa
    NCBITaxa().update_taxonomy_database()

## Usage & Examples 
### Basic usage
    MetaDirectMS1 -full_uniprot_fasta /path/to/uniprot_fasta_file.fasta -mzml_folder /path/to/folder/with/mzml_files -ms1searchpy /path/to/ms1searchpy -biosaur2 /path/to/biosaur2
Though it is highly recommended to specify `-outdir`, by default `MetaDirectMS1` creates folder `current_date_MetaDirectMS1` in the input files folder.
### Config file
    MetaDirectMS1 -cfg /path/to/config_file

OR

    MetaDirectMS1 -cfg /path/to/config_file -cfg_category my_settings
Each category of config file is independent set of parameters and only one is read from file. By default it uses DEFAULT category. Example of config file presented in the project folder on github.
### Modes
`MetaDirectMS1` has `-mode` option that allows to run algorithm starting with every stage, rewriting results or continuing previous analysis. However, generally, to rewrite starting with `3) Blind search` it needs intermediate results from both previous stages, to unite protein database in `4) Precise search` individual `.fasta` files are needed and so on.

So if in any case analysis stopped before sucessfull ending, `-mode 0` option with exact same inputs and output folder should help to continue workflow from last completed stage. By default `-mode 1` is using to rewrite any old results in output folder.

### Tooltips
There is an important point in the analysis where user's attention and interpretation is needed. By default MetaDirectMS1 is running all its stages one-by-one, however it is helpfull to check the size of the automatically generated precise search protein database in the log file or manually after the blind search is done. It is recommended to keep the protein database size lower than 200-300 thousands proteins for the search engine to work efficently. To control the database size the file `outdir/results/blind_identified_proteins_GROUP.tsv` contains the column `include in combined fasta` (where `GROUP` is a taxonomy level for the analysis, `OX` by default). Accordingly, after the blind search stage is done it is recommended to check the `blind_identified_proteins_GROUP.tsv` and manually edit the column `include in combined fasta` excluding poorly detected taxonomy groups (or choosing the groups of interest if any prior information about the samples is availiable) to avoid extra wide search space.

### Known issues
Running MetaDirectMS1 in Jupyter notebook an error "ValueError: Key backend: 'module://matplotlib_inline.backend_inline' is not a valid value for backend; supported values ..." could happen. This is specific to Jupyter and caused by conflict of the virtual environments and is solved by installing matplotlib-inline inside the notebook.

## Workflow stages & Output description

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

## Full cmd options description
options:
  `-h`, `--help`            show this help message and exit
  
  `-logs` [{DEBUG,INFO,WARNING,ERROR,CRITICAL,NOTSET,debug,info,warning,error,critical,nonset}]
                        level of logging, (DEBUG, INFO, WARNING, ERROR,
                        CRITICAL) (default: INFO)
                        
  `-log_path` [LOG_PATH]  Path to logging file. By default it is in the outdir.
                        (default: )
                        
  `-full_uniprot_fasta` [FULL_UNIPROT_FASTA]
                        Path to full bacterial SwissProt+TrEMBL .fasta file to
                        search organisms in (default: )
                        
  `-uniprot_folder` [UNIPROT_FOLDER]
                        Path to folder to store .fasta files splited by
                        taxonomic identifiers from uniprot (default: )
                        
  `-sprot_folder` [SPROT_FOLDER]
                        Path to folder to store .fasta files splited by
                        taxonomic identifiers from swissprot (default: )
                        
  `-input_files` [INPUT_FILES ...]
                        input .raw, .mzML or features.tsv files (default: )
                        
  `-raw_folder` [RAW_FOLDER]
                        Input directory with .raw files (default: )
                        
  `-outdir` [OUTDIR]      Output directory (default: )
  
  `-mzml_folder` [MZML_FOLDER]
                        Directory to search or to store .mzML files after
                        convertation (default: )
                        
  `-feature_folder` [FEATURE_FOLDER]
                        Path to folder where to store and search for biosaur2
                        files that ends with features.tsv, default:
                        outdir/features (default: )
                        
  `-db_parsing_folder` [DB_PARSING_FOLDER]
                        Path to folder to store results of parsing initial
                        .fasta database, that could help reproduce analysis
                        faster (default: outdir/db_parsing_folder) (default: )
                        
  `-cfg` [CFG]            Path to config file (default: )
  
  `-biosaur2` [BIOSAUR2]  path to Biosaur2 (default: )
  
  `-mono` [MONO]          path to mono to use ThermoRawParser (default: )
  
  `-thermorawparser` [THERMORAWPARSER]
                        path to ThermoRawParser.exe to use ThermoRawParser
                        (default: )
                        
  `-ms1searchpy` [MS1SEARCHPY]
                        path to ms1searchpy (default: )
                        
  `-precise_search_fasta` [PRECISE_SEARCH_FASTA]
                        path to custom .fasta file for precise search. By
                        default fasta will be created automatically. (default:
                        )
                        
  `-sample_file` [SAMPLE_FILE]
                        path to file with grouping of input files into
                        comparison groups. Tab-separated, one pair file_name
                        (without extension) - group_name by line. (default: )
                        
  `-cfg_category` [CFG_CATEGORY]
                        Name of category in config file to use (default:
                        DEFAULT) (default: DEFAULT)
                        
  `-decoy_prefix` [DECOY_PREFIX]
                        String added to the protein name to showcase that it
                        is a decoy (default: DECOY) (default: DECOY)
                        
  `-sprot_suf` [SPROT_SUF]
                        String added to the taxid of the organism's .fasta
                        file to showcase swissprot database (default: _sp)
                        (default: _sp)
                        
  `-uniprot_suf` [UNIPROT_SUF]
                        String added to the taxid of the organism's .fasta
                        file to showcase uniprot database (default: _un)
                        (default: _un)
                        
  `-min_pept_length` [MIN_PEPT_LENGTH]
                        Minimal peptide lenght to use in cleaving 10%
                        database, lowering the value is not recommended
                        (default: 9) (default: 9)
                        
  `-missed_cleavages` [MISSED_CLEAVAGES]
                        Number of missed cleavages for cleaving peptide in
                        both searches. (default: 0) (default: 0)
                        
  `-cleavage_rule` [CLEAVAGE_RULE]
                        cleavage rule in quotes!. X!Tandem style for cleavage
                        rules: "[RK]|{P}" for trypsin, "[X]|[D]" for asp-n or
                        "[RK]|{P},[K]|[X]" for mix of trypsin and lys-c
                        (default: [RK]|{P})
                        
  `-fmods` [FMODS]        fixed modifications. Use "[" and "]" for N-term and
                        C-term amino acids. in
                        psiname1@aminoacid1,psiname2@aminoacid2 format
                        (default: Carbamidomethyl@C)
                        
  `-fmods_legend` [FMODS_LEGEND]
                        PSI Names for extra fixed modifications. Oxidation,
                        Carbamidomethyl and TMT6plex are stored by default in
                        source code. in psiname1@monomass1,psiname2@monomass2
                        format (default: )
                        
  `-mass_accuracy` [MASS_ACCURACY]
                        Mass accuracy in ppm for blind search. (default: 4)
                        (default: 4)
                        
  `-mz_for_mass_accuracy` [MZ_FOR_MASS_ACCURACY]
                        Approximate maximum m/z value to use with
                        mass_accuracy. (default: 1000) (default: 1000)
                        
  `-allowed_ranks` [ALLOWED_RANKS]
                        Allowed taxonomy categories between <...> and <...>,
                        default="strain,subspecies,forma
                        specialis,isolate,serotype,serogroup,no rank"
                        (default: strain,subspecies,forma
                        specialis,isolate,serotype,serogroup,no rank)
                        
  `-exclude_sp_uncul` [{0,1}]
                        Exclude sp. and uncultured. organisms (0 - no
                        exclusion, 1 - exclude) (default: 0)
                        
  `-min_prot` [MIN_PROT]  Minimal number of proteins in leader's fasta to write
                        in 10% organisms fasta. (default: 200) (default: 200)
                        
  `-taxid_group` [TAXID_GROUP]
                        taxid level to filter organisms found in blind search,
                        availible options: 'OX', 'genus', 'family', 'order',
                        'class', 'phylum', 'kingdom', 'domain' (default: 'OX')
                        (default: OX)
                        
  `-taxid_presence_thr` [TAXID_PRESENCE_THR]
                        Percentage threshold to filter taxid with low
                        abundance in sample from united .fasta file for
                        precise search (default: 0.02)
                        
  `-score_threshold` [SCORE_THRESHOLD]
                        Minimal number of matched proteins in blind search to
                        report taxid (default: 4) (default: 4)
                        
  `-num_top_spec` [NUM_TOP_SPEC]
                        Number of taxids to report in blind search. Only top N
                        taxid by number of proteins found are reported.
                        (default: 10) (default: 10)
                        
  `-generate_figures` [{0,1}]
                        Generate figures default: 1 (default: 1)
                        
  `-mode` [{0,1,2,3,4,5}]
                        Mode for MetaDirectMS1 to work: 0 - try to continue
                        existing analysis in selected outdir without rewriting
                        anything, 1 - run all stages of analysis overwriting
                        results in outdir, 2 - overwrite all stages except
                        initial fasta parcing, 3 - overwrite all stages except
                        initial fasta parcing and feature generation (start
                        with blind search), 4 - overwrite precise search and
                        quantitation, 5 - overwrite quantitation (default: 1)
                        
  `-bio2_args` [BIO2_ARGS]
                        String of additional arguments to submit into Biosaur2
                        (in command line the string should be in double
                        quotes: '" "', in cfg file in single quotes) except:
                        -o; default: "" (default: )
                        
  `-ms1searchpy_args` [MS1SEARCHPY_ARGS]
                        String of additional arguments to submit into
                        ms1searchpy (in command line the string should be in
                        double quotes: '" "', in cfg file in single quotes)
                        except: -d, -deeplc, -e, -ad, -prefix, -ml, -ts, -o;
                        default: "" (default: )
## Contacts
Valeriy Postoenko - v.i.postoenko@gmail.com

Elizaveta Kazakova - kazakovaem@gmail.com

Mark Ivanov - markmipt@gmail.com

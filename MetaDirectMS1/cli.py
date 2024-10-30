import argparse
import sys
from os import path, makedirs
import logging
from . import utils, workflow

logger = logging.getLogger()


def run():
    argparser = argparse.ArgumentParser(
        description = 'Ultrafast metaproteomics for quantitative assessment of strain isolates and microbiomes',
        epilog = '''
    ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@')

    argparser.add_argument('-logs', nargs='?', help='level of logging, (DEBUG, INFO, WARNING, ERROR, CRITICAL)', type=str, default='INFO', const='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'NOTSET', 'debug', 'info', 'warning', 'error', 'critical', 'nonset'])
    argparser.add_argument('-log_path', nargs='?', help='path to logging file', type=str, default='./logs.txt', const='./logs.txt')
    # parser.add_argument('-example_cfg', nargs='?', help='Path to create example config .ini file or not if not stated', type=str, default='', const='')
    argparser.add_argument('-full_uniprot_fasta', nargs='?', help='Path to full bacterial SwissProt+TrEMBL .fasta file to search organisms in', type=str, default='', )
    argparser.add_argument('-uniprot_folder', nargs='?', help='Path to folder to store .fasta files splited by taxonomic identifiers from uniprot', type=str, default='', )
    argparser.add_argument('-sprot_folder', nargs='?', help='Path to folder to store .fasta files splited by taxonomic identifiers from swissprot', type=str, default='', )

    argparser.add_argument('-input_files', nargs='*', help='input .raw, .mzML or features.tsv files', type=str, default='',)
    argparser.add_argument('-raw_folder', nargs='?', help='Input directory with .raw files', type=str, default='', const='')
    argparser.add_argument('-outdir', nargs='?', help='Output directory', type=str, default='', const='')
    argparser.add_argument('-mzml_folder', nargs='?', help='Directory to search or to store .mzML files after convertation', type=str, default='', const='')
    argparser.add_argument('-feature_folder', nargs='?', help='Path to folder where to store and search for biosaur2 files that ends with features.tsv, default: outdir/features', type=str, default='', const='')
    argparser.add_argument('-db_parsing_folder', nargs='?', help='Path to folder to store results of parsing initial .fasta database, that could help reproduce analysis faster (default: outdir/db_parsing_folder)', type=str, default='', const='')
    argparser.add_argument('-cfg', nargs='?', help='Path to config file', type=str, default='', )
    # argparser.add_argument('-10per_fasta', nargs='?', help='Path to save .fasta with 10% proteins from original full database', type=str, default='', )
    argparser.add_argument('-biosaur2', nargs='?', help='path to Biosaur2', type=str, default='', const='')
    argparser.add_argument('-mono', nargs='?', help='path to mono to use ThermoRawParser', type=str, default='', const='')
    argparser.add_argument('-thermorawparser', nargs='?', help='path to ThermoRawParser.exe to use ThermoRawParser', type=str, default='', const='')
    argparser.add_argument('-ms1searchpy', nargs='?', help='path to ms1searchpy', type=str, default='', const='')
    # argparser.add_argument('-deeplc', nargs='?', help='path to deepLC for ms1searchpy to use', type=str, default='', const='')
    argparser.add_argument('-precise_search_fasta', nargs='?', help='path to custom .fasta file for precise search. By default fasta will be created automatically.', type=str, default='', const='')
    argparser.add_argument('-sample_file', nargs='?', help='path to file with grouping of input files into comparison groups. Tab-separated, one pair file_name (without extension) - group_name by line.', type=str, default='', const='')

    argparser.add_argument('-cfg_category', nargs='?', help='Name of category in config file to use (default: DEFAULT)', type=str, default='DEFAULT', const='DEFAULT')
    argparser.add_argument('-decoy_prefix', nargs='?', help='String added to the protein name to showcase that it is a decoy (default: DECOY)', type=str, default='DECOY', const='DECOY')
    argparser.add_argument('-sprot_suf', nargs='?', help='String added to the taxid of the organism\'s .fasta file to showcase swissprot database (default: _sp)', type=str, default='_sp', const='_sp')
    argparser.add_argument('-uniprot_suf', nargs='?', help='String added to the taxid of the organism\'s .fasta file to showcase uniprot database (default: _un)', type=str, default='_un', const='_un') 
    argparser.add_argument('-min_pept_length', nargs='?', help='Minimal peptide lenght to use in cleaving 10% database, lowering the value is not recommended (default: 9)', type=int, default=9, const=9)
    argparser.add_argument('-missed_cleavages', nargs='?', help='Number of missed cleavages for cleaving peptide in both searches. (default: 0)', type=int, default=0, const=0)
    argparser.add_argument('-cleavage_rule', nargs='?', help='cleavage rule in quotes!. X!Tandem style for cleavage rules: "[RK]|{P}" for trypsin, "[X]|[D]" for asp-n or "[RK]|{P},[K]|[X]" for mix of trypsin and lys-c', default='[RK]|{P}', const='[RK]|{P}')
    # argparser.add_argument('-cmin', nargs='?', help='Min precursor charge', default=2, type=int, const=3)
    # argparser.add_argument('-cmax', nargs='?', help='Max precursor charge', default=3, type=int, const=3)
    argparser.add_argument('-fmods', nargs='?', help='fixed modifications. Use "[" and "]" for N-term and C-term amino acids. in psiname1@aminoacid1,psiname2@aminoacid2 format', default='Carbamidomethyl@C', const='Carbamidomethyl@C')
    argparser.add_argument('-fmods_legend', nargs='?', help='PSI Names for extra fixed modifications. Oxidation, Carbamidomethyl and TMT6plex are stored by default in source code. in psiname1@monomass1,psiname2@monomass2 format', default='', const='')
    argparser.add_argument('-mass_accuracy', nargs='?', help='Mass accuracy in ppm for blind search. (default: 4)', default=4, type=float, const=4)
    argparser.add_argument('-mz_for_mass_accuracy', nargs='?', help='Approximate maximum m/z value to use with mass_accuracy. (default: 1000)', default=1000, type=float, const=1000)
    argparser.add_argument('-allowed_ranks', nargs='?', help='Allowed taxonomy categories between <...> and <...>, default="strain,subspecies,forma specialis,isolate,serotype,serogroup,no rank" ', type=str, default='strain,subspecies,forma specialis,isolate,serotype,serogroup,no rank', const='strain,subspecies,forma specialis,isolate,serotype,serogroup,no rank')
    argparser.add_argument('-min_prot', nargs='?', help="Minimal number of proteins in leader's fasta to write in 10% organisms fasta. (default: 200)", default=200, type=float, const=200)
    argparser.add_argument('-taxid_group', nargs='?', help="taxid level to filter organisms found in blind search, availible options: OX, genus, family, species. (default: OX)", default='OX', type=str, const='OX')
    argparser.add_argument('-taxid_presence_thr', nargs='?', help='Percentage threshold to filter taxid with low abundance in sample from united .fasta file for precise search', default=0.02, type=float, const=0.02)
    argparser.add_argument('-score_threshold', nargs='?', help='Minimal number of matched proteins in blind search to report taxid (default: 4)', default=4, type=float, const=4)
    argparser.add_argument('-num_top_spec', nargs='?', help='Number of taxids to report in blind search. Only top N taxid by number of proteins found are reported. (default: 15)', default=15, type=float, const=15)
    ##############
    argparser.add_argument('-generate_figures', nargs='?', help='Generate figures default: 1', default=1, type=int, choices=[0, 1,], const=1)
    argparser.add_argument('-mode', nargs='?', help='Mode for MetaDirectMS1 to work: 0 - try to continue existing analisys in selected outdir without rewriting anything, 1 - run all stages of analisys overwriting results in outdir, 2 - overwrite all stages except initial fasta parcing, 3 - overwrite all stages except initial fasta parcing and feature generation (start with blind search), 4 - overwrite precise search and following steps, 5 - overwrite quantitation', default=1, type=int, choices=[0, 1, 2, 3, 4, 5], const=1)
    ##############
    
    console_config = vars(argparser.parse_args()) 
    console_keys = [x[1:] for x in sys.argv if x.startswith('-')] 
    default_config = vars(argparser.parse_args([]))
    additional_config = {}

    if console_config['cfg'] :
        if path.exists(console_config['cfg']) :
            if console_config['cfg_category'] :
                s, additional_keys = utils.read_cfg(console_config['cfg'] , console_config['cfg_category'] )
            else :
                s, additional_keys = utils.read_cfg( console_config['cfg'], default_config['cfg_category'] )
            additional_config = vars(argparser.parse_args(s))
        else :
            logging.warning('Path to config file does not exist')

    args = default_config
    if additional_config :
        for k in additional_keys :
            args.update({k: additional_config[k]})
    for k in console_keys :
        args.update({k: console_config[k]})
    
    loglevel = args['logs'].upper()
    numeric_level = getattr(logging, loglevel, None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    logger.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s')
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args['log_path'] :
        args['log_path'] = path.abspath(path.normpath(args['log_path']))
        log_directory = path.dirname(args['log_path'])
        makedirs(log_directory, exist_ok=True)
        fh = logging.FileHandler(args['log_path'], mode='w', encoding='utf-8')
        fh.setLevel(numeric_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    logging.getLogger('matplotlib').setLevel(logging.ERROR)

    logger.info('MetaDirectMS1 started')

    # if args['example_cfg'] :
    #     p = os.path.abspath(os.path.normpath(args['example_cfg']))
    #     if os.path.exists(os.path.dirname(p)) or os.path.exists(p) :
    #         if os.path.exists(p) :
    #             logger.info('Example cfg would be overwrited')
    #         utils.write_example_cfg(args['example_cfg'], default_config)
    #         logger.info('Example cfg created')
    #         return 0
    #     else :
    #         logger.warning('Invalid path for example cfg creation. Directory does not exist')
    #         return 1

    workflow.process_files(args)
    
    logger.info('MetaDirectMS1 has completed the workflow.')

if __name__ == '__main__':
    run()
import argparse
import sys
import os
import logging
from . import utils, workflow

logger = logging.getLogger()


def run():
    parser = argparse.ArgumentParser(
        description = 'Ultrafast metaproteomics for quantitative assessment of strain isolates and microbiomes',
        epilog = '''
    ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@')

    parser.add_argument('-logs', nargs='?', help='level of logging, (DEBUG, INFO, WARNING, ERROR, CRITICAL)', type=str, default='INFO', const='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'NOTSET', 'debug', 'info', 'warning', 'error', 'critical', 'nonset'])
    parser.add_argument('-log_path', nargs='?', help='path to logging file', type=str, default='./logs.txt', const='./logs.txt')
    # parser.add_argument('-example_cfg', nargs='?', help='Path to create example config .ini file or not if not stated', type=str, default='', const='')
    parser.add_argument('-full_uniprot_fasta', nargs='?', help='Path to full bacterial SwissProt+TrEMBL .fasta file to search organisms in', type=str, default='', )
    parser.add_argument('-uniprot_folder', nargs='?', help='Path to folder to store .fasta files splited by taxonomic identifiers from uniprot', type=str, default='', )
    parser.add_argument('-swissprot_folder', nargs='?', help='Path to folder to store .fasta files splited by taxonomic identifiers from swissprot', type=str, default='', )
    
    parser.add_argument('-input_files', nargs='?', help='input .raw or .mzML files to search, separated by whitespace', type=str, default='',)
    parser.add_argument('-raw_folder', nargs='?', help='Input directory with .raw files', type=str, default='', const='')
    parser.add_argument('-outdir', nargs='?', help='Output directory', type=str, default='', const='')
    parser.add_argument('-mzML_folder', nargs='?', help='Directory to search or to store .mzML files after convertation', type=str, default='', const='')
    parser.add_argument('-feature_folder', nargs='?', help='Path to folder where to store and search for biosaur2 features .tsv files', type=str, default='', const='')

    parser.add_argument('-decoy_prefix', nargs='?', help='String added to the protein name to showcase that it is a decoy (default: DECOY)', type=str, default='DECOY', const='DECOY')
    parser.add_argument('-sprot_suf', nargs='?', help='String added to the taxid of the organism\'s .fasta file to showcase swissprot database (default: _sp)', type=str, default='_sp', const='_sp')
    parser.add_argument('-uniprot_suf', nargs='?', help='String added to the taxid of the organism\'s .fasta file to showcase uniprot database (default: _sp)', type=str, default='_un', const='_un') 
    parser.add_argument('-min_pept_length', nargs='?', help='Minimal peptide lenght to use in cleaving 10% database, lowering the value is not recommended (default: 9)', type=int, default=9, const=9)
    parser.add_argument('-missed_cleavages', nargs='?', help='Number of missed cleavages for cleaving peptide in both searches (default: 0)', type=int, default=0, const=0)
    parser.add_argument('-cleavage_rule', help='cleavage rule in quotes!. X!Tandem style for cleavage rules: "[RK]|{P}" for trypsin, "[X]|[D]" for asp-n or "[RK]|{P},[K]|[X]" for mix of trypsin and lys-c', default='[RK]|{P}')
    parser.add_argument('-cmin', help='Min precursor charge', default=2, type=int, const=3)
    parser.add_argument('-cmax', help='Max precursor charge', default=3, type=int, const=3)
    parser.add_argument('-mass_accuracy', help='Mass accuracy in ppm for initial search', default=4, type=float, const=4)
    parser.add_argument('-mz_for_mass_accuracy', help='Approximate maximum m/z value to use with mass_accuracy', default=1000, type=float, const=1000)
    parser.add_argument('-allowed_ranks', nargs='?', help='Allowed taxonomy categories between <...> and <...>, default="strain,subspecies,forma specialis,isolate,serotype,serogroup,no rank" ', type=str, default='strain,subspecies,forma specialis,isolate,serotype,serogroup,no rank',)
    ##############
    parser.add_argument('-score_threshold', help='', default=4, type=float, const=4)
    parser.add_argument('-rewrite_individual_taxid_fasta', help='', default=0, type=int, const=0, choice=[0, 1])
    parser.add_argument('-rewrite_10per_fasta', help='', default=0, type=int, const=0, choice=[0, 1])
    parser.add_argument('-rewrite_blind_search', help='', default=0, type=int, const=0, choice=[0, 1])
    # self.aa_mass = args['aa_mass'] # см. парсинг в ms1searchpy
    
    # может в единый dict их засунуть?
    # self.path_to_uniprot_taxid_set = args['path_to_uniprot_taxid_set']
    # self.path_to_sprot_taxid_set = args['path_to_sprot_taxid_set']
    # self.path_to_len_fasta_uniprot = args['path_to_len_fasta_uniprot']
    # self.path_to_len_fasta_sprot = args['path_to_len_fasta_sprot']
    # self.path_to_species_descendants = args['path_to_species_descendants']
    # self.path_to_leaders_uniprot = args['path_to_leaders_uniprot']
    # self.path_to_leaders_sprot = args['path_to_leaders_sprot']
    # self.path_to_10_perc_fasta = args['path_to_10_perc_fasta']
    # self.path_to_output_prot_set = args['path_to_output_prot_set']
    # self.path_to_output_specmap_id = args['path_to_output_specmap_id']
    # self.path_to_output_cnt_to_spec = args['path_to_output_cnt_to_spec']
    # self.path_to_out_strain_statistics = args['path_to_out_strain_statistics']
    # self.path_to_top15_leaders_fasta = args['path_to_top15_leaders_fasta']
    ##############

    
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
        args['log_path'] = os.path.abspath(os.path.normpath(args['log_path']))
        log_directory = os.path.dirname(args['log_path'])
        os.makedirs(log_directory, exist_ok=True)
        fh = logging.FileHandler(args['log_path'], mode='w', encoding='utf-8')
        fh.setLevel(numeric_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    # logging.getLogger('matplotlib').setLevel(logging.ERROR)

    logger.info('Started')

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

if __name__ == '__main__':
    run()
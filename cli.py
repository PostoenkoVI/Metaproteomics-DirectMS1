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
    parser.add_argument('-uniprot_fold', nargs='*', help='Path to folder to store .fasta files splited by taxonomic identifiers from uniprot', type=str, default='', )
    parser.add_argument('-swissprot_fold', nargs='*', help='Path to folder to store .fasta files splited by taxonomic identifiers from swissprot', type=str, default='', )
    
    parser.add_argument('-outdir', nargs='?', help='name of directory to store results', type=str, default='', const='')
    # parser.add_argument('-feature_folder', nargs='?', help='directory to store features', type=str, default='', const='')

    parser.add_argument('-decoy_prefix', nargs='?', help='String added to the protein name to showcase that it is a decoy (default: DECOY)', type=str, default='DECOY', const='DECOY')

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
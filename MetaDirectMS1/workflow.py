from os import path, listdir, makedirs
import pickle
import logging
import pandas as pd
from ete3 import NCBITaxa
from pyteomics import fasta
import datetime
from .utils import log_subprocess_output, feature_generation, mzml_generation, call_ms1searchpy, call_DirectMS1quantmulti, call_ms1groups, call_ThermoRawFileParser, call_Biosaur2, get_aa_mass_with_fixed_mods, load, save, optimize_md, calc_sf_all, prepare_df, noisygaus, plot_identification_hist, plot_tax_barplot, plot_tax_boxplot, fc_plot, unite_fasta, Fasta_manipulations

logger = logging.getLogger(__name__)


def log_parsing(log_path:str, stage_sep:str='Stage', stage_id:str='') :
    lines = []
    with open(log_path, mode='r') as logfile :
        if stage_id == 'All' :
            for line in logfile.readlines() :
                lines.append(line)
        else :
            pattern1 = stage_sep + ' ' + stage_id
            pattern2 = stage_sep
            flag = False
            for line in logfile.readlines() :
                if pattern1 in line :
                    # print(line)
                    flag = True
                elif pattern2 in line and flag :
                    # print(line)
                    flag = False
                    break
                if flag :
                    lines.append(line)

    counter = {'WARNING':0, 'ERROR':0, 'CRITICAL':0}
    for line in lines :
        for key in counter.keys() :
            if key in line :
                counter[key] += 1
    return counter


def get_outdir(args) :
    outdir = args['outdir']
    if not outdir :
        date = str(datetime.datetime.today()).split()[0].replace('-', '_')
        drname = '_'.join(['metadirectms1', date])
        
        if args['input_files'] :
            f0 = args['input_files'].split(' ')[0]
            indir = path.dirname(f0)
        elif args['feature_folder'] :
            indir = args['feature_folder']
        elif args['mzml_folder'] :
            indir = args['mzml_folder']
        elif args['raw_folder'] :
            indir = args['raw_folder']
        else :
            indir = ''
        
        if indir :
            return path.join(indir, drname)
        else :
            return ''
    else :
        return outdir


def initialisation(args) :
    all_paths = dict()
    all_paths['outdir'] = args['outdir']
    rewrite = dict()
    # input
    input_type = False
    input_files = []
    input_dir = ''
    samples = []
    if args['input_files'] :
        
        if type(args['input_files']) == str :
            input_files = args['input_files'].split(' ')
        elif type(args['input_files']) == list :
            input_files = args['input_files']
        input_type = input_files[0].split('.')[-1].lower()
        if input_type == 'tsv' :
            input_type = 'feature'
        all_paths[input_type] = {}
    else :
        if args['feature_folder'] and path.exists(path.abspath(args['feature_folder'])) :
            logger.debug('Searching for features')
            input_type = 'feature'
            all_paths[input_type] = {}
            folder = path.abspath(args['feature_folder'])
            input_files = [path.join(folder, input_file) for input_file in listdir(folder) if input_file.lower().endswith('features.tsv')]
            if (len(input_files) == 0) :
                logger.info('No .tsv files found in %s', args['feature_folder'])
        if (len(input_files) == 0) and args['mzml_folder'] and path.exists(path.abspath(args['mzml_folder'])) :
            input_type = 'mzml'
            all_paths[input_type] = {}
            folder = path.abspath(args['mzml_folder'])
            input_files = [path.join(folder, input_file) for input_file in listdir(folder) if input_file.lower().endswith('.mzml')]
            if (len(input_files) == 0) :
                logger.info('No .mzml files found in %s', args['mzml_folder'])
        if (len(input_files) == 0) and args['raw_folder'] and path.exists(path.abspath(args['raw_folder'])) :
            input_type = 'raw'
            all_paths[input_type] = {}
            folder = path.abspath(args['raw_folder'])
            input_files = [path.join(folder, input_file) for input_file in listdir(folder) if input_file.lower().endswith('.raw')]
            if (len(input_files) == 0) :
                logger.warning('No .raw files found in %s', args['raw_folder'])
        if (len(input_files) == 0) and path.isdir(args['outdir']) :
            if path.isdir(path.join(args['outdir'], 'features')) :
                logger.debug('Searching for features %s', path.join(args['outdir'], 'features'))
                input_type = 'feature'
                all_paths[input_type] = {}
                folder = path.abspath(path.join(args['outdir'], 'features'))
                input_files = [path.join(folder, input_file) for input_file in listdir(folder) if input_file.lower().endswith('features.tsv')]
                if (len(input_files) == 0) :
                    logger.info('No .features.tsv files found in %s', path.join(args['outdir'], 'features'))
            if (len(input_files) == 0) and path.isdir(path.join(args['outdir'], 'mzml')) :
                input_type = 'mzml'
                all_paths[input_type] = {}
                folder = path.abspath(path.join(args['outdir'], 'mzml'))
                input_files = [path.join(folder, input_file) for input_file in listdir(folder) if input_file.lower().endswith('.mzml')]
                if (len(input_files) == 0) :
                    logger.info('No .mzml files found in %s', path.join(args['outdir'], 'features'))
            if (len(input_files) == 0) and path.isdir(path.join(args['outdir'], 'raw')) :
                input_type = 'raw'
                all_paths[input_type] = {}
                folder = path.abspath(path.join(args['outdir'], 'raw'))
                input_files = [path.join(folder, input_file) for input_file in listdir(folder) if input_file.lower().endswith('.raw')]
                if (len(input_files) == 0) :
                    logger.warning('No .raw files found in %s', path.join(args['outdir'], 'raw'))
    
    if len(input_files) > 0 :
        input_dir = path.dirname(input_files[0])
    else :
        logger.critical('No input files added. Either input files should be added through "-input_files" or at least one of the folders ("raw_folder", "mzml_folder", "feature_folder") should contain them.')
        return {}, {}, {}
    
    input_type_dct = {
        'feature' : 2,
        'mzml' : 1,
        'raw' : 0,
    }
    if input_type == 'feature' :
        all_paths['mzml_folder'] = ''
        all_paths['raw_folder'] = ''
        all_paths['feature_folder'] = input_dir
    elif input_type == 'mzml' :
        all_paths['raw_folder'] = ''
        all_paths['mzml_folder'] = input_dir
        if not args['feature_folder'] :
            all_paths['feature_folder'] = path.join(all_paths['outdir'], 'features')
        else :
            all_paths['feature_folder'] = args['feature_folder']
    elif input_type == 'raw' :
        all_paths['raw_folder'] = input_dir
        if not args['mzml_folder'] :
            all_paths['mzml_folder'] = path.join(all_paths['outdir'], 'mzml')
        else :
            all_paths['mzml_folder'] = args['mzml_folder']
        if not args['feature_folder'] :
            all_paths['feature_folder'] = path.join(all_paths['outdir'], 'features')
        else :
            all_paths['feature_folder'] = args['feature_folder']
    
    for t in ['raw', 'mzml', 'feature'] :
        if t != input_type :
            all_paths[t] = {}
    for input_file in input_files :
        if input_type != 'feature' :
            sample = path.basename(input_file.strip()).rsplit('.', maxsplit=1)[0]
            samples.append(sample)
        else :
            sample = path.basename(input_file.strip()).rsplit('.', maxsplit=2)[0]
            samples.append(sample)
        all_paths[input_type][sample] = input_file.strip()        
    if input_type == 'feature' :
        for sample in samples :
            all_paths['raw'][sample] = ''
            all_paths['mzml'][sample] = ''
    elif input_type == 'mzml' :
        for sample in samples :
            all_paths['raw'][sample] = ''
            all_paths['feature'][sample] = path.join(all_paths['feature_folder'], sample+'.features.tsv')
    elif input_type == 'raw' :
        for sample in samples :
            all_paths['mzml'][sample] = path.join(all_paths['mzml_folder'], sample+'.mzML')
            all_paths['feature'][sample] = path.join(all_paths['feature_folder'], sample+'.features.tsv')
    if all_paths['mzml_folder'] :
        makedirs(all_paths['mzml_folder'], exist_ok=True)
    if all_paths['feature_folder'] :
        makedirs(all_paths['feature_folder'], exist_ok=True)
    
    # structure
    if args['full_uniprot_fasta'] and path.isfile(path.abspath(args['full_uniprot_fasta'])) :
        all_paths['full_uniprot_fasta'] = args['full_uniprot_fasta']
    else :
        all_paths['full_uniprot_fasta'] = ''
        logger.warning('Path to full uniprot fasta is not stated or does not exists. It is critical error if database parsing stage wasn\'t completed earlier.')
    
    if args['db_parsing_folder'] :
        all_paths['db_parsing_folder'] = path.abspath(args['db_parsing_folder'])
    else :
        all_paths['db_parsing_folder'] = path.join(all_paths['outdir'], 'db_parsing_folder')
    makedirs(all_paths['db_parsing_folder'], exist_ok=True)
    all_paths['temp_folder'] = path.join(all_paths['db_parsing_folder'], 'tmp')
    makedirs(all_paths['temp_folder'], exist_ok=True)
    
    if not args['uniprot_folder'] :
        all_paths['uniprot_folder'] = path.join(all_paths['db_parsing_folder'], 'uniprot_fasta_fold')
    else :
        all_paths['uniprot_folder'] = args['uniprot_folder']
    if not args['sprot_folder'] :
        all_paths['sprot_folder'] = path.join(all_paths['db_parsing_folder'], 'sprot_fasta_fold')
    else :
        all_paths['sprot_folder'] = args['sprot_folder']
    makedirs(all_paths['uniprot_folder'], exist_ok=True)
    makedirs(all_paths['sprot_folder'], exist_ok=True)
    
    all_paths['10per_fasta'] = path.join(all_paths['db_parsing_folder'], '10per_fasta.fasta')
    
    # if not args['feature_folder'] :
    #     all_paths['feature_folder'] = path.join(all_paths['outdir'], 'features')
    # else :
    #     all_paths['feature_folder'] = args['feature_folder']
    # makedirs(all_paths['feature_folder'], exist_ok=True)
    
    all_paths['blind_search'] = path.join(all_paths['outdir'], 'blind_search')
    makedirs(all_paths['blind_search'], exist_ok=True)
    
    all_paths['top_leaders_search1'] = {}
    for sample in samples :
        all_paths['top_leaders_search1'][sample] = path.join(all_paths['blind_search'], sample+'_search1.fasta')

    all_paths['precise_search'] = path.join(all_paths['outdir'], 'precise_search')
    makedirs(all_paths['precise_search'], exist_ok=True)
    
    if args['precise_search_fasta'] :
        if path.isfile(args['precise_search_fasta']) :
            all_paths['search2_fasta'] = args['precise_search_fasta']
        else :
            logger.info('Path to .fasta file for precise search (precise_search_fasta) does not exists. Fasta will be generated automatically.')
            all_paths['search2_fasta'] = path.join(all_paths['precise_search'], 'precise_search.fasta')
    else :
        all_paths['search2_fasta'] = path.join(all_paths['precise_search'], 'precise_search.fasta')

    
    all_paths['quantitation'] = path.join(all_paths['outdir'], 'quantitation')
    makedirs(all_paths['quantitation'], exist_ok=True)
    
    all_paths['results'] = path.join(all_paths['outdir'], 'results')
    makedirs(all_paths['results'], exist_ok=True)
    
    if path.isfile(args['sample_file']) :
        all_paths['sample_file'] = path.abspath(args['sample_file'])
    else :
        logger.warning('Sample file does not exists: %s\tSearching in default location: %s', args['sample_file'], path.join(all_paths['quantitation'], 'sample_file.tsv'))
        all_paths['sample_file'] = path.join(all_paths['quantitation'], 'sample_file.tsv')
    all_paths['uniprot_taxid_set'] = path.join(all_paths['temp_folder'], 'uniprot_taxid_set.pickle')
    all_paths['sprot_taxid_set'] = path.join(all_paths['temp_folder'], 'sprot_taxid_set.pickle')
    all_paths['len_fasta_uniprot'] = path.join(all_paths['temp_folder'], 'len_fasta_uniprot.pickle')
    all_paths['len_fasta_sprot'] = path.join(all_paths['temp_folder'], 'len_fasta_sprot.pickle')
    all_paths['species_descendants'] = path.join(all_paths['temp_folder'], 'species_descendants.pickle')
    all_paths['leaders_uniprot'] = path.join(all_paths['temp_folder'], 'leaders_uniprot.pickle')
    all_paths['leaders_sprot'] = path.join(all_paths['temp_folder'], 'leaders_sprot.pickle')
    all_paths['exclude_names'] = path.join(all_paths['temp_folder'], 'exclude_names.pickle')
    all_paths['prot_set'] = path.join(all_paths['temp_folder'], 'prot_set.pickle')
    all_paths['specmap_id'] = path.join(all_paths['temp_folder'], 'specmap_id.pickle')
    all_paths['cnt_to_spec'] = path.join(all_paths['temp_folder'], 'cnt_to_spec.pickle')
    all_paths['protsN'] = path.join(all_paths['temp_folder'], 'protsN.pickle')
    all_paths['accurate_mz_map'] = path.join(all_paths['temp_folder'], 'accurate_mz_map.pickle')
    all_paths['spec_map_id_reversed'] = path.join(all_paths['temp_folder'], 'spec_map_id_reversed.pickle')
        
    groups = list(set(['OX', args['taxid_group'] ]))
    
    a = []
    b = [all_paths['mzml'][sample] for sample in samples]
    c = [all_paths['feature'][sample] for sample in samples]
    
    stage_result_dct = {
        'input':{
            'mzml' : a if input_type_dct[input_type] > 0 else b, 
            'feature' : c if input_type_dct[input_type] < 2 else a,
        }, 
        'db_parsing' : {
            'taxid_set' : [all_paths['uniprot_taxid_set'], all_paths['sprot_taxid_set']],
            'len_fasta' : [all_paths['len_fasta_uniprot'], all_paths['len_fasta_sprot']],
            'species_descendants' : [all_paths['species_descendants']],
            'leaders_uniprot' : [all_paths['leaders_uniprot']],
            'leaders_sprot' : [all_paths['leaders_sprot']],
            'exclude_names' : [all_paths['exclude_names']],
            '10per_fasta' : [all_paths['10per_fasta']],
        },
        'blind_search' : {
            'prot_set' : [all_paths['prot_set'], all_paths['specmap_id'], all_paths['cnt_to_spec'], 
                          all_paths['protsN'], all_paths['accurate_mz_map'], all_paths['spec_map_id_reversed']],
            'blind_search' : [all_paths['top_leaders_search1'][sample] for sample in samples],
            'ms1search' : [path.join(all_paths['blind_search'], sample+'.features_PFMs_ML.tsv') for sample in samples],
            'ms1groups' : [path.join(all_paths['blind_search'], sample+'_'+group+'.tsv') for sample in samples for group in groups],
            'blind_combined_results' : [path.join(all_paths['results'], 'blind_identified_proteins_'+group+'.tsv') for group in groups],
        },
        'precise_search' : {
            'search2_fasta' : [all_paths['search2_fasta']],
            # 'precise_search' : [path.join(all_paths['precise_search'], sample+'_'+group+'.tsv') for sample in samples for group in groups],
            'precise_search' : [path.join(all_paths['precise_search'], sample+'.features_PFMs_ML.tsv') for sample in samples],
            'ms1groups' : [path.join(all_paths['precise_search'], sample+'_'+group+'.tsv') for sample in samples for group in groups],
            'precise_combined_results' : [path.join(all_paths['results'], 'precise_identified_proteins_'+group+'.tsv') for group in groups],
        },
        'quantitation' : {
            'quantitation' : [path.join(all_paths['quantitation'], 'DQmulti_proteins_LFQ.tsv'),],
        },
    }
    # rewrite
    rewrite_dict = {
        'input' : {
            'mzml' : True,
            'feature' : True,
        },
        'db_parsing' : {
            'taxid_set' : True,
            'len_fasta' : True,
            'species_descendants' : True,
            'leaders_uniprot' : True,
            'leaders_sprot' : True,
            'exclude_names' : True,
            '10per_fasta' : True,
        },
        'blind_search' : {
            'prot_set' : True,
            'blind_search' : True,
            'ms1search' : True,
            'ms1groups' : True,
            'blind_combined_results' : True,
        },
        'precise_search' : {
            'search2_fasta' : True,
            'precise_search' : True,
            'ms1groups' : True,
            'precise_combined_results' : True,
        },
        'quantitation' : {
            'quantitation' : True,
        },
    }
    
    stages = ['db_parsing', 'blind_search', 'precise_search', 'quantitation', 'input']
    if args['mode'] == 0 :
        for key in stages :
            for k in rewrite_dict[key].keys() :                    
                if all([path.isfile(file) for file in stage_result_dct[key][k]]) :
                    rewrite_dict[key][k] = False
                else :
                    rewrite_dict[key][k] = True
        logger.debug('stage_result_dct = '+str(stage_result_dct))
        logger.debug('rewrite_dict = '+str(rewrite_dict))
    elif args['mode'] == 1 :
        for key in stages :
            for k in rewrite_dict[key].keys() :
                rewrite_dict[key][k] = True
    elif args['mode'] == 2 :
        exclude = ['db_parsing']
        for key in stages :
            if key in exclude : 
                for k in rewrite_dict[key].keys() :
                    rewrite_dict[key][k] = False
    elif args['mode'] == 3 :
        exclude = ['input', 'db_parsing']
        for key in stages :
            if key in exclude : 
                for k in rewrite_dict[key].keys() :
                    rewrite_dict[key][k] = False
    elif args['mode'] == 4 :
        exclude = ['input', 'db_parsing', 'blind_search']
        for key in stages :
            if key in exclude : 
                for k in rewrite_dict[key].keys() :
                    rewrite_dict[key][k] = False
    elif args['mode'] == 5 :
        exclude = ['input', 'db_parsing', 'blind_search', 'precise_search']
        for key in stages :
            if key in exclude : 
                for k in rewrite_dict[key].keys() :
                    rewrite_dict[key][k] = False
    
    if args['mode'] in [1, 2] :
        key = 'input'
        for k in rewrite_dict[key].keys() :
            rewrite_dict[key][k] = False if input_type_dct[input_type] >= input_type_dct[k] else True
        

    aa_mass, aa_to_psi = get_aa_mass_with_fixed_mods(args['fmods'], args['fmods_legend'])
        
    return all_paths, rewrite_dict, aa_mass


def process_files(args) :
    
    all_paths, rewrite_dict, aa_mass = initialisation(args)
    if not all_paths :
        return 1
    # logger.debug(all_paths)

    fmanip = Fasta_manipulations(all_paths, args)
    allowed_ranks = list(map(str.strip, args['allowed_ranks'].split(',') ) )
    if rewrite_dict['db_parsing']['taxid_set'] :
        logger.info('Generating %s.', 'uniprot_taxid_set')
        fmanip.prepare_uniprot_taxid_set(dump=True)
    else :
        try :
            logger.info('Trying to load %s.', 'uniprot_taxid_set')
            fmanip.uniprot_taxid_set = load(all_paths['uniprot_taxid_set'])
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'uniprot_taxid_set', all_paths['uniprot_taxid_set'])
    if rewrite_dict['db_parsing']['taxid_set'] :
        logger.info('Generating %s.', 'sprot_taxid_set')
        fmanip.prepare_sprot_taxid_set(dump=True)
    else :
        try :
            logger.info('Trying to load %s.', 'sprot_taxid_set')
            fmanip.sprot_taxid_set = load(all_paths['sprot_taxid_set'])
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'sprot_taxid_set', all_paths['sprot_taxid_set'])

    if rewrite_dict['db_parsing']['len_fasta'] :
        logger.info('Generating %s.', 'len_fasta_uniprot')
        fmanip.calc_lens_for_uniprot_taxid_set(dump=True)
    else :
        try :
            logger.info('Trying to load %s.', 'len_fasta_uniprot')
            fmanip.len_fasta_uniprot = load(all_paths['len_fasta_uniprot'])
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'len_fasta_uniprot', all_paths['len_fasta_uniprot'])
    if rewrite_dict['db_parsing']['len_fasta'] :
        logger.info('Generating %s.', 'len_fasta_sprot')
        fmanip.calc_lens_for_sprot_taxid_set(dump=True)
    else :
        try :
            logger.info('Trying to load %s.', 'len_fasta_sprot')
            fmanip.len_fasta_sprot = load(all_paths['len_fasta_sprot'])
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'len_fasta_sprot', all_paths['len_fasta_sprot'])
    if rewrite_dict['db_parsing']['species_descendants'] :
        logger.info('Generating %s.', 'species_descendants')
        fmanip.get_descendents_dict(allowed_ranks, dump=True)
    else :
        try :
            logger.info('Trying to load %s.', 'species_descendants')
            fmanip.species_descendants = load(all_paths['species_descendants'])
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'species_descendants', all_paths['species_descendants'])
    if rewrite_dict['db_parsing']['leaders_uniprot'] :
        logger.info('Generating %s.', 'leaders_uniprot')
        fmanip.get_leaders_uniprot(dump=True)
    else :
        try :
            logger.info('Trying to load %s.', 'leaders_uniprot')
            fmanip.species_leader_uniprot = load(all_paths['leaders_uniprot'])
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'leaders_uniprot', all_paths['leaders_uniprot'])
    if rewrite_dict['db_parsing']['leaders_sprot'] :
        logger.info('Generating %s.', 'leaders_sprot')
        fmanip.get_leaders_sprot(dump=True)
        fmanip.leaders = set(fmanip.species_leader_sprot.values()).union(set(fmanip.species_leader_uniprot.values()))
    else :
        try :
            logger.info('Trying to load %s.', 'leaders_sprot')
            fmanip.species_leader_sprot = load(all_paths['leaders_sprot'])
            fmanip.leaders = set(fmanip.species_leader_sprot.values()).union(set(fmanip.species_leader_uniprot.values()))
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'leaders_sprot', all_paths['leaders_sprot'])
    if rewrite_dict['db_parsing']['exclude_names'] :
        logger.info('Generating %s.', 'exclude_names')
        fmanip.exclude_wrong(dump=True)
    else :
        try :
            logger.info('Trying to load %s.', 'exclude_names')
            fmanip.exclude_names = load(all_paths['exclude_names'])
        except :
            logger.warning('Unable to load %s. Try mode 0 or 1. %s', 'exclude_names', all_paths['exclude_names'])
    
    if rewrite_dict['db_parsing']['10per_fasta'] :
        taxid_count = fmanip.write_10_perc_fasta(treshold=args['min_prot'])
        logger.info('Number of found taxid in 10%% fasta: %d', taxid_count)
    
    if rewrite_dict['input']['mzml'] :
        mzml_generation(all_paths['raw'], all_paths['mzml_folder'], path_to_mono=args['mono'], path_to_parser=args['thermorawparser'], str_of_other_args='')
    
    if rewrite_dict['input']['feature'] :
        feature_generation(all_paths['mzml'], all_paths['feature_folder'], path_to_fd=args['biosaur2'], str_of_other_args=args['bio2_args'])
    
    if rewrite_dict['blind_search']['prot_set'] :
        if path.isfile(all_paths['prot_set']) and path.isfile(all_paths['specmap_id']) and path.isfile(all_paths['cnt_to_spec']) :
            logger.info('Trying to load %s.', 'prot_set')
            fmanip.prot_sets = load(all_paths['prot_set'])
            fmanip.spec_map_id = load(all_paths['specmap_id'])
            fmanip.cnt_to_spec = load(all_paths['cnt_to_spec'])
        else :
            logger.info('Generating prot_set.')
            fmanip.prepare_protein_set(dump=True)
        fmanip.reverse_spec_map_id()
        fmanip.mz_map()
    
    if rewrite_dict['blind_search']['blind_search'] :
        logger.info('Starting blind search peptide matching.' )
        for sample, feature_path in all_paths['feature'].items() :
            df = prepare_df(feature_path)
            out_stat_path = path.join(all_paths['blind_search'], sample+'_strain_statistics.tsv')
            exitscore = fmanip.blind_search(df, path_to_out_fasta=all_paths['top_leaders_search1'][sample], path_to_out_strain_statistics=out_stat_path, exclude_sp_uncul=args['exclude_sp_uncul'])
            if exitscore != 0 :
                logger.critical('Something went wrong during blind search for file %s', feature_path)
                return 1
        logger.info('Blind search peptide matching is finished.' )
    
    if args['taxid_group'] != 'OX' :
        groups = ['OX', args['taxid_group'] ]
    else :
        groups = ['OX']
    logger.debug(groups)
    if rewrite_dict['blind_search']['ms1search'] :
        logger.info('Starting MS1 search for blind search.' )
        for sample, feature_path in all_paths['feature'].items() :
            call_ms1searchpy(args['ms1searchpy'], 
                             feature_path, 
                             all_paths['top_leaders_search1'][sample], 
                             outdir=all_paths['blind_search'],
                             cleavage_rule=args['cleavage_rule'], 
                             decoy_prefix=args['decoy_prefix'],
                             str_of_other_args=args['ms1searchpy_args'],
                             shuffle=True
                            )
        logger.info('MS1 search for blind search is finished.' )
    if rewrite_dict['blind_search']['ms1groups'] :
        logger.info('Starting to recalculate fdr for blind search.' )
        for sample, feature_path in all_paths['feature'].items() :
            for group in groups :
                PFMs_ML = path.join(all_paths['blind_search'], path.basename(all_paths['feature'][sample]).replace('.tsv', '_PFMs_ML.tsv'))
                # print(PFMs_ML)
                call_ms1groups(args['ms1searchpy'].replace('ms1searchpy', 'ms1groups'), 
                               PFMs_ML,
                               fasta_path=all_paths['top_leaders_search1'][sample].replace('.fasta', '_shuffled.fasta'), 
                               out=path.join(all_paths['blind_search'], sample+'_'), 
                               group=group, 
                               fdr=5, 
                               nproc=4, 
                               decoy_prefix=args['decoy_prefix'],
                              )
        logger.info('Recalculation of fdr for blind search is finished.' )
    # output tables merging
    if rewrite_dict['blind_search']['blind_combined_results'] :
        logger.info('Starting combining results for blind search.' )
        for group in groups :
            logger.debug('Combining blind search results for group %s', group)
            i = 0
            for sample, feature_path in all_paths['feature'].items() :
                p = path.join(all_paths['blind_search'], sample+'_'+group+'.tsv')
                tdf = pd.read_csv(p, sep='\t')
                tdf.rename(columns={'proteins':sample}, inplace=True)
                if i == 0 :
                    merge_df = tdf.copy()
                    i += 1
                else :
                    merge_df = merge_df.merge(tdf, on=['taxid', 'group'], how='outer')
            names_dct = NCBITaxa().get_taxid_translator(merge_df['taxid'].to_list())
            merge_df['name'] = merge_df['taxid'].apply(lambda x: names_dct[x])
            meaningfull_cols = ['group', 'taxid', 'name', 'mean identified proteins', 'median identified proteins',]
            samples = [col for col in all_paths['feature'].keys()]
            merge_df.fillna(0, inplace=True)
            merge_df['mean identified proteins'] = merge_df[samples].mean(axis=1)
            merge_df['median identified proteins'] = merge_df[samples].median(axis=1)
            if group == args['taxid_group'] :
                for sample in all_paths['feature'].keys() :
                    merge_df['needed_'+sample] = merge_df[sample]/merge_df[sample].sum() >= args['taxid_presence_thr']
                
                merge_df['include in combined fasta'] = merge_df[['needed_'+sample for sample in all_paths['feature'].keys()]].apply(any, axis=1)
                meaningfull_cols.append('include in combined fasta')
                merge_df.drop(columns=['needed_'+sample for sample in all_paths['feature'].keys()], inplace=True)

            if group == 'OX' :
                merge_df['len_fasta'] = merge_df['taxid'].apply(lambda x: fmanip.len_fasta_uniprot[x])
                meaningfull_cols += ['len_fasta']
            else :
                oxdf = pd.read_csv(path.join(all_paths['results'], 'blind_identified_proteins_OX.tsv'), sep='\t')
                dct = dict((group_taxid, []) for group_taxid in merge_df['taxid'])
                for group_taxid in merge_df['taxid'] :
                    for ox in oxdf['taxid'].values :
                        if group_taxid in NCBITaxa().get_lineage(ox) :
                            dct[group_taxid].append(fmanip.len_fasta_uniprot[ox])
                merge_df['len_fasta_by_ox'] = merge_df['taxid'].apply(lambda x: dct[x])
                merge_df['len_fasta_sum'] = merge_df['len_fasta_by_ox'].apply(sum)
                meaningfull_cols += ['len_fasta_by_ox', 'len_fasta_sum']
            cols = meaningfull_cols + samples
            merge_df = merge_df[cols].sort_values(by='mean identified proteins', ascending=False)
            merged_path = path.join(all_paths['results'], 'blind_identified_proteins_'+group+'.tsv')
            merge_df.to_csv(merged_path, sep='\t', index=False)
            merged_path = path.join(all_paths['blind_search'], 'blind_identified_proteins_'+group+'.tsv')
            merge_df.to_csv(merged_path, sep='\t', index=False)
        if args['generate_figures'] :
            logger.info('Generating figures for blind search combined results.' )
            p = path.join(all_paths['results'], 'blind_identified_proteins_OX.tsv')
            plot_identification_hist(p, search='blind')
            for group in groups :
                p = path.join(all_paths['results'], 'blind_identified_proteins_'+group+'.tsv')
                plot_tax_barplot(p, group=group, search='blind', ascending=True)
                plot_tax_boxplot(p, group=group, search='blind', ascending=True)
        logger.info('Combining results for blind search is finished.' )
        
    if rewrite_dict['precise_search']['precise_search'] :
        logger.info('Starting precise search.' )
        if rewrite_dict['precise_search']['search2_fasta'] and not path.isfile(all_paths['search2_fasta']) :
            logger.info('United .fasta database would be created: %s', all_paths['search2_fasta'])
            group = args['taxid_group']
            identification_table = path.join(all_paths['results'], 'blind_identified_proteins_'+group+'.tsv')
            unite_fasta(identification_table, 
                        all_paths['search2_fasta'],
                        all_paths['uniprot_folder'], 
                        threshold=args['taxid_presence_thr'], 
                        uniprot_suf=args['uniprot_suf'],
                       )
            f = fasta.read(all_paths['search2_fasta'])
            fasta_len = 0
            for _, _ in f :
                fasta_len += 1
            logger.info('Length of the united fasta: %d', fasta_len)
        
        shuffle_fasta = True
        for sample, feature_path in all_paths['feature'].items() :
            if shuffle_fasta :
                call_ms1searchpy(args['ms1searchpy'], 
                                 feature_path, 
                                 all_paths['search2_fasta'], 
                                 outdir=all_paths['precise_search'],
                                 cleavage_rule=args['cleavage_rule'], 
                                 decoy_prefix=args['decoy_prefix'],
                                 str_of_other_args=args['ms1searchpy_args'],
                                 shuffle=shuffle_fasta,
                                )
                shuffle_fasta = False
            else :
                call_ms1searchpy(args['ms1searchpy'], 
                                 feature_path, 
                                 all_paths['search2_fasta'].replace('.fasta', '_shuffled.fasta'), 
                                 outdir=all_paths['precise_search'],
                                 cleavage_rule=args['cleavage_rule'], 
                                 decoy_prefix=args['decoy_prefix'],
                                 str_of_other_args=args['ms1searchpy_args'],
                                )
        logger.info('Precise search is finished.' )
    if rewrite_dict['precise_search']['ms1groups'] :
        logger.info('Starting to recalculate fdr for precise search.' )
        for sample, feature_path in all_paths['feature'].items() :
            for group in groups :
                PFMs_ML = path.join(all_paths['precise_search'], path.basename(all_paths['feature'][sample]).replace('.tsv', '_PFMs_ML.tsv'))
                # print(PFMs_ML)
                call_ms1groups(args['ms1searchpy'].replace('ms1searchpy', 'ms1groups'), 
                               PFMs_ML,
                               fasta_path=all_paths['search2_fasta'].replace('.fasta', '_shuffled.fasta'), 
                               out=path.join(all_paths['precise_search'], sample+'_'), 
                               group=group, 
                               fdr=5, 
                               nproc=4, 
                               decoy_prefix=args['decoy_prefix'],
                              )
        logger.info('Recalculation of fdr for precise search is finished.' )
    # output tables merging
    if rewrite_dict['precise_search']['precise_combined_results'] :
        logger.info('Starting combining results for precise search.' )
        for group in groups :
            i = 0
            for sample, feature_path in all_paths['feature'].items() :
                p = path.join(all_paths['precise_search'], sample+'_'+group+'.tsv')
                tdf = pd.read_csv(p, sep='\t')
                tdf.rename(columns={'proteins':sample}, inplace=True)
                if i == 0 :
                    merge_df = tdf.copy()
                    i += 1
                else :
                    merge_df = merge_df.merge(tdf, on=['taxid', 'group'], how='outer')
            names_dct = NCBITaxa().get_taxid_translator(merge_df['taxid'].to_list())
            merge_df['name'] = merge_df['taxid'].apply(lambda x: names_dct[x])
            meaningfull_cols = ['group', 'taxid', 'name', 'mean identified proteins', 'median identified proteins',]
            samples = [col for col in all_paths['feature'].keys()]
            merge_df.fillna(0, inplace=True)
            merge_df['mean identified proteins'] = merge_df[samples].mean(axis=1)
            merge_df['median identified proteins'] = merge_df[samples].median(axis=1)
            cols = meaningfull_cols + samples
            merge_df = merge_df[cols].sort_values(by='mean identified proteins', ascending=False)
            
            merged_path = path.join(all_paths['results'], 'precise_identified_proteins_'+group+'.tsv')
            merge_df.to_csv(merged_path, sep='\t', index=False)
            merged_path = path.join(all_paths['precise_search'], 'precise_identified_proteins_'+group+'.tsv')
            merge_df.to_csv(merged_path, sep='\t', index=False)
        if args['generate_figures'] :
            logger.info('Generating figures for precise search combined results.' )
            p = path.join(all_paths['results'], 'precise_identified_proteins_OX.tsv')
            plot_identification_hist(p, search='precise')
            for group in groups :
                p = path.join(all_paths['results'], 'precise_identified_proteins_'+group+'.tsv')
                plot_tax_barplot(p, group=group, search='precise', ascending=True)
                plot_tax_boxplot(p, group=group, search='precise', ascending=True)
        logger.info('Combining results for precise search is finished.' )

    if rewrite_dict['quantitation']['quantitation'] :
        if path.isfile(all_paths['sample_file']) :
            logger.info('Starting quantitation stage.' )
            call_DirectMS1quantmulti(args['ms1searchpy'].replace('ms1searchpy', 'directms1quantmulti'),
                                     all_paths['precise_search'],
                                     fasta_path=all_paths['search2_fasta'].replace('.fasta', '_shuffled.fasta'),
                                     samples=all_paths['sample_file'],
                                     out_folder=all_paths['quantitation'],
                                     # min_signif_for_pept=999,
                                    )
            logger.info('Quantitation stage is finished.' )
        else :
            logger.warning('Path to sample file does not exists: %s', all_paths['sample_file'])
        if args['generate_figures'] :
            for file in listdir(all_paths['quantitation']) :
                if file.endswith('_quant_peptides.tsv') :
                    logger.info('Generating figures for file %s', file)
                    p = path.join(all_paths['quantitation'], file)
                    fc_plot(p, fasta=all_paths['search2_fasta'])
    else :
        logger.info('Quantitation stage is skipped.' )
    
    return 0
            
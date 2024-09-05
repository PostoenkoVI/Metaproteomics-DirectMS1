from os import path, listdir, makedirs
import pickle
import logging
from .utils import log_subprocess_output, feature_generation, mzml_generation, call_ms1searchpy, call_ThermoRawFileParser, call_Biosaur2, get_aa_mass_with_fixed_mods, load, save, optimize_md, calc_sf_all, prepare_df, noisygaus, Fasta_manipulations

logger = logging.getLogger(__name__)


def initialisation(args) :
    all_paths = dict()
    rewrite = dict()
    # input
    input_type = False
    input_files = []
    input_dir = ''
    samples = []
    if args['input_files'] :
        input_type = args['input_files'].split(' ')[0].split('.')[-1].lower()
        if input_type == 'tsv' :
            input_type = 'feature'
        all_paths[input_type] = {}
        input_files = args['input_files'].split(' ')
    else :
        if args['feature_folder'] and path.exists(path.abspath(args['feature_folder'])) :
            input_type = 'feature'
            all_paths[input_type] = {}
            folder = path.abspath(args['feature_folder'])
            input_files = [path.join(folder, input_file) for input_file in listdir(folder) if input_file.lower().endswith('_features.tsv')]
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
    
    if len(input_files) > 0 :
        input_dir = path.dirname(input_files[0])
    else :
        logger.critical('No input files added. Either input files should be added through "-input_files" or at least one of the folders ("raw_folder", "mzml_folder", "feature_folder") should contain them.')
        return {}, {}
            
    if args['outdir'] :
        all_paths['outdir'] = args['outdir']
        makedirs(all_paths['outdir'], exist_ok=True)
    else :
        all_paths['outdir'] = input_dir
        logger.warning('Path to output directory is not specified. Using input files directory instead.')
            
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
        sample = path.basename(input_file.strip()).split('.')[0]
        samples.append(sample)            
        all_paths[input_type][sample] = input_file.strip()
    if input_type == 'feature' :
        for sample in samples :
            all_paths['raw'][sample] = ''
            all_paths['mzml'][sample] = ''
    elif input_type == 'mzml' :
        for sample in samples :
            all_paths['raw'][sample] = ''
            all_paths['feature'][sample] = path.join(all_paths['feature_folder'], sample+'_features.tsv')
    elif input_type == 'raw' :
        for sample in samples :
            all_paths['mzml'][sample] = path.join(all_paths['mzml_folder'], sample+'.mzML')
            all_paths['feature'][sample] = path.join(all_paths['feature_folder'], sample+'_features.tsv')
    if all_paths['mzml_folder'] :
        makedirs(all_paths['mzml_folder'], exist_ok=True)
    if all_paths['feature_folder'] :
        makedirs(all_paths['feature_folder'], exist_ok=True)
    
    # structure
    if args['full_uniprot_fasta'] and path.exists(path.abspath(args['full_uniprot_fasta'])) :
        all_paths['full_uniprot_fasta'] = args['full_uniprot_fasta']
    else :
        logger.warning('Path to full uniprot fasta is not stated or does not exists.')
        
    if not args['uniprot_folder'] :
        all_paths['uniprot_folder'] = path.join(all_paths['outdir'], 'uniprot_fasta_fold')
    else :
        all_paths['uniprot_folder'] = args['uniprot_folder']
    if not args['sprot_folder'] :
        all_paths['sprot_folder'] = path.join(all_paths['outdir'], 'sprot_fasta_fold')
    else :
        all_paths['sprot_folder'] = args['sprot_folder']
    makedirs(all_paths['uniprot_folder'], exist_ok=True)
    makedirs(all_paths['sprot_folder'], exist_ok=True)
    
    if not args['feature_folder'] :
        all_paths['feature_folder'] = path.join(all_paths['outdir'], 'features')
    else :
        all_paths['feature_folder'] = args['feature_folder']
    makedirs(all_paths['feature_folder'], exist_ok=True)
    
    if not args['10per_fasta'] :
        all_paths['10per_fasta'] = path.join(all_paths['outdir'], '10per_fasta.fasta')
    else :
        all_paths['10per_fasta'] = args['10per_fasta']
    
    all_paths['top_leaders_fasta'] = {}
    
    if not args['temp_folder'] :
        all_paths['temp_folder'] = path.join(all_paths['outdir'], 'tmp')
    else :
        all_paths['temp_folder'] = args['temp_folder']
    makedirs(all_paths['temp_folder'], exist_ok=True)
    
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
    
    # rewrite
    rewrite_dict = {}
    if args['rewrite_individual_taxid_fasta'] :
        rewrite_dict['taxid_set'] = True
        rewrite_dict['len_fasta'] = True
    else :
        rewrite_dict['taxid_set'] = False
        rewrite_dict['len_fasta'] = False
    if args['rewrite_10per_fasta'] :
        rewrite_dict['species_descendants'] = True
        rewrite_dict['leaders_uniprot'] = True
        rewrite_dict['leaders_sprot'] = True
        rewrite_dict['exclude_names'] = True
        rewrite_dict['10per_fasta'] = True
        rewrite_dict['prot_set'] = True
    else :
        rewrite_dict['species_descendants'] = False
        rewrite_dict['leaders_uniprot'] = False
        rewrite_dict['leaders_sprot'] = False
        rewrite_dict['exclude_names'] = False
        rewrite_dict['10per_fasta'] = False
        rewrite_dict['prot_set'] = False
    if args['rewrite_blind_search'] :
        rewrite_dict['blind_search'] = True
    else :
        rewrite_dict['blind_search'] = False
    if args['rewrite_mzml'] :
        rewrite_dict['mzml'] = True
    else :
        rewrite_dict['mzml'] = False
    if args['rewrite_features'] :
        rewrite_dict['features'] = True
    else :
        rewrite_dict['features'] = False
    if args['rewrite_ms1search'] :
        rewrite_dict['ms1search'] = True
    else :
        rewrite_dict['ms1search'] = False
    aa_mass, aa_to_psi = get_aa_mass_with_fixed_mods(args['fmods'], args['fmods_legend'])
        
    return all_paths, rewrite_dict, aa_mass
    
    
def process_files(args) :
    
    all_paths, rewrite_dict, aa_mass = initialisation(args)
    logger.debug(all_paths)
    # parse_uniprot()
    fmanip = Fasta_manipulations(all_paths, args)
    
    if path.exists(all_paths['uniprot_taxid_set']) and (not rewrite_dict['taxid_set']) :
        fmanip.uniprot_taxid_set = load(all_paths['uniprot_taxid_set'])
    else :
        fmanip.prepare_uniprot_taxid_set(dump=True)
    if path.exists(all_paths['sprot_taxid_set']) and (not rewrite_dict['taxid_set']) :
        fmanip.sprot_taxid_set = load(all_paths['sprot_taxid_set'])
    else :
        fmanip.prepare_sprot_taxid_set(dump=True)
        
    if path.exists(all_paths['len_fasta_uniprot']) and (not rewrite_dict['len_fasta']) :
        fmanip.len_fasta_uniprot = load(all_paths['len_fasta_uniprot'])
    else :
        fmanip.calc_lens_for_uniprot_taxid_set(dump=True)
    if path.exists(all_paths['len_fasta_sprot']) and (not rewrite_dict['len_fasta']) :
        fmanip.len_fasta_sprot = load(all_paths['len_fasta_sprot'])
    else :
        fmanip.calc_lens_for_sprot_taxid_set(dump=True)
    
    allowed_ranks = list(map(str.strip, args['allowed_ranks'].split(',') ) )
    if path.exists(all_paths['species_descendants']) and (not rewrite_dict['species_descendants']) :
        fmanip.species_descendants = load(all_paths['species_descendants'])
    else :
        fmanip.get_descendents_dict(allowed_ranks, dump=True)

    if path.exists(all_paths['leaders_uniprot']) and (not rewrite_dict['leaders_uniprot']) :
        fmanip.species_leader_uniprot = load(all_paths['leaders_uniprot'])
    else :
        fmanip.get_leaders_uniprot(dump=True)
    if path.exists(all_paths['leaders_sprot']) and (not rewrite_dict['leaders_sprot']) :
        fmanip.species_leader_sprot = load(all_paths['leaders_sprot'])
    else :
        fmanip.get_leaders_sprot(dump=True)
    fmanip.leaders = set(fmanip.species_leader_sprot.values()).union(set(fmanip.species_leader_uniprot.values()))
    
    if path.exists(all_paths['exclude_names']) and (not rewrite_dict['exclude_names']) :
        fmanip.exclude_names = load(all_paths['exclude_names'])
    else :
        fmanip.exclude_wrong(dump=True)
    
    if rewrite_dict['10per_fasta'] or (not path.exists(all_paths['10per_fasta'])) :
        taxid_count = fmanip.write_10_perc_fasta(treshold=args['min_prot'])
        logger.debug(taxid_count)
    
    if rewrite_dict['mzml'] :
        mzml_generation(all_paths['raw'], all_paths['mzml_folder'], path_to_mono=args['mono'], path_to_parser=args['thermorawparser'], str_of_other_args='')
    
    if rewrite_dict['features'] :
        feature_generation(all_paths['mzml'], all_paths['feature_folder'], path_to_fd=args['biosaur2'], str_of_other_args='')
    
    if path.exists(all_paths['prot_set']) and path.exists(all_paths['specmap_id']) and path.exists(all_paths['cnt_to_spec']) and (not rewrite_dict['prot_set']) :
        fmanip.prot_sets = load(all_paths['prot_set'])
        fmanip.spec_map_id = load(all_paths['specmap_id'])
        fmanip.cnt_to_spec = load(all_paths['cnt_to_spec'])
    else :
        fmanip.prepare_protein_set(dump=True)
    fmanip.reverse_spec_map_id()
    fmanip.mz_map()
    
    for sample, feature_path in all_paths['feature'].items() :
        if rewrite_dict['blind_search'] or (not path.exists(path.join(all_paths['outdir'], sample+'_top_organisms.fasta') )) :
            df = prepare_df(feature_path)
            out_stat_path = path.join(all_paths['outdir'], sample+'_strain_statistics.tsv')
            out_fasta = path.join(all_paths['outdir'], sample+'_top_organisms.fasta')
            all_paths['top_leaders_fasta'][sample] = out_fasta
            fmanip.blind_search(df, path_to_out_fasta=out_fasta, path_to_out_strain_statistics=out_stat_path)
            
    if rewrite_dict['ms1search'] :
        for sample, feature_path in all_paths['feature'].items() :
            call_ms1searchpy(args['ms1searchpy'], 
                             feature_path, 
                             all_paths['top_leaders_fasta'][sample], 
                             path_to_deeplc=args['deeplc'],
                             cleavage_rule=args['cleavage_rule'], 
                             decoy_prefix=args['decoy_prefix'],
                             str_of_other_args=''
                            )
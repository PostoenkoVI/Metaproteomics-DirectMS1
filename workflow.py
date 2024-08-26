from os import path, listdir
import pickle
import logging
from .utils import parse_fasta_for_organisms, calc_lens_for_taxid_set, calc_lens_for_taxid_set, get_descendents_dict, get_leaders, exclude_wrong

# переписать под словари путей на входе
def feature_generation(mzml_paths, feature_folder, path_to_fd='', str_of_other_args='') :
    for mzml_path in mzml_paths :
        outpath = path.join(feature_folder, path.basename(mzml_path).replace('.mzML', '_features.tsv') )
        call_Biosaur2(path_to_fd, mzml_path, outpath, str_of_other_args=str_of_other_args) :
            
# переписать под словари путей на входе
def mzml_generation(raw_paths, mzml_folder, path_to_parser='', str_of_other_args='') :
    for raw_path in raw_paths :
        outpath = path.join(feature_folder, path.basename(raw_path).replace('.raw', '.mzML') )
        call_ThermoRawFileParser(path_to_parser, raw_path, outpath, str_of_other_args=str_of_other_args)
        
        
def initialisation(args) :
    all_paths = dict()
    rewrite = dict()
    
    # input
    input_type=False
    input_dir = ''
    samples = []
    if args['input_files'] :
        input_type = args['input_files'].split(' ')[0].split('.')[-1].lower()
        all_paths[input_type] = {}
        input_files = args['input_files'].split(' ')
    elif args['raw_folder'] :
        input_type = 'raw'
        all_paths[input_type] = {}
        folder = path.abspath(args['raw_folder'])
        input_files = [input_file for input_file in os.listdir(folder) if input_file.lower().endswith('.raw')]
    elif args['mzml_folder'] :
        input_type = 'mzml'
        all_paths[input_type] = {}
        folder = path.abspath(args['mzml_folder'])
        input_files = [input_file for input_file in os.listdir(folder) if input_file.lower().endswith('.mzml')]
    
    if input_files[0] :
        input_dir = os.path.dirname(input_files[0])
    else :
        logger.critical('No input files added. At least one of the arguments "input_files", "raw_folder", "mzml_folder" should be filled correctly.')
        return {}, {}
        
    for input_file in input_files :
        sample = path.basename(input_file.strip()).split('.')[0]
        samples.append(sample)            
        all_paths[input_type][sample] = input_file.strip()
    if input_type = 'mzml' :
        for sample in samples :
            all_paths['raw'][sample] = ''
    
    # structure
    if args['full_uniprot_fasta'] and path.exists(path.abspath(args['full_uniprot_fasta'])) :
        all_paths['full_uniprot_fasta'] = args['full_uniprot_fasta']
    else :
        logger.warning('Path to full uniprot fasta is not stated or does not exists.')

    if args['outdir'] :
        all_paths['outdir'] = args['outdir']
        os.makedirs(all_paths['outdir'], exists_ok=True)
    else :
        all_paths['outdir'] = input_dir
        logger.warning('Path to output directory is not specified. Using input files directory instead.')
        
    if not args['uniprot_folder'] :
        all_paths['uniprot_folder'] = path.join(all_paths['outdir'], 'uniprot_fasta_fold')
    else :
        all_paths['uniprot_folder'] = args['uniprot_folder']
    if not args['sprot_folder'] :
        all_paths['sprot_folder'] = path.join(all_paths['outdir'], 'sprot_fasta_fold')
    else :
        all_paths['sprot_folder'] = args['sprot_folder']
    os.makedirs(all_paths['uniprot_folder'], exists_ok=True)
    os.makedirs(all_paths['sprot_folder'], exists_ok=True)

    if not args['feature_folder'] :
        all_paths['feature_folder'] = path.join(all_paths['outdir'], 'features')
    else :
        all_paths['feature_folder'] = args['feature_folder']
    os.makedirs(all_paths['feature_folder'], exists_ok=True)
    
    path_to_top15_leaders_fasta
    
    if not args['temp_folder'] :
        all_paths['temp_folder'] = path.join(all_paths['outdir'], 'tmp')
    else :
        all_paths['temp_folder'] = args['temp_folder']
    os.makedirs(all_paths['temp_folder'], exists_ok=True)
    
    all_paths['uniprot_taxid_set'] = path.join(all_paths['temp_folder'], 'uniprot_taxid_set.pickle')
    all_paths['sprot_taxid_set'] = path.join(all_paths['temp_folder'], 'sprot_taxid_set.pickle')
    all_paths['len_fasta_uniprot'] = path.join(all_paths['temp_folder'], 'len_fasta_uniprot.pickle')
    all_paths['len_fasta_sprot'] = path.join(all_paths['temp_folder'], 'len_fasta_sprot.pickle')
    all_paths['species_descendants'] = path.join(all_paths['temp_folder'], 'species_descendants.pickle')
    all_paths['leaders_uniprot'] = path.join(all_paths['temp_folder'], 'leaders_uniprot.pickle')a
    all_paths['leaders_sprot'] = path.join(all_paths['temp_folder'], 'leaders_sprot.pickle')
    all_paths['prot_set'] = path.join(all_paths['temp_folder'], 'prot_set.pickle')
    all_paths['specmap_id'] = path.join(all_paths['temp_folder'], 'specmap_id.pickle')
    all_paths['cnt_to_spec'] = path.join(all_paths['temp_folder'], 'cnt_to_spec.pickle')
    all_paths['protsN'] = path.join(all_paths['temp_folder'], 'protsN.pickle')
    all_paths['accurate_mz_map'] = path.join(all_paths['temp_folder'], 'accurate_mz_map.pickle')
    all_paths['spec_map_id_reversed'] = path.join(all_paths['temp_folder'], 'spec_map_id_reversed.pickle')
    
    # rewrite
    rewrite_dict = {}
    
    
    return all_paths, rewrite    
    
    
def process_files(args) :
    
    all_paths, rewrite_dict = initialisation(args)
    
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
    
    if rewrite_10per_fasta or (not path.exists(all_paths['10_perc_fasta'])) :
        taxid_count = fmanip.write_10_perc_fasta(treshold=args['min_prot'])
        logger.debug(taxid_count)
    
    if write_mzml :
        mzml_generation(all_paths['raw'], all_paths['mzml_folder'], path_to_parser=all_paths['ThermoRawParser'], str_of_other_args='')
    
    if generate_features :
        feature_generation(all_paths['mzml'], all_paths['feature_folder'], path_to_fd=all_paths['Biosaur2'], str_of_other_args='')
    
    if rewrite_blind_search or (not path.exists(fmanip.path_to_top15_leaders_fasta)) :
        for feature_path, out_stat_path in zip(feature_paths, out_stat_paths) :
            df = pd.read_csv(feature_path, sep='\t')
            fmanip.blind_search(df, path_to_out_strain_statistics=out_stat_path)
        
    
    ms1searchpy()
    directMS1()
    generate_images()
from os import path, listdir
import pickle
import logging
from .utils import parse_fasta_for_organisms, calc_lens_for_taxid_set, calc_lens_for_taxid_set, get_descendents_dict, get_leaders, exclude_wrong


def feature_generation(mzml_paths, feature_folder, path_to_fd='', str_of_other_args='') :
    for mzml_path in mzml_paths :
        outpath = path.join(feature_folder, path.basename(mzml_path).replace('.mzML', '_features.tsv') )
        call_Biosaur2(path_to_fd, mzml_path, outpath, str_of_other_args=str_of_other_args) :
            
            
def mzml_generation(raw_paths, mzml_folder, path_to_parser='', str_of_other_args='') :
    for raw_path in raw_paths :
        outpath = path.join(feature_folder, path.basename(raw_path).replace('.raw', '.mzML') )
        call_ThermoRawFileParser(path_to_parser, raw_path, outpath, str_of_other_args=str_of_other_args)
        
        
def read_paths(args) :
    all_paths = dict()
    rewrite = dict()
    
    # input
    
    input_type=False
    samples = []
    if args['input_files'] :
        input_type = args['input_files'].split(' ')[0].split('.')[-1].lower()
        all_paths[input_type] = {}
        for input_file in args['input_files'].split(' ') :
            sample = path.basename(input_file.strip()).split('.')[0]
            samples.append(sample)            
            all_paths[input_type][sample] = input_file.strip()
        
    elif args['raw_folder'] :
        input_type = 'raw'
        all_paths[input_type] = {}
        folder = path.abspath(args['raw_folder'])
        for input_file in os.listdir(folder) :
            if input_file.lower().endswith('.raw') :
                sample = input_file.split('.')[0]
                samples.append(sample)
                all_paths[input_type][sample] = input_file
            
    elif args['mzml_folder'] :
        input_type = 'mzml'
        all_paths[input_type] = {}
        folder = path.abspath(args['mzml_folder'])
        for input_file in os.listdir(folder) :
            if input_file.lower().endswith('.mzml') :
                sample = input_file.split('.')[0]
                samples.append(sample)
                all_paths[input_type][sample] = input_file
    else :
        logger.critical('No input files added. At least one of the arguments "input_files", "raw_folder", "mzml_folder" should be filled.')
        return -1, -1
        
    if input_type = 'mzml' :
        for sample in samples :
            all_paths['raw'][sample] = ''
        
    
    return all_paths, rewrite    
    
    
def process_files(args) :
    
    , rewrite_dict = read_paths(args)
    
    # parse_uniprot()
    fmanip = Fasta_manipulations(args)
    
    if path.exists(fmanip.path_to_uniprot_taxid_set) and (not rewrite_dict['taxid_set']) :
        fmanip.uniprot_taxid_set = load(fmanip.path_to_uniprot_taxid_set)
    else :
        fmanip.prepare_uniprot_taxid_set(dump=True)
    if path.exists(fmanip.path_to_sprot_taxid_set) and (not rewrite_dict['taxid_set']) :
        fmanip.sprot_taxid_set = load(fmanip.path_to_sprot_taxid_set)
    else :
        fmanip.prepare_sprot_taxid_set(dump=True)
        
    if path.exists(fmanip.path_to_len_fasta_uniprot) and (not rewrite_dict['len_fasta']) :
        fmanip.len_fasta_uniprot = load(fmanip.path_to_len_fasta_uniprot)
    else :
        fmanip.calc_lens_for_uniprot_taxid_set(dump=True)
    if path.exists(fmanip.path_to_len_fasta_sprot) and (not rewrite_dict['len_fasta']) :
        fmanip.len_fasta_sprot = load(fmanip.path_to_len_fasta_sprot)
    else :
        fmanip.calc_lens_for_sprot_taxid_set(dump=True)
    
    if path.exists(fmanip.path_to_species_descendants) and (not rewrite_dict['species_descendants']) :
        fmanip.species_descendants = load(fmanip.path_to_species_descendants)
    else :
        fmanip.get_descendents_dict(allowed_ranks, dump=True)

    if path.exists(fmanip.path_to_leaders_uniprot) and (not rewrite_dict['leaders_uniprot']) :
        fmanip.species_leader_uniprot = load(fmanip.path_to_leaders_uniprot)
    else :
        fmanip.get_leaders_uniprot(dump=True)
    if path.exists(fmanip.path_to_leaders_sprot) and (not rewrite_dict['leaders_sprot']) :
        fmanip.species_leader_sprot = load(fmanip.path_to_leaders_sprot)
    else :
        fmanip.get_leaders_sprot(dump=True)
    fmanip.leaders = set(fmanip.species_leader_sprot.values()).union(set(fmanip.species_leader_uniprot.values()))
    
    if path.exists(fmanip.path_to_exclude_names) and (not rewrite_dict['exclude_names']) :
        fmanip.exclude_names = load(fmanip.path_to_exclude_names)
    else :
        fmanip.exclude_wrong(dump=True)
    taxid_count = fmanip.write_10_perc_fasta(treshold=200)
    
    mzml_generation()
    feature_generation()
    
    blind_search()
    
    ms1searchpy()
    directMS1()
    generate_images()
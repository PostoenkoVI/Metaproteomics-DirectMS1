from os import path
import pickle
import logging
from .utils import parse_fasta_for_organisms, calc_lens_for_taxid_set, calc_lens_for_taxid_set, get_descendents_dict, get_leaders, exclude_wrong

def forming_dicts(
    path_to_uniprot:str='/home/fasta/uniprot_bacteria.fasta',
    path_to_uniprot_dbs:str='/home/kae-13-1/fasta/bacts_bases_uniprot/',
    path_to_swissprot_dbs:str='/home/kae-13-1/fasta/bacts_bases_sprot/',
    
    path_to_len_fasta_uniprot:str='', 
    path_to_len_fasta_sprot:str='', 
    path_to_species_descendants:str='',
    path_to_leaders_uniprot:str='',
    path_to_leaders_sprot:str='',
    path_to_exclude_names:str='',
    
    allowed_ranks=('strain', 'subspecies', 'forma specialis', 'isolate', 'serotype', 'serogroup', 'no rank'),
    rewrite:bool=False,
) :
    
    if not path.exists(path_to_uniprot_dbs) :
        os.makedirs(path_to_uniprot_dbs, exist_ok=True)
    
    if not path.exists(path_to_swissprot_dbs) :
        os.makedirs(path_to_swissprot_dbs, exist_ok=True)
    
    uniprot_taxid_set, swissprot_taxid_set = parse_fasta_for_organisms(path_to_uniprot, 
                                                                       path_to_uniprot_dbs, 
                                                                       path_to_swissprot_dbs)
    
    if (not path.exists(path_to_len_fasta_uniprot)) or (rewrite) :
        len_fasta_uniprot = calc_lens_for_taxid_set(uniprot_taxid_set, 
                                                    path_to_uniprot_dbs, 
                                                    path_to_dump=path_to_len_fasta_uniprot)
    else :
        with open(path_to_len_fasta_uniprot, 'rb') as f :
            len_fasta_uniprot = pickle.load(f) 
    
    if (not path.exists(path_to_len_fasta_sprot)) or (rewrite) :
        len_fasta_sprot = calc_lens_for_taxid_set(swissprot_taxid_set, 
                                                  path_to_swissprot_dbs, 
                                                  path_to_dump=path_to_len_fasta_sprot, 
                                                  swissprot=True)
    else :
        with open(path_to_len_fasta_sprot, 'rb') as f :
            len_fasta_sprot = pickle.load(f) 
    
    if (not path.exists(path_to_species_descendants)) or (rewrite) :
        species_descendants = get_descendents_dict(uniprot_taxid_set, 
                                                   allowed_ranks, 
                                                   path_to_species_descendants=path_to_species_descendants) 
    else :
        with open(path_to_species_descendants, 'rb') as f :
            species_descendants = pickle.load(f)

    if (not path.exists(path_to_leaders_uniprot)) or (not path.exists(path_to_leaders_sprot)) or rewrite :
        species_leader_uniprot, species_leader_sprot = get_leaders(species_descendants, 
                                                                   swissprot_taxid_set, 
                                                                   len_fasta_uniprot, 
                                                                   len_fasta_sprot, 
                                                                   path_to_leaders_uniprot=path_to_leaders_uniprot, 
                                                                   path_to_leaders_sprot=path_to_leaders_sprot)
    else :
        with open(path_to_leaders_uniprot, 'rb') as f :
            species_leader_uniprot = pickle.load(f)
        with open(path_to_leaders_sprot, 'rb') as f :
            species_leader_sprot = pickle.load(f)
    
    if (not path.exists(path_to_exclude_names)) or (rewrite) :
        leaders = set(species_leader_sprot.values()).union(set(species_leader_uniprot.values()))
        exclude_names = exclude_wrong(leaders, path_to_exclude_names=path_to_exclude_names)
    else :
        with open(path_to_exclude_names, 'rb') as f :
            exclude_names = pickle.load(f)
            
    return species_leader_uniprot, species_leader_sprot, len_fasta_uniprot, len_fasta_sprot, exclude_names

###
# NOT WORKING YET
###
def blind_search(species_leader_uniprot, 
                 species_leader_sprot, 
                 len_fasta_uniprot, 
                 len_fasta_sprot, 
                 exclude_names,
                 
                 outfasta_path:str='',
                ) :
    
    feature_tsv_path = '/home/lab006/data/VGNKI/2024-04-12_microLC/features/B2P1_ms1_3ug.features.tsv'
    df1 = prepare_df(feature_tsv_path, charge_range=(2, 3), nIsotopes=3, mz_step=mz_step)
    df1.head()
    
    cnt, top_5_k, md_ar1, id_ar1 = get_matches(df1, [-mass_accuracy, mass_accuracy], accurate_mz_map, protsN, score_threshold=4.0)
    
    md_ar2 = []
    for z1, z2 in zip(md_ar1, id_ar1):
        if cnt_to_spec[z2-1] in top_5_k:
            md_ar2.append(z1)

    shift, sigma = optimize_md(md_ar2)
    print('Optimized mass shift and sigma: ', shift, sigma)

    custom_range_mass_accuracy = [shift-2*sigma, shift+2*sigma]
    
    cnt, _, md_ar1, id_ar1 = get_matches(df1, custom_range_mass_accuracy, accurate_mz_map, protsN, score_threshold=4.0)
    
    spec_map_id_reversed = dict()
    for k, v in spec_map_id.items():
        spec_map_id_reversed[v] = k
    
    top_species_names = get_top_species_names(cnt, 
                                              spec_map_id_reversed, 
                                              len_fasta_uniprot, 
                                              exclude_names, 
                                              number_of_top_proteins=15, 
                                              max_len_uniprot=220000
                                             )
    
    exitscore = write_top_species_fasta(top_species_names, 
                                        path_to_fasta_dir_sprot='/home/lerost/Metaproteomics-DirectMS1/Metaproteomics-DirectMS1_close_rep/testing_folder/test1/sprot_fasta_fold/', 
                                        path_to_fasta_dir_uniprot='/home/lerost/Metaproteomics-DirectMS1/Metaproteomics-DirectMS1_close_rep/testing_folder/test1/uniprot_fasta_fold/', 
                                        outfasta_path='/home/lerost/Metaproteomics-DirectMS1/Metaproteomics-DirectMS1_close_rep/testing_folder/test1/top_species.fasta', 
                                        out_strain_statistics_path='/home/lerost/Metaproteomics-DirectMS1/Metaproteomics-DirectMS1_close_rep/testing_folder/test1/B2P1_ms1_3ug.strain_statistics.csv',
                                       ) 
    print(exitscore)

def process_files(args) :
    ncbi = NCBITaxa()
    species_leader_uniprot, species_leader_sprot, len_fasta_uniprot, len_fasta_sprot, exclude_names = forming_dicts(
            path_to_uniprot='/home/fasta/uniprot_bacteria.fasta',
            path_to_uniprot_dbs='/home/kae-13-1/fasta/bacts_bases_uniprot/',
            path_to_swissprot_dbs='/home/kae-13-1/fasta/bacts_bases_sprot/',

            path_to_len_fasta_uniprot='',
            path_to_len_fasta_sprot='',
            path_to_species_descendants='',
            path_to_leaders_uniprot='',
            path_to_leaders_sprot='',
            path_to_exclude_names='',
        
            allowed_ranks=('strain', 'subspecies', 'forma specialis', 'isolate', 'serotype', 'serogroup', 'no rank'),
            rewrite=False,
        )
    
    if 
        exitscore = write_10_perc_fasta(outfasta_path='', 
                                        path_to_fasta_dir_uniprot='', 
                                        len_uniprot=len_fasta_uniprot, 
                                        len_sprot=len_fasta_sprot, 
                                        species_leader_uniprot=species_leader_uniprot,
                                        species_leader_sprot=species_leader_sprot, 
                                        treshold=200
                                       )
        logger.debug('write_10_perc_fasta : {}', exitscore)
    else :
        pass
    
    aa_mass = prepare_aa_mass()
    
    prot_sets,spec_map_id,cnt_to_spec = prepare_protein_set(
                                                            fasta_path='', 
                                                            decoy_prefix='DECOY_', 
                                                            cleave_rule='[RK]', 
                                                            min_length=9, 
                                                            missed_cleavages=0, 
                                                            aa_mass=aa_mass,
                                                            path_to_output_prot_set=path_to_output_prot_set, 
                                                            path_to_output_specmap_id=path_to_output_specmap_id, 
                                                            path_to_output_cnt_to_spec=path_to_output_cnt_to_spec,
                                                           )
    
    
    
    
    
    return 0
import os
import pickle
from collections import Counter, defaultdict
import operator
from ete3 import NCBITaxa
import pandas as pd
import numpy as np
from pyteomics import fasta, mass, parser, cmass, auxiliary as aux
from scipy.stats import binom
from copy import copy
from os import path, listdir, environ
import matplotlib.pyplot as plt
import random


logger = logging.getLogger(__name__)


class Fasta_manipulations :
    def __init__(self, args:dict) :
        self.path_to_uniprot = args['path_to_uniprot']
        self.path_to_uniprot_dbs = args['path_to_uniprot_dbs']
        self.path_to_sprot_dbs = args['path_to_sprot_dbs']
    # maybe add here possible load for intermediate dicts and sets if args['rewrite'] False (by default)
        
        self.sprot_suf = args['sprot_suf']
        self.uniprot_suf = args['uniprot_suf']
        self.aa_mass = args['aa_mass']
        self.min_pept_length = args['min_pept_length']
        self.missed_cleavages = args['missed_cleavages']
        self.cleave_rule = args['cleave_rule']
        self.mass_accuracy = args['mass_accuracy'] # (in ppm)
        self.mz_for_mass_accuracy = args['mz_for_mass_accuracy'] # (approximate max mz value)
        self.mz_step = self.mass_accuracy * 1e-6 * self.mz_for_mass_accuracy
        
        self.path_to_len_fasta_uniprot = args['path_to_len_fasta_uniprot']
        self.path_to_len_fasta_sprot = args['path_to_len_fasta_sprot']
        self.path_to_species_descendants = args['path_to_species_descendants']
        self.path_to_leaders_uniprot = args['path_to_leaders_uniprot']
        self.path_to_leaders_sprot = args['path_to_leaders_sprot']
        self.path_to_10_perc_fasta = args['path_to_10_perc_fasta']
        self.path_to_output_prot_set = args['path_to_output_prot_set']
        self.path_to_output_specmap_id = args['path_to_output_specmap_id']
        self.path_to_output_cnt_to_spec = args['path_to_output_cnt_to_spec']
        self.path_to_output_prot_set = args['path_to_output_prot_set']
        self.path_to_out_strain_statistics = args['path_to_out_strain_statistics']
        self.path_to_top15_leaders_fasta = args['path_to_top15_leaders_fasta']
    
    def prepare_uniprot_taxid_set(self) :
        uniprot_taxid_set = set()
        for p in fasta.read(path_to_uniprot):
            spec_i = p[0].split('OX=')[-1].split(' ')[0]
            fasta.write([(p[0], p[1])], output = path.join(self.path_to_uniprot_dbs, '{}'+self.uniprot_suf+'.fasta'.format(spec_i)),
                            file_mode = 'a')
            if spec_i not in uniprot_taxid_set:
                uniprot_taxid_set.update([int(spec_i)])
        self.uniprot_taxid_set = uniprot_taxid_set

    def prepare_sprot_taxid_set(self) :
        sprot_taxid_set = set()
        for p in fasta.read(path_to_uniprot):
            if p[0].startswith('sp'):
                spec_i = p[0].split('OX=')[-1].split(' ')[0]
                # fasta.write([(p[0], p[1])], output = path.join(path_to_swissprot_dbs, '{}.fasta'.format(spec_i)),
                fasta.write([(p[0], p[1])], output = path.join(self.path_to_sprot_dbs, '{}'+self.sprot_suf+'.fasta'.format(spec_i)),
                                file_mode = 'a')
                if spec_i not in sprot_taxid_set:
                    sprot_taxid_set.update([int(spec_i)])
        self.sprot_taxid_set = sprot_taxid_set

    def calc_lens_for_sprot_taxid_set(self, dump:bool=False) :
        taxid_lens_dict = {}
        for i in self.sprot_taxid_set :
            filename = '{}'+self.sprot_suf+'.fasta'.format(i)
            file = path.join(self.path_to_sprot_dbs, filename)
            n = sum(1 for _ in fasta.read(file))
            taxid_lens_dict[i] = n
        if dump :
            if path.exists(path.abspath(self.path_to_len_fasta_sprot)) :
            with open(self.path_to_len_fasta_sprot, 'wb') as f :
                pickle.dump(
                    taxid_lens_dict, 
                    f, 
                    protocol=pickle.HIGHEST_PROTOCOL
                )
        self.len_fasta_sprot = taxid_lens_dict

    def calc_lens_for_uniprot_taxid_set(self, dump:bool=False) :
        taxid_lens_dict = {}
        for i in self.uniprot_taxid_set :
            filename = '{}'+self.uniprot_suf+'.fasta'.format(i)
            file = path.join(self.path_to_uniprot_dbs, filename)
            n = sum(1 for _ in fasta.read(file))
            taxid_lens_dict[i] = n
        if dump :
            if path.exists(path.abspath(self.path_to_len_fasta_uniprot)) :
            with open(self.path_to_len_fasta_uniprot, 'wb') as f :
                pickle.dump(
                    taxid_lens_dict, 
                    f, 
                    protocol=pickle.HIGHEST_PROTOCOL
                )
        self.len_fasta_uniprot = taxid_lens_dict

    def get_descendents_dict(self, allowed_ranks:set, dump:bool=False) :
        ###
        # Описание!!!
        ###
        species_descendants = defaultdict(set)
        used = set()

        for i in taxid_set:
            if i not in used:
                rank = NCBITaxa().get_rank([i])
                if rank:
                    if rank[i] == 'species':
                        descendants = set(NCBITaxa().get_descendant_taxa(i) + [i])
        #                 descendants = [j for j in descendants if j in taxid_set]
                        species_descendants[i].update(descendants.intersection(taxid_set))
                        used.update(species_descendants[i])
                    elif rank[i] in allowed_ranks:
                        lineage = NCBITaxa().get_lineage(i)
                        ranks = NCBITaxa().get_rank(lineage)
                        species = [k for k in ranks.keys() if ranks[k] == 'species'][0]
                        descendants = set(NCBITaxa().get_descendant_taxa(species) + [species])
        #                 descendants = [j for j in descendants if j in taxid_set]
                        species_descendants[species].update(descendants.intersection(taxid_set))
                        species_descendants[species].add(i)
                        used.update(species_descendants[species])
        used = set()
        if dump :
            if path.exists(path.abspath(self.path_to_species_descendants)) :
                with open(self.path_to_species_descendants, 'wb') as f :
                    pickle.dump(
                        species_descendants, 
                        f, 
                        protocol=pickle.HIGHEST_PROTOCOL
                )
        self.species_descendants = species_descendants
    
    def get_leaders_uniprot(self, dump:bool=False) :
        ###
        # Описание!!!
        ###
        species_leader_uniprot = {}
        for i in self.species_descendants.keys() :
            strains = self.species_descendants[i]
            lens = {j:self.len_fasta_uniprot[j] for j in strains}
            if len(lens) != 0 :
                lead = max(lens.items(), key=operator.itemgetter(1))[0]
                species_leader_uniprot[i] = lead
        if dump :
            if path.exists(path.abspath(self.path_to_leaders_uniprot)) :
                with open(self.path_to_leaders_uniprot, 'wb') as f :
                    pickle.dump(
                        species_leader_uniprot, 
                        f, 
                        protocol=pickle.HIGHEST_PROTOCOL
                    )
        self.species_leader_uniprot = species_leader_uniprot
        
    def get_leaders_sprot(self, dump:bool=False) :
        ###
        # Описание!!!
        ###
        species_leader_sprot = {}
        for i in self.species_descendants.keys() :
            strains = self.species_descendants[i]
            strains_sp = [i for i in strains if i in sprot_taxid_set]
            lens_sp = {j:self.len_fasta_sprot[j] for j in strains_sp}
            if len(lens_sp) != 0 :
                lead = max(lens_sp.items(), key=operator.itemgetter(1))[0]
                species_leader_sprot[i] = lead
        if dump :
            if path.exists(path.abspath(self.path_to_leaders_sprot)) :
                with open(self.path_to_leaders_sprot, 'wb') as f :
                    pickle.dump(
                        species_leader_sprot, 
                        f, 
                        protocol=pickle.HIGHEST_PROTOCOL
                    )
        self.species_leader_sprot = species_leader_sprot
        
    def exclude_wrong(self, dump:bool=False) :
        ###
        # Описание!!!
        ###
        i = 0
        exclude_names = set()
        leaders = set(self.species_leader_sprot.values()).union(set(self.species_leader_uniprot.values()))

        for k in leaders:
            name = list(NCBITaxa().get_taxid_translator([k]).values())[0]
            if 'sp.' in name:
                if name.split(' ')[1] == 'sp.':
                    exclude_names.update([k])
                    i+=1
            if name.startswith('uncultured'):
                exclude_names.update([k])
        if dump :
            if path.exists(path.abspath(self.path_to_exclude_names)) :
                with open(self.path_to_exclude_names, 'wb') as f :
                    pickle.dump(
                        exclude_names, 
                        f, 
                        protocol=pickle.HIGHEST_PROTOCOL
                    )
        self.leaders = leaders
        self.exclude_names = exclude_names
        
    def write_10_perc_fasta(self, treshold:int=200) :
        random_dict = {}
        for k in self.len_fasta_uniprot.keys():
            random_dict[k] = 2000 / self.len_fasta_uniprot[k]

        outf = open(self.path_to_10_perc_fasta, 'w')
        outf.close()
        c = 0
        for f in self.leaders :
            prots = [] 
            if f in self.species_leader_sprot.values():
                n_prot = self.len_fasta_sprot[f]
            else: 
                n_prot = self.len_fasta_sprot[f]
            if n_prot >= treshold :
                for p in fasta.read(path.join(self.path_to_uniprot_dbs, str(f) + self.uniprot_suf + '.fasta')):
                    if p[0][:2] == 'sp' or random_dict[f] >= random.random():
                        prots.append(p)
                with open(self.path_to_10_perc_fasta, 'a') as f :
                    fasta.write(prots, output=f)
                    c += 1
        return c
    
    def prepare_protein_set(self, dump=False, ) :        
        prot_sets = defaultdict(list)
        spec_map_id = dict()
        cnt_to_spec = list()
        spec_map_id_max = 0

        cnt = 0
        for p in fasta.read(self.path_to_10_perc_fasta) :

            if self.decoy_prefix not in p[0]:
                cnt += 1
                if cnt % 1000000 == 0:
                    logger.debug('Proteins scored: '+str(cnt))

                spec = p[0].split('OX=')[-1].split(' ')[0]
                if spec not in spec_map_id:
                    spec_map_id_max += 1
                    spec_map_id[spec] = spec_map_id_max

                spec_id = spec_map_id[spec]
                cnt_to_spec.append(spec_id)

                peptides = parser.cleave(p[1], self.cleave_rule, self.missed_cleavages, min_length=self.min_pept_length)
                mz_list = []

                dont_use_fast_valid = parser.fast_valid(p[1])
                for pep in peptides:
                    plen = len(pep) 
                    if plen <= 15:
                        if dont_use_fast_valid or parser.fast_valid(pep) :
                            mz_list.append(cmass.fast_mass(pep, aa_mass=self.aa_mass))
                for mz in set(mz_list):
                    prot_sets[mz].append(cnt)
        self.prot_sets = prot_sets
        self.spec_map_id = spec_map_id
        self.cnt_to_spec = cnt_to_spec
        
        if dump :
            if path.exists(path.abspath(self.path_to_output_prot_set) :
                with open(self.path_to_output_prot_set, 'wb') as f :
                    pickle.dump(prot_sets, 
                                f, protocol=pickle.HIGHEST_PROTOCOL)
                           
            if path.exists(path.abspath(self.path_to_output_specmap_id) :
                with open(self.path_to_output_specmap_id, 'wb') as f :
                    pickle.dump(prot_sets, 
                                f, protocol=pickle.HIGHEST_PROTOCOL)
                           
            if path.exists(path.abspath(self.path_to_output_cnt_to_spec) :
                with open(self.path_to_output_cnt_to_spec, 'wb') as f :
                    pickle.dump(prot_sets, 
                                f, protocol=pickle.HIGHEST_PROTOCOL)
                           
    def mz_map(self) :
        protsN = Counter()
        accurate_mz_map = defaultdict(list)

        for v, prots in self.prot_sets.items():
            v_int = int(v/mz_step)
            accurate_mz_map[v_int].append(v)
            protsN.update(prots)
        self.protsN = protsN
        self.accurate_mz_map = accurate_mz_map

    def reverse_spec_map_id(self) :
        spec_map_id_reversed = dict()
        for k, v in self.spec_map_id.items():
            spec_map_id_reversed[v] = k
        self.spec_map_id_reversed = spec_map_id_reversed
    
    def blind_search(self, df:pandas.DataFrame, ) :
        # self.score_threshold
        cnt, top_5_k, md_ar1, id_ar1 = self.get_matches(df, [-self.mass_accuracy, self.mass_accuracy], score_threshold=self.score_threshold)
        md_ar2 = []
        for z1, z2 in zip(md_ar1, id_ar1):
            if self.cnt_to_spec[z2-1] in top_5_k:
                md_ar2.append(z1)
        
        shift, sigma = optimize_md(md_ar2, bin_width=0.1) ######################## bin_width не менять?
        custom_range_mass_accuracy = [shift-2*sigma, shift+2*sigma]
        cnt, _, md_ar1, id_ar1 = self.get_matches(df, custom_range_mass_accuracy, score_threshold=self.score_threshold)
        
        number_of_top_proteins = self.number_of_top_proteins
        top_100_species_names = set()
        for k, v in cnt.most_common():
            if len(top_100_species_names) < number_of_top_proteins:
                k_orig = self.spec_map_id_reversed[k]
                if self.len_fasta_uniprot[int(k_orig)] < 220000:
                    if not int(k_orig) in self.exclude_names: 
                        top_100_species_names.add(k_orig) 
                        orig_name = ncbi.get_taxid_translator([k_orig,])[int(k_orig)]
                        # print(k, k_orig, orig_name, v)
        prots = []
        cnt = 0    
        report = pd.DataFrame()
        for leader in top_100_species_names :
            SP = 0
            UN = 0
            sp_filename = str(leader)+self.sprot_suf+'.fasta'
            if filename in listdir(self.path_to_sprot_dbs):
                for p in fasta.read(path.join(self.path_to_sprot_dbs, filename)):
                    SP+=1
            un_filename = str(leader)+self.uniprot_suf+'.fasta'
            for p in fasta.read(path.join(self.path_to_uniprot_dbs, unip_filename)):
                prots.append(p)
                UN+=1
            report = pd.concat([report, pd.DataFrame.from_dict({'ID':[leader],
                                                               'Sprot':[SP],
                                                               'Uniprot':[UN]})])

        if self.path_to_out_strain_statistics :
            report.to_csv(self.path_to_out_strain_statistics, index=False)
        random.shuffle(prots)
        with open(self.path_to_top15_leaders_fasta, 'w') as f :
            fasta.write(prots, output=f)
                           
    def get_matches(self, 
                    df1:pd.DataFrame,
                    custom_range_mass_accuracy:list,
                    score_threshold:float=self.score_threshold
                    ) :

        prots_spc = Counter()
        md_ar1 = []
        id_ar1 = []

        nmasses = df1['massCalib'].values
        charges = df1['charge'].values
        nmasses_int = df1['massCalib_int'].values
        nmasses_int_dict = defaultdict(set)
        for idx, nm in enumerate(nmasses_int) :
            nmasses_int_dict[nm].add(idx)
            nmasses_int_dict[nm-1].add(idx)
            nmasses_int_dict[nm+1].add(idx)

        mz_acc_checked = set()

        for mz_int in self.accurate_mz_map :
            if mz_int in nmasses_int_dict :
                for mz_acc in accurate_mz_map[mz_int] :
                    if mz_acc not in mz_acc_checked :
                        mz_acc_checked.add(mz_acc)
                        for idx_nm in nmasses_int_dict[mz_int] :
                            nm_val = nmasses[idx_nm]
                            mass_diff_ppm = (mz_acc - nm_val) / mz_acc * 1e6
    #                         if abs(mass_diff_ppm) <= mass_accuracy:
                            if custom_range_mass_accuracy[0] <= mass_diff_ppm <= custom_range_mass_accuracy[1] :
    #                                 md_ar1.append(mass_diff_ppm)
                                prots_spc.update(self.prot_sets[mz_acc])
                                for prot in self.prot_sets[mz_acc] :
                                    md_ar1.append(mass_diff_ppm)
                                    id_ar1.append(prot)
                                break

        top100decoy_score = [prots_spc.get(dprot, 0) for dprot in self.protsN]
        top100decoy_N = [val for key, val in self.protsN.items()]
        p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
        logger.debug('p=%s', (np.mean(top100decoy_score) / np.mean(top100decoy_N)))

        names_arr = np.array(list(prots_spc.keys()))
        v_arr = np.array([prots_spc[k] for k in names_arr])
        n_arr = np.array([self.protsN[k] for k in names_arr])

        prots_spc2 = dict()
        all_pvals = calc_sf_all(v_arr, n_arr, p)
        for idx, k in enumerate(names_arr):
            prots_spc2[k] = all_pvals[idx]

        cnt = Counter()
        for k, v in prots_spc2.items():
            if v >= self.score_threshold:
                sp_id = cnt_to_spec[k-1]
                cnt[sp_id] += 1

        top_5_k = set()
        for k, v in cnt.most_common():
            if len(top_5_k) < 5:
                top_5_k.add(k)

        return cnt, top_5_k, md_ar1, id_ar1


def save(obj, file ) :
    if file :
        with open(file, 'wb') as f :
            pickle.dump(obj, 
                        f, 
                        protocol=pickle.HIGHEST_PROTOCOL)
        return 0
    else :
        return 1


def load(file) :
    with open(file, 'rb') as f:
        obj = pickle.load(f)
    return obj


def calc_sf_all(v:int, n:int, p:float) :
    ###
    # calculates survival function for binomial discrete random variable
    # v <= n + 1
    # 0 <= p <= 1
    ###
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[np.isinf(sf_values)] = max(sf_values[~np.isinf(sf_values)])
    return sf_values


def get_matches(
    df1:pd.DataFrame,
    custom_range_mass_accuracy:list,
    accurate_mz_map:dict,
    protsN:collections.Counter,
    score_threshold=4.0
    ) :
    
    prots_spc = Counter()
    md_ar1 = []
    id_ar1 = []

    nmasses = df1['massCalib'].values
    charges = df1['charge'].values
    nmasses_int = df1['massCalib_int'].values
    nmasses_int_dict = defaultdict(set)
    for idx, nm in enumerate(nmasses_int) :
        nmasses_int_dict[nm].add(idx)
        nmasses_int_dict[nm-1].add(idx)
        nmasses_int_dict[nm+1].add(idx)

    mz_acc_checked = set()

    for mz_int in accurate_mz_map :
        if mz_int in nmasses_int_dict :
            for mz_acc in accurate_mz_map[mz_int] :
                if mz_acc not in mz_acc_checked :
                    mz_acc_checked.add(mz_acc)
                    for idx_nm in nmasses_int_dict[mz_int] :
                        nm_val = nmasses[idx_nm]
                        mass_diff_ppm = (mz_acc - nm_val) / mz_acc * 1e6
#                         if abs(mass_diff_ppm) <= mass_accuracy:
                        if custom_range_mass_accuracy[0] <= mass_diff_ppm <= custom_range_mass_accuracy[1] :
#                                 md_ar1.append(mass_diff_ppm)
                            prots_spc.update(prot_sets[mz_acc])
                            for prot in prot_sets[mz_acc] :
                                md_ar1.append(mass_diff_ppm)
                                id_ar1.append(prot)
                            break

    top100decoy_score = [prots_spc.get(dprot, 0) for dprot in protsN]
    top100decoy_N = [val for key, val in protsN.items()]
    p = np.mean(top100decoy_score) / np.mean(top100decoy_N)
    logger.debug('p=%s', (np.mean(top100decoy_score) / np.mean(top100decoy_N)))

    names_arr = np.array(list(prots_spc.keys()))
    v_arr = np.array([prots_spc[k] for k in names_arr])
    n_arr = np.array([protsN[k] for k in names_arr])

    prots_spc2 = dict()
    all_pvals = calc_sf_all(v_arr, n_arr, p)
    for idx, k in enumerate(names_arr):
        prots_spc2[k] = all_pvals[idx]

    cnt = Counter()
    for k, v in prots_spc2.items():
        if v >= score_threshold:
            sp_id = cnt_to_spec[k-1]
            cnt[sp_id] += 1

    top_5_k = set()
    for k, v in cnt.most_common():
        if len(top_5_k) < 5:
            top_5_k.add(k)

    return cnt, top_5_k, md_ar1, id_ar1


def call_Biosaur2(path_to_fd, mzml_path, outpath, str_of_other_args='') :
    wrong_add_args = [
                      '-o',
                     ]
    if str_of_other_args :
        for wrong_add_arg in wrong_add_args :
            if str_of_other_args.find(wrong_add_arg)>=0 :
                msg = 'Ignoring {} parameter from additional args to Biosaur2'.format(wrong_add_arg)
                logger.warning(msg)
                str_of_other_args = str_of_other_args.split(wrong_add_arg)[0] + ' '.join(str_of_other_args.split(wrong_add_arg)[-1].strip().split()[1:])
        other_args = [x.strip() for x in str_of_other_args.split(' ')]
    else :
        other_args = []
    final_args = [path_to_fd, mzml_path, '-o', outpath, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore


def call_ThermoRawFileParser(path_to_parser, raw_file, outpath, str_of_other_args='') :
    wrong_add_args = ['-L', '--msLevel=', 
                      '-o', '--output_directory=',
                      '-d', '--input_directory=', 
                      '-b', '--output=',
                     ]
    if str_of_other_args :
        for wrong_add_arg in wrong_add_args :
            if str_of_other_args.find(wrong_add_arg)>=0 :
                msg = 'Ignoring {} parameter from additional args to ThermoRawFileParser'.format(wrong_add_arg)
                logger.warning(msg)
                str_of_other_args = str_of_other_args.split(wrong_add_arg)[0] + ' '.join(str_of_other_args.split(wrong_add_arg)[-1].strip().split()[1:])
        other_args = [x.strip() for x in str_of_other_args.split(' ')]
    else :
        other_args = []
    other_args += ['-L', 1]
    final_args = [path_to_fd, '-i', raw_file, '-o', outpath, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore


def prepare_df(feature_tsv_path:str, charge_range:tuple=(2, 3), nIsotopes:int=3, mz_step:float=0.004) :
    df1 = pd.read_table(feature_tsv_path)
    df1 = df1[df1['nIsotopes'] >= nIsotopes]
    df1 = df1[df1['charge'] >= charge_range[0]]
    df1 = df1[df1['charge'] <= charge_range[1]]
    df1['massCalib_int'] = df1['massCalib'].apply(lambda x: int(x/mz_step))
    return df1


def noisygaus(x, a, x0, sigma, b) :
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b
    
    
def optimize_md(md_ar2:np.ndarray, bin_width=0.1) :
    bbins = np.arange(np.min(md_ar2), np.max(md_ar2), bin_width)
    H2, b2 = np.histogram(md_ar2, bins=bbins)
    m, mi, s = np.max(H2), b2[np.argmax(H2)], (np.max(md_ar2) - np.min(md_ar2))/6
    noise = np.min(H2)
    logger.debug('p0: {}', [m, mi, s, noise])
    popt, pcov = curve_fit(noisygaus, b2[1:], H2, p0=[m, mi, s, noise])
    shift, sigma = popt[1], abs(popt[2])
    logger.debug('Optimized mass shift and sigma: {}, {}', shift, sigma)
    return shift, sigma


def get_top_species_names(cnt:dict, 
                          spec_map_id_reversed:dict, 
                          len_uniprot:dict, 
                          exclude_names:set, 
                          number_of_top_proteins:int=15, 
                          max_len_uniprot:int=220000
                         ) :
    top_species_names = set()
    for k, v in cnt.most_common():
        if len(top_species_names) < number_of_top_proteins:
            k_orig = spec_map_id_reversed[k]
            if len_uniprot[int(k_orig)] < max_len_uniprot :
                if not int(k_orig) in exclude_names: 
                    top_species_names.add(k_orig) # OR k_orig???
                    orig_name = ncbi.get_taxid_translator([k_orig,])[int(k_orig)]
                    logger.debug(k, k_orig, orig_name, v)
    return top_species_names


def write_top_species_fasta(top_species_names:set, 
                            path_to_fasta_dir_sprot:str, 
                            path_to_fasta_dir_uniprot:str, 
                            outfasta_path:str, 
                            out_strain_statistics_path:str=''
                           ) :
    prots = []
    report = pd.DataFrame()
    for leader in top_species_names :
        SP = 0
        UN = 0

        if str(leader)+'.fasta' in listdir(path_to_fasta_dir_sprot):
            for p in fasta.read(path.join(path_to_fasta_dir_sprot, str(leader)+'_sp.fasta')):
                SP+=1
        for p in fasta.read(path.join(path_to_fasta_dir_uniprot, str(leader)+'.fasta')):
            prots.append(p)
            UN+=1
        report = pd.concat([report, pd.DataFrame.from_dict({'ID':[leader],
                                                           'Sprot':[SP],
                                                           'Uniprot':[UN]})])

    if out_strain_statistics_path :
        report.to_csv(out_strain_statistics_path, index=False)
    random.shuffle(prots)
    with open(out_strain_statistics_path, 'w') as f :
        fasta.write(prots, output=f)
    return 0


def call_ms1searchpy(path_to_ms1searchpy:str, feat_path:str, fasta_path:str, cleavage_rule:str='', decoy_prefix:str='') :
    wrong_add_args = ['-d', '-deeplc', 
                      '-e', '-ad', 
                      '-prefix', '-ml',
                      '-ts', 
                     ]
    if str_of_other_args :
        for wrong_add_arg in wrong_add_args :
            if str_of_other_args.find(wrong_add_arg)>=0 :
                msg = 'Ignoring {} parameter from additional args to ms1searchpy'.format(wrong_add_arg)
                logger.info(msg)
                str_of_other_args = str_of_other_args.split(wrong_add_arg)[0] + ' '.join(str_of_other_args.split(wrong_add_arg)[-1].strip().split()[1:])
        other_args = [x.strip() for x in str_of_other_args.split(' ')]
    else :
        other_args = []
    other_args = ['-d', fasta_path, '-e', cleavage_rule, '-prefix', decoy_prefix, '-ad', 1, '-deeplc', 1, '-ml', 1, '-ts', 2] + other_args
    final_args = [path_to_ms1searchpy, feat_path, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore


def log_subprocess_output(pipe):
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        logger.info('From subprocess: %s', line.decode(sys.stdout.encoding).rstrip('\n'))


def parse_fasta_for_organisms(path_to_uniprot:str, path_to_uniprot_dbs:str, path_to_swissprot_dbs:str) :
    ###
    # Описание!!!
    ###
    uniprot_taxid_set = set()
    swissprot_taxid_set = set()
    for p in fasta.read(path_to_uniprot):
        spec_i = p[0].split('OX=')[-1].split(' ')[0]
        fasta.write([(p[0], p[1])], output = path.join(path_to_uniprot_dbs, '{}.fasta'.format(spec_i)),
                        file_mode = 'a')
        if spec_i not in uniprot_taxid_set:
            uniprot_taxid_set.update([int(spec_i)])

        if p[0].startswith('sp'):
            # fasta.write([(p[0], p[1])], output = path.join(path_to_swissprot_dbs, '{}.fasta'.format(spec_i)),
            fasta.write([(p[0], p[1])], output = path.join(path_to_swissprot_dbs, '{}_sp.fasta'.format(spec_i)),
                            file_mode = 'a')
            if spec_i not in swissprot_taxid_set:
                swissprot_taxid_set.update([int(spec_i)])
        
    return uniprot_taxid_set, swissprot_taxid_set


def calc_lens_for_taxid_set(taxid_set:set, path_to_dbs:str, path_to_dump='', swissprot=False) :
    ###
    # Описание!!!
    ###
    taxid_lens_dict = {}
    for i in taxid_set :
        if swissprot :
            filename = '{}_sp.fasta'.format(i)
        else :
            filename = '{}.fasta'.format(i)
        file = path.join(path_to_dbs, filename)
        # logger.debug(file)

        n = sum(1 for _ in fasta.read(file))
        taxid_lens_dict[i] = n
    
    if path_to_dump :
        with open(path_to_dump, 'wb') as f :
            pickle.dump(
                taxid_lens_dict, 
                f, 
                protocol=pickle.HIGHEST_PROTOCOL
            )
    return taxid_lens_dict


def get_descendents_dict(taxid_set:set, allowed_ranks:set, path_to_species_descendants='') :
    ###
    # Описание!!!
    ###
    species_descendants = defaultdict(set)
    used = set()

    for i in taxid_set:
        if i not in used:
            rank = NCBITaxa().get_rank([i])
            if rank:
                if rank[i] == 'species':
                    descendants = set(NCBITaxa().get_descendant_taxa(i) + [i])
    #                 descendants = [j for j in descendants if j in taxid_set]
                    species_descendants[i].update(descendants.intersection(taxid_set))
                    used.update(species_descendants[i])
                elif rank[i] in allowed_ranks:
                    lineage = NCBITaxa().get_lineage(i)
                    ranks = NCBITaxa().get_rank(lineage)
                    species = [k for k in ranks.keys() if ranks[k] == 'species'][0]
                    descendants = set(NCBITaxa().get_descendant_taxa(species) + [species])
    #                 descendants = [j for j in descendants if j in taxid_set]
                    species_descendants[species].update(descendants.intersection(taxid_set))
                    species_descendants[species].add(i)
                    used.update(species_descendants[species])
    
    if path_to_species_descendants :
        with open(path_to_species_descendants, 'wb') as f :
            pickle.dump(
                species_descendants, 
                f, 
                protocol=pickle.HIGHEST_PROTOCOL
            )
    return species_descendants


def get_leaders(species_descendants:dict, sprot_taxid_set:set, len_fasta_uniprot:dict, len_fasta_sprot:dict, path_to_leaders_uniprot='', path_to_leaders_sprot='') :
    ###
    # Описание!!!
    ###
    species_leader_sprot = {}
    species_leader_uniprot = {}
    for i in species_descendants.keys() :
        strains = species_descendants[i]
        strains_sp = [i for i in strains if i in sprot_taxid_set]
        lens = {j:len_fasta_uniprot[j] for j in strains}
        lens_sp = {j:len_fasta_sprot[j] for j in strains_sp}
        if len(lens) != 0 :
            lead = max(lens.items(), key=operator.itemgetter(1))[0]
            species_leader_uniprot[i] = lead
        if len(lens_sp) != 0 :
            lead = max(lens_sp.items(), key=operator.itemgetter(1))[0]
            species_leader_sprot[i] = lead
            
    if path_to_leaders_uniprot :
        with open(path_to_leaders_uniprot, 'wb') as f :
            pickle.dump(
                species_leader_uniprot, 
                f, 
                protocol=pickle.HIGHEST_PROTOCOL
            )
    if path_to_leaders_sprot :
        with open(path_to_leaders_sprot, 'wb') as f :
            pickle.dump(
                species_leader_sprot, 
                f, 
                protocol=pickle.HIGHEST_PROTOCOL
            )
    return species_leader_uniprot, species_leader_sprot


def exclude_wrong(leaders:set, path_to_exclude_names='') :
    ###
    # Описание!!!
    ###
    i = 0
    exclude_names = set()
    # leaders = set(species_leader_sprot.values()).union(set(species_leader_uniprot.values()))

    for k in leaders:
        name = list(NCBITaxa().get_taxid_translator([k]).values())[0]
        if 'sp.' in name:
            if name.split(' ')[1] == 'sp.':
                exclude_names.update([k])
                i+=1
        if name.startswith('uncultured'):
            exclude_names.update([k])

    if path_to_exclude_names :
        with open(path_to_exclude_names, 'wb') as f :
            pickle.dump(
                exclude_names, 
                f, 
                protocol=pickle.HIGHEST_PROTOCOL
            )
    return exclude_names
                           
                           
def write_10_perc_fasta(outfasta_path:str='', 
                        path_to_fasta_dir_uniprot:str='', 
                        len_uniprot:dict={}, 
                        len_sprot:dict={}, 
                        species_leader_sprot:set=set(), 
                        species_leader_uniprot:set=set(), 
                        treshold:int=200
                       ) :
    random_dict = {}
    for k in len_uniprot.keys():
        random_dict[k] = 2000 / len_uniprot[k]

    outf = open(outfasta_path, 'w')
    outf.close()

    leaders = set(species_leader_sprot.values()).union(set(species_leader_uniprot.values()))
    for f in leaders:
        prots = [] 
        if f in species_leader_sprot.values():
            n_prot = len_sprot[f]
        else: 
            n_prot = len_uniprot[f]
        if n_prot >= treshold :
            for p in fasta.read(path.join(path_to_fasta_dir_uniprot, str(f) + '.fasta')):
                if p[0][:2] == 'sp' or random_dict[f] >= random.random():
                    prots.append(p)
            with open(outfasta_path, 'a') as f :
                fasta.write(prots, output=f)
    return 0
            
            
def prepare_protein_set(fasta_path:str, 
                        decoy_prefix:str='DECOY_', 
                        cleave_rule:str='[RK]', 
                        min_length:int=9, 
                        missed_cleavages:int=0, 
                        aa_mass:dict=copy(mass.std_aa_mass),
                        path_to_output_prot_set:str='', 
                        path_to_output_specmap_id:str='', 
                        path_to_output_cnt_to_spec:str='',
                        rewrite=True
                       ) :
    # Add fixed modifications here
    # aa_mass = copy(mass.std_aa_mass)
    # aa_mass['C'] += 57.021464
    # aa_mass
    
    if path.exists(path_to_output_prot_set) and path.exists(path_to_output_specmap_id) and path.exists(path_to_output_cnt_to_spec) and (not rewrite) :
        with open(path_to_output_prot_set, 'wb') as f :
            prot_sets = pickle.load(f)
        with open(path_to_output_spec_map_id, 'wb') as f :
            spec_map_id = pickle.load(f)
        with open(path_to_output_cnt_to_spec, 'wb') as f :
            cnt_to_spec = pickle.load(f)
        return prot_sets, spec_map_id, cnt_to_spec
    
    else :
        prot_sets = defaultdict(list)
        spec_map_id = dict()
        cnt_to_spec = list()
        spec_map_id_max = 0

        cnt = 0
        for p in fasta.read(fasta_path):   

            if decoy_prefix not in p[0]:
                cnt += 1
                if cnt % 1000000 == 0:
                    logger.debug('Proteins scored: '+str(cnt))

                spec = p[0].split('OX=')[-1].split(' ')[0]
                if spec not in spec_map_id:
                    spec_map_id_max += 1
                    spec_map_id[spec] = spec_map_id_max

                spec_id = spec_map_id[spec]
                cnt_to_spec.append(spec_id)

                peptides = parser.cleave(p[1], cleave_rule, missed_cleavages, min_length=min_length)
                mz_list = []

                dont_use_fast_valid = parser.fast_valid(p[1])
                for pep in peptides:
                    plen = len(pep) 
                    if plen <= 15:
                        if dont_use_fast_valid or parser.fast_valid(pep) :
                            mz_list.append(cmass.fast_mass(pep, aa_mass=aa_mass))
                for mz in set(mz_list):
                    prot_sets[mz].append(cnt)

        if path_to_output_prot_set :
            with open(path_to_output_prot_set, 'wb') as f :
                pickle.dump(prot_sets, 
                            f, protocol=pickle.HIGHEST_PROTOCOL)

        if path_to_output_specmap_id :
            with open(path_to_output_specmap_id, 'wb') as f :
                pickle.dump(spec_map_id, 
                            f, protocol=pickle.HIGHEST_PROTOCOL)

        if path_to_output_cnt_to_spec :
            with open(path_to_output_cnt_to_spec, 'wb') as f :
                pickle.dump(cnt_to_spec, 
                            f, protocol=pickle.HIGHEST_PROTOCOL)

        return prot_sets, spec_map_id, cnt_to_spec


def mz_map(
    prot_sets:dict
    ) :
    protsN = Counter()
    accurate_mz_map = defaultdict(list)

    for v, prots in prot_sets.items():
        v_int = int(v/mz_step)
        accurate_mz_map[v_int].append(v)
        protsN.update(prots)
    return protsN, accurate_mz_map
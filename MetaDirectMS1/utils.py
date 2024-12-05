import pickle
from collections import Counter, defaultdict
import operator
from ete3 import NCBITaxa
import pandas as pd
import numpy as np
from pyteomics import fasta, mass, parser, cmass, auxiliary as aux
from scipy.stats import binom
from scipy.optimize import curve_fit
from copy import copy, deepcopy
from os import path, listdir

import matplotlib.pyplot as plt
import random
import logging
import subprocess
import sys


logger = logging.getLogger(__name__)


class Fasta_manipulations :
    
    def __init__(self, all_paths:dict, args:dict) :
        self.path_to_uniprot = all_paths['full_uniprot_fasta']
        self.path_to_uniprot_dbs = all_paths['uniprot_folder']
        self.path_to_sprot_dbs = all_paths['sprot_folder']
    # maybe add here possible load for intermediate dicts and sets if args['rewrite'] False (by default)
        
        self.sprot_suf = args['sprot_suf']
        self.uniprot_suf = args['uniprot_suf']
        self.decoy_prefix = args['decoy_prefix']
        aa_mass, aa_to_psi = get_aa_mass_with_fixed_mods(args['fmods'], args['fmods_legend'])
        self.aa_mass = aa_mass
        self.min_pept_length = args['min_pept_length']
        self.missed_cleavages = args['missed_cleavages']
        self.cleave_rule = args['cleavage_rule']
        self.mass_accuracy = args['mass_accuracy'] # (in ppm)
        self.mz_for_mass_accuracy = args['mz_for_mass_accuracy'] # (approximate max mz value)
        self.mz_step = self.mass_accuracy * 1e-6 * self.mz_for_mass_accuracy
        self.score_threshold = args['score_threshold']
        self.num_top_spec = args['num_top_spec']
        
        self.path_to_uniprot_taxid_set = all_paths['uniprot_taxid_set']
        self.path_to_sprot_taxid_set = all_paths['sprot_taxid_set'] 
        self.path_to_len_fasta_uniprot = all_paths['len_fasta_uniprot']
        self.path_to_len_fasta_sprot = all_paths['len_fasta_sprot']
        self.path_to_species_descendants = all_paths['species_descendants']
        self.path_to_leaders_uniprot = all_paths['leaders_uniprot']
        self.path_to_leaders_sprot = all_paths['leaders_sprot']
        self.path_to_exclude_names = all_paths['exclude_names']
        self.path_to_10_perc_fasta = all_paths['10per_fasta']
        self.path_to_output_prot_set = all_paths['prot_set']
        self.path_to_output_specmap_id = all_paths['specmap_id']
        self.path_to_output_cnt_to_spec = all_paths['cnt_to_spec']
        self.path_to_protsN = all_paths['protsN']
        self.path_to_accurate_mz_map = all_paths['accurate_mz_map']
        self.path_to_spec_map_id_reversed = all_paths['spec_map_id_reversed']
    
    def prepare_uniprot_taxid_set(self, dump:bool=False) :
        uniprot_taxid_set = set()
        for p in fasta.read(self.path_to_uniprot):
            spec_i = p[0].split('OX=')[-1].split(' ')[0]
            if int(spec_i) not in uniprot_taxid_set:
                uniprot_taxid_set.update([int(spec_i)])
                fasta.write([(p[0], p[1])], output = path.join(self.path_to_uniprot_dbs, '{}{}.fasta'.format(spec_i, self.uniprot_suf)),
                                file_mode = 'w')
            else :
                fasta.write([(p[0], p[1])], output = path.join(self.path_to_uniprot_dbs, '{}{}.fasta'.format(spec_i, self.uniprot_suf)),
                                file_mode = 'a')
        if dump :
            with open(self.path_to_uniprot_taxid_set, 'wb') as f :
                pickle.dump(
                    uniprot_taxid_set, 
                    f, 
                    protocol=pickle.HIGHEST_PROTOCOL
                )
        self.uniprot_taxid_set = uniprot_taxid_set

    def prepare_sprot_taxid_set(self, dump:bool=False) :
        sprot_taxid_set = set()
        for p in fasta.read(self.path_to_uniprot):
            if p[0].startswith('sp'):
                spec_i = p[0].split('OX=')[-1].split(' ')[0]
                if int(spec_i) not in sprot_taxid_set:
                    sprot_taxid_set.update([int(spec_i)])
                # fasta.write([(p[0], p[1])], output = path.join(path_to_swissprot_dbs, '{}.fasta'.format(spec_i)),
                    fasta.write([(p[0], p[1])], output = path.join(self.path_to_sprot_dbs, '{}{}.fasta'.format(spec_i, self.sprot_suf)),
                                    file_mode = 'w')
                else :
                    fasta.write([(p[0], p[1])], output = path.join(self.path_to_sprot_dbs, '{}{}.fasta'.format(spec_i, self.sprot_suf)),
                                    file_mode = 'a')
                
        if dump :
            with open(self.path_to_sprot_taxid_set, 'wb') as f :
                pickle.dump(
                    sprot_taxid_set, 
                    f, 
                    protocol=pickle.HIGHEST_PROTOCOL
                )
        self.sprot_taxid_set = sprot_taxid_set
        

    def calc_lens_for_sprot_taxid_set(self, dump:bool=False) :
        taxid_lens_dict = {}
        for i in self.sprot_taxid_set :
            filename = '{}{}.fasta'.format(i, self.sprot_suf)
            file = path.join(self.path_to_sprot_dbs, filename)
            n = sum(1 for _ in fasta.read(file))
            taxid_lens_dict[i] = n
        if dump :
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
            filename = '{}{}.fasta'.format(i, self.uniprot_suf)
            file = path.join(self.path_to_uniprot_dbs, filename)
            n = sum(1 for _ in fasta.read(file))
            taxid_lens_dict[i] = n
        if dump :
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

        for i in self.uniprot_taxid_set:
            if i not in used:
                rank = NCBITaxa().get_rank([i])
                if rank:
                    if rank[i] == 'species':
                        descendants = set(NCBITaxa().get_descendant_taxa(i) + [i])
                        species_descendants[i].update(descendants.intersection(self.uniprot_taxid_set))
                        used.update(species_descendants[i])
                    elif rank[i] in allowed_ranks:
                        lineage = NCBITaxa().get_lineage(i)
                        ranks = NCBITaxa().get_rank(lineage)
                        species = [k for k in ranks.keys() if ranks[k] == 'species'][0]
                        descendants = set(NCBITaxa().get_descendant_taxa(species) + [species])
                        species_descendants[species].update(descendants.intersection(self.uniprot_taxid_set))
                        species_descendants[species].add(i)
                        used.update(species_descendants[species])
        used = set()
        if dump :
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
            strains_sp = [i for i in strains if i in self.sprot_taxid_set]
            lens_sp = {j:self.len_fasta_sprot[j] for j in strains_sp}
            if len(lens_sp) != 0 :
                lead = max(lens_sp.items(), key=operator.itemgetter(1))[0]
                species_leader_sprot[i] = lead
        if dump :
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
                n_prot = self.len_fasta_uniprot[f]
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
            with open(self.path_to_output_prot_set, 'wb') as f :
                pickle.dump(prot_sets, 
                            f, protocol=pickle.HIGHEST_PROTOCOL)

            with open(self.path_to_output_specmap_id, 'wb') as f :
                pickle.dump(spec_map_id, 
                            f, protocol=pickle.HIGHEST_PROTOCOL)

            with open(self.path_to_output_cnt_to_spec, 'wb') as f :
                pickle.dump(cnt_to_spec, 
                            f, protocol=pickle.HIGHEST_PROTOCOL)
                           
    def mz_map(self) :
        protsN = Counter()
        accurate_mz_map = defaultdict(list)

        for v, prots in self.prot_sets.items():
            v_int = int(v/self.mz_step)
            accurate_mz_map[v_int].append(v)
            protsN.update(prots)
        self.protsN = protsN
        self.accurate_mz_map = accurate_mz_map

    def reverse_spec_map_id(self) :
        spec_map_id_reversed = dict()
        for k, v in self.spec_map_id.items():
            spec_map_id_reversed[v] = k
        self.spec_map_id_reversed = spec_map_id_reversed
    
    def blind_search(self, df:pd.DataFrame, path_to_out_fasta='', path_to_out_strain_statistics='') :
#         using under hood 
#         self.mass_accuracy
#         self.score_threshold
#         self.cnt_to_spec
#         self.num_top_spec
#         self.spec_map_id_reversed
#         self.len_fasta_uniprot
#         self.exclude_names
#         self.sprot_suf
#         self.path_to_sprot_dbs
#         self.uniprot_suf
#         self.path_to_uniprot_dbs
#         self.accurate_mz_map
#         self.prot_sets
#         self.protsN
        
        cnt, top_5_k, md_ar1, id_ar1 = self.get_matches(df, [-self.mass_accuracy, self.mass_accuracy], score_threshold=self.score_threshold)
        md_ar2 = []
        for z1, z2 in zip(md_ar1, id_ar1):
            if self.cnt_to_spec[z2-1] in top_5_k:
                md_ar2.append(z1)
        if len(md_ar2) == 0 :
            logger.debug('Something went wrong during blind search')
            return 1
         
        shift, sigma = optimize_md(md_ar2, bin_width=0.1) ######################## bin_width не менять?
        custom_range_mass_accuracy = [shift-2*sigma, shift+2*sigma]
        cnt, _, md_ar1, id_ar1 = self.get_matches(df, custom_range_mass_accuracy, score_threshold=self.score_threshold)
        
        top_100_species_names = set()
        for k, v in cnt.most_common():
            if len(top_100_species_names) < self.num_top_spec:
                k_orig = self.spec_map_id_reversed[k]
                if self.len_fasta_uniprot[int(k_orig)] < 220000:
                    if not int(k_orig) in self.exclude_names: 
                        top_100_species_names.add(k_orig) 
                        orig_name = NCBITaxa().get_taxid_translator([k_orig,])[int(k_orig)]
                        # print(k, k_orig, orig_name, v)
        prots = []
        cnt = 0    
        report = pd.DataFrame()
        for leader in top_100_species_names :
            SP = 0
            UN = 0
            sp_filename = str(leader)+self.sprot_suf+'.fasta'
            if sp_filename in listdir(self.path_to_sprot_dbs):
                for p in fasta.read(path.join(self.path_to_sprot_dbs, sp_filename)):
                    SP+=1
            un_filename = str(leader)+self.uniprot_suf+'.fasta'
            for p in fasta.read(path.join(self.path_to_uniprot_dbs, un_filename)):
                prots.append(p)
                UN+=1
            report = pd.concat([report, pd.DataFrame.from_dict({'ID':[leader],
                                                               'Sprot':[SP],
                                                               'Uniprot':[UN]})])

        if path_to_out_strain_statistics :
            report.to_csv(path_to_out_strain_statistics, index=False)
        random.shuffle(prots)
        with open(path_to_out_fasta, 'w') as f :
            fasta.write(prots, output=f)
        return 0
                           
    def get_matches(self, 
                    df1:pd.DataFrame,
                    custom_range_mass_accuracy:list,
                    score_threshold:float=4.0
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
                for mz_acc in self.accurate_mz_map[mz_int] :
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
                sp_id = self.cnt_to_spec[k-1]
                cnt[sp_id] += 1

        top_5_k = set()
        for k, v in cnt.most_common():
            if len(top_5_k) < 5:
                top_5_k.add(k)

        return cnt, top_5_k, md_ar1, id_ar1


def noisygaus(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b


def prepare_df(feature_tsv_path:str, charge_range:tuple=(2, 3), nIsotopes:int=3, mz_step:float=0.004) :
    df1 = pd.read_table(feature_tsv_path)
    df1 = df1[df1['nIsotopes'] >= nIsotopes]
    df1 = df1[df1['charge'] >= charge_range[0]]
    df1 = df1[df1['charge'] <= charge_range[1]]
    df1['massCalib_int'] = df1['massCalib'].apply(lambda x: int(x/mz_step))
    return df1


def calc_sf_all(v:int, n:int, p:float) :
    ###
    # calculates survival function for binomial discrete random variable
    # v <= n + 1
    # 0 <= p <= 1
    ###
    sf_values = -np.log10(binom.sf(v-1, n, p))
    sf_values[np.isinf(sf_values)] = max(sf_values[~np.isinf(sf_values)])
    return sf_values


def optimize_md(md_ar2:np.ndarray, bin_width=0.1) :
    bbins = np.arange(np.min(md_ar2), np.max(md_ar2), bin_width)
    H2, b2 = np.histogram(md_ar2, bins=bbins)
    m, mi, s = np.max(H2), b2[np.argmax(H2)], (np.max(md_ar2) - np.min(md_ar2))/6
    noise = np.min(H2)
    p0 = (m, mi, s, noise)
    logger.debug('p0: %s', p0 )
    popt, pcov = curve_fit(noisygaus, b2[1:], H2, p0=p0)
    shift, sigma = popt[1], abs(popt[2])
    logger.debug('Optimized mass shift and sigma: %s, %s', shift, sigma)
    return shift, sigma


def save(obj, file ) :
    if file :
        with open(file, 'wb') as f :
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        return 0
    else :
        return 1


def load(file) :
    with open(file, 'rb') as f:
        obj = pickle.load(f)
    return obj


def get_aa_mass_with_fixed_mods(
        fmods, 
        fmods_legend, 
        mods_custom_dict = {
            'Oxidation': 15.994915,
            'Carbamidomethyl': 57.021464,
            'TMT6plex': 229.162932,
            }
    ):

    if fmods_legend:
        for mod in fmods_legend.split(','):
            psiname, m = mod.split('@')
            mods_custom_dict[psiname] = float(m)

    aa_mass = deepcopy(mass.std_aa_mass)
    aa_to_psi = dict()

    mass_h2o = mass.calculate_mass('H2O')
    for k in list(aa_mass.keys()):
        aa_mass[k] = round(mass.calculate_mass(sequence=k) - mass_h2o, 7)

    if fmods:
        for mod in fmods.split(','):
            psiname, aa = mod.split('@')
            if psiname not in mods_custom_dict:
                logger.error('PSI Name for modification %s is missing in the modification legend' % (psiname, ))
                raise Exception('Exception: missing PSI Name for modification')
            if aa == '[':
                aa_mass['Nterm'] = float(mods_custom_dict[psiname])#float(m)
                aa_to_psi['Nterm'] = psiname
            elif aa == ']':
                aa_mass['Cterm'] = float(mods_custom_dict[psiname])#float(m)
                aa_to_psi['Cterm'] = psiname
            else:
                aa_mass[aa] += float(mods_custom_dict[psiname])#float(m)
                aa_to_psi[aa] = psiname

    logger.debug(aa_mass)

    return aa_mass, aa_to_psi


def call_Biosaur2(path_to_bio2:str, mzml_path:str, outpath:str, str_of_other_args:str='') :
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
    final_args = [path_to_bio2, mzml_path, '-o', outpath, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    logger.debug(final_args)
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore


def call_ThermoRawFileParser(path_to_mono:str, path_to_parser:str, raw_file:str, outdir:str, str_of_other_args:str='') :
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
    other_args += ['-L', '1']
    final_args = [path_to_mono, path_to_parser, '-i', raw_file, '-o', outdir, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    logger.debug(final_args)
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore


def call_ms1searchpy(path_to_ms1searchpy:str, feat_path:str, fasta_path:str, outdir:str='', cleavage_rule:str='', decoy_prefix:str='', str_of_other_args:str='') :
    wrong_add_args = ['-d', '-deeplc', 
                      '-e', '-ad', 
                      '-prefix', '-ml',
                      '-ts', '-o'
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
    other_args = ['-o', outdir, '-d', fasta_path, '-e', cleavage_rule, '-prefix', decoy_prefix, '-ad', '1', '-deeplc', '1', '-ml', '1', '-ts', '2'] + other_args
    final_args = [path_to_ms1searchpy, feat_path, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    logger.debug(final_args)
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore            
            
def call_ms1groups(path_to_ms1searchpy:str, PFM_ML:str, fasta_path:str, out:str='', group:str='genus', fdr:float=5, nproc:int=4, decoy_prefix:str='',) :
    
    other_args = ['-d', fasta_path, 
                  '-prefix', decoy_prefix, 
                  '-fdr', str(fdr), 
                  '-out', out, 
                  '-groups', group, 
                  '-nproc', str(nproc)]
    final_args = [path_to_ms1searchpy, PFM_ML, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    logger.debug(final_args)
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout :
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore

def call_DirectMS1quantmulti(path_to_ms1searchpy:str, pdir:str, fasta_path:str, samples:str, 
                             out:str='DQmulti', 
                             out_folder:str='',
                             pep_min_non_missing_samples:float=0.5, 
                             min_signif_for_pept:int=1, 
                             decoy_prefix:str='DECOY_',
                            ) :
    other_args = ['-d', fasta_path,
                  '-pdir', pdir,
                  '-samples', samples,
                  '-prefix', decoy_prefix, 
                  '-pep_min_non_missing_samples', str(pep_min_non_missing_samples), 
                  '-min_signif_for_pept', str(min_signif_for_pept),
                  '-out', out,
                  '-out_folder', out_folder,
                  ]
    final_args = [path_to_ms1searchpy, ] + other_args
    final_args = list(filter(lambda x: False if x=='' else True, final_args))
    logger.debug(final_args)
    process = subprocess.Popen(final_args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    with process.stdout :
        log_subprocess_output(process.stdout)
    exitscore = process.wait()
    return exitscore


def feature_generation(mzml_paths:dict, feature_folder:str, path_to_fd:str='', str_of_other_args:str='') :
    for sample, mzml_path in mzml_paths.items() :
        outpath = path.join(feature_folder, '.'.join([path.basename(mzml_path).rsplit('.', maxsplit=1)[0], 'features.tsv']) )
        # logger.debug(path_to_fd, mzml_path, outpath, str_of_other_args)
        call_Biosaur2(path_to_fd, mzml_path, outpath, str_of_other_args=str_of_other_args)
            

def mzml_generation(raw_paths:dict, mzml_folder:str, path_to_mono:str='', path_to_parser:str='', str_of_other_args:str='') :
    for sample, raw_path in raw_paths.items() :
        outdir = mzml_folder
        # logger.debug(path_to_mono, path_to_parser, raw_path, outdir, str_of_other_args)
        call_ThermoRawFileParser(path_to_mono, path_to_parser, raw_path, outdir, str_of_other_args=str_of_other_args)


def log_subprocess_output(pipe):
    for line in iter(pipe.readline, b''): # b'\n'-separated lines
        logger.info('From subprocess: %s', line.decode(sys.stdout.encoding).rstrip('\n'))
        

def unite_fasta(identification_table:str, out_fasta:str, uniprot_folder:str, threshold:float=0.02, uniprot_suf:str='_un', sep='\t') :
    group = identification_table.rsplit('_')[-1].split('.')[0]
    id_df = pd.read_csv(identification_table, sep=sep)
    exclude = ['name', 'group', 'taxid', 'include in combined fasta', 'len_fasta_by_ox', 'len_fasta_sum']
    samples = [col for col in id_df.columns if not col in exclude]
    
    needed_group_taxids = set()
    for sample in samples :
        t = set(id_df[ id_df['include in combined fasta'] ]['taxid'].values)
        needed_group_taxids.update(t)
    
    needed_taxids = set()
    for sample in samples :
        id_df_ox = pd.read_csv(identification_table.replace(group, 'OX'), sep=sep)
        t = list(id_df_ox['taxid'].values)
        for i in t :
            lineage = NCBITaxa().get_lineage(i)
            if any( [group_taxid in lineage for group_taxid in needed_group_taxids] ) :
                needed_taxids.update([i])

    with open(out_fasta, mode='w') as fout :
        for tid in needed_taxids :
            f = fasta.read(path.join(uniprot_folder, str(tid)+uniprot_suf+'.fasta'))
            fasta.write(f, output=fout)
            
            
def plot_tax_barplot(path_to_file:str, group:str='OX', search:str='blind', ascending:bool=True) : #save:bool=False) :
    # path_to_file = '/home/kae-13-1/bact_VGNKI_Apr2024/target_group_by_genus.tsv'
    plt.ioff()
    df = pd.read_csv(path_to_file, sep = '\t')
    excl = ['group', 'taxid', 'name', 'include in the combined fasta']
    file_cols = [col for col in df.columns if not col in excl]
    for file in file_cols :
        df = df.sort_values(by = file, ascending = ascending)
        y = df[file].values
        x = df['name'].values
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        colors = ['#4472c4ff','#133054ff']
        bars = ax.barh(x, y, color = colors)
        ax.bar_label(bars, padding = 5)
        ax.grid(axis='x')
        max_y = round(max(y)+100, -2)
        ax.set_xlim((0, max_y))
        # ax.set_xticks(fontsize = 14)
        ax.tick_params(axis='x', which='major', labelsize=14)
        ax.set_xlabel('# proteins', fontweight='bold')
        ax.set_title(file + ' ' + group)
        if 'include in the combined fasta' in df.columns :
            if ascending :
                h = len(df) - df['include in the combined fasta'].sum() - 0.5
            else :
                h = df['include in the combined fasta'].sum() + 0.5
            ax.text(max_y*0.3, h+0.05, 'border to include taxid in combined fasta')
            ax.plot((0, max_y), (h, h) , color='r')
        savepath = path.join(path.dirname(path_to_file), '_'.join([search, file, group])+'.png')
        fig.savefig(savepath, dpi = 300, transparent = False, bbox_inches = 'tight')
        
        
def plot_identification_hist(path_to_file:str, search:str='blind' ) : #save:bool=False) :
    # path_to_file = '/home/kae-13-1/bact_VGNKI_Apr2024/target_group_by_genus.tsv'
    plt.ioff()
    df = pd.read_csv(path_to_file, sep = '\t')
    excl = ['group', 'taxid', 'name', 'include in the combined fasta']
    file_cols = [col for col in df.columns if not col in excl]
    x = np.arange(len(file_cols))
    y = [df[file].sum() for file in file_cols]
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    colors = ['#4472c4ff','#133054ff']
    ax.bar(x, y, color = colors)
    # ax.bar_label(bars, padding = 5)
    # ax.grid(axis='x')
    max_y = round(max(y)+100, -2)
    # ax.set_xlim((0, max_y))
    ax.set_xticks(x, labels=file_cols, fontsize = 14, rotation=30)
    ax.set_ylabel('# proteins', fontweight='bold')
    ax.set_title('Identified proteins per file in {} search'.format(search))
    savepath = path.join(path.dirname(path_to_file), '_'.join([search, 'identified_proteins'])+'.png')
    fig.savefig(savepath, dpi = 300, transparent = False, bbox_inches = 'tight')
    
    
def read_cfg(file, category) :
    with open(file, 'r') as ff :
        f = ff.read()

    start = int(f.find('['+category+']')) + 2 + len(category)
    end = min(f.find('\n\n', start), len(f))

    cfg_string = f[start:end]
    while '#' in cfg_string :
        l = cfg_string.find('#')
        r = cfg_string.find('\n', l) + 1
        cfg_string = cfg_string[:l] + cfg_string[r:]

    lst_of_strings = cfg_string.lstrip().split('\n')
    final = []
    keys = []
    for el in lst_of_strings :
        if el :
            key_value = el.split(' = ')
            if len(key_value) > 1 :
                key = key_value[0]
                value = key_value[1]
            else :
                key = key_value[0]
                value = None
            keys.append(key)
            key = '-' + key
            final.append(key)

            if value :
                if value.startswith('\"') or value.startswith("\'") :
                    final.append(value)
                else :
                    vals = value.split()
                    for v in vals :
                        final.append(v)
    return final, keys


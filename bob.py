#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Usage:
        bob.py -p <FILE> -g <FILE> [-m <INT> -t <FILE> -o <STR>]


Options:
        -h --help                               show this
        -p, --phy <FILE>                        Phylip alignment file of all transcripts
        -m, --max_4d <INT>                      Maximum number of 4Xdegenerate sites. Zero means all. (restricts total metrics) [default: 0]
        -g, --groups <FILE>                     Orthogroups 
        -t, --target_ogs <FILE>                 List of orthogroups to analyse (restricts total metrics)
        -o, --outprefix <STRING>                Outprefix to use for output
        
Example commands:

    - Single copy orthologues: ./bob.py -p pieris_brassicae.phy -o singlecopy -g Orthogroups.sylvanus_fixed.txt -t SingleCopyOrthogroups.txt
    - Restrict by max4d count: ./bob.py -p pieris_brassicae.phy -m 614 -o 614 -g ../pi_by_locus/Orthogroups.sylvanus_fixed.txt

"""

from __future__ import division
from docopt import docopt
from collections import Counter, defaultdict
import sys
import os
from tqdm import tqdm
from timeit import default_timer as timer

gen_code_dict = { 'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L', \
        'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S', \
        'tat': 'Y', 'tac': 'Y', 'taa': 'X', 'tag': 'X', \
        'tgt': 'C', 'tgc': 'C', 'tga': 'X', 'tgg': 'W', \
        'ctt': 'L', 'ctc': 'L', 'cta': 'L', 'ctg': 'L', \
        'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', \
        'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q', \
        'cgt': 'R', 'cgc': 'R', 'cga': 'R', 'cgg': 'R', \
        'att': 'I', 'atc': 'I', 'ata': 'I', 'atg': 'M', \
        'act': 'T', 'acc': 'T', 'aca': 'T', 'acg': 'T', \
        'aat': 'N', 'aac': 'N', 'aaa': 'K', 'aag': 'K', \
        'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R', \
        'gtt': 'V', 'gtc': 'V', 'gta': 'V', 'gtg': 'V', \
        'gct': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A', \
        'gat': 'D', 'gac': 'D', 'gaa': 'E', 'gag': 'E', \
        'ggt': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G'} 

# ['ctg', 'ctg', 'ctg', 'ttg']
# [2, 2, 2, None]

#dict to look up degeneracy 
degen_dict = {'ttt': '002', 'ttc': '002', 'tta': '202', 'ttg': '202', \
        'tct': '004', 'tcc': '004', 'tca': '004', 'tcg': '004', \
        'tat': '002', 'tac': '002', 'taa': '022', 'tag': '002', \
        'tgt': '002', 'tgc': '002', 'tga': '020', 'tgg': '000', \
        'ctt': '004', 'ctc': '004', 'cta': '204', 'ctg': '204', \
        'cct': '004', 'ccc': '004', 'cca': '004', 'ccg': '004', \
        'cat': '002', 'cac': '002', 'caa': '002', 'cag': '002', \
        'cgt': '004', 'cgc': '004', 'cga': '204', 'cgg': '204', \
        'att': '003', 'atc': '003', 'ata': '003', 'atg': '000', \
        'act': '004', 'acc': '004', 'aca': '004', 'acg': '004', \
        'aat': '002', 'aac': '002', 'aaa': '002', 'aag': '002', \
        'agt': '002', 'agc': '002', 'aga': '202', 'agg': '202', \
        'gtt': '004', 'gtc': '004', 'gta': '004', 'gtg': '004', \
        'gct': '004', 'gcc': '004', 'gca': '004', 'gcg': '004', \
        'gat': '002', 'gac': '002', 'gaa': '002', 'gag': '002', \
        'ggt': '004', 'ggc': '004', 'gga': '004', 'ggg': '004'}


def compute_pi(hetA, hetB, hetAB, fixed, invariant):
    total_sites = sum([hetA, hetB, hetAB, fixed, invariant])
    pi_A, pi_B, pi_between, pi_total = 0.0, 0.0, 0.0, 0.0 
    if total_sites:
        pi_A = float("%.8f" % ((hetA + hetAB) / total_sites))
        pi_B = float("%.8f" % ((hetB + hetAB) / total_sites))
        pi_between = float("%.8f" % ((((hetA + hetB + hetAB) / 2) + fixed) / total_sites))
        pi_total = float("%.8f" % (((4.0 / 6) * pi_between) + ((1.0 / 6) * pi_A) + ((1.0 / 6) * pi_B)))
    else:
        pi_A = 'NA'
        pi_B = 'NA'
        pi_between = 'NA'
        pi_total = 'NA'
    return pi_A, pi_B, pi_between, pi_total, total_sites

def get_degeneracies(codons, variant_pos):
    return [int(degen_dict.get(codon)[variant_pos]) for codon in codons]

def get_genotypes(codons, variant_pos):
    return [codon[variant_pos] for codon in codons]

def get_AA(codons):
    return [gen_code_dict.get(codon, None) for codon in codons]

def get_variant_pos(codons):
    variant_pos = []
    for codon in codons:
        for idx, base in enumerate(codon):
            if not base == codons[0][idx]:
                variant_pos.append(idx)
    return list(set(variant_pos))

fourXdeg_genotype = {
        'tct': 't', 'tcc': 'c', 'tca': 'a', 'tcg': 'g', \
        'ctt': 't', 'ctc': 'c', 'cta': 'a', 'ctg': 'g', \
        'cct': 't', 'ccc': 'c', 'cca': 'a', 'ccg': 'g', \
        'cgt': 't', 'cgc': 'c', 'cga': 'a', 'cgg': 'g', \
        'act': 't', 'acc': 'c', 'aca': 'a', 'acg': 'g', \
        'gtt': 't', 'gtc': 'c', 'gta': 'a', 'gtg': 'g', \
        'gct': 't', 'gcc': 'c', 'gca': 'a', 'gcg': 'g', \
        'ggt': 't', 'ggc': 'c', 'gga': 'a', 'ggg': 'g'}

def get_fourXdeg_genotypes(codons):
    return [fourXdeg_genotype.get(codon, None) for codon in codons]


def get_zygosity(geno1, geno2):
    if geno1 == geno2:
        return 'HOM'
    else:
        return 'HET'

class DataObj(object):
    def __init__(self):
        self.alnObjs = []
        self.sites = Counter()
        self.codons = Counter()
        self.species = None
        self.samples = []
        self.metrics = {
            '0_pi_A' : 0,
            '0_pi_B' : 0,
            '0_pi_between' : 0,
            '0_pi_total' : 0,
            '4_pi_A' : 0,
            '4_pi_B' : 0,
            '4_pi_between' : 0,
            '4_pi_total' : 0
            }

    def __str__(self):
        return "\t".join([str(x) for x in [self.species, 'all_loci', self.codons['proper'], self.codons['multivariant'], self.codons['multidegenerate'], self.codons['valid'], self.sites['0_all'], self.sites['0_invariant'], self.sites['0_hetA'], self.sites['0_hetB'], self.sites['0_hetAB'], self.sites['0_fixed'], self.metrics['0_pi_A'], self.metrics['0_pi_B'], self.metrics['0_pi_between'], self.metrics['0_pi_total'], self.sites['4_all'], self.sites['4_invariant'], self.sites['4_hetA'], self.sites['4_hetB'], self.sites['4_hetAB'], self.sites['4_fixed'], self.metrics['4_pi_A'], self.metrics['4_pi_B'], self.metrics['4_pi_between'], self.metrics['4_pi_total']]])

    def parse_groups(self, groups_f, target_groups_f):
        target_ogs = defaultdict(lambda: False)
        if target_groups_f:
            with open(target_groups_f) as target_groups_fh:
                for line in target_groups_fh:
                    target_ogs[line.rstrip("\n")] = True
        else:
            target_ogs = defaultdict(lambda: True)
        group_id_by_protein_id = {}
        with open(groups_f) as groups_fh:
            for line in groups_fh:
                temp = line.rstrip("\n").split()
                group_id = temp[0].rstrip(":")
                if target_ogs[group_id]:
                    for i in range(1, len(temp)):
                        protein_id = "_".join(temp[i].split("_")[0:4])
                        group_id_by_protein_id[protein_id] = group_id
        return group_id_by_protein_id

    def parse_phy(self, phy_f):
        alnObj = AlnObj()
        with open(phy_f) as phy_fh:
            for line in phy_fh:
                temp = line.rstrip("\n").split()
                if temp: 
                    if temp[0] == '4':  # header
                        if alnObj.locus:
                            self.add_alnObj(alnObj)
                            alnObj = AlnObj()
                        alnObj.aln_length = int(temp[1])
                    else: # seq
                        alnObj.add_seq(temp)
        self.add_alnObj(alnObj)

    def add_alnObj(self, alnObj):
        if not self.species:
            self.species = alnObj.species
            self.samples = alnObj.samples
        else:
            if not self.species == alnObj.species:
                sys.exit("[X] Inconsistent species: %s != %s" % (self.species, species))  
            if not self.samples == alnObj.samples:
                sys.exit("[X] Inconsistent samples: %s != %s" % (self.samples, samples))
        self.alnObjs.append(alnObj)

    def analyse(self, max_4d, target_groups_f, groups_f):
        group_id_by_protein_id = self.parse_groups(groups_f, target_groups_f)
        for alnObj in tqdm(self.alnObjs, total=len(self.alnObjs), desc="[%] ", ncols=80):
            alnObj.group = group_id_by_protein_id.get(alnObj.protein, 'NA')
            alnObj.count_sites(max_4d)
            alnObj.compute_transcript_pi()
            if not alnObj.group == 'NA':
                self.sites += Counter(alnObj.sites)
                self.codons += Counter(alnObj.codons)
        self.metrics['0_pi_A'], self.metrics['0_pi_B'], self.metrics['0_pi_between'], self.metrics['0_pi_total'], self.sites['0_all'] = compute_pi(self.sites['0_hetA'], self.sites['0_hetB'], self.sites['0_hetAB'], self.sites['0_fixed'], self.sites['0_invariant'])
        self.metrics['4_pi_A'], self.metrics['4_pi_B'], self.metrics['4_pi_between'], self.metrics['4_pi_total'], self.sites['4_all'] = compute_pi(self.sites['4_hetA'], self.sites['4_hetB'], self.sites['4_hetAB'], self.sites['4_fixed'], self.sites['4_invariant'])

    def write_results(self, outprefix):
        loci_header = "\t".join(['species', 'locus', 'protein', 'group', 'proper_codons', 'multivariant_codons', 'multidegenerate_codons', 'valid_codons', '0_all', '0_invariant', '0_hetA', '0_hetB', '0_hetAB', '0_fixed', '0_pi_%s' % (self.samples[0]), '0_pi_%s' % (self.samples[1]), '0_pi_between', '0_pi_total', '4_all', '4_invariant', '4_hetA', '4_hetB', '4_hetAB', '4_fixed', '4_pi_%s' % (self.samples[0]), '4_pi_%s' % (self.samples[1]), '4_pi_between', '4_pi_total'])
        loci_output = []
        loci_output.append(loci_header)
        for alnObj in self.alnObjs:
            loci_output.append(str(alnObj))
        loci_fn = "loci.metrics.txt"
        if outprefix:
            loci_fn = outprefix + "." + loci_fn
        with open(loci_fn, 'w') as loci_fh:
            loci_fh.write("\n".join(loci_output) + "\n")
        total_header = "\t".join(['species', 'locus', 'proper_codons', 'multivariant_codons', 'multidegenerate_codons', 'valid_codons', '0_all', '0_invariant', '0_hetA', '0_hetB', '0_hetAB', '0_fixed', '0_pi_%s' % (self.samples[0]), '0_pi_%s' % (self.samples[1]), '0_pi_between', '0_pi_total', '4_all', '4_invariant', '4_hetA', '4_hetB', '4_hetAB', '4_fixed', '4_pi_%s' % (self.samples[0]), '4_pi_%s' % (self.samples[1]), '4_pi_between', '4_pi_total'])
        total_output = []
        total_output.append(total_header)
        total_output.append(str(self))
        total_fn = "total.metrics.txt"
        if outprefix:
            total_fn = outprefix + "." + total_fn
        with open(total_fn, 'w') as total_fh:
            total_fh.write("\n".join(total_output) + "\n")

class AlnObj(object):
    def __init__(self):
        self.locus = None
        self.protein = None
        self.species = None
        self.group = None
        self.aln_length = 0
        self.samples = [] # order matters (A, B)
        self.sequence_by_haplotype_by_sample = {}
        self.sites = {
            '0_all' : 0,
            '0_invariant' : 0,
            '0_hetA' : 0,
            '0_hetB' : 0,
            '0_hetAB' : 0,
            '0_fixed' : 0,
            '4_all' : 0,
            '4_invariant' : 0,
            '4_hetA' : 0,
            '4_hetB' : 0,
            '4_hetAB' : 0,
            '4_fixed' : 0,
        }
        self.codons = {
            'proper' : 0,
            'multivariant': 0,
            'multidegenerate': 0,
            'valid' : 0
        }
        self.metrics = {
            '0_pi_A' : 0,
            '0_pi_B' : 0,
            '0_pi_between' : 0,
            '0_pi_total' : 0,
            '4_pi_A' : 0,
            '4_pi_B' : 0,
            '4_pi_between' : 0,
            '4_pi_total' : 0
        }

    def __str__(self):
        return "\t".join([str(x) for x in [self.species, self.locus, self.protein, self.group, self.codons['proper'], self.codons['multivariant'], self.codons['multidegenerate'], self.codons['valid'], self.sites['0_all'], self.sites['0_invariant'], self.sites['0_hetA'], self.sites['0_hetB'], self.sites['0_hetAB'], self.sites['0_fixed'], self.metrics['0_pi_A'], self.metrics['0_pi_B'], self.metrics['0_pi_between'], self.metrics['0_pi_total'], self.sites['4_all'], self.sites['4_invariant'], self.sites['4_hetA'], self.sites['4_hetB'], self.sites['4_hetAB'], self.sites['4_fixed'], self.metrics['4_pi_A'], self.metrics['4_pi_B'], self.metrics['4_pi_between'], self.metrics['4_pi_total']]])

    def add_seq(self, temp):
        header = temp[0].split(".")
        seq = temp[1]
        species = header[0]
        sample = header[1]
        haplotype = header[2]
        locus = header[3]
        self.check_consistency(species, locus)
        if not sample in self.sequence_by_haplotype_by_sample:
            self.samples.append(sample)
            self.sequence_by_haplotype_by_sample[sample] = {}
        self.sequence_by_haplotype_by_sample[sample][haplotype] = seq.lower()

    def check_consistency(self, species, locus):
        if not self.species and not self.locus:
            self.species = species
            self.locus = locus
            self.protein = ".".join([species.split("_")[1], locus.split("(")[0]])
        else:
            if not self.species == species:
                sys.exit("[X] Inconsistent species: %s != %s" % (self.species, species))  
            if not self.locus == locus:
                sys.exit("[X] Inconsistent locus: %s != %s" % (self.locus, locus))

    def yield_codons(self):
        for start in range(0, self.aln_length, 3):
            codons = []
            for sample in self.samples:
                for haplotype, sequence in sorted(self.sequence_by_haplotype_by_sample[sample].items()):
                    codons.append(sequence[start:start + 3])
            yield codons

    def count_sites(self, max_4d):
        '''
        Only mutations at a single 4D/0D site per codon are considered    
        Only mutations between 0D sites are considered
        '''
        #print self.locus
        for codons in self.yield_codons():
            
            if max_4d:
                if self.sites['fourX'] >= max_4d:
                    break
            if not "-" in set("".join(codons)):
                self.codons['proper'] += 1
                variant_positions = get_variant_pos(codons)
                aminoacids = get_AA(codons)
                if not variant_positions: # invariant codon
                    self.codons['valid'] += 1
                    #print "[+] Invariant: %s [%s]" % (codons, aminoacids)
                    for idx in range(0, 3):
                        degeneracy = get_degeneracies(codons, idx)[0]  # taking first, since all the same
                        key = "%s_invariant" % degeneracy
                        if key in self.sites:
                            self.sites[key] += 1
                else:  # variant codon 
                    if not len(variant_positions) == 1: # if more than one variant, ignore ... but record
                        self.codons['multivariant'] += 1
                        #print "[+] MultiVariant: %s [%s] var_pos=%s" % (codons, aminoacids, variant_positions)
                    else: # only single variant                    
                        variant_position = variant_positions[0]
                        degeneracies = get_degeneracies(codons, variant_position)
                        if not len(set(degeneracies)) == 1: # multi degeneracy variant 
                            self.codons['multidegenerate'] += 1
                            #print "[+] MultiDegenerate: %s [%s] var_pos=%s deg=%s" % (codons, aminoacids, variant_positions, degeneracies)
                        else:
                            self.codons['valid'] += 1
                            genotypes = get_genotypes(codons, variant_position)
                            if len(set(genotypes)) == 2:
                                degeneracy = degeneracies[0]
                                A1, A2, B1, B2 = genotypes
                                if degeneracy in set([0, 4]):
                                    #print "[+] Variant: %s %s var_pos=%s deg=%s" % (codons, genotypes, variant_position, degeneracy)
                                    A_zygosity = get_zygosity(A1, A2)
                                    B_zygosity = get_zygosity(B1, B2)
                                    key = None
                                    if A_zygosity == 'HOM' and B_zygosity == 'HOM':
                                        #print "A1=%s A2=%s B1=%s B2=%s Deg=%s [fixed]" % (A1, A2, B1, B2, degeneracy)
                                        key = "%s_fixed" % degeneracy
                                    elif A_zygosity == 'HET' and B_zygosity == 'HET':
                                        #print "A1=%s A2=%s B1=%s B2=%s Deg=%s [hetAB]" % (A1, A2, B1, B2, degeneracy)
                                        key = "%s_hetAB" % degeneracy
                                    elif A_zygosity == 'HET' and B_zygosity == 'HOM':
                                        #print "A1=%s A2=%s B1=%s B2=%s Deg=%s [hetA]" % (A1, A2, B1, B2, degeneracy)
                                        key = "%s_hetA" % degeneracy
                                    elif A_zygosity == 'HOM' and B_zygosity == 'HET':
                                        #print "A1=%s A2=%s B1=%s B2=%s Deg=%s [hetB]" % (A1, A2, B1, B2, degeneracy)
                                        key = "%s_hetB" % degeneracy
                                    else:
                                        sys.exit("[X] A1=%s A2=%s B1=%s B2=%s [???]" % (A1, A2, B1, B2, degeneracy))
                                    self.sites[key] += 1

    def compute_transcript_pi(self):
        self.metrics['0_pi_A'], self.metrics['0_pi_B'], self.metrics['0_pi_between'], self.metrics['0_pi_total'], self.sites['0_all'] = compute_pi(self.sites['0_hetA'], self.sites['0_hetB'], self.sites['0_hetAB'], self.sites['0_fixed'], self.sites['0_invariant'])
        self.metrics['4_pi_A'], self.metrics['4_pi_B'], self.metrics['4_pi_between'], self.metrics['4_pi_total'], self.sites['4_all'] = compute_pi(self.sites['4_hetA'], self.sites['4_hetB'], self.sites['4_hetAB'], self.sites['4_fixed'], self.sites['4_invariant'])

if __name__ == "__main__":
    __version__ = 0.2 

    try:
        start_time = timer()
        args = docopt(__doc__)
        # print "[V] bob.py v%s" % __version__
        dataObj = DataObj()
        target_groups_f = args['--target_ogs']
        groups_f = args['--groups']
        print "[+] Reading alignments ..."
        dataObj.parse_phy(args['--phy'])
        print "[+] Analysing alignments ..."
        dataObj.analyse(int(args['--max_4d']), target_groups_f, groups_f)
        print "[+] Writing results ..."
        dataObj.write_results(args['--outprefix'])
        print "[+] Total runtime: %.3fs" % (timer() - start_time)
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)

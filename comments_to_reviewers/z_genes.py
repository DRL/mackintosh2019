#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
usage: z_genes.py [-h] -d DIR [-m FILE] [-f FILE] [-u FILE] [--4pi] [-p STR]

A script which outputs useful data concerning whether gene/transcripts/loci
are Z-linked given heterozygosity in males and females

optional arguments:
  -h, --help                show this help message and exit
  -d, --dir DIR             Directory containing *.phy.SCO.loci.metrics.txt
  -m, --males FILE          list of male IDs
  -f, --females FILE        list of female IDs
  -u, --unknowns FILE       list of unknown IDs
  -p, --prefix STR          Output prefix [default: z_genes]
  --4pi                     use 4_pi values (default is using 0_pi)
'''

'''
Installation:
1. # create conda enviroment with dependencies
conda env create -f=z_genes.conda.yaml -n zgenes
'''


from collections import defaultdict
from timeit import default_timer as timer
from sys import stderr, exit
from docopt import docopt
from os.path import splitext, join
from os import listdir

def get_fs(directory, suffix):
    fs = []
    for f in sorted(listdir(directory)):
        if f.endswith(suffix):
            fs.append(join(directory, f))
    return fs

def write_lines(lines, f):
    with open(f, 'w') as fh:
        fh.write("\n".join(lines))

def get_sex(args):
    fs = [args['--males'], args['--females'], args['--unknowns']]
    sex_by_sample_id = defaultdict(lambda: 'NA') # default is not X
    sample_ids_by_sex = defaultdict(list)
    for sex, f in zip(SEXES, fs):
        with open(f, 'r') as fh:
            for line in fh:
                sample_id = line.rstrip()
                sex_by_sample_id[sample_id] = sex
                sample_ids_by_sex[sex].append(sample_id)
    return (sex_by_sample_id, sample_ids_by_sex)

def parse_sco_loci_metrics(args, sex_by_sample_id):
    print("[+] Processing *.phy.SCO.loci.metrics.txt ...")
    fs = get_fs(directory=args['--dir'], suffix='.phy.SCO.loci.metrics.txt')
    orthogroups = OrthogroupsObj()
    speciesObjs = []
    for f in fs:
        species_id = splitext(f)[0]
        speciesObj = SpeciesObj(species_id)
        with open(f, 'r') as fh:
            header_col = fh.readline().rstrip().split()
            # Sample A
            sample_id_A = header_col[14].lstrip("0_pi")
            sex_A = sex_by_sample_id[sample_id_A]
            sampleObj_A = SampleObj(sample_id_A, species_id, sex_A)
            speciesObj.add_sampleObj(sampleObj_A)
            # Sample B
            sample_id_B = header_col[15].lstrip("0_pi")
            sex_B = sex_by_sample_id[sample_id_B]
            sampleObj_B = SampleObj(sample_id_B, species_id, sex_B)
            speciesObj.add_sampleObj(sampleObj_B)
            # pi
            for line in fh:
                col = line.rstrip().split()
                orthogroup_id = col[3]
                if not orthogroup_id == 'NA':
                    # orthogroup
                    orthogroups.add_sample_ids_for_orthogroup(sample_id_A, sample_id_B, orthogroup_id)
                    # pi
                    pi_by_sample_id = {}
                    if args['--4pi']:
                        pi_by_sample_id[sample_id_A] = col[24]
                        pi_by_sample_id[sample_id_B] = col[25]
                    else:
                        pi_by_sample_id[sample_id_A] = col[14]
                        pi_by_sample_id[sample_id_B] = col[15]
                    speciesObj.add_orthogroup(orthogroup_id, pi_by_sample_id)
        speciesObjs.append(speciesObj)
    return (orthogroups, speciesObjs)

def write_genotypes_by_species(args, speciesObjs):
    lines = []
    lines.append("\t".join(['species', 'sex', 'orthogroups', 'missing', 'HOM_M', 'HOM_F', 'HOM_X', 'HOM_M_%', 'HOM_F_%', 'HOM_X_%']))
    for speciesObj in speciesObjs:
        values_by_label = speciesObj.get_values_by_label()
        line = []
        line.append(speciesObj.species_id)
        line.append("".join(sorted(speciesObj.sex_by_sample_id.values())))
        line.append(len(speciesObj.orthogroup_ids))
        line.append(len(speciesObj.orthogroup_ids) - values_by_label.get('total', 0))
        line.append(values_by_label.get('HOM_M', 0))
        line.append(values_by_label.get('HOM_F', 0))
        line.append(values_by_label.get('HOM_X', 0))
        line.append(values_by_label.get('HOM_M_%', 0.0))
        line.append(values_by_label.get('HOM_F_%', 0.0))
        line.append(values_by_label.get('HOM_X_%', 0.0))
        lines.append("\t".join([str(l) for l in line]))
    outfile = "%s.genotype_count.0pi.tsv" % args['--prefix'] 
    if args['--4pi']:
        outfile = "%s.genotype_count.4pi.tsv" % args['--prefix']
    write_lines(lines, outfile)
    return outfile

def write_orthogroup_homs(args, orthogroups, speciesObjs):
    for orthogroup in orthogroups
class SampleObj(object):
    def __init__(self, sample_id, species_id, sex):
        self.sample_id = sample_id
        self.species_id = species_id
        self.genotype_by_orthogroup_id = defaultdict(lambda: 'NA')
        self.sex = sex

    def add_genotype_by_orthogroup_id_from_pi(self, pi, orthogroup):
        if pi == '0.0':
            self.genotype_by_orthogroup_id[orthogroup] = 'HOM'
        else:
            self.genotype_by_orthogroup_id[orthogroup] = 'HET'

class SpeciesObj(object):
    def __init__(self, species_id):
        self.species_id = species_id
        self.sampleObj_by_sample_id = {}
        self.orthogroup_ids = []
        self.sex_by_sample_id = {}

    def add_sampleObj(self, sampleObj):
        self.sampleObj_by_sample_id[sampleObj.sample_id] = sampleObj
        self.sex_by_sample_id[sampleObj.sample_id] = sampleObj.sex

    def add_orthogroup(self, orthogroup_id, pi_by_sample_id):
        self.orthogroup_ids.append(orthogroup_id)
        for sample_id, pi in pi_by_sample_id.items():
            self.sampleObj_by_sample_id[sample_id].add_genotype_by_orthogroup_id_from_pi(pi, orthogroup_id)

    def get_values_by_label(self):
        values_by_label = {}
        for sample_id, sex in self.sex_by_sample_id.items():
            label = "HOM_%s" % sex
            label_perc = "HOM_%s_%%" % sex
            sampleObj = self.sampleObj_by_sample_id[sample_id]
            values_by_label[label] = len([genotype for genotype in sampleObj.genotype_by_orthogroup_id.values() if genotype == 'HOM'])
            values_by_label[label_perc] = "{:.2}".format(values_by_label[label] / len(self.orthogroup_ids))
            values_by_label['total'] = len(sampleObj.genotype_by_orthogroup_id.values())
        return values_by_label

class OrthogroupsObj(object):
    def __init__(self):
        self.orthogroup_ids = set()
        self.orthogroup_ids_by_sample_ids = defaultdict(list)

    def add_sample_ids_for_orthogroup(self, sample_id_A, sample_id_B, orthogroup_id):
        self.orthogroup_ids.add(orthogroup_id)
        self.orthogroup_ids_by_sample_ids[sample_id_A].append(orthogroup_id)
        self.orthogroup_ids_by_sample_ids[sample_id_B].append(orthogroup_id)

    def __len__(self):
        return len(self.orthogroup_ids)

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        sex_by_sample_id, sample_ids_by_sex = get_sex(args)
        print("[+] Parsed lists ...\n%s" % ("\n".join(["\t [%s] %s sample IDs" % (sex, len(sample_ids)) for sex, sample_ids in sample_ids_by_sex.items()])))
        orthogroups, speciesObjs = parse_sco_loci_metrics(args, sex_by_sample_id)
        print("[+] Parsed data for %s orthogroups ..." % (len(orthogroups))) # should be 1314
        print("[+] Parsed data for %s species ..." % (len(speciesObjs))) # should be 38
        genotype_count_f = write_genotypes_by_species(args, speciesObjs)
        print("[+] Written %s ..." % genotype_count_f) 
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    SEXES = ['M', 'F', 'X']
    main()
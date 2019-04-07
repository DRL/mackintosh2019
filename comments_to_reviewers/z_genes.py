#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
usage: z_genes.py [-h] -d DIR -y DIR -m FILE [-p STR]

A script which outputs useful data concerning whether gene/transcripts/loci
are Z-linked given heterozygosity in males and females

optional arguments:
  -h, --help                show this help message and exit
  -d, --sco DIR             Directory containing *.phy.SCO.loci.metrics.txt
  -y, --phy DIR             Directory containing *.phy
  -m, --metadata FILE       File with metadata 
  -p, --prefix STR          Output prefix [default: z_genes]
'''

'''
TASKS:
    - Transdecoder output
        - annotation for all orthogroups
    - Do HOM_M/HOM_F for all transcripts
    - read in alignments
        - length, numdiff
        - num hetM, hetF
        => Identity : proxy for pairwise heterozygosity

- Males     = ZZ
- Females   = ZW (W has no genes)

# 'euchlodes sylvanus' : MM
# 'cionympha arcana' : MM
'''


from collections import defaultdict
from timeit import default_timer as timer
from sys import stderr, exit
from docopt import docopt
from os.path import basename, join
from os import listdir
from tqdm import tqdm

def get_fs(directory, suffix):
    fs = []
    for f in sorted(listdir(directory)):
        if f.endswith(suffix):
            fs.append(join(directory, f))
    return fs

def write_lines(lines, f):
    with open(f, 'w') as fh:
        fh.write("\n".join(lines))


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


#def parse_phy_by_species(phy_dir):
#    alnObj = AlnObj()
#        with open(phy_f) as phy_fh:
#            for line in phy_fh:
#                temp = line.rstrip("\n").split()
#                if temp: 
#                    if temp[0] == '4':  # header
#                        if alnObj.locus:
#                            self.add_alnObj(alnObj)
#                            alnObj = AlnObj()
#                        alnObj.aln_length = int(temp[1])
#                    else: # seq
#                        alnObj.add_seq(temp)
#        self.add_alnObj(alnObj)
#
#    def add_alnObj(self, alnObj):
#        if not self.species:
#            self.species = alnObj.species
#            self.samples = alnObj.samples
#        else:
#            if not self.species == alnObj.species:
#                sys.exit("[X] Inconsistent species: %s != %s" % (self.species, species))  
#            if not self.samples == alnObj.samples:
#                sys.exit("[X] Inconsistent samples: %s != %s" % (self.samples, samples))
#        self.alnObjs.append(alnObj)

#########################################################################

def parse_metadata(metadata_f):
    metadataObj = MetadataObj()
    with open(metadata_f, 'r') as fh:
        header = fh.readline().rstrip().split(",")
        for line in fh:
            col = line.rstrip().split(",")
            metadataObj.add_data(*col)
    return metadataObj

class SeqdataObj(object):
    def __init__(self):
        self.transcriptObj_by_transcript_id_by_species_id = defaultdict(dict)

    def add_transcriptObj(self, transcriptObj):
        self.transcriptObj_by_transcript_id_by_species_id[transcriptObj.species_id][transcriptObj.locus_id] = transcriptObj

    def get_transcriptObjs(self, species_id):
        return self.transcriptObj_by_transcript_id_by_species_id[species_id].values()

def parse_sco_loci_metrics(directory, metadataObj):
    print("[+] Processing *.phy.SCO.loci.metrics.txt ...")
    fs = get_fs(directory=directory, suffix='.phy.SCO.loci.metrics.txt')
    seqdataObj = SeqdataObj()
    for f in tqdm(fs, total=len(fs), desc="[%] ", ncols=200):
        species_id = basename(f).split(".")[0]
        with open(f, 'r') as fh:
            header_col = fh.readline().rstrip().split()
            sample_id_A = header_col[14].lstrip("0_pi")
            sample_id_B = header_col[15].lstrip("0_pi")
            for line in fh:
                col = line.rstrip().split()
                locus_id = col[1]
                orthogroup_id = col[3] if not col[3] == 'NA' else None
                transcriptObj = TranscriptObj(locus_id, species_id, orthogroup_id)
                transcriptObj.add_genotype('zero_d', sample_id_A, col[14])
                transcriptObj.add_genotype('zero_d', sample_id_B, col[15])
                transcriptObj.add_genotype('four_d', sample_id_A, col[24])
                transcriptObj.add_genotype('four_d', sample_id_B, col[25])
                seqdataObj.add_transcriptObj(transcriptObj)
    return seqdataObj

def parse_phy(directory, seqdataObj, metadataObj):
    print("[+] Processing *.phy ...")
    fs = get_fs(directory=directory, suffix='.phy')
    for f in tqdm(fs, total=len(fs), desc="[%] ", ncols=200):
        species_id = basename(f).split(".")[0]
        with open(f, 'r') as fh:
            for line in fh:
                if line.startswith(species_id):
                    line1 = line.split()
                    line2 = fh.readline().split()
                    sample_id = line1[0].split(".")[1]
                    locus_id = line1[0].split(".")[3]
                    genotype = 'HOM' if line1[1] == line2[1] else 'HET'
                    transcriptObj = seqdataObj.transcriptObj_by_transcript_id_by_species_id[species_id][locus_id]
                    transcriptObj.infer_length_from_phy(line1[1])
                    transcriptObj.add_genotype('phy', sample_id, genotype)
    return seqdataObj

class TranscriptObj(object):
    def __init__(self, locus_id, species_id, orthogroup_id):
        self.locus_id = locus_id
        self.species_id = species_id
        self.orthogroup_id = orthogroup_id
        self.four_d_genotype_by_sample_id = {}
        self.zero_d_genotype_by_sample_id = {}
        self.all_genotype_by_sample_id = {}
        self.length = 0

    def infer_length_from_phy(self, seq):
        self.length = len(seq.replace("-",""))

    def add_genotype(self, category, sample_id, genotype):
        if category == 'phy':
            self.all_genotype_by_sample_id[sample_id] = genotype
        else:
            if category == 'four_d':
                self.four_d_genotype_by_sample_id[sample_id] = 'HOM' if genotype == '0.0' else 'HET'
            elif category == 'zero_d':
                self.zero_d_genotype_by_sample_id[sample_id] = 'HOM' if genotype == '0.0' else 'HET'

class MetadataObj(object):
    def __init__(self):
        self.sex_by_sample_id = {}
        self.species_id_by_sample_id = {}
        self.sexed_by_species_id = defaultdict(bool)
        self.sample_ids_by_sex = defaultdict(list)
        self.sample_ids_by_species_id = defaultdict(list)

    def add_data(self, species_id, sample_id, sex):
        self.sex_by_sample_id[sample_id] = sex
        self.species_id_by_sample_id[sample_id] = species_id
        self.sample_ids_by_sex[sex].append(sample_id)
        self.sample_ids_by_species_id[species_id].append(sample_id)
        if len(self.sample_ids_by_species_id[species_id]) == 2:
            sexes = set([self.sex_by_sample_id[sample_id] for sample_id in self.sample_ids_by_species_id[species_id]])
            if sexes == set(['F', 'M']):
                self.sexed_by_species_id[species_id] = True

def analyse_transcriptObjs(seqdataObj, metadataObj):
    species_ids_by_composition_by_orthogroup_id = defaultdict(lambda: defaultdict(list))
    count_by_composition_by_species_ids = defaultdict(lambda: defaultdict(int))
    for species_id, sexed in tqdm(metadataObj.sexed_by_species_id.items(), total=len(metadataObj.sexed_by_species_id), desc="[%] ", ncols=200):
        if sexed:
            transcriptObjs = seqdataObj.get_transcriptObjs(species_id)
            for transcriptObj in transcriptObjs:
                pattern = [genotype for sample_id, genotype in transcriptObj.all_genotype_by_sample_id.items()]
                if pattern[0] == pattern[1]:
                    composition = '%s_%s' % (pattern[0], 'both')
                    if transcriptObj.orthogroup_id:
                        species_ids_by_composition_by_orthogroup_id[transcriptObj.orthogroup_id][composition].append(species_id)
                    count_by_composition_by_species_ids[species_id][composition] += 1
                else:
                    for sample_id, genotype in transcriptObj.all_genotype_by_sample_id.items():
                        sex = metadataObj.sex_by_sample_id[sample_id]
                        composition = '%s_%s' % (genotype, sex)
                        if transcriptObj.orthogroup_id:
                            species_ids_by_composition_by_orthogroup_id[transcriptObj.orthogroup_id][composition].append(species_id)
                        count_by_composition_by_species_ids[species_id][composition] += 1

    header_int = "\t".join(['orthogroup_id', 'species', 'HOM_both', 'HET_both', 'HOM_M', 'HOM_F'])
    lines_int = [header_int]
    header_str = "\t".join(['orthogroup_id', 'HOM_both', 'HET_both', 'HOM_M', 'HOM_F'])
    lines_str = [header_str]
    for orthogroup_id in species_ids_by_composition_by_orthogroup_id:
        list_of_species = []
        for composition, species_ids in species_ids_by_composition_by_orthogroup_id[orthogroup_id].items():
            list_of_species.extend(species_ids)
        total = len(set(list_of_species))
        line_int = [orthogroup_id, total]
        line_str = [orthogroup_id]
        for composition in ['HOM_both', 'HET_both', 'HOM_M', 'HOM_F']:
            line_int.append("%.2f" % (len(species_ids_by_composition_by_orthogroup_id[orthogroup_id].get(composition, [])) / total))
            line_str.append(",".join(species_ids_by_composition_by_orthogroup_id[orthogroup_id].get(composition, ['NA'])))
        lines_int.append("\t".join([str(l) for l in line_int]))
        lines_str.append("\t".join(line_str))
    write_lines(lines_int, "orthogroups.genotype_count.tsv")
    write_lines(lines_str, "orthogroups.genotype_species.tsv")

    header = "\t".join(['species_id', 'transcripts', 'HOM_both', 'HET_both', 'HOM_M', 'HOM_F'])
    lines = [header]
    for species_id in count_by_composition_by_species_ids:
        total = sum([count for composition, count in count_by_composition_by_species_ids[species_id].items()])
        line = [species_id, total]
        for composition in ['HOM_both', 'HET_both', 'HOM_M', 'HOM_F']:
            line.append("%.2f" % (count_by_composition_by_species_ids[species_id].get(composition, 0) / total))
        lines.append("\t".join([str(l) for l in line]))
    write_lines(lines, "species.composition_proportion.tsv")

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        metadataObj = parse_metadata(args['--metadata'])
        seqdataObj = parse_sco_loci_metrics(args['--sco'], metadataObj)
        seqdataObj = parse_phy(args['--phy'], seqdataObj, metadataObj)
        analyse_transcriptObjs(seqdataObj, metadataObj)
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    SEXES = ['M', 'F', 'X']
    main()
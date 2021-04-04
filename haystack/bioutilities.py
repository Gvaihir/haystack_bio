import math
import os
import subprocess
import numpy as np
from bx.intervals.intersection import Intersecter, Interval
from typing import List


try:
    import cPickle as cp
except:
    import pickle as cp

try:
    import pysam
except:
    pass


def smooth(raw_values: List, window: str='hanning', window_len: int=200) -> List:
    allowed_window_types = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
    if window not in allowed_window_types:
        raise ValueError(f"smooth(): provided window {window}. Must be one of {allowed_window_types}")
    reversed_first_window = raw_values[window_len - 1:0:-1]
    reversed_last_window = raw_values[-1:-window_len:-1]
    concat_window_parts = np.r_[reversed_first_window, raw_values, reversed_last_window]
    if window == 'flat':
        window = np.ones(window_len, 'd')
    else:
        window = eval('np.' + window + '(window_len)')
    convolution = np.convolve(window / window.sum(), concat_window_parts, mode='valid')
    return convolution[int(window_len/2):int(-window_len/2+1)]


def mask_lower(nucleotide: str) -> str:
    return nucleotide if nucleotide.isupper() else 'N'

def chunks(l, n):    return [l[i:i + n] for i in range(0, len(l), n)]


class constant_minus_one_dict(dict):
    def __missing__(self, key):
        return -1


class constant_n_dict(dict):
    def __missing__(self, key):
        return 'n'


nt2int = constant_minus_one_dict({'a': 0, 'c': 1, 'g': 2, 't': 3})
int2nt = constant_n_dict({0: 'a', 1: 'c', 2: 'g', 3: 't'})
nt_complement = dict({'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'})


def read_sequence_from_fasta(fin, bpstart, bpend, line_length=50.0):
    bpstart = bpstart - 1

    fin.seek(0)
    fin.readline()  # read the first line; the pointer is at the second line

    nbp = bpend - bpstart
    offset = int(
        bpstart + math.floor(bpstart / line_length))  # assuming each line contains 50 characters; add 1 offset per line

    if offset > 0:
        fin.seek(int(offset), 1)

    seq = fin.read(nbp + int(math.floor(nbp / line_length)) + 1)
    seq = seq.replace('\n', '')

    if len(seq) < nbp:
        print('Coordinate out of range:', bpstart, bpend)

    return seq[0:nbp].lower()


def fasta_iterator(handle):
    while True:
        line = handle.readline()
        if line == "": return
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError("Fasta files should start with '>'")
        else:
            descr = line[1:].rstrip()
            id = descr.split()[0]
            name = id

        lines = []
        line = handle.readline()
        while True:
            if not line: break
            if line[0] == ">": break
            lines.append(line.rstrip().replace(" ", "").replace("\r", ""))
            line = handle.readline()

        yield "".join(lines)

        if not line: return  # End



# intersecter
class Coordinates_Intersecter:
    def __init__(self, coordinates):
        self.interval_tree = dict()

        for c in coordinates:
            if c.chr_id not in self.interval_tree:
                self.interval_tree[c.chr_id] = Intersecter()

            self.interval_tree[c.chr_id].add_interval(Interval(c.bpstart, c.bpend, c))

    def find_intersections(self, c):
        coords_in_common = list()
        if self.interval_tree.has_key(c.chr_id):
            coords_hits = self.interval_tree[c.chr_id].find(c.bpstart - 1, c.bpend + 1)

            for coord_hit in coords_hits:
                coords_in_common.append(coord_hit.value)

        return coords_in_common











def build_motif_in_seq_matrix(bed_filename, genome_directory, meme_motifs_filename, bg_filename, genome_mm=True,
                              temp_directory='./', mask_repetitive=False, p_value=1.e-4, check_only_presence=False):
    print('Loading coordinates  from bed')
    target_coords = Coordinate.bed_to_coordinates(bed_filename)

    print('Initialize Genome')
    if genome_mm:
        genome = Genome(genome_directory)
    else:
        genome = Genome(genome_directory)

    print('Initilize Fimo and load motifs')
    fimo = Fimo(meme_motifs_filename, bg_filename, temp_directory=temp_directory, p_value=p_value)

    print('Initialize the matrix')
    motifs_in_sequences_matrix = np.zeros((len(target_coords), len(fimo.motif_names)))

    for idx_seq, c in enumerate(target_coords):
        seq = genome.extract_sequence(c, mask_repetitive)
        print(idx_seq, len(target_coords))
        if check_only_presence:
            motifs_in_sequences_matrix[idx_seq, fimo.extract_motifs(seq, report_mode='indexes_set')] = 1
        else:
            motifs_in_sequences_matrix[idx_seq, :] += fimo.extract_motifs(seq, report_mode='fq_array')

    return motifs_in_sequences_matrix, fimo.motif_names, fimo.motif_ids


def build_motif_profile(target_coords, genome, meme_motifs_filename, bg_filename, genome_mm=True, temp_directory='./',
                        mask_repetitive=False, p_value=1.e-4, check_only_presence=False):
    # print 'Initilize Fimo and load motifs'
    fimo = Fimo(meme_motifs_filename, bg_filename, temp_directory=temp_directory, p_value=p_value)

    # print 'Allocate memory'
    if check_only_presence:
        motifs_in_sequences_profile = np.zeros(len(fimo.motif_names))
    else:
        motifs_in_sequences_profile = {'fq': np.zeros((len(target_coords), len(fimo.motif_names))),
                                       'presence': np.zeros(len(fimo.motif_names))}
    for idx_seq, c in enumerate(target_coords):
        seq = genome.extract_sequence(c, mask_repetitive)
        print(idx_seq, len(target_coords))
        if check_only_presence:
            motifs_in_sequences_profile[fimo.extract_motifs(seq, report_mode='indexes_set')] += 1
        else:
            motifs_in_sequences = fimo.extract_motifs(seq, report_mode='fq_and_presence')
            # print motifs_in_sequences,motifs_in_sequences['presence'],motifs_in_sequences['fq']
            motifs_in_sequences_profile['presence'][list(motifs_in_sequences['presence'])] += 1
            motifs_in_sequences_profile['fq'][idx_seq, :] += motifs_in_sequences['fq']

    return motifs_in_sequences_profile, fimo.motif_names, fimo.motif_ids


'''
given a set of coordinates in a bed file and a bam/sam file calculate the profile matrix
'''


def calculate_profile_matrix_bed_bam(bed_filename, sam_filename, window_size=5000, resolution=50, fragment_length=200,
                                     use_strand=False, return_coordinates=False):
    cs = Coordinate.bed_to_coordinates(bed_filename)
    cs = Coordinate.coordinates_of_intervals_around_center(cs, window_size)
    samfile = pysam.Samfile(sam_filename)
    n_bins = (window_size / resolution) + 1

    profile_matrix = np.zeros((len(cs), n_bins))

    for idx_c, c in enumerate(cs):
        n_start = max(0, c.bpcenter - window_size / 2)
        n_end = c.bpcenter + window_size / 2

        for rd in samfile.fetch(c.chr_id, n_start, n_end):

            if rd.is_reverse:
                rd_st = rd.positions[-1] - fragment_length - 1
                rd_end = rd.positions[-1]
            else:
                rd_st = rd.positions[0]
                rd_end = max(rd.positions[-1], rd.positions[1] + fragment_length)

            bin_idx_st = max(0, 1 + (rd_st - n_start) / resolution)
            bin_idx_en = min(n_bins, 1 + (rd_end - n_start) / resolution)

            # print rd.positions[1],rd.positions[-1],rd_st, rd_end,rd_end-rd_st,rd_st-n_start,rd_end-n_start,bin_idx_st,bin_idx_en

            profile_matrix[idx_c, bin_idx_st:bin_idx_en] += 1

        if use_strand and c.strand == '-':
            profile_matrix[idx_c, :] = profile_matrix[idx_c, ::-1]

    if return_coordinates:
        return profile_matrix, samfile.mapped, cs
    else:
        return profile_matrix, samfile.mapped


def extract_bg_from_bed(bed_filename, genome_directory, bg_filename, genome_mm=True):
    acgt_fq = {'a': 0.0, 'c': 0.0, 'g': 0.0, 't': 0.0}
    total = 0

    print('Loading coordinates  from bed')
    target_coords = Coordinate.bed_to_coordinates(bed_filename)

    print('Initialize Genome')
    if genome_mm:
        genome = Genome(genome_directory)
    else:
        genome = Genome(genome_directory)

    for idx_seq, c in enumerate(target_coords):
        seq = genome.extract_sequence(c)
        for nt in ['a', 'c', 't', 'g']:
            acgt_fq[nt] += seq.count(nt)

    total = sum(acgt_fq.values())

    for nt in ['a', 'c', 't', 'g']:
        acgt_fq[nt] /= total

    print(acgt_fq)
    with open(bg_filename, 'w+') as out_file:
        for nt in ['a', 'c', 't', 'g']:
            out_file.write('%s\t%1.4f\n' % (nt, acgt_fq[nt]))


# hgWiggle wrapper
def read_from_wig(c, wig_path, wig_mask='.phastCons44way.hg18.compiled', only_average=False):
    position = c.chr_id + ':' + str(c.bpstart) + '-' + str(c.bpend)
    wig_file = os.path.join(wig_path, c.chr_id + wig_mask)

    if only_average:
        command = ' '.join(['hgWiggle', '-position=' + position, wig_file, "-doStats | sed  -e '1,3d' | cut -f4,10"])
    else:
        command = ' '.join(['hgWiggle', '-position=' + position, wig_file, " -rawDataOut "])

    wig_process = subprocess.Popen(command, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output = wig_process.communicate()[0]

    if only_average:
        return tuple(output.split())
    else:
        values = output.split()
        return len(values), values





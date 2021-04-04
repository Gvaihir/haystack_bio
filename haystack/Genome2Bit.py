class Genome2Bit:

    def __init__(self, genome_2bit_file, verbose=False):
        self.genome = TwoBitFile(open(genome_2bit_file, 'rb'))
        self.chr_len = dict()
        self.verbose = verbose

        for chr_id in self.genome.keys():
            self.chr_len[chr_id] = len(self.genome[chr_id])
        if verbose:
            print('Genome initializated')

    def extract_sequence(self, coordinate, mask_repetitive=False):
        if mask_repetitive:
            seq = ''.join(
                [mask(c) for c in self.genome[coordinate.chr_id][coordinate.bpstart - 1:coordinate.bpend]]).lower()
        else:
            seq = self.genome[coordinate.chr_id][coordinate.bpstart - 1:coordinate.bpend].lower()

        if coordinate.strand == '-':
            return Sequence.reverse_complement(seq)
        else:
            return seq

    def estimate_background(self):
        counting = {'a': .0, 'c': .0, 'g': .0, 't': .0}
        all = 0.0

        for chr_id in self.genome.keys():
            if self.verbose:
                start_time = time.time()
                print('Counting on:', chr_id)

            for nt in counting.keys():
                count_nt = self.genome[chr_id][:].lower().count(nt)
                counting[nt] += count_nt
                all += count_nt

            print('elapsed:', time.time() - start_time)

        if self.verbose:
            print(counting)

        for nt in counting.keys():
            counting[nt] /= all

        return counting

    def write_meme_background(self, filename):
        counting = self.estimate_background()
        with open(filename, 'w+') as outfile:
            for nt in counting.keys():
                outfile.write('%s\t%2.4f\n' % (nt, counting[nt]))

    def write_chr_len(self, filename):
        with open(filename, 'w+') as outfile:
            for chr_id in self.genome.keys():
                outfile.write('%s\t%s\n' % (chr_id, self.chr_len[chr_id]))


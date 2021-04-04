class Genome:

    def __init__(self, genome_directory, release='ND', verbose=False):
        self.chr = dict()
        self.genome_directory = genome_directory
        self.release = release
        self.chr_len = dict()
        self.verbose = verbose

        for infile in glob.glob(os.path.join(genome_directory, '*.fa')):
            mm_filename = infile.replace('.fa', '.mm')

            filename = infile.replace(genome_directory, '').replace('.fa', '')
            chr_id = os.path.basename(infile).replace('.fa', '')

            if not os.path.isfile(mm_filename):
                if verbose:
                    print(
                        'Missing:' + chr_id + ' generating memory mapped file (This is necessary only the first time) \n')
                with open(infile) as fi:
                    with open(os.path.join(genome_directory, chr_id + '.mm'), 'w+') as fo:
                        # skip header
                        fi.readline()
                        for line in fi:
                            fo.write(line.strip())
                    if verbose:
                        print('Memory mapped file generated for:', chr_id)

            with open(mm_filename, 'r') as f:
                self.chr[chr_id] = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                self.chr_len[chr_id] = len(self.chr[chr_id])
                if verbose:
                    print("Chromosome:%s Read" % chr_id)

        if verbose:
            print('Genome initializated')

    def extract_sequence(self, coordinate, mask_repetitive=False):
        if mask_repetitive:
            seq = ''.join(
                [mask(c) for c in self.chr[coordinate.chr_id][coordinate.bpstart - 1:coordinate.bpend]]).lower()
        else:
            seq = self.chr[coordinate.chr_id][coordinate.bpstart - 1:coordinate.bpend].lower()

        if coordinate.strand == '-':
            return Sequence.reverse_complement(seq).lower()
        else:
            return seq

    def estimate_background(self):
        counting = {'a': .0, 'c': .0, 'g': .0, 't': .0}
        all = 0.0
        for chr_id in self.chr.keys():
            if self.verbose:
                print('Counting on:', chr_id)

            for nt in counting.keys():
                count_nt = self.chr[chr_id][:].lower().count(nt)
                counting[nt] += count_nt
                all += count_nt

        if self.verbose:
            print(counting)

        for nt in counting.keys():
            counting[nt] /= all

        return counting

class Sequence:

    def __init__(self, seq=''):
        self.seq = string.lower(seq)

    def set_bg_model(self, ACGT_probabilities):
        self.bg_model = rv_discrete(name='bg', values=([0, 1, 2, 3], ACGT_probabilities))

    def __str__(self):
        return self.seq

    def __len__(self):
        return len(self.seq)

    # seq_complement
    def _reverse_complement(self):
        return "".join([nt_complement[c] for c in self.seq[-1::-1]])

    reverse_complement = property(_reverse_complement)

    @classmethod
    def reverse_complement(self, seq):
        return "".join([nt_complement[c] for c in seq[-1::-1]])

    @classmethod
    def generate_random(cls, n,
                        bgmodel=rv_discrete(name='bg', values=([0, 1, 2, 3], [0.2955, 0.2045, 0.2045, 0.2955]))):
        int_seq = cls.bg_model.rvs(size=n)
        return ''.join([int2nt[c] for c in int_seq])


class Seq_set(dict):
    def add_name(self, name):
        self.name = name

    def add_seq(self, seq, cord):
        self[cord] = seq

    def remove_seq_from_cord(self, cord):
        if cord in self.keys():
            del self[cord]
        else:
            print('sequence not found')


class Genome:
    def __init__(self, genome_directory, release='ND', verbose=False):
        self.chr = dict()
        self.genome_directory = genome_directory
        self.release = release
        self.chr_len = dict()
        self.verbose = verbose

        for infile in glob.glob(os.path.join(genome_directory, '*.fa')):
            try:

                chr_id = os.path.basename(infile).replace('.fa', '')
                self.chr[chr_id] = open(infile, 'r')

                self.chr_len[chr_id] = 0

                self.chr[chr_id].readline()
                for line in self.chr[chr_id]:
                    self.chr_len[chr_id] += len(line.strip())

                if verbose:
                    print('Read:' + infile)
            except:
                if verbose:
                    print('Error, not loaded:', infile)

        if verbose:
            print('Genome initializated')

    def estimate_background(self):
        counting = {'a': .0, 'c': .0, 'g': .0, 't': .0}
        all = 0.0

        for chr_id in self.chr.keys():
            if self.verbose:
                print('Counting on:', chr_id)

            self.chr[chr_id].seek(0)
            self.chr[chr_id].readline()

            for line in self.chr[chr_id]:
                for nt in counting.keys():
                    count_nt = line.lower().count(nt)
                    counting[nt] += count_nt
                    all += count_nt

        if self.verbose:
            print(counting)

        for nt in counting.keys():
            counting[nt] /= all

        return counting

    def extract_sequence(self, coordinate, mask_repetitive=False, line_length=50.0):
        if not self.chr.has_key(coordinate.chr_id):
            if self.verbose:
                print("Warning: chromosome %s not present in the genome" % coordinate.chr_id)
        else:

            bpstart = coordinate.bpstart - 1
            bpend = coordinate.bpend

            self.chr[coordinate.chr_id].seek(0)
            self.chr[coordinate.chr_id].readline()

            nbp = bpend - bpstart
            offset = int(bpstart + math.floor(bpstart / line_length)) - 1

            if offset > 0:
                self.chr[coordinate.chr_id].seek(offset, 1)

            seq = self.chr[coordinate.chr_id].read(nbp + int(math.floor(nbp / line_length)) + 1)
            seq = seq.replace('\n', '')

            if len(seq) < nbp:
                if self.verbose:
                    print('Warning: coordinate out of range:', bpstart, bpend)

            if mask_repetitive:
                seq = ''.join([mask(c) for c in seq[0:nbp]]).lower()
            else:
                seq = seq[0:nbp].lower()

            if coordinate.strand == '-':
                return Sequence.reverse_complement(seq).lower()
            else:
                return seq

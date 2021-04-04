class Ngram:
    def __init__(self, ngram_length=4, alphabet_size=4):
        # print 'Ngram initialization'
        # build useful dictionary
        self.alphabet_size = alphabet_size
        self.ngram_length = ngram_length

        all_ngram = map(self.int2ngram, range(self.alphabet_size ** self.ngram_length),
                        [self.ngram_length] * (self.alphabet_size ** self.ngram_length),
                        [self.alphabet_size] * (self.alphabet_size ** self.ngram_length))
        self.ngram2index = dict()
        self.index2ngram = dict()
        self.non_redundant_ngram2index = dict()
        self.non_redundant_index2ngram = dict()
        self.ngram_rev_complement = dict()
        self.number_of_ngrams = len(all_ngram)

        for idx, ngram in enumerate(all_ngram):
            self.ngram2index[ngram] = idx
            self.ngram_rev_complement[ngram] = Sequence.reverse_complement(ngram)
            self.index2ngram[idx] = ngram

        for idx, ngram in enumerate([all_ngram[i] for i in self.non_redundant_idx()]):
            self.non_redundant_ngram2index[ngram] = idx
            self.non_redundant_index2ngram[idx] = ngram

        self.number_of_non_redundant_ngrams = len(self.non_redundant_ngram2index)

    def int2ngram(self, idx, ngram_length, alphabet_size):
        l = []
        for _ in range(ngram_length):
            l.append(int2nt[idx % alphabet_size])
            idx /= alphabet_size

        return "".join(l)[-1::-1]

    def non_redundant_idx(self):
        ngram_taken = set()
        non_redundant_idxs = []
        for idx in range(self.number_of_ngrams):
            n = self.index2ngram[idx]
            if (self.ngram2index[n] in ngram_taken) or (self.ngram2index[self.ngram_rev_complement[n]] in ngram_taken):
                pass
            else:
                non_redundant_idxs.append(idx)
                ngram_taken.add(self.ngram2index[n])
                ngram_taken.add(self.ngram2index[self.ngram_rev_complement[n]])

        return tuple(non_redundant_idxs)

    def build_ngram_fq_vector(self, seq):

        ngram_vector = zeros(self.number_of_ngrams)

        for i in xrange(len(seq) - self.ngram_length + 1):
            try:
                ngram_vector[self.ngram2index[seq[i:i + self.ngram_length]]] += 1
            except:
                pass
            # ngram_vector[self.ngram2index[self.ngram_rev_complement[seq[i:i+self.ngram_length]]]]+=1

        return ngram_vector

    def build_ngram_fq_vector_non_redundant(self, seq):

        ngram_vector = zeros(self.number_of_non_redundant_ngrams)

        for i in xrange(len(seq) - self.ngram_length + 1):

            try:
                ngram_vector[self.non_redundant_ngram2index[seq[i:i + self.ngram_length]]] += 1
            except:
                pass
            try:
                ngram_vector[
                    self.non_redundant_ngram2index[self.ngram_rev_complement[seq[i:i + self.ngram_length]]]] += 1
            except:
                pass

        return ngram_vector

    def build_ngrams_fq_matrix(self, seq_set, non_redundant=True):
        if non_redundant:
            ngram_matrix = zeros((len(seq_set), self.number_of_non_redundant_ngrams))
            for idx_seq, seq in enumerate(seq_set):
                ngram_matrix[idx_seq, :] = self.build_ngram_fq_vector_non_redundant(seq)
        else:
            ngram_matrix = zeros((len(seq_set), self.number_of_ngrams))
            for idx_seq, seq in enumerate(seq_set):
                ngram_matrix[idx_seq, :] = self.build_ngram_fq_vector(seq)

        return ngram_matrix, np.array([len(seq) for seq in seq_set], ndmin=2)

    def build_ngram_profile_vector_non_redundant(self, seq):

        profile = zeros(len(seq) - self.ngram_length + 1)

        for i in xrange(len(seq) - self.ngram_length + 1):

            try:
                profile[i] = self.non_redundant_ngram2index[seq[i:i + self.ngram_length]]
            except:
                pass
            try:
                profile[i] = self.non_redundant_ngram2index[self.ngram_rev_complement[seq[i:i + self.ngram_length]]]
            except:
                pass

        return profile

    def build_ngrams_profile_matrix_non_redundant(self, seq_set):
        profile_matrix = zeros((len(seq_set), len(seq_set[0]) - self.ngram_length + 1))
        for idx_seq, seq in enumerate(seq_set):
            profile_matrix[idx_seq, :] = self.build_ngram_profile_vector_non_redundant(seq)

        return profile_matrix

    def save_to_file(self, filename=None):
        if not filename:
            filename = 'ng' + str(self.ngram_length)
        with open(filename, 'wb+') as outfile:
            print('saving...')
            cp.dump(self, outfile, 2)
            print('done')

    @classmethod
    def load_from_file(cls, filename):
        with open(filename, 'rb') as infile:
            return cp.load(infile)

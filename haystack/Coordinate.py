class Coordinate(object):

    def __init__(self, chr_id, bpstart, bpend, name=None, score=None, strand=None):
        self.chr_id = chr_id
        self.bpstart = bpstart
        self.bpend = bpend
        self.name = name
        self.score = score
        self.strand = strand

    def chr_id2_ord(self):
        id = self.chr_id.replace('chr', '')
        try:
            id = int(id)
        except:
            cs = 0
            for c in id:
                cs += ord(c)
            id = cs

        return id



    def bpcenter(self):
        return (self.bpstart + self.bpend) / 2

    def __eq__(self, other):
        return (self.chr_id == other.chr_id) & (self.bpstart == other.bpstart) & (self.bpend == other.bpend)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if self.chr_id2_ord < other.chr_id2_ord:
            return True
        elif self.chr_id2_ord > other.chr_id2_ord:
            return False
        elif self.bpstart < other.bpstart:
            return True
        elif self.bpstart > other.bpstart:
            return False
        else:
            return False

    def __gt__(self, other):
        return not (self.__lt__(other) or self.__eq__(other))

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)

    def __hash__(self):
        return hash((self.chr_id, self.bpstart, self.bpend))

    def __str__(self):
        return self.chr_id + ':' + str(self.bpstart) + '-' + str(self.bpend) + (
            ' ' + self.name if self.name else '') + (' ' + str(self.score) if self.score else '') + (
                   ' ' + self.strand if self.strand else '')

    def __repr__(self):
        return self.chr_id + ':' + str(self.bpstart) + '-' + str(self.bpend) + (
            ' ' + self.name if self.name else '') + (' ' + str(self.score) if self.score else '') + (
                   ' ' + self.strand if self.strand else '')

    def __len__(self):
        return self.bpend - self.bpstart

    def __and__(self, other):
        if self.bpend < other.bpstart or other.bpend < self.bpstart or self.chr_id != other.chr_id:
            return None
        else:
            return Coordinate(self.chr_id, max(self.bpstart, other.bpstart), min(self.bpend, other.bpend))

    def upstream(self, offset=2000):
        return Coordinate(self.chr_id, self.bpend + 1, self.bpend + offset)

    def downstream(self, offset=2000):
        return Coordinate(self.chr_id, max(0, self.bpstart - offset), self.bpstart - 1)

    @classmethod
    def bed_to_coordinates(cls, bed_filename, header_lines=0, cl_chr_id=0, cl_bpstart=1, cl_bpend=2, cl_name=3,
                           cl_score=4, cl_strand=5):
        with open(bed_filename, 'r') as infile:
            coordinates = list()

            for _ in range(header_lines):
                infile.readline()

            for line in infile:
                line = line.strip()
                try:
                    coord = line.split('\t')
                    chr_id = coord[cl_chr_id]
                    bpstart = coord[cl_bpstart]
                    bpend = coord[cl_bpend]
                    try:
                        name = coord[cl_name]
                    except:
                        name = 'ND'
                    try:
                        score = float(coord[cl_score])
                    except:
                        score = 0
                    try:
                        strand = coord[cl_strand]
                    except:
                        strand = 'ND'
                    coordinates.append(Coordinate(str(chr_id), int(bpstart), int(bpend), name, score, strand))
                except:
                    print('Skipping line:', line)

            return coordinates

    @classmethod
    def bed_to_coordinates_dict(cls, bed_file, header_lines=0, cl_chr_id=0, cl_bpstart=1, cl_bpend=2, cl_name=3,
                                cl_score=4, cl_strand=5):
        with open(bed_file, 'r') as infile:
            coordinates = dict()

            for _ in range(header_lines):
                infile.readline()

            for idx, line in enumerate(infile):
                coord = line.split()
                chr_id = coord[cl_chr_id]
                bpstart = coord[cl_bpstart]
                bpend = coord[cl_bpend]
                try:
                    name = coord[cl_name]
                except:
                    name = 'ND'
                try:
                    score = coord[cl_score]
                except:
                    score = 'ND'
                try:
                    strand = coord[cl_strand]
                except:
                    strand = 'ND'
                c = Coordinate(str(chr_id), int(bpstart), int(bpend), name, score, strand)
                print(c)
                coordinates[c] = idx + 1

            return coordinates

    @classmethod
    def coordinates_to_bed(cls, coordinates, bed_file, minimal_format=False):
        with open(bed_file, 'w+') as outfile:
            for c in coordinates:
                if minimal_format:
                    outfile.write('%s\t%d\t%d\n' % (c.chr_id, c.bpstart, c.bpend))
                else:
                    outfile.write('%s\t%d\t%d\t%s\t%f\t%s\n' % (
                    c.chr_id, c.bpstart, c.bpend, c.name if c.name else 'ND', c.score if c.score else 0,
                    c.strand if c.strand else '+'))

    @classmethod
    def coordinates_to_nscore_format(cls, coordinates, bed_file):
        with open(bed_file, 'w+') as outfile:
            for c in coordinates:
                outfile.write('%s\t%d\n' % (c.chr_id, c.bpcenter))

    @classmethod
    def coordinates_from_interval(cls, chr_id, interval):
        return Coordinate(chr_id, interval.start, interval.end, )

    @classmethod
    def coordinates_to_fasta(cls, coordinates, fasta_file, genome, chars_per_line=50, mask_repetitive=False):
        with open(fasta_file, 'w+') as outfile:
            for c in coordinates:
                seq = genome.extract_sequence(c, mask_repetitive)
                outfile.write('>' + str(c) + '\n' + '\n'.join(chunks(seq, chars_per_line)) + '\n')

    @classmethod
    def calculate_intersection(cls, coords1, coords2, build_matrix=False):

        coords_in_common = set()
        intersection_indexes = set()
        if build_matrix:
            intersection_matrix = sparse_matrix((len(coords1), len(coords2)))
        coord_to_row_index = dict()
        row_index = 0
        interval_tree = dict()

        # Build the interval tree on the first set of coordinates
        for c in coords1:
            if c.chr_id not in interval_tree:
                interval_tree[c.chr_id] = Intersecter()

            interval_tree[c.chr_id].add_interval(Interval(c.bpstart, c.bpend))
            coord_to_row_index[c] = row_index
            row_index += 1

        # Calculating the intersection
        # for each coordinates on the second set check intersection and fill the matrix
        for cl_index, c in enumerate(coords2):
            if interval_tree.has_key(c.chr_id):
                coords_hits = interval_tree[c.chr_id].find(c.bpstart, c.bpend)
                # coords_in_common+=coords_hits
                for coord_hit in coords_hits:
                    c_to_add = Coordinate.coordinates_from_interval(c.chr_id, coord_hit)
                    coords_in_common.add(c_to_add)
                    row_index = coord_to_row_index[c_to_add]
                    intersection_indexes.add(row_index)

                    if build_matrix:
                        intersection_matrix[row_index, cl_index] += 1

        if build_matrix:
            return list(coords_in_common), list(intersection_indexes), intersection_matrix
        else:
            return list(coords_in_common), list(intersection_indexes)

    bpcenter = property(bpcenter)
    chr_id2_ord = property(chr_id2_ord)

    @classmethod
    def coordinates_of_intervals_around_center(cls, coords, window_size):
        half_window = window_size / 2
        return [Coordinate(c.chr_id, c.bpcenter - half_window, c.bpcenter + half_window, strand=c.strand, name=c.name,
                           score=c.score) for c in coords]
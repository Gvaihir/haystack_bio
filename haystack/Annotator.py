from Coordinate import Coordinate
import numpy as np

class Annotator(object):

    def __init__(self, input_filename, annotations_filenames, annotation_names=None):

        self.annotations_filenames = annotations_filenames
        self.input_filename = input_filename

        # print(annotation_names)

        if annotation_names is None:
            self.annotation_names = annotations_filenames
        else:
            self.annotation_names = annotation_names

        assert len(annotations_filenames) == len(annotation_names)

        # associate to each name a prime number for the multiple annotation trick..
        self.annotation_names_to_prime = {name: prime for (name, prime) in
                                          zip(self.annotation_names, self._primes(len(self.annotation_names)))}
        print(self.annotation_names_to_prime)

        self.input_coordinates = Coordinate.bed_to_coordinates(input_filename)

    def annotate(self):
        self.annotation_track = np.ones(len(self.input_coordinates), dtype=np.int)

        self.interval_tree = dict()
        self.coord_to_row_index = dict()
        self.row_index = 0

        # Build the interval Tree
        for c in self.input_coordinates:
            if c.chr_id not in self.interval_tree:
                self.interval_tree[c.chr_id] = Intersecter()

            self.interval_tree[c.chr_id].add_interval(Interval(c.bpstart, c.bpend))
            self.coord_to_row_index[c] = self.row_index
            self.row_index += 1

        for idx, bed_filename in enumerate(self.annotations_filenames):
            coordinates = Coordinate.bed_to_coordinates(bed_filename)

            prime_number = self.annotation_names_to_prime[self.annotation_names[idx]]
            for idx_intersection in self._intersection_indexes(coordinates):
                self.annotation_track[idx_intersection] *= prime_number

    def coordinates_by_annotation_name(self, annotation_name):
        if annotation_name == 'out':
            return [self.input_coordinates[idx] for idx, value in enumerate(self.annotation_track) if value == 1]
        else:
            prime_number = self.annotation_names_to_prime[annotation_name]
            return [self.input_coordinates[idx] for idx, value in enumerate(self.annotation_track) if
                    (value % prime_number) == 0]


    def save_annotation_track_bed(self, filename):
        with open(filename, 'w+') as outfile:

            for idx, c in enumerate(self.input_coordinates):
                annotated = False

                line = '%s\t%s\t%d\t' % (c.chr_id, c.bpstart, c.bpend,)
                for name in self.annotation_names_to_prime.keys():
                    if self.annotation_track[idx] % self.annotation_names_to_prime[name] == 0:
                        line += name + ','
                        annotated = True

                line = line[:-1]
                outfile.write(line + '\n')

        print('Annotation track saved to:', filename)

    def _intersection_indexes(self, coordinates):
        intersection_indexes = set()

        for cl_index, c in enumerate(coordinates):
            if self.interval_tree.has_key(c.chr_id):

                coords_hits = self.interval_tree[c.chr_id].find(c.bpstart, c.bpend)
                for coord_hit in coords_hits:
                    intersection_indexes.add(
                        self.coord_to_row_index[Coordinate.coordinates_from_interval(c.chr_id, coord_hit)])

        return intersection_indexes

    def _gen_primes(self):

        D = {}
        q = 2

        while True:
            if q not in D:
                yield q
                D[q * q] = [q]
            else:

                for p in D[q]:
                    D.setdefault(p + q, []).append(p)
                del D[q]

            q += 1

    def _primes(self, n):
        primes = self._gen_primes()
        return [next(primes) for i in range(n)]

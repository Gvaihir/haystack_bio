from Coordinate import Coordinate

class Gene(object):
    regions = [8000, 2000, 0, 1000]

    def __init__(self, access, name, coordinate, regions=None, exons=None, introns=None) -> None:
        self.access = access
        self.c = coordinate
        self.name = name
        self.exons = exons
        self.introns = introns

        if regions:
            self.regions = regions

    def __str__(self):
        return self.name + '_' + self.access + '_' + str(self.c)

    def __repr__(self):
        return self.name + '_' + self.access + '_' + str(self.c)

    def tss(self):
        if self.c.strand == '-':
            return self.c.bpend
        else:
            return self.c.bpstart

    def tes(self):
        if self.c.strand == '+':
            return self.c.bpend
        else:
            return self.c.bpstart

    def end(self):
        if self.c.strand == '+':
            return self.c.bpend
        else:
            return self.c.bpstart

    def distal_c(self):
        if self.c.strand == '+':
            return Coordinate(self.c.chr_id, self.tss - self.regions[0], self.tss - self.regions[1] - 1, strand='+')
        else:
            return Coordinate(self.c.chr_id, self.tss + self.regions[1] + 1, self.tss + self.regions[0], strand='-')

    def promoter_c(self):
        if self.c.strand == '+':
            return Coordinate(self.c.chr_id, self.tss - self.regions[1], self.tss + self.regions[2] - 1, strand='+')
        else:
            return Coordinate(self.c.chr_id, self.tss - self.regions[2] + 1, self.tss + self.regions[1], strand='-')

    def intra_c(self):
        if self.c.strand == '+':
            return Coordinate(self.c.chr_id, self.tss + self.regions[2], self.end + self.regions[3], strand='+')
        else:
            return Coordinate(self.c.chr_id, self.end - self.regions[3], self.tss - self.regions[2], strand='-')

    def full_c(self):
        if self.c.strand == '+':
            return Coordinate(self.c.chr_id, self.tss - self.regions[0], self.end + self.regions[3], strand='+')
        else:
            return Coordinate(self.c.chr_id, self.end - self.regions[3], self.tss + self.regions[0], strand='-')

    
    def load_from_annotation(self, gene_annotation_file, load_exons_introns_info=False, header_lines=1, regions=regions,
                             format='txt'):
        genes_list = []
        with open(gene_annotation_file, 'r') as genes_file:

            for _ in range(header_lines):
                genes_file.readline()

            if format == 'txt':
                idx_chr_id = 2
                idx_bpstart = 4
                idx_bpend = 5
                idx_strand = 3
                idx_access = 1
                idx_name = -4
                idx_exon_starts = 9
                idx_exon_ends = 10
            else:
                idx_chr_id = 0
                idx_bpstart = 1
                idx_bpend = 2
                idx_strand = 5
                idx_access = 3
                idx_name = 3
                idx_exon_starts = 10
                idx_exon_ends = 11

            for gene_line in genes_file:
                fields = gene_line.split('\t')
                chr_id = fields[idx_chr_id]
                bpstart = int(fields[idx_bpstart])
                bpend = int(fields[idx_bpend])
                strand = fields[idx_strand]
                access = fields[idx_access]
                name = fields[idx_name]
                c = Coordinate(chr_id, bpstart, bpend, strand=strand)

                if load_exons_introns_info:
                    exon_starts = map(int, fields[idx_exon_starts].split(',')[:-1])
                    exon_ends = map(int, fields[idx_exon_ends].split(',')[:-1])

                    exons = [Coordinate(chr_id, bpstart, bpend, strand=strand) for bpstart, bpend in
                             zip(exon_starts, exon_ends)]
                    introns = [Coordinate(chr_id, bpstart + 1, bpend - 1, strand=strand) for bpstart, bpend in
                               zip(exon_ends[:-1], exon_starts[1:])]
                    genes_list.append(Gene(access, name, c, exons=exons, introns=introns, regions=regions))
                else:
                    genes_list.append(Gene(access, name, c, regions=regions))

        return genes_list

    
    def exons_from_annotations(self, gene_annotation_file, header_lines=1):
        exons_list = []

        with open(gene_annotation_file, 'r') as genes_file:

            for _ in range(header_lines):
                genes_file.readline()

            for gene_line in genes_file:
                fields = gene_line.split('\t')
                chr_id = fields[2]
                strand = fields[3]
                access = fields[1]
                exon_starts = map(int, fields[9].split(',')[:-1])
                exon_ends = map(int, fields[10].split(',')[:-1])
                exons = [Coordinate(chr_id, bpstart, bpend, strand=strand, name=access) for bpstart, bpend in
                         zip(exon_starts, exon_ends)]
                exons_list += exons

            return exons_list

    
    def introns_from_annotations(self, gene_annotation_file, header_lines=1):
        introns_list = []

        with open(gene_annotation_file, 'r') as genes_file:

            for _ in range(header_lines):
                genes_file.readline()

            for gene_line in genes_file:
                fields = gene_line.split('\t')
                chr_id = fields[2]
                strand = fields[3]
                access = fields[1]
                exon_starts = map(int, fields[9].split(',')[:-1])
                exon_ends = map(int, fields[10].split(',')[:-1])
                introns = [Coordinate(chr_id, bpstart - 1, bpend, strand=strand, name=access) for bpstart, bpend in
                           zip(exon_ends[:-1], exon_starts[1:])]
                introns_list += introns

            return introns_list

    
    def genes_coordinates_from_annotations(self, gene_annotation_file, header_lines=1):
        genes_coordinates = []
        with open(gene_annotation_file, 'r') as genes_file:

            for _ in range(header_lines):
                genes_file.readline()

            for gene_line in genes_file:
                fields = gene_line.split('\t')
                chr_id = fields[2]
                bpstart = int(fields[4])
                bpend = int(fields[5])
                strand = fields[3]
                access = fields[1]
                name = fields[-4]

                genes_coordinates.append(Coordinate(chr_id, bpstart, bpend, strand=strand, name=access))

            return genes_coordinates

    
    def promoter_coordinates_from_annotations(self, gene_annotation_file, load_exons_introns_info=False, header_lines=1,
                                              promoter_region=None):
        if promoter_region:
            return [g.promoter_c for g in
                    self.load_from_annotation(gene_annotation_file, load_exons_introns_info=False, header_lines=1,
                                             regions=[self.regions[0], promoter_region[0], promoter_region[1],
                                                      self.regions[3]])]
        else:
            return [g.promoter_c for g in
                    self.load_from_annotation(gene_annotation_file, load_exons_introns_info=False, header_lines=1)]

    
    def write_genomic_regions_bed(self, gene_annotation_file, genome_name='', minimal_format=True, promoter_region=None):

        exons = self.exons_from_annotations(gene_annotation_file)
        Coordinate.coordinates_to_bed(exons, genome_name + 'exons.bed', minimal_format=minimal_format)
        del (exons)

        introns = self.introns_from_annotations(gene_annotation_file)
        Coordinate.coordinates_to_bed(introns, genome_name + 'introns.bed', minimal_format=minimal_format)
        del (introns)

        promoters = self.promoter_coordinates_from_annotations(gene_annotation_file, promoter_region=promoter_region)
        Coordinate.coordinates_to_bed(promoters, genome_name + 'promoters.bed', minimal_format=True)
        del (promoters)

    tss = property(tss)
    tes = property(tes)
    distal_c = property(distal_c)
    promoter_c = property(promoter_c)
    intra_c = property(intra_c)
    full_c = property(full_c)

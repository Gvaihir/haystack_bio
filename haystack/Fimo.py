class Fimo:
    def __init__(self, meme_motifs_filename, bg_filename, p_value=1.e-4, temp_directory='./'):
        # be aware that they have changed the command line interface recently!
        self.fimo_command = 'fimo --text --thresh ' + str(
            p_value) + '  --bgfile ' + bg_filename + ' ' + meme_motifs_filename
        self.temp_directory = temp_directory

        with open(meme_motifs_filename) as infile:
            self.motif_id_to_name = dict()
            self.motif_id_to_index = dict()
            self.motif_names = []
            self.motif_name_to_index = dict()
            self.motif_ids = []
            motif_index = 0
            for line in infile:
                try:
                    if 'MOTIF' in line:

                        # in meme format sometime you don't have the name and id but only one string for both
                        # like in jolma 2013...
                        motif_id = line.split()[1]

                        try:
                            motif_name = line.split()[2]
                        except:
                            motif_name = motif_id

                        self.motif_id_to_name[motif_id] = motif_name
                        self.motif_id_to_index[motif_id] = motif_index
                        self.motif_name_to_index[motif_name] = motif_name
                        self.motif_ids.append(motif_id)
                        self.motif_names.append(motif_name)
                        motif_index += 1
                except:
                    print('problem with this line:', line)

    def extract_motifs(self, seq, report_mode='full'):
        if report_mode == 'indexes_set':
            motifs_in_sequence = set()
        elif report_mode == 'fq_array':
            motifs_in_sequence = np.zeros(len(self.motif_names))
        elif report_mode == 'full':
            motifs_in_sequence = list()
        elif report_mode == 'fq_and_presence':
            motifs_in_sequence = {'presence': set(), 'fq': np.zeros(len(self.motif_names))}
        else:
            raise Exception('report_mode not recognized')

        with tempfile.NamedTemporaryFile('w+', dir=self.temp_directory, delete=False) as tmp_file:
            tmp_file.write(''.join(['>S\n', seq, '\n']))
            tmp_filename = tmp_file.name
            tmp_file.close()

            fimo_process = subprocess.Popen(self.fimo_command + ' ' + tmp_filename, stdin=None, stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE, shell=True)
            output = fimo_process.communicate()[0]
            fimo_process.wait()

            os.remove(tmp_filename)

            lines = output.split('\n')
            lines = lines[1:]
            for line in lines:
                if line:
                    fields = line.split('\t')
                    motif_id = fields[0]
                    motif_name = self.motif_id_to_name[motif_id]

                    if report_mode == 'full':
                        c_start = int(fields[2])
                        c_end = int(fields[3])
                        strand = fields[4]
                        score = float(fields[5])
                        p_value = float(fields[6])
                        length = len(fields[8])

                        motifs_in_sequence.append(
                            {'id': motif_id, 'name': motif_name, 'start': c_start, 'end': c_end, 'strand': strand,
                             'score': score, 'p_value': p_value, 'length': length})
                    elif report_mode == 'indexes_set':
                        motifs_in_sequence.add(self.motif_name_to_index[motif_name])

                    elif report_mode == 'fq_array':
                        motifs_in_sequence[self.motif_name_to_index[motif_name]] += 1

                    elif report_mode == 'fq_and_presence':
                        motifs_in_sequence['presence'].add(self.motif_name_to_index[motif_name])
                        motifs_in_sequence['fq'][self.motif_name_to_index[motif_name]] += 1

            return motifs_in_sequence if report_mode == 'fq_array' or 'fq_and_presence' else list(motifs_in_sequence)

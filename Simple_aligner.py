for_rev_read_files = ('A2_S1_L001_R1_001.fastq', 'A2_S1_L001_R2_001.fastq')
ref_file = './S1.fasta'

for record in SeqIO.parse(ref_file, "fasta"):
    seq_len = len(record.seq)
    # create string of O's to insert sequence matches
    string_base = 'O' * seq_len
    read_list = []
    insert_list = []
    ref_sequence_parsed = str(record.seq)
    circular_ref_seq = str(ref_sequence_parsed.upper() + ref_sequence_parsed.upper())
    # iterate through both R1 and R2 files
    for read_file in for_rev_read_files:
        for seq_string in SeqIO.parse(read_file, "fastq"):
            score_list = seq_string.letter_annotations["phred_quality"]
            score = sum(score_list) / float(len(score_list))
            if score < 20:
                continue
            read_sec_string = str(seq_string.seq)
            # ensure case is same
            if str(read_sec_string.upper()) in ref_sequence_parsed.upper():
                index_loc_start = int(ref_sequence_parsed.upper().index(read_sec_string.upper()))
                index_loc_end = int(index_loc_start) + int(len(seq_string.seq))
                read_string_base = str(string_base[0:index_loc_start] + read_sec_string + string_base[index_loc_end:])
                read_string_base_list = list(read_string_base)
                read_list.append(read_string_base_list)
                exact_match_count += 1

            elif str(read_sec_string.upper()) in circular_ref_seq:
                # ensure the start location of match is before the repeat of the sequence
                if int(circular_ref_seq.index(read_sec_string.upper())) < seq_len:
                    index_loc_start = int(circular_ref_seq.upper().index(read_sec_string.upper()))
                    index_loc_end = int(index_loc_start) + int(seq_len)
                    start_read_slice = seq_len - index_loc_start
                    end_read_slice = len(read_sec_string) - start_read_slice
                    start = read_sec_string[:start_read_slice]
                    end = read_sec_string[start_read_slice:]
                    read_string_base = end + string_base[end_read_slice:index_loc_start] + start
                    read_string_base_list = list(read_string_base)
                    read_list.append(read_string_base_list)
                    circ_match_count += 0

            else:
                # Allows for 1 SNP per read, 1 deletion
                fuzzy_matches = find_near_matches(read_sec_string, ref_sequence_parsed, max_l_dist=1)
                if fuzzy_matches:
                    if (int(fuzzy_matches[0][1]) - int(fuzzy_matches[0][0])) == len(read_sec_string):
                        fuzzy_start_loc = fuzzy_matches[0][0]
                        fuzzy_end_loc = fuzzy_matches[0][1]
                        read_string_base = str(
                            string_base[0:fuzzy_start_loc] + read_sec_string + string_base[fuzzy_end_loc:])
                        read_string_base_list = list(read_string_base)
                        read_list.append(read_string_base_list)
                        snp_match_count += 1

                    elif (int(fuzzy_matches[0][1]) - int(fuzzy_matches[0][0])) >= len(read_sec_string):
                        mismatch = ([a == b for (a_i, a) in enumerate(ref_sequence_parsed[fuzzy_matches[0][0]:fuzzy_matches[0][1]])
                                     for(b_i, b) in enumerate(read_sec_string.upper()) if a_i == b_i].index(False))
                        fuzzy_start_loc = fuzzy_matches[0][0]
                        fuzzy_end_loc = fuzzy_matches[0][1]
                        read_string_base = str(string_base[0:fuzzy_start_loc] + read_sec_string[0:mismatch] + 'O' +
                                               read_sec_string[(mismatch):fuzzy_end_loc] + string_base[fuzzy_end_loc:])
                        read_string_base_list = list(read_string_base)
                        read_list.append(read_string_base_list)
                        del_match_count += 1

                    elif (int(fuzzy_matches[0][1]) - int(fuzzy_matches[0][0])) <= len(read_sec_string):
                        insert_item_list = []
                        mismatch = (
                        [a == b for (a_i, a) in enumerate(ref_sequence_parsed[fuzzy_matches[0][0]:fuzzy_matches[0][1]]) for
                         (b_i, b) in enumerate(read_sec_string.upper()) if a_i == b_i].index(False))
                        fuzzy_start_loc = fuzzy_matches[0][0]
                        fuzzy_end_loc = fuzzy_matches[0][1]
                        without_insert_read = read_sec_string[:mismatch] + read_sec_string[mismatch + 1:]
                        loc_match_key = fuzzy_start_loc + mismatch
                        read_string_base = str(
                            string_base[0:fuzzy_start_loc] + without_insert_read + string_base[fuzzy_end_loc:])
                        read_string_base_list = list(read_string_base)
                        read_list.append(read_string_base_list)
                        insert_list.append((loc_match_key, read_sec_string[mismatch]))
                        ins_match_count += 1

                    else:
                        pass

    read_array = np.array(read_list)
    read_array_trans = np.transpose(read_array)
    made_sequence = ''
    coverage = []
    for array in read_array_trans:
        # find top two highest counts to account for the 'O's
        array_counts = (Counter(array).most_common(2))
        # needed for fuzzy string matching
        non_zero_bases = [i for i in array_counts if i[0] != 'O']
        try:
            if non_zero_bases[0][1] > 1:
                made_sequence += non_zero_bases[0][0]
                # find coverage at each position
                coverage.append(non_zero_bases[0][1])
        except IndexError:
            pass
    insert_overall_count = Counter(insert_list)
    counter_dict = dict(insert_overall_count)
    desc_counter_dict = dict(sorted(counter_dict.items(), key=lambda pair: pair[0], reverse=True))
    inserts_sequence = ''
    for key, value in counter_dict.items():
        if value >= 2:
            inserts_sequence = made_sequence[:key[0]] + key[1] + made_sequence[key[0]:]
            coverage.insert(key[0] - 1, value)
    print(inserts_sequence)
    print(coverage)

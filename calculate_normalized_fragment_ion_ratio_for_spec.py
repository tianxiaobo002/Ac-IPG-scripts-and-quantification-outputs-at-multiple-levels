def normalize_ratio_for_spec (work_path,fragment_ions_matched_intens):

    ###Based on the theroretical ratio, calculating the normalize numeber
    desired_ratio = input ('please indicattes the mixed sample ratio')
    normalize_number = 0
    for each_ratio_number in desired_ratio:
        normalize_number = normalize_number + int(each_ratio_number)

    ###For 3-plex labeled fragment ion, if any channel lost the signal, the one will be discarded. Here, defining the minimum number of macthed 4-plex labeled ion.
    ###If the spectrum has less matched 4-plex labeled, the spectrum will be marked as 'not counted'
    minimum_number_of_matched_ion_for_calculate_ratio = int (input ('please indicattes minimum number of matched ion for single spectrum'))
    file_path = work_path
    gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time = {}
    
    for each_retention_time in fragment_ions_matched_intens.keys():
        all_ion_name = []
        gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time] = []
        ion_name_match_three_plex_intensity = {}
        collect_three_intensities = []
        sum_three_intensities = []
        
        for each_frag_ion in fragment_ions_matched_intens [each_retention_time]:

            ### To exclude the y1 ions,the '-H2O' and '-NH3' ions. If the original 3-plex labeled fragment ions are needed, just skip these two row.
            if each_frag_ion [0][0:3] == 'y1_' or '-H2O' in each_frag_ion [0] or '-NH3' in each_frag_ion [0]:
                ###print (each_retention_time,each_frag_ion,'discarded')
                continue
            if each_frag_ion [0][0:3] == 'y1_' or '-H2O' in each_frag_ion [0] or '-NH3' in each_frag_ion [0]:
                print (each_retention_time,each_frag_ion,'Check point')
                
            if each_frag_ion [0][-5:] == '13C+0':

                all_ion_name.append (each_frag_ion [0][0:-6])
                ion_name_match_three_plex_intensity [each_frag_ion [0][0:-6]] =  []
                
                collect_three_intensities.append (float (each_frag_ion [3]))
                ion_name_match_three_plex_intensity[each_frag_ion [0][0:-6]].append (each_frag_ion)
                
            elif each_frag_ion [0][-5:] == '13C+1':
                collect_three_intensities.append (float (each_frag_ion [3]))
                ion_name_match_three_plex_intensity[each_frag_ion [0][0:-6]].append (each_frag_ion)
            elif each_frag_ion [0][-5:] == '13C+2':
                collect_three_intensities.append (float (each_frag_ion [3]))
                ion_name_match_three_plex_intensity[each_frag_ion [0][0:-6]].append (each_frag_ion)

            ### 3-plex labeling was used in the Ac-IPG tag 
                ###elif each_frag_ion [0][-5:] == '13C+3':
                
                ###collect_four_intensities.append (float (each_frag_ion [3]))
                ###ion_name_match_four_plex_intensity[each_frag_ion [0][0:-6]].append (each_frag_ion)
                
                ###'''The b ions and y ions has the complementary ratio. For keep the final ratio of y ions consist ratio of b ions, reversing the intensity order.''' Attention!
                ### Only y_ions were multiple labelled, so we don't need to reverse the ratio
                ###if each_frag_ion [0][0:1] == 'y':
                    ###collect_four_intensities.reverse ()

                ion_name_match_three_plex_intensity [each_frag_ion [0][0:-6]].insert (0,collect_three_intensities)
               
                ###The '13C+3' is used as the end of every fragment ion. At the end of every fragment ion, the total intensity of that single fragment ion will be caculated.
                three_channels_total_intensity = 0
                for each_intensity in collect_three_intensities:
                    three_channels_total_intensity = three_channels_total_intensity + float (each_intensity)
                sum_three_intensities.append (three_channels_total_intensity)
                ion_name_match_three_plex_intensity [each_frag_ion [0][0:-6]].insert (1,sum_three_intensities)

                ###The calculated sum_four_intensities will be used to produce normalize number.
                ###For example, for 1515 theroretical ratio the sum_four_intensities/(1+5+1+5) is the normalize nuber which represent the intensity of '1'.
                
                normalized_three_intensities_ratio = []
                uncounted_channel_signal_loss = ['uncounted_due_to_channel_signal_loss',0,0]
                if 0 in collect_three_intensities:
                    ion_name_match_three_plex_intensity [each_frag_ion [0][0:-6]].insert (0,uncounted_channel_signal_loss)
                else:
                    for each_intensity in collect_three_intensities:
                        normalized_three_intensities_ratio.append (round (float (each_intensity)/(three_channels_total_intensity/normalize_number),2))
                    ion_name_match_three_plex_intensity [each_frag_ion [0][0:-6]].insert (0,normalized_three_intensities_ratio)
                
                collect_three_intensities = []
                sum_three_intensities = []
                                                 
            elif each_frag_ion [0][0:4] == 'spec':

                gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].append (ion_name_match_three_plex_intensity)
                
                gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].append (each_frag_ion)

        ###To determine which ions should be used for calculating ratio on the sacn level.If the precursor ion selection was bothered by 13C contribution, only the b1,2 and y1,2 will be used.

        number_of_match_ion = 0
        counted_matched_ions = []
        sum_intensity_spec_ratio = [0,0,0]
        normalized_spec_ratio = [0,0,0]
        
        if fragment_ions_matched_intens [each_retention_time][-1:][0][1][-1:][0][1] == 1:
        
            for each_kind_of_fragment_ion in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0].keys():
                ###if len (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0][each_kind_of_fragment_ion][0]) == 3:
                if 'uncounted_due_to_channel_signal_loss' not in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0][each_kind_of_fragment_ion][0]:
                    number_of_match_ion = number_of_match_ion + 1
                    counted_matched_ions.append (each_kind_of_fragment_ion)
                    sum_intensity_spec_ratio [0] = sum_intensity_spec_ratio [0] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][0])
                    sum_intensity_spec_ratio [1] = sum_intensity_spec_ratio [1] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][1])
                    sum_intensity_spec_ratio [2] = sum_intensity_spec_ratio [2] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][2])
                    ###sum_intensity_spec_ratio [3] = sum_intensity_spec_ratio [3] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][3])
        elif fragment_ions_matched_intens [each_retention_time][-1:][0][1][-1:][0][1] == 0:
            for each_kind_of_fragment_ion in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0].keys():
                ####if len (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0][each_kind_of_fragment_ion][0]) == 4 and int (each_kind_of_fragment_ion.split ('_')[0][1:]) < 3:
                if 'uncounted_due_to_channel_signal_loss' not in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0][each_kind_of_fragment_ion][0] and int (each_kind_of_fragment_ion.split ('_')[0][1:]) < 3:
                    number_of_match_ion = number_of_match_ion + 1
                    counted_matched_ions.append (each_kind_of_fragment_ion)

                    sum_intensity_spec_ratio [0] = sum_intensity_spec_ratio [0] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][0])
                    sum_intensity_spec_ratio [1] = sum_intensity_spec_ratio [1] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][1])
                    sum_intensity_spec_ratio [2] = sum_intensity_spec_ratio [2] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][2])
                    ###sum_intensity_spec_ratio [3] = sum_intensity_spec_ratio [3] + float (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time][0] [each_kind_of_fragment_ion][1][3])
            
        sum_intensity_spec_ratio [0] = round (sum_intensity_spec_ratio [0],4)
        sum_intensity_spec_ratio [1] = round (sum_intensity_spec_ratio [1],4)
        sum_intensity_spec_ratio [2] = round (sum_intensity_spec_ratio [2],4)
        ###sum_intensity_spec_ratio [3] = round (sum_intensity_spec_ratio [3],4)
        
        if number_of_match_ion >= minimum_number_of_matched_ion_for_calculate_ratio:
                
            spec_all_channel_total_intensity = 0
            for spec_each_channel_total_intensity in sum_intensity_spec_ratio:
                spec_all_channel_total_intensity = spec_all_channel_total_intensity + float (spec_each_channel_total_intensity)
            normalized_spec_ratio [0] = round (sum_intensity_spec_ratio [0]/(spec_all_channel_total_intensity/normalize_number),2)
            normalized_spec_ratio [1] = round (sum_intensity_spec_ratio [1]/(spec_all_channel_total_intensity/normalize_number),2)
            normalized_spec_ratio [2] = round (sum_intensity_spec_ratio [2]/(spec_all_channel_total_intensity/normalize_number),2)
            ###normalized_spec_ratio [3] = round (sum_intensity_spec_ratio [3]/(spec_all_channel_total_intensity/normalize_number),2)
            gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].insert (0,normalized_spec_ratio)
            gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].insert (1,sum_intensity_spec_ratio)
            gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].insert (2,counted_matched_ions)

        else:
            uncounted_too_less_matched_ion = ['uncounted_due_to_too_less_matched_ion',0,0]
            no_ion_was_used = ['no_ion_was_used']
            gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].insert (0,uncounted_too_less_matched_ion)
            gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].insert (1,sum_intensity_spec_ratio)
            gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time [each_retention_time].insert (2,counted_matched_ions)
                    
    output_file = open (file_path + 'Dict_all_fragment_ion_add_matched_intens_normalized_ratio.txt', 'w')
    output_file.write(str(gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time))
    output_file.close()
    
    ############################################
    import csv

    write_by_rows = [['retention_time','ion_name','normalized_ratio_frag_ion_level','every_channel_intensity__frag_ion_level','total_intensity_frag_ion_level','fragment_ion_mass','peptide', 'protein']]
    
    for each_retention_time in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time.keys():
        ###print (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time])
        
        new_line = [str (each_retention_time)]
        for each_fragment_ion_frag_ion_level in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3].keys():
            new_line.append (str (each_fragment_ion_frag_ion_level))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][0]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][1]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][2][0]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][3][2]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][6][1]))
            write_by_rows.append (new_line)

            new_line = [str (each_retention_time)]
       
    
    out_putfile = open (file_path + 'Fragment_ions_normalized_ratio_for_quik_check.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing Fragment_ions_normalized_ratio_for_quik_check csv file')
    out_put_csv_file_writing.writerows(write_by_rows)
    output_file.close()
    

    ############################################
    import csv

    ###write_by_rows = [['retention_time','ion_name','nor_ratio_frag_ion_level_plus0','nor_ratio_frag_ion_level_plus1','nor_ratio_frag_ion_level_plus2','nor_ratio_frag_ion_level_plus3','every_channel_intensity__frag_ion_level','total_intensity_frag_ion_level','fragment_ion_mass','peptide', 'protein']]
    write_by_rows = [['retention_time','ion_name','nor_ratio_frag_ion_level_plus0','nor_ratio_frag_ion_level_plus1','nor_ratio_frag_ion_level_plus2','every_channel_intensity__frag_ion_level','total_intensity_frag_ion_level','fragment_ion_mass','peptide', 'protein']]

    for each_retention_time in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time.keys():
        
        new_line = [str (each_retention_time)]
        for each_fragment_ion_frag_ion_level in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3].keys():
            
            new_line.append (str (each_fragment_ion_frag_ion_level))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][0][0]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][0][1]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][0][2]))
            ###new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][0][3]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][1]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][2][0]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][3][each_fragment_ion_frag_ion_level][3][2]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]))
            new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][6][1]))
            write_by_rows.append (new_line)

            new_line = [str (each_retention_time)]
       
    
    out_putfile = open (file_path + 'Fragment_ions_normalized_ratio_for_figure.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing Fragment_ions_normalized_ratio_for_figure csv file')
    out_put_csv_file_writing.writerows(write_by_rows)
    output_file.close()
    
    ############################################
    import csv

    write_by_rows = [['retention_time','normalized_ratio','total_intensity_on_spec_level','counted_frag_ions','peptide', 'protein']]
    how_many_retention_time_in_retention_time_csv = 0
    for each_retention_time in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time.keys():
        
        new_line = [str (each_retention_time)]
        new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][0]))
        new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1]))
        new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][2]))
        new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]))
        new_line.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][6][1]))
        write_by_rows.append (new_line)
           
    out_putfile = open (file_path + 'Every_spec_normalized_ratio_for_quik_check.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing csv file')
    out_put_csv_file_writing.writerows(write_by_rows)
    output_file.close()

    ##############################################
    write_by_rows_for_figure = [['retention_time','nor_ratio_Gly_plus0_Pro','nor_ratio_Gly_plus1_Pro','nor_ratio_Gly_plus2_Pro','intensity_Gly_plus0_Pro','intensity_Gly_plus1_Pro','intensity_Gly_plus2_Pro','total_intensity','counted_frag_ions','peptide', 'protein']]
    
    for each_retention_time in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time.keys():
    
        new_line_figure = [str (each_retention_time)]
        
        if 0 not in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][0]:
            
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][0][0]))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][0][1]))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][0][2]))
            ###new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][0][3]))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1][0]))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1][1]))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1][2]))
            ###new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1][3]))
            new_line_figure.append (str (sum (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1])))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][2]))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]))
            new_line_figure.append (str (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][6][1]))
        write_by_rows_for_figure.append (new_line_figure)
           
    out_putfile = open (file_path + 'Every_spec_normalized_ratio_for_figure_on_spec_level.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing csv file')
    out_put_csv_file_writing.writerows(write_by_rows_for_figure)
    output_file.close()
    
    ##############################################
    top_n_spectrum_for_peptide_ratio = int(input ('please indicates top N spectrum used for peptide ratio'))
    all_intensity_collection_peptide_level = {}
    all_peptides_collection = []
    
    for each_retention_time in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time.keys():

        ### the spectrum which less than 3 matched y-ion needs to be excluded

        if 'uncounted_due_to_too_less_matched_ion'in gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][0]:
            continue
        
        ### the next step is to creat the intensities collection grouped by peptide sequence, the first item is the all 3-plex intensities collection which will be used to select top n
        ### the second item is the records of the specific spectrum, the third is the peptide information

        if gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1] not in all_peptides_collection:
            all_peptides_collection.append (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1])
            all_intensity_collection_peptide_level [gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]] = []
            all_intensity_collection_peptide_level [gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]].append([gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1]])
            all_intensity_collection_peptide_level [gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]].append ({each_retention_time:gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1]})
            all_intensity_collection_peptide_level [gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]].append (gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][-1:])
        else:
            all_intensity_collection_peptide_level [gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]][0].append(gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1])
            all_intensity_collection_peptide_level [gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][4][1][1][1]][1][each_retention_time] = gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time[each_retention_time][1]

    ### the next is to find out the top n spectra 
    import heapq
    for each_peptide_seq in all_intensity_collection_peptide_level.keys():
        all_intensity_collection_peptide_level [each_peptide_seq][0] =heapq.nlargest(top_n_spectrum_for_peptide_ratio,all_intensity_collection_peptide_level [each_peptide_seq][0])
    
    
    for each_peptide_seq in all_intensity_collection_peptide_level.keys():
        
        top_n_spectra_sum_intensity = [0,0,0]
        for each_top_selected_intensity in all_intensity_collection_peptide_level [each_peptide_seq][0]:
            
            top_n_spectra_sum_intensity [0]= round((top_n_spectra_sum_intensity [0] + each_top_selected_intensity [0]),4)
            top_n_spectra_sum_intensity [1]= round((top_n_spectra_sum_intensity [1] + each_top_selected_intensity [1]),4)
            top_n_spectra_sum_intensity [2]= round((top_n_spectra_sum_intensity [2] + each_top_selected_intensity [2]),4)
            ###top_n_spectra_sum_intensity [3]= round((top_n_spectra_sum_intensity [3] + each_top_selected_intensity [3]),4)

        if sum(top_n_spectra_sum_intensity) != 0:
            normalized_ratio_on_peptide_level = [round ((top_n_spectra_sum_intensity [0]/(sum(top_n_spectra_sum_intensity)/normalize_number)),2),round ((top_n_spectra_sum_intensity [1]/(sum(top_n_spectra_sum_intensity)/normalize_number)),2),round ((top_n_spectra_sum_intensity [2]/(sum(top_n_spectra_sum_intensity)/normalize_number)),2)]
        else:
            normalized_ratio_on_peptide_level = [0,0,0]
               
        all_intensity_collection_peptide_level [each_peptide_seq].insert (0,top_n_spectra_sum_intensity)
        all_intensity_collection_peptide_level [each_peptide_seq].insert (0,normalized_ratio_on_peptide_level)

    write_by_rows = [['peptide','normalized_ratio_peptide_level','total_intensity_every_channel_peptide_level','total_intensity_all_channels_peptide_level']]

    for each_peptide_seq in all_intensity_collection_peptide_level.keys():
        
        if 0 in all_intensity_collection_peptide_level[each_peptide_seq][0]:
            continue
    
        new_line = [str (each_peptide_seq)]
        new_line.append (str (all_intensity_collection_peptide_level[each_peptide_seq][0]))
        new_line.append (str (all_intensity_collection_peptide_level[each_peptide_seq][1]))
        new_line.append (str (sum (all_intensity_collection_peptide_level[each_peptide_seq][1])))
        
        
        write_by_rows.append (new_line)
       
    
    out_putfile = open (file_path + 'Every_peptide_normalized_ratio_for_quik_check.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing csv file')
    out_put_csv_file_writing.writerows(write_by_rows)
    output_file.close()


    write_by_rows = [['peptide','normalized_ratio_peptide_level_Gly_plus0_Pro','normalized_ratio_peptide_level_Gly_plus1_Pro','normalized_ratio_peptide_level_Gly_plus2_Pro','total_intensity_every_channel_peptide_level','total_intensity_all_channels_peptide_level']]

    for each_peptide_seq in all_intensity_collection_peptide_level.keys():
        if 0 in all_intensity_collection_peptide_level[each_peptide_seq][0]:
            continue
    
        new_line = [str (each_peptide_seq)]
        new_line.append (str (all_intensity_collection_peptide_level[each_peptide_seq][0][0]))
        new_line.append (str (all_intensity_collection_peptide_level[each_peptide_seq][0][1]))
        new_line.append (str (all_intensity_collection_peptide_level[each_peptide_seq][0][2]))
        ###new_line.append (str (all_intensity_collection_peptide_level[each_peptide_seq][0][3]))
        new_line.append (str (all_intensity_collection_peptide_level[each_peptide_seq][1]))
        new_line.append (str (sum (all_intensity_collection_peptide_level[each_peptide_seq][1])))
        
        
        write_by_rows.append (new_line)
       
    
    out_putfile = open (file_path + 'Every_peptide_normalized_ratio_for_figure.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing csv file')
    out_put_csv_file_writing.writerows(write_by_rows)
    output_file.close()
    #####################################################

    top_n_peptides_for_protein_ratio = int (input ('please indicates top N peptides used for protein ratio'))
    all_intensity_collection_protein_level = {}
    all_protein_accession_collection = []
    for each_peptide_sequence in all_intensity_collection_peptide_level.keys():
        
        temp_related_proteins = all_intensity_collection_peptide_level[each_peptide_sequence][4][0][1][6][1].split(':')
        for each_temp_related_protein in temp_related_proteins:
            if each_temp_related_protein not in all_protein_accession_collection:
                all_protein_accession_collection.append (each_temp_related_protein)
                all_intensity_collection_protein_level [each_temp_related_protein] = []
                all_intensity_collection_protein_level [each_temp_related_protein].append([all_intensity_collection_peptide_level[each_peptide_sequence][1]])
                all_intensity_collection_protein_level [each_temp_related_protein].append ([{each_peptide_sequence:all_intensity_collection_peptide_level[each_peptide_sequence]}])
            else:
                all_intensity_collection_protein_level [each_temp_related_protein][0].append(all_intensity_collection_peptide_level[each_peptide_sequence][1])
                all_intensity_collection_protein_level [each_temp_related_protein][1].append ({each_peptide_sequence:all_intensity_collection_peptide_level[each_peptide_sequence]})

    ### the next is to find out the top n spectra 
    import heapq
    for each_selected_protein in all_intensity_collection_protein_level.keys():
        all_intensity_collection_protein_level [each_selected_protein][0] =heapq.nlargest(top_n_peptides_for_protein_ratio,all_intensity_collection_protein_level [each_selected_protein][0])
    
    
    for each_selected_protein in all_intensity_collection_protein_level.keys():
        
        top_n_peptide_sum_intensity = [0,0,0]
        for each_top_selected_peptide_intensity in all_intensity_collection_protein_level [each_selected_protein][0]:
            
            top_n_peptide_sum_intensity [0]= round((top_n_peptide_sum_intensity [0] + each_top_selected_peptide_intensity [0]),4)
            top_n_peptide_sum_intensity [1]= round((top_n_peptide_sum_intensity [1] + each_top_selected_peptide_intensity [1]),4)
            top_n_peptide_sum_intensity [2]= round((top_n_peptide_sum_intensity [2] + each_top_selected_peptide_intensity [2]),4)
            ###top_n_peptide_sum_intensity [3]= round((top_n_peptide_sum_intensity [3] + each_top_selected_peptide_intensity [3]),4)
            
        if sum(top_n_peptide_sum_intensity) != 0:
            normalized_ratio_on_protein_level = [round ((top_n_peptide_sum_intensity [0]/(sum(top_n_peptide_sum_intensity)/normalize_number)),2),round ((top_n_peptide_sum_intensity [1]/(sum(top_n_peptide_sum_intensity)/normalize_number)),2),round ((top_n_peptide_sum_intensity [2]/(sum(top_n_peptide_sum_intensity)/normalize_number)),2)]
        else:
            normalized_ratio_on_protein_level = [0,0,0]
               
        all_intensity_collection_protein_level [each_selected_protein].insert (0,top_n_peptide_sum_intensity)
        all_intensity_collection_protein_level [each_selected_protein].insert (0,normalized_ratio_on_protein_level)
    #######################################################
    write_by_rows = [['protein','normalized_ratio_protein_level','total_intensity_every_channel_protein_level','total_intensity_all_channels_protein_level']]

    for each_protein_accession in all_intensity_collection_protein_level.keys():
        if 0 in all_intensity_collection_protein_level[each_protein_accession][0]:
            continue
    
        new_line = [str (each_protein_accession)]
        new_line.append (str (all_intensity_collection_protein_level[each_protein_accession][0]))
        new_line.append (str (all_intensity_collection_protein_level[each_protein_accession][1]))
        new_line.append (str (sum (all_intensity_collection_protein_level[each_protein_accession][1])))
        
        write_by_rows.append (new_line)
        
    out_putfile = open (file_path + 'Every_protein_normalized_ratio_for_quik_check.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing Every_protein_normalized_ratio_for_quik_check csv file')
    out_put_csv_file_writing.writerows(write_by_rows)
    output_file.close()
    #######################################################
    
    write_by_rows = [['protein','normalized_ratio_protein_level_Gly_plus0_Pro','normalized_ratio_protein_level_Gly_plus1_Pro','normalized_ratio_protein_level_Gly_plus2_Pro','total_intensity_every_channel_protein_level','total_intensity_all_channels_protein_level']]

    for each_protein_accession in all_intensity_collection_protein_level.keys():

        if 0 in all_intensity_collection_protein_level[each_protein_accession][0]:
            continue
    
        new_line = [str (each_protein_accession)]
        new_line.append (str (all_intensity_collection_protein_level[each_protein_accession][0][0]))
        new_line.append (str (all_intensity_collection_protein_level[each_protein_accession][0][1]))
        new_line.append (str (all_intensity_collection_protein_level[each_protein_accession][0][2]))
        ###new_line.append (str (all_intensity_collection_protein_level[each_protein_accession][0][3]))
        new_line.append (str (all_intensity_collection_protein_level[each_protein_accession][1]))
        new_line.append (str (sum (all_intensity_collection_protein_level[each_protein_accession][1])))
        
        write_by_rows.append (new_line)
        
    out_putfile = open (file_path + 'Every_protein_normalized_ratio_for_figure.csv', 'w', newline = '')
    #creating a output file. the newline=''is to avoid the blank row in the produced csv 
    out_put_csv_file_writing = csv.writer(out_putfile)
    print ('preparing Every_protein_normalized_ratio_for_figure csv file')
    out_put_csv_file_writing.writerows(write_by_rows)
    output_file.close()
    #######################################################
    
    
    
    return gather_ion_name_intense_ratio_sum_intensities_for_all_retention_time
    


'''test_dict = {6302: [['b1_13C+0', 'S', 349.0705, 1902503.625], ['b1_13C+1', 'S', 350.0739, 1562874.25], ['b1_13C+2', 'S', 351.0773, 753759.8125], ['b1_13C+3', 'S', 352.0807, 1368160.25], ['b2_13C+0', 'SE', 478.1131, 923492.0625], ['b2_13C+1', 'SE', 479.1165, 1314714.25], ['b2_13C+2', 'SE', 480.1199, 995121.625], ['b2_13C+3', 'SE', 481.1233, 845429.625], ['b3_13C+0', 'SEI', 591.1945, '0'], ['b3_13C+1', 'SEI', 592.1979, '0'], ['b3_13C+2', 'SEI', 593.2013, '0'], ['b3_13C+3', 'SEI', 594.2047, '0'], ['b4_13C+0', 'SEIA', 662.2316, '0'], ['b4_13C+1', 'SEIA', 663.235, '0'], ['b4_13C+2', 'SEIA', 664.2384, '0'], ['b4_13C+3', 'SEIA', 665.2418, '0'], ['b5_13C+0', 'SEIAH', 799.2905, '0'], ['b5_13C+1', 'SEIAH', 800.2939, '0'], ['b5_13C+2', 'SEIAH', 801.2973, '0'], ['b5_13C+3', 'SEIAH', 802.3007, '0'], ['b6_13C+0', 'SEIAHR', 955.3916, 240325.2031], ['b6_13C+1', 'SEIAHR', 956.395, 265952.5312], ['b6_13C+2', 'SEIAHR', 957.3984, '0'], ['b6_13C+3', 'SEIAHR', 958.4018, 308347.6562], ['b7_13C+0', 'SEIAHRF', 1102.46, '0'], ['b7_13C+1', 'SEIAHRF', 1103.4634, '0'], ['b7_13C+2', 'SEIAHRF', 1104.4668, '0'], ['b7_13C+3', 'SEIAHRF', 1105.4702, '0'], ['b1_-H2O_13C+0', 'S', 331.0599, 1554115.75], ['b1_-H2O_13C+1', 'S', 332.0633, 1185224.5], ['b1_-H2O_13C+2', 'S', 333.0667, 1092276.0], ['b1_-H2O_13C+3', 'S', 334.0701, 1082826.0], ['b1_-NH3_13C+0', 'S', 332.044, 1185224.5], ['b1_-NH3_13C+1', 'S', 333.0474, 1092276.0], ['b1_-NH3_13C+2', 'S', 334.0508, 1082826.0], ['b1_-NH3_13C+3', 'S', 335.0542, '0'], ['b2_-H2O_13C+0', 'SE', 460.1025, 1563130.625], ['b2_-H2O_13C+1', 'SE', 461.1059, 2031038.375], ['b2_-H2O_13C+2', 'SE', 462.1093, 1456113.5], ['b2_-H2O_13C+3', 'SE', 463.1127, 1470833.125], ['b2_-NH3_13C+0', 'SE', 461.0866, 2031038.375], ['b2_-NH3_13C+1', 'SE', 462.09, 1456113.5], ['b2_-NH3_13C+2', 'SE', 463.0934, 1470833.125], ['b2_-NH3_13C+3', 'SE', 464.0968, '0'], ['b3_-H2O_13C+0', 'SEI', 573.1839, '0'], ['b3_-H2O_13C+1', 'SEI', 574.1873, '0'], ['b3_-H2O_13C+2', 'SEI', 575.1907, '0'], ['b3_-H2O_13C+3', 'SEI', 576.1941, '0'], ['b3_-NH3_13C+0', 'SEI', 574.168, '0'], ['b3_-NH3_13C+1', 'SEI', 575.1714, '0'], ['b3_-NH3_13C+2', 'SEI', 576.1748, '0'], ['b3_-NH3_13C+3', 'SEI', 577.1782, '0'], ['b4_-H2O_13C+0', 'SEIA', 644.221, '0'], ['b4_-H2O_13C+1', 'SEIA', 645.2244, '0'], ['b4_-H2O_13C+2', 'SEIA', 646.2278, '0'], ['b4_-H2O_13C+3', 'SEIA', 647.2312, '0'], ['b4_-NH3_13C+0', 'SEIA', 645.2051, '0'], ['b4_-NH3_13C+1', 'SEIA', 646.2085, '0'], ['b4_-NH3_13C+2', 'SEIA', 647.2119, '0'], ['b4_-NH3_13C+3', 'SEIA', 648.2153, '0'], ['b5_-H2O_13C+0', 'SEIAH', 781.2799, '0'], ['b5_-H2O_13C+1', 'SEIAH', 782.2833, '0'], ['b5_-H2O_13C+2', 'SEIAH', 783.2867, '0'], ['b5_-H2O_13C+3', 'SEIAH', 784.2901, '0'], ['b5_-NH3_13C+0', 'SEIAH', 782.264, '0'], ['b5_-NH3_13C+1', 'SEIAH', 783.2674, '0'], ['b5_-NH3_13C+2', 'SEIAH', 784.2708, '0'], ['b5_-NH3_13C+3', 'SEIAH', 785.2742, '0'], ['b6_-H2O_13C+0', 'SEIAHR', 937.381, '0'], ['b6_-H2O_13C+1', 'SEIAHR', 938.3844, '0'], ['b6_-H2O_13C+2', 'SEIAHR', 939.3878, '0'], ['b6_-H2O_13C+3', 'SEIAHR', 940.3912, '0'], ['b6_-NH3_13C+0', 'SEIAHR', 938.3651, '0'], ['b6_-NH3_13C+1', 'SEIAHR', 939.3685, '0'], ['b6_-NH3_13C+2', 'SEIAHR', 940.3719, '0'], ['b6_-NH3_13C+3', 'SEIAHR', 941.3753, '0'], ['b7_-H2O_13C+0', 'SEIAHRF', 1084.4494, '0'], ['b7_-H2O_13C+1', 'SEIAHRF', 1085.4528, '0'], ['b7_-H2O_13C+2', 'SEIAHRF', 1086.4562, '0'], ['b7_-H2O_13C+3', 'SEIAHRF', 1087.4596, '0'], ['b7_-NH3_13C+0', 'SEIAHRF', 1085.4335, '0'], ['b7_-NH3_13C+1', 'SEIAHRF', 1086.4369, '0'], ['b7_-NH3_13C+2', 'SEIAHRF', 1087.4403, '0'], ['b7_-NH3_13C+3', 'SEIAHRF', 1088.4437, '0'], ['y1_13C+0', 'K', 260.161, 589575.25], ['y1_13C+1', 'K', 261.1644, 304855.3125], ['y1_13C+2', 'K', 262.1678, 558615.875], ['y1_13C+3', 'K', 263.1712, 303316.7188], ['y2_13C+0', 'KF', 407.2294, 624678.125], ['y2_13C+1', 'KF', 408.2328, 349911.875], ['y2_13C+2', 'KF', 409.2362, 677629.625], ['y2_13C+3', 'KF', 410.2396, 849589.8125], ['y3_13C+0', 'KFR', 563.3305, 1533017.5], ['y3_13C+1', 'KFR', 564.3339, 1507670.75], ['y3_13C+2', 'KFR', 565.3373, 1972881.875], ['y3_13C+3', 'KFR', 566.3407, 1339827.5], ['y4_13C+0', 'KFRH', 700.3894, 3417394.5], ['y4_13C+1', 'KFRH', 701.3928, 2721766.25], ['y4_13C+2', 'KFRH', 702.3962, 3689109.0], ['y4_13C+3', 'KFRH', 703.3996, 3383925.5], ['y5_13C+0', 'KFRHA', 771.4265, 5251915.0], ['y5_13C+1', 'KFRHA', 772.4299, 5554512.5], ['y5_13C+2', 'KFRHA', 773.4333, 6313523.0], ['y5_13C+3', 'KFRHA', 774.4367, 6475699.5], ['y6_13C+0', 'KFRHAI', 884.5079, 1271257.375], ['y6_13C+1', 'KFRHAI', 885.5113, 1481627.875], ['y6_13C+2', 'KFRHAI', 886.5147, 2040247.875], ['y6_13C+3', 'KFRHAI', 887.5181, 2084363.25], ['y7_13C+0', 'KFRHAIE', 1013.5505, 949708.1875], ['y7_13C+1', 'KFRHAIE', 1014.5539, 946727.9375], ['y7_13C+2', 'KFRHAIE', 1015.5573, 835570.9375], ['y7_13C+3', 'KFRHAIE', 1016.5607, 1060332.625], ['y1_-H2O_13C+0', 'K', 242.1504, '0'], ['y1_-H2O_13C+1', 'K', 243.1538, 274019.5938], ['y1_-H2O_13C+2', 'K', 244.1572, '0'], ['y1_-H2O_13C+3', 'K', 245.1606, '0'], ['y1_-NH3_13C+0', 'K', 243.1345, 274019.5938], ['y1_-NH3_13C+1', 'K', 244.1379, '0'], ['y1_-NH3_13C+2', 'K', 245.1413, '0'], ['y1_-NH3_13C+3', 'K', 246.1447, '0'], ['y2_-H2O_13C+0', 'KF', 389.2188, '0'], ['y2_-H2O_13C+1', 'KF', 390.2222, '0'], ['y2_-H2O_13C+2', 'KF', 391.2256, '0'], ['y2_-H2O_13C+3', 'KF', 392.229, '0'], ['y2_-NH3_13C+0', 'KF', 390.2029, '0'], ['y2_-NH3_13C+1', 'KF', 391.2063, '0'], ['y2_-NH3_13C+2', 'KF', 392.2097, '0'], ['y2_-NH3_13C+3', 'KF', 393.2131, '0'], ['y3_-H2O_13C+0', 'KFR', 545.3199, '0'], ['y3_-H2O_13C+1', 'KFR', 546.3233, '0'], ['y3_-H2O_13C+2', 'KFR', 547.3267, '0'], ['y3_-H2O_13C+3', 'KFR', 548.3301, '0'], ['y3_-NH3_13C+0', 'KFR', 546.304, '0'], ['y3_-NH3_13C+1', 'KFR', 547.3074, '0'], ['y3_-NH3_13C+2', 'KFR', 548.3108, '0'], ['y3_-NH3_13C+3', 'KFR', 549.3142, '0'], ['y4_-H2O_13C+0', 'KFRH', 682.3788, '0'], ['y4_-H2O_13C+1', 'KFRH', 683.3822, 173730.7969], ['y4_-H2O_13C+2', 'KFRH', 684.3856, '0'], ['y4_-H2O_13C+3', 'KFRH', 685.389, '0'], ['y4_-NH3_13C+0', 'KFRH', 683.3629, 173730.7969], ['y4_-NH3_13C+1', 'KFRH', 684.3663, '0'], ['y4_-NH3_13C+2', 'KFRH', 685.3697, '0'], ['y4_-NH3_13C+3', 'KFRH', 686.3731, '0'], ['y4_2charge_13C+0', 'KFRH', 351.2025, 184259.4219], ['y4_2charge_13C+1', 'KFRH', 352.2059, 897756.875], ['y4_2charge_13C+2', 'KFRH', 353.2093, '0'], ['y4_2charge_13C+3', 'KFRH', 354.2127, '0'], ['y5_-H2O_13C+0', 'KFRHA', 753.4159, '0'], ['y5_-H2O_13C+1', 'KFRHA', 754.4193, '0'], ['y5_-H2O_13C+2', 'KFRHA', 755.4227, 286368.1562], ['y5_-H2O_13C+3', 'KFRHA', 756.4261, '0'], ['y5_-NH3_13C+0', 'KFRHA', 754.4, '0'], ['y5_-NH3_13C+1', 'KFRHA', 755.4034, 286368.1562], ['y5_-NH3_13C+2', 'KFRHA', 756.4068, '0'], ['y5_-NH3_13C+3', 'KFRHA', 757.4102, '0'], ['y5_2charge_13C+0', 'KFRHA', 386.721, 527228.9375], ['y5_2charge_13C+1', 'KFRHA', 387.7244, '0'], ['y5_2charge_13C+2', 'KFRHA', 388.7278, '0'], ['y5_2charge_13C+3', 'KFRHA', 389.7312, '0'], ['y6_-H2O_13C+0', 'KFRHAI', 866.4973, '0'], ['y6_-H2O_13C+1', 'KFRHAI', 867.5007, '0'], ['y6_-H2O_13C+2', 'KFRHAI', 868.5041, '0'], ['y6_-H2O_13C+3', 'KFRHAI', 869.5075, '0'], ['y6_-NH3_13C+0', 'KFRHAI', 867.4814, '0'], ['y6_-NH3_13C+1', 'KFRHAI', 868.4848, '0'], ['y6_-NH3_13C+2', 'KFRHAI', 869.4882, '0'], ['y6_-NH3_13C+3', 'KFRHAI', 870.4916, '0'], ['y6_2charge_13C+0', 'KFRHAI', 443.2618, 311066.1875], ['y6_2charge_13C+1', 'KFRHAI', 444.2652, 626713.6875], ['y6_2charge_13C+2', 'KFRHAI', 445.2686, '0'], ['y6_2charge_13C+3', 'KFRHAI', 446.272, '0'], ['y7_-H2O_13C+0', 'KFRHAIE', 995.5399, '0'], ['y7_-H2O_13C+1', 'KFRHAIE', 996.5433, '0'], ['y7_-H2O_13C+2', 'KFRHAIE', 997.5467, '0'], ['y7_-H2O_13C+3', 'KFRHAIE', 998.5501, 211476.4062], ['y7_-NH3_13C+0', 'KFRHAIE', 996.524, '0'], ['y7_-NH3_13C+1', 'KFRHAIE', 997.5274, '0'], ['y7_-NH3_13C+2', 'KFRHAIE', 998.5308, 211476.4062], ['y7_-NH3_13C+3', 'KFRHAIE', 999.5342, '0'], ['y7_2charge_13C+0', 'KFRHAIE', 507.7831, 353702.7188], ['y7_2charge_13C+1', 'KFRHAIE', 508.7865, 187330.3594], ['y7_2charge_13C+2', 'KFRHAIE', 509.7899, '0'], ['y7_2charge_13C+3', 'KFRHAIE', 510.7933, '0'], ['spectrum_information', [['m_over_z', 682.8141], ['retention_time', 1807.2], ['peptide', 'SEIAHRFK'], ['mass', 1363.6183], ['charge_states', 2], ['logP', 36.75], ['length', 8], ['ppm', -3.4], ['area', '9.4865E8'], ['protein_accession', 'P02769'], ['determining_precursor_selection', 1]]]]}
test_dict_large_peptide ={'1734.27354':[['b1', 'S', 116.0711, '0'], ['b2', 'SE', 245.1137, 11881.1], ['b3', 'SEI', 358.1951, 5838.0], ['b4', 'SEIA', 429.2322, 1736.5], ['b5', 'SEIAH', 566.2911, '0'], ['b6', 'SEIAHR', 722.3922, '0'], ['b7', 'SEIAHRF', 869.4606, '0'], ['b1_-H2O', 'S', 98.0605, '0'], ['b1_-NH3', 'S', 99.0446, '0'], ['b2_-H2O', 'SE', 227.1031, '0'], ['b2_-NH3', 'SE', 228.0872, '0'], ['b3_-H2O', 'SEI', 340.1845, '0'], ['b3_-NH3', 'SEI', 341.1686, '0'], ['b4_-H2O', 'SEIA', 411.2216, '0'], ['b4_-NH3', 'SEIA', 412.2057, '0'], ['b5_-H2O', 'SEIAH', 548.2805, '0'], ['b5_-NH3', 'SEIAH', 549.2646, '0'], ['b6_-H2O', 'SEIAHR', 704.3816, '0'], ['b6_-NH3', 'SEIAHR', 705.3657, '0'], ['b7_-H2O', 'SEIAHRF', 851.45, '0'], ['b7_-NH3', 'SEIAHRF', 852.4341, '0'], ['y1_13C+0', 'K', 301.1875, '0'], ['y1_13C+1', 'K', 302.1909, '0'], ['y1_13C+2', 'K', 303.1943, '0'], ['y2_13C+0', 'FK', 448.2559, '0'], ['y2_13C+1', 'FK', 449.2593, '0'], ['y2_13C+2', 'FK', 450.2627, '0'], ['y3_13C+0', 'RFK', 604.357, '0'], ['y3_13C+1', 'RFK', 605.3604, '0'], ['y3_13C+2', 'RFK', 606.3638, '0'], ['y4_13C+0', 'HRFK', 741.4159, 2923.6], ['y4_13C+1', 'HRFK', 742.4193, 3320.7], ['y4_13C+2', 'HRFK', 743.4227, 2604.6], ['y5_13C+0', 'AHRFK', 812.453, 2678.6], ['y5_13C+1', 'AHRFK', 813.4564, 3183.4], ['y5_13C+2', 'AHRFK', 814.4598, 7245.4], ['y6_13C+0', 'IAHRFK', 925.5344, 3611.7], ['y6_13C+1', 'IAHRFK', 926.5378, 3840.5], ['y6_13C+2', 'IAHRFK', 927.5412, 3007.2], ['y7_13C+0', 'EIAHRFK', 1054.577, 10002.2], ['y7_13C+1', 'EIAHRFK', 1055.5804, 12959.2], ['y7_13C+2', 'EIAHRFK', 1056.5838, 10919.7], ['y8_13C+0', 'SEIAHRFK', 1169.6403, '0'], ['y8_13C+1', 'SEIAHRFK', 1170.6437, '0'], ['y8_13C+2', 'SEIAHRFK', 1171.6471, '0'], ['y1_-H2O_13C+0', 'K', 283.1769, '0'], ['y1_-H2O_13C+1', 'K', 284.1803, '0'], ['y1_-H2O_13C+2', 'K', 285.1837, '0'], ['y1_-NH3_13C+0', 'K', 284.161, '0'], ['y1_-NH3_13C+1', 'K', 285.1644, '0'], ['y1_-NH3_13C+2', 'K', 286.1678, '0'], ['y2_-H2O_13C+0', 'FK', 430.2453, '0'], ['y2_-H2O_13C+1', 'FK', 431.2487, '0'], ['y2_-H2O_13C+2', 'FK', 432.2521, '0'], ['y2_-NH3_13C+0', 'FK', 431.2294, '0'], ['y2_-NH3_13C+1', 'FK', 432.2328, '0'], ['y2_-NH3_13C+2', 'FK', 433.2362, '0'], ['y3_-H2O_13C+0', 'RFK', 586.3464, 14544.6], ['y3_-H2O_13C+1', 'RFK', 587.3498, '0'], ['y3_-H2O_13C+2', 'RFK', 588.3532, '0'], ['y3_-NH3_13C+0', 'RFK', 587.3305, '0'], ['y3_-NH3_13C+1', 'RFK', 588.3339, '0'], ['y3_-NH3_13C+2', 'RFK', 589.3373, '0'], ['y4_-H2O_13C+0', 'HRFK', 723.4053, '0'], ['y4_-H2O_13C+1', 'HRFK', 724.4087, '0'], ['y4_-H2O_13C+2', 'HRFK', 725.4121, '0'], ['y4_-NH3_13C+0', 'HRFK', 724.3894, '0'], ['y4_-NH3_13C+1', 'HRFK', 725.3928, '0'], ['y4_-NH3_13C+2', 'HRFK', 726.3962, '0'], ['y4_2charge_13C+0', 'HRFK', 371.2118, '0'], ['y4_2charge_13C+1', 'HRFK', 371.7135, '0'], ['y4_2charge_13C+2', 'HRFK', 372.2152, '0'], ['y5_-H2O_13C+0', 'AHRFK', 794.4424, '0'], ['y5_-H2O_13C+1', 'AHRFK', 795.4458, '0'], ['y5_-H2O_13C+2', 'AHRFK', 796.4492, '0'], ['y5_-NH3_13C+0', 'AHRFK', 795.4265, '0'], ['y5_-NH3_13C+1', 'AHRFK', 796.4299, '0'], ['y5_-NH3_13C+2', 'AHRFK', 797.4333, '0'], ['y5_2charge_13C+0', 'AHRFK', 406.7304, '0'], ['y5_2charge_13C+1', 'AHRFK', 407.2321, '0'], ['y5_2charge_13C+2', 'AHRFK', 407.7338, '0'], ['y6_-H2O_13C+0', 'IAHRFK', 907.5238, '0'], ['y6_-H2O_13C+1', 'IAHRFK', 908.5272, '0'], ['y6_-H2O_13C+2', 'IAHRFK', 909.5306, '0'], ['y6_-NH3_13C+0', 'IAHRFK', 908.5079, '0'], ['y6_-NH3_13C+1', 'IAHRFK', 909.5113, '0'], ['y6_-NH3_13C+2', 'IAHRFK', 910.5147, '0'], ['y6_2charge_13C+0', 'IAHRFK', 463.2711, '0'], ['y6_2charge_13C+1', 'IAHRFK', 463.7728, '0'], ['y6_2charge_13C+2', 'IAHRFK', 464.2745, '0'], ['y7_-H2O_13C+0', 'EIAHRFK', 1036.5664, '0'], ['y7_-H2O_13C+1', 'EIAHRFK', 1037.5698, '0'], ['y7_-H2O_13C+2', 'EIAHRFK', 1038.5732, '0'], ['y7_-NH3_13C+0', 'EIAHRFK', 1037.5505, '0'], ['y7_-NH3_13C+1', 'EIAHRFK', 1038.5539, '0'], ['y7_-NH3_13C+2', 'EIAHRFK', 1039.5573, '0'], ['y7_2charge_13C+0', 'EIAHRFK', 527.7924, '0'], ['y7_2charge_13C+1', 'EIAHRFK', 528.2941, '0'], ['y7_2charge_13C+2', 'EIAHRFK', 528.7958, '0'], ['y8_-H2O_13C+0', 'SEIAHRFK', 1151.6297, '0'], ['y8_-H2O_13C+1', 'SEIAHRFK', 1152.6331, '0'], ['y8_-H2O_13C+2', 'SEIAHRFK', 1153.6365, '0'], ['y8_-NH3_13C+0', 'SEIAHRFK', 1152.6138, '0'], ['y8_-NH3_13C+1', 'SEIAHRFK', 1153.6172, '0'], ['y8_-NH3_13C+2', 'SEIAHRFK', 1154.6206, '0'], ['y8_2charge_13C+0', 'SEIAHRFK', 585.3241, 13170.9], ['y8_2charge_13C+1', 'SEIAHRFK', 585.8258, 11758.4], ['y8_2charge_13C+2', 'SEIAHRFK', 586.3275, 14544.6], ['spectrum_information', [['m_over_z', 663.874], ['peptide', 'SEIAHRFK'], ['mass', 1325.7366], ['charge_states', '2+'], ['localization_confidence', 'Ac-Ile(+0)-Pro-Gly(+2) (8: Very Confident)'], ['ppm', -2.3783], ['protein_accession', 'P02769'], ['determining_precursor_selection', 1]]]]}
working_path = 'G:\\Desktop\\NAPI-tag\\Scripts for Ac-IPG tag\\'

simpilified_test_dict = normalize_ratio_for_spec (working_path,test_dict_large_peptide)'''

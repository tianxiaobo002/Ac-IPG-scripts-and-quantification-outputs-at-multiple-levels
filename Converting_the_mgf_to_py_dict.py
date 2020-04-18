### the mgf file was produced by the RawConverter
def converting_raw_mgf_to_dict (raw_mgf_file_route_and_file_name):

    file_route = raw_mgf_file_route_and_file_name [0]
    file_name = raw_mgf_file_route_and_file_name [1]

    raw_mgf_file = open (file_route + file_name, 'r')

    raw_mgf_MS2_spectra_list = {}

    scan_number = ''
    peak_list = []
    spectrum_information = []

    for row in raw_mgf_file:
        
        if row == 'BEGIN IONS\n':
            continue
        elif row.split('=')[0] == 'TITLE':
            continue
        elif row == '\n':
            continue
        elif row.split (' ')[0] == 'END':
            ###spectrum_information.append(['retention_seconds', retention_seconds])
            spectrum_information.append(['scan_number', scan_number])

            spectrum_information.append(['pep_precursor', pep_precursor])
            spectrum_information.append(['charge_states', charge_states])
            peak_list.extend([spectrum_information])                            
            raw_mgf_MS2_spectra_list [retention_seconds] = peak_list
            peak_list = []
            scan_number = ''
            spectrum_information = []
            continue
        elif row.split('=')[0] == 'SCANS':
            scan_number =  int (row.split('SCANS=')[1])
            '''print (scan_number)'''
                         
        elif row.split('=')[0] == 'RTINSECONDS':
            retention_seconds = row.split ('=')[1][:-1]
        
        elif row.split('=')[0] == 'PEPMASS':
            pep_precursor = round(float((row.split('=')[1]).split (' ')[0]),4)
        
        elif row.split('=')[0] == 'CHARGE':
            charge_states = row.split ('=')[1][:-1]
        
        elif type (float(row.split (' ')[0])) == float :

            '''print ('input_m/z')'''

            obs_m_to_z = round(float(row.split(' ')[0]), 4)
            obs_inten = round(float(row.split(' ')[1][:-1]), 4)
            peak_list.append([obs_m_to_z, obs_inten])
    raw_mgf_file.close()

    

    output_file = open (file_route + 'raw_mgf_to_dict.txt', 'w')
    output_file.write(str(raw_mgf_MS2_spectra_list))
    output_file.close()

    return raw_mgf_MS2_spectra_list

'''file_path ='G:\\Desktop\\NAPI-tag\\Scripts for Ac-IPG tag\\'
test_file_1 = 'Ac-IPG-BSA-1-1-1.mgf'



raw_mgf_MS2_spectra_list = converting_raw_mgf_to_dict ([file_path,test_file_1])'''


###def extract_PSMs_and_return_all_theoretial_ions (all_unique_peptides,file_route_and_file_name):
def extract_PSMs_and_return_all_theoretial_ions (file_route_and_file_name):

    import csv
#importing csv module
    import re
#importing re module for removing the numbers in row of .csv file

    space = ''
    cvs_every_spectrum = {}
    retention_time_collection = []
    
    file_route = file_route_and_file_name [0]
    file_name = file_route_and_file_name [1]

    with open(file_route + file_name,encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
#calling the csv reader and read the target csv file
    
        for row in reader:
            spectrum_information = []
            retention_time = row['RT']
            print (retention_time)
            if retention_time == '':
                break
            else:
                retention_time = row['RT']
                m_over_z = round (float (row['m/z']),4)
                ###retention_time = round (float (row['RT'])*60,4)
                ###peptide = space.join (re.findall('[A-Z]+',str (row['Peptide'])))
                peptide = row['Sequence']
                mass = round (float (row['Theoretical Mass']),4)
                charge_states = row['Measured Charge']
                ###logP = round (float (row['-10lgP']),2)
                localization_confidence = row['Localization Confidence']
                ###length = int (row['Length'])
                ppm = round (float (row['Precursor m/z Error [ppm]']),4)
                ###area = row['Area']
                protein_accession = row['Protein(s)']

                ###if peptide in all_unique_peptides and retention_time not in retention_time_collection:
                    
                pep_seq = peptide
                import ion_produceor_Ac_IPG_3channels
                all_theoretial_ions = ion_produceor_Ac_IPG_3channels.combine_b_y_theoretical_ions(pep_seq)
                #print (all_theoretial_ions)
                        
                spectrum_information.append (['m_over_z',m_over_z])
                ###spectrum_information.append (['retention_time',retention_time])
                spectrum_information.append (['peptide',peptide])
                spectrum_information.append (['mass',mass])
                spectrum_information.append (['charge_states',charge_states])
                spectrum_information.append (['localization_confidence',localization_confidence])
                ###spectrum_information.append (['length',length])
                spectrum_information.append (['ppm',ppm])
                ###spectrum_information.append (['area',area])
                spectrum_information.append (['protein_accession',protein_accession])
                all_theoretial_ions.append (['spectrum_information',spectrum_information])
                cvs_every_spectrum [retention_time] = all_theoretial_ions
                retention_time_collection.append (retention_time)

    output_file = open (file_route + 'Dict_all_fragment_ions.txt', 'w')
    output_file.write(str(cvs_every_spectrum))
    output_file.close()

    return cvs_every_spectrum

'''file_path ='G:\\Desktop\\NAPI-tag\\Scripts for Ac-IPG tag\\'
test_file = 'Default PSM Report.csv'

all_uni_peps = extract_PSMs_and_return_all_theoretial_ions ([file_path,test_file])'''



        

       

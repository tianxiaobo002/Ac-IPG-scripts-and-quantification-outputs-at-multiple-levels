def b_ion_type_name_mass_channels_name_mass_producor (pep_seq):
    single_aa_mass = {"G":57.0215,"A":71.0371,"L":113.0814,"I":113.0814,"V":99.0684,"F":147.0684,"M":131.0405,"P":97.0528,"D":115.0269,"E":129.0426,
                  "N":114.0429,"Q":128.0586,"S":87.0320,"T":101.0477,"R":156.1011,"K":128.0950,"H":137.0589,"W":186.0793,"Y":163.0633,"C":103.0092}
    modification = {"H_hydrogen":1.0078,"H_proton":1.0078,"OH_hydroxy":17.0027,"Oxidation_Met":15.9949,"Carbamidomethylation":57.0215,"Ac-Cys-Ma at N":262.0385,
                "Ac-Ala at C-K":114.0555,"13C+1-Ac-Cys-Ma at N":263.0419,"13C+1-Ac-Ala at C-K":115.0589,"13C+2-Ac-Cys-Ma at N":264.0453,"13C+2-Ac-Ala at C-K":116.0623,
                "13C+3-Ac-Cys-Ma at N":265.0487,"13C+3-Ac-Ala at C-K":117.0657, "H2O":18.0106,"NH3":17.0265,"Pro-Gly+0":155.0820,"Pro-Gly+1":156.0854,"Pro-Gly+2":157.0888,"Dimeth_N-termi":29.0391}
    
    all_b_ions = []
    ion_type_name_seq_mass = []
    
    for i in range (1,len (pep_seq)):
        ion_type_name_temp = 'b'+str (i)
        ion_type_seq_temp = pep_seq [0:i]
        ion_type_name_seq_mass.append ([ion_type_name_temp,ion_type_seq_temp])
        
    for j in range (0,len (ion_type_name_seq_mass)):
        mass_cal_temp = 0
        for each_amino_acid in ion_type_name_seq_mass[j][1]:
            mass_cal_temp = mass_cal_temp + single_aa_mass [each_amino_acid]
            if each_amino_acid == 'C':
                mass_cal_temp = mass_cal_temp + modification ['Carbamidomethylation']###because Cys was fixed modified
        mass_cal_temp = mass_cal_temp + modification ['Dimeth_N-termi']
        ion_type_name_seq_mass [j].append (mass_cal_temp)
        
    ion_type_name_seq_mass_Met_oxida = []
    for j in range (0,len (ion_type_name_seq_mass)):
        Met_counter = ion_type_name_seq_mass[j][1].count ('M')
        if Met_counter != 0:
            for oxidation_number in range (1,(Met_counter+1)):
                ion_type_name_temp = ion_type_name_seq_mass[j][0] + '_'+ str (oxidation_number) +'Met_oxida'
                mass_cal_temp = ion_type_name_seq_mass[j][2] + (15.9949*oxidation_number)
                ion_type_name_seq_mass_Met_oxida.append ([ion_type_name_temp,ion_type_name_seq_mass[j][1],mass_cal_temp])
    ion_type_name_seq_mass.extend (ion_type_name_seq_mass_Met_oxida)
                
    ion_type_name_seq_mass_loss_H2O_NH3_2charge = []            
    for k in range (0,len (ion_type_name_seq_mass)):
        ion_type_name_temp_loss_H2O = ion_type_name_seq_mass[k][0] + '_'+ '-H2O'
        mass_cal_temp_loss_H2O = ion_type_name_seq_mass[k][2] - 18.0106
        ion_type_name_temp_loss_NH3 = ion_type_name_seq_mass[k][0] + '_'+ '-NH3'
        mass_cal_temp_loss_NH3 = ion_type_name_seq_mass[k][2] - 17.0265
        ion_type_name_seq_mass_loss_H2O_NH3_2charge .append ([ion_type_name_temp_loss_H2O,ion_type_name_seq_mass[k][1],mass_cal_temp_loss_H2O])
        ion_type_name_seq_mass_loss_H2O_NH3_2charge .append ([ion_type_name_temp_loss_NH3,ion_type_name_seq_mass[k][1],mass_cal_temp_loss_NH3])

        basic_amino_acid_counter = ion_type_name_seq_mass[k][1].count ('K')+ion_type_name_seq_mass[k][1].count ('R')+ion_type_name_seq_mass[k][1].count ('H')
        if basic_amino_acid_counter > 2:
            ion_type_name_temp_2charge = ion_type_name_seq_mass[k][0] + '_' + '2charge'
            mass_cal_temp_2charge = (ion_type_name_seq_mass[k][2]+1.0078)/2
            ion_type_name_seq_mass_loss_H2O_NH3_2charge .append ([ion_type_name_temp_2charge,ion_type_name_seq_mass[k][1],mass_cal_temp_2charge])
    ion_type_name_seq_mass.extend (ion_type_name_seq_mass_loss_H2O_NH3_2charge)

    '''ion_type_name_seq_mass_modifi_loss_all_channels = []###This was made for multiplex labelled b_ions. But, for Ac-IPG tag, b_ions from all channels are same.
    for l in range (0,len (ion_type_name_seq_mass)):
        for m in range (0,4):
            ion_type_name_temp_channel = ion_type_name_seq_mass[l][0] + '_' + '13C+' + str (m)
            mass_cal_temp_channel = round (ion_type_name_seq_mass[l][2] + (1.0034*m),4)
            ion_type_name_seq_mass_modifi_loss_all_channels.append ([ion_type_name_temp_channel,ion_type_name_seq_mass[l][1],mass_cal_temp_channel])'''
    
    for each_mass in ion_type_name_seq_mass:
        each_mass [2] = round (each_mass [2],4)
        
    all_b_ions = ion_type_name_seq_mass
    return all_b_ions

def y_ion_type_name_mass_channels_name_mass_producor (original_pep_seq):
    single_aa_mass = {"G":57.0215,"A":71.0371,"L":113.0814,"I":113.0814,"V":99.0684,"F":147.0684,"M":131.0405,"P":97.0528,"D":115.0269,"E":129.0426,
                  "N":114.0429,"Q":128.0586,"S":87.0320,"T":101.0477,"R":156.1011,"K":128.0950,"H":137.0589,"W":186.0793,"Y":163.0633,"C":103.0092}
    modification = {"H_hydrogen":1.0078,"H_proton":1.0078,"OH_hydroxy":17.0027,"Oxidation_Met":15.9949,"Carbamidomethylation":57.0215,"Ac-Cys-Ma at N":262.0385,
                "Ac-Ala at C-K":114.0555,"13C+1-Ac-Cys-Ma at N":263.0419,"13C+1-Ac-Ala at C-K":115.0589,"13C+2-Ac-Cys-Ma at N":264.0453,"13C+2-Ac-Ala at C-K":116.0623,
                "13C+3-Ac-Cys-Ma at N":265.0487,"13C+3-Ac-Ala at C-K":117.0657, "H2O":18.0106,"NH3":17.0265,"Pro-Gly+0":155.0820,"Pro-Gly+1":156.0854,"Pro-Gly+2":157.0888,"Dimeth_N-termi":29.0391}
    pep_seq_to_list = [] 
    for each_amico_acid in original_pep_seq:
        pep_seq_to_list.append(each_amico_acid)
    pep_seq_to_list.reverse ()
    pep_seq = ''.join(pep_seq_to_list)
    
    all_y_ions = []
    ion_type_name_seq_mass = []
    for i in range (1,len (pep_seq)+1):###plus 1 is for the quantification ion which was produced by the neutral lose of Ac-Ile 
        
        ion_type_name_temp = 'y'+str (i)
        ion_type_seq_temp = original_pep_seq [(0-i):]
        ion_type_name_seq_mass.append ([ion_type_name_temp,ion_type_seq_temp])
        
    for j in range (0,len (ion_type_name_seq_mass)):
        mass_cal_temp = 0
        for each_amino_acid in ion_type_name_seq_mass[j][1]:
            mass_cal_temp = mass_cal_temp + single_aa_mass [each_amino_acid]
            if each_amino_acid == 'C':
                mass_cal_temp = mass_cal_temp + modification ['Carbamidomethylation']
        if ion_type_name_seq_mass[j][1] == original_pep_seq:
            mass_cal_temp = mass_cal_temp + modification ['OH_hydroxy'] + modification ['Pro-Gly+0'] + modification ['Dimeth_N-termi']
        else:
            mass_cal_temp = mass_cal_temp + modification ['OH_hydroxy'] + modification ['H_proton'] + modification ['Pro-Gly+0']
        ion_type_name_seq_mass [j].append (mass_cal_temp)
        
    ion_type_name_seq_mass_Met_oxida = []
    for j in range (0,len (ion_type_name_seq_mass)):
        Met_counter = ion_type_name_seq_mass[j][1].count ('M')
        if Met_counter != 0:
            for oxidation_number in range (1,(Met_counter+1)):
                ion_type_name_temp = ion_type_name_seq_mass[j][0] + '_'+ str (oxidation_number) +'Met_oxida'
                mass_cal_temp = ion_type_name_seq_mass[j][2] + (15.9949*oxidation_number)
                ion_type_name_seq_mass_Met_oxida.append ([ion_type_name_temp,ion_type_name_seq_mass[j][1],mass_cal_temp])
    ion_type_name_seq_mass.extend (ion_type_name_seq_mass_Met_oxida)
                
    ion_type_name_seq_mass_loss_H2O_NH3_2charge = []            
    for k in range (0,len (ion_type_name_seq_mass)):
        ion_type_name_temp_loss_H2O = ion_type_name_seq_mass[k][0] + '_'+ '-H2O'
        mass_cal_temp_loss_H2O = ion_type_name_seq_mass[k][2] - 18.0106
        ion_type_name_temp_loss_NH3 = ion_type_name_seq_mass[k][0] + '_'+ '-NH3'
        mass_cal_temp_loss_NH3 = ion_type_name_seq_mass[k][2] - 17.0265
        ion_type_name_seq_mass_loss_H2O_NH3_2charge .append ([ion_type_name_temp_loss_H2O,ion_type_name_seq_mass[k][1],mass_cal_temp_loss_H2O])
        ion_type_name_seq_mass_loss_H2O_NH3_2charge .append ([ion_type_name_temp_loss_NH3,ion_type_name_seq_mass[k][1],mass_cal_temp_loss_NH3])

        basic_amino_acid_counter = ion_type_name_seq_mass[k][1].count ('K')+ion_type_name_seq_mass[k][1].count ('R')+ion_type_name_seq_mass[k][1].count ('H')
        if basic_amino_acid_counter > 2 or ion_type_name_seq_mass[k][1] == original_pep_seq:
            ion_type_name_temp_2charge = ion_type_name_seq_mass[k][0] + '_' + '2charge'
            mass_cal_temp_2charge = (ion_type_name_seq_mass[k][2]+1.0078)/2
            ion_type_name_seq_mass_loss_H2O_NH3_2charge .append ([ion_type_name_temp_2charge,ion_type_name_seq_mass[k][1],mass_cal_temp_2charge])
        
    ion_type_name_seq_mass.extend (ion_type_name_seq_mass_loss_H2O_NH3_2charge)

    ion_type_name_seq_mass_modifi_loss_all_channels = []
    for l in range (0,len (ion_type_name_seq_mass)):

        ###print(ion_type_name_seq_mass[l][0][-7:])
        
        if ion_type_name_seq_mass[l][0][-7:] == '2charge':
            charge_status = 2
        else:
            charge_status = 1
            
        for m in range (0,3):###3 represent the 3-plex
            ion_type_name_temp_channel = ion_type_name_seq_mass[l][0] + '_' + '13C+' + str (m)
            mass_cal_temp_channel = round (ion_type_name_seq_mass[l][2] + ((1.0034/charge_status)*m),4)
            ion_type_name_seq_mass_modifi_loss_all_channels.append ([ion_type_name_temp_channel,ion_type_name_seq_mass[l][1],mass_cal_temp_channel])
    all_y_ions = ion_type_name_seq_mass_modifi_loss_all_channels
    return all_y_ions
        
'''def combine_b_y_theoreticals (b_theoretial_ions,y_theoretial_ions):
    theoretial_ions = []
    theoretial_ions.extend (b_theoretial_ions)
    theoretial_ions.extend (y_theoretial_ions)
    return theoretial_ions'''
def combine_b_y_theoretical_ions (unique_peptide):
    theoretial_ions = []
    b_theoretial_ions = b_ion_type_name_mass_channels_name_mass_producor (unique_peptide)
    y_theoretial_ions = y_ion_type_name_mass_channels_name_mass_producor (unique_peptide)
    
    theoretial_ions.extend (b_theoretial_ions)
    theoretial_ions.extend (y_theoretial_ions)
    
    return theoretial_ions



'''###sequence = 'GTDWLANK'
###sequence = 'GWLYRAK'
sequence = 'SEIAHRFK'
###sequence = 'YNGVFQECCQAEDK'

fragment_ions = combine_b_y_theoretical_ions (sequence)

print (fragment_ions)'''





'''import Extract_unique_peptides_return_uni_pep_list

file_path ='G:\\Desktop\\IPTL\\4-plex labeling quantification scripts\\Under writing\\'
test_file_1 = 'protein-peptides_Ac-Cys-plus2.csv'
test_file_2 = 'protein-peptides_Ac-Cys-plus3.csv'
test_file_3 = 'protein-peptides_Ac-Cys-plus0.csv'
test_file_4 = 'protein-peptides_Ac-Cys-plus1.csv'

unique_peptides = Extract_unique_peptides_return_uni_pep_list.extract_and_return_uni_pep_list ([file_path,test_file_1],[file_path,test_file_2],[file_path,test_file_3],[file_path,test_file_4])'''

import indentified_PSM_csv_to_all_information_list

file_path ='G:\\Desktop\\NAPI-tag\\Ac-IPG tag manuscript preparing\\BSA-Yeast analysis\\'
test_file = 'Default PSM Report with non-validated matches.csv'

all_thero_fragment_ions = indentified_PSM_csv_to_all_information_list.extract_PSMs_and_return_all_theoretial_ions ([file_path,test_file])

import Converting_the_mgf_to_py_dict

file_path = 'G:\\Desktop\\NAPI-tag\\Ac-IPG tag manuscript preparing\\BSA-Yeast analysis\\'
test_file_1 = 'Ac-IPG-BSA1510_yeast_125.mgf'

raw_mgf_to_dict = Converting_the_mgf_to_py_dict.converting_raw_mgf_to_dict ([file_path,test_file_1])

import matching_theo_frag_ion_with_raw_mgf_peak_list

"""file_path ='G:\\Desktop\\NAPI-tag\\Ac-IPG tag manuscript preparing\\Analysis of yeast 111\\'"""

all_thero_fragment_ions_added_measured_intens = matching_theo_frag_ion_with_raw_mgf_peak_list.adding_measured_intensity_to_thero_frag_ions (file_path,all_thero_fragment_ions,raw_mgf_to_dict)

import calculate_normalized_fragment_ion_ratio_for_spec

"""file_path ='G:\\Desktop\\NAPI-tag\\Ac-IPG tag manuscript preparing\\Analysis of yeast 111\\'"""

all_spec_normalized_ratio_and_thero_fragment_ions_added_measured_intens = calculate_normalized_fragment_ion_ratio_for_spec.normalize_ratio_for_spec (file_path,all_thero_fragment_ions_added_measured_intens)



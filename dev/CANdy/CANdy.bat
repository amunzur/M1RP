
SET alpha=0.01
SET cohort=M1RP


python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\normal_samples_cn\create_normal_gene_distribution.py" 
ECHO Normal samples up to date
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\calculate_p_value.py" %cohort% %alpha%
ECHO Adjusted copy number calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\calculate_min_tc.py" %cohort% %alpha%
ECHO Minimum tumor content calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\calculate_copy_number.py" %cohort% %alpha%
ECHO Copy number call calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\visualize_cn_stats.py" %cohort%
ECHO Visualized per sample copy number plots
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\visualize_patient_cn.py" %cohort%
ECHO Visualized per patient copy number plots

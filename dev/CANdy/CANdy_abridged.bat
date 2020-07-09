
SET alpha=0.005
SET cohort=M1B
SET abridged=True


python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\normal_samples_cn\create_normal_gene_distribution.py" 
ECHO Normal samples up to date
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\calculate_p_value.py" %cohort% %alpha% %abridged%
ECHO Adjusted copy number calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\calculate_min_tc.py" %cohort% %alpha% %abridged%
ECHO Minimum tumor content calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\Targeted_sequencing_error_correction\calculate_copy_number_abridged.py" %cohort% %alpha%
ECHO Copy number call calculated

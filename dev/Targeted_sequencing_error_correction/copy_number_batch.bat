
SET alpha=0.001

python "G:\Andy Murtha\Ghent\M1RP\dev\Targeted_sequencing_error_correction\cn_add_error.py"
ECHO Adjusted copy number calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\Targeted_sequencing_error_correction\calculate_min_tc.py" %alpha%
ECHO Minimum tumor content calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\Targeted_sequencing_error_correction\add_copy_number.py" %alpha%
ECHO Copy number call calculated
python "G:\Andy Murtha\Ghent\M1RP\dev\Targeted_sequencing_error_correction\visualize_cn_stats.py"
ECHO Visualized per sample copy number plots
python "G:\Andy Murtha\Ghent\M1RP\dev\Targeted_sequencing_error_correction\visualize_patient_cn.py"
ECHO Visualized per patient copy number plots



set cohort=M1B
set min_reads=75
set lr_max="0.3"

python "C:\Users\amurtha\Dropbox\Ghent M1 2019\sandbox\mutations\betastasis\melt_betastasis.py" "C:\Users\amurtha\Dropbox\Ghent M1 2019\sandbox\mutations\betastasis\%cohort%_betastasis_all.xlsx" "C:\Users\amurtha\Dropbox\Ghent M1 2019\sandbox\mutations\%cohort%_mutations.tsv" 

python "G:\Andy Murtha\Ghent\M1RP\dev\tumor_fraction_calcuation\tumor_fraction_calculation.py" %cohort% %min_reads% %lr_max%

python "G:\Andy Murtha\Ghent\M1RP\dev\tumor_fraction_calcuation\uncalled_mutation_check.py" %cohort% %min_reads% %lr_max%

python "G:\Andy Murtha\Ghent\M1RP\dev\SNP_analysis\snp_tc_calculator.py" %cohort%
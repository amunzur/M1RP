
SET cohort=M1RP
SET nascent=False
SET min_reads=30
SET max_lr=0.3

echo "Update mutations"
python "C:\Users\amurtha\Dropbox\Ghent M1 2019\sandbox\mutations\betastasis\melt_betastasis.py" %cohort%
echo "Update independent mutation"
python "C:\Users\amurtha\Dropbox\Ghent M1 2019\sandbox\mutations\betastasis\melt_betastasis_inclDependent.py" %cohort%
echo "Update CNA"
python "G:\Andy Murtha\Ghent\M1RP\dev\CANdy\normal_samples_cn\cn_to_tsv.py" 

echo "Mutation TC calls"
python "G:\Andy Murtha\Ghent\M1RP\dev\tumor_fraction_calcuation\mut_tc_calculation.py" %cohort% %min_reads% %max_lr%
echo "SNP TC calls"
python "G:\Andy Murtha\Ghent\M1RP\dev\SNP_analysis\snp_tc_calculator.py" %cohort% %nascent%
echo "Create final file and plot"
python "G:\Andy Murtha\Ghent\M1RP\dev\tumor_fraction_calcuation\snpTC_mutTC_byGroup.py" 
echo "Done"
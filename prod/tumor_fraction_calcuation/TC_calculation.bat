
SET cohort=M1RP
SET nascent=False
SET min_reads=30
SET max_lr=0.3

echo "Update mutations"
python "C:\Users\amurtha\Dropbox\Ghent M1 2019\Mar2021_datafreeze\mutations\betastasis\melt_betastasis.py" %cohort%
echo "Update independent mutation"
python "C:\Users\amurtha\Dropbox\Ghent M1 2019\Mar2021_datafreeze\mutations\betastasis\melt_betastasis_inclDependent.py" %cohort%
echo "Update SNPs"
python "G:\Andy Murtha\Ghent\M1RP\prod\SNP_analysis\melt_snps.py"
echo "SNP TC calls"
python "G:\Andy Murtha\Ghent\M1RP\prod\SNP_analysis\snp_tc_calculator.py" %cohort% %nascent%

echo "Mutation TC calls"
python "G:\Andy Murtha\Ghent\M1RP\prod\tumor_fraction_calcuation\mut_tc_calculation.py" %cohort% %min_reads% %max_lr%

echo "Create final file and plot"
python "G:\Andy Murtha\Ghent\M1RP\prod\tumor_fraction_calcuation\snpTC_mutTC_byGroup.py" 
echo "Create final file and plot"
python "G:\Andy Murtha\Ghent\M1RP\prod\tumor_fraction_calcuation\TC_bland_altman.py" 
echo "Done"
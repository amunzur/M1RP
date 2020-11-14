cd /groups/wyattgrp/eritch/projects/ghent_m1rp/wxs/phylowgs
conda activate phylowgs
for pt in *;
do
	cd /groups/wyattgrp/eritch/projects/ghent_m1rp/wxs/phylowgs/$pt;
	mkdir test_results;
	cd test_results;
	python /home/amurtha/anaconda3/envs/phylowgs/share/phylowgs/write_results.py $pt ../chains/trees.zip $pt.summ.json.gz $pt.muts.json.gz $pt.mutass.zip;
	cd ..;
	mv test_results/* /home/amurtha/anaconda3/envs/phylowgs/share/phylowgs/witness/data/test_results;
	cd /home/amurtha/anaconda3/envs/phylowgs/share/phylowgs/witness
	gunzip data/*/*.gz;
	python index_data.py;
done;
python2 -m SimpleHTTPServer

cd /p/keles/fandingzhou/volumeA/software
for j in 2 3 4 6
do
java -Xmx2g -jar juicer_tools_1.22.01.jar pre -r 40000 /p/keles/fandingzhou/volumeA/DCI/juicebox/GM12878_rep${j}_chr1.binPairs.chr1 /p/keles/fandingzhou/volumeA/DCI/juicebox/GM12878_rep$j.hic hg19
done

j=2
for i in 2 3 4 6
do
#for i in 2 3 4 6
#do
java -Xmx2g -jar juicer_tools_1.22.01.jar pre -r 40000 /p/keles/fandingzhou/volumeA/DCI/juicebox/GM12878_strategy4_FC${i}_rep${j}_chr1.binPairs.chr1 /p/keles/fandingzhou/volumeA/DCI/juicebox/GM12878_strategy4_FC${i}_rep$j.hic hg19
#done
done

for i in 3 4 6
do
scp fzhou49@hodor01.stat.wisc.edu:/p/keles/fandingzhou/volumeA/DCI/juicebox/\{GM12878_strategy5_FC${i}_rep2.hic,GM12878_strategy5_FC${i}_rep3.hic,GM12878_strategy5_FC${i}_rep4.hic,GM12878_strategy5_FC${i}_rep6.hic\} .
done

cd /p/keles/fandingzhou/volumeA/DCI/juicebox
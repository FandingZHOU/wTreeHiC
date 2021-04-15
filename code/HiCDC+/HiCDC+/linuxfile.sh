for j in {11..20};
do
cd seed$j
mkdir 4vs4
cd 4vs4
mkdir filter_as_HiCDCPlus  filter_as_diffHic  filter_as_multiHiCcompare
cd ..
done

maxsize=100
for seed in {5..10};
do
for FC in 2 3 4 6;
do
OnTAD2.sh ${FC} 1000 ${seed}
done
done



for k in {11..20}
do
mkdir seed${k}
cd seed${k}
mkdir 4vs4 2vs2
for h in 4vs4 2vs2
do
cd $h
mkdir filter_as_HiCDCPlus  filter_as_diffHic  filter_as_multiHiCcompare
for j in filter_as_HiCDCPlus  filter_as_diffHic  filter_as_multiHiCcompare;
do
cd $j
mkdir  HiCDCPlus  multiHiCcompare  Trueset diffHic
for i in HiCDCPlus  multiHiCcompare diffHic
do
cd $i
mkdir grouptable  originaltable  weighttable
cd ..
done
cd ..
done
cd ..
done
cd ..
done

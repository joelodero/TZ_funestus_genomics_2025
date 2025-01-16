ADMIXPATH=~/software/admixture_linux-1.3.0/admixture
PLINKPATH=/home/dennist/lstm_scratch/funestus/plink/2RL.50000.0.05.0.0.EA.bed

for SEED in 12345 58364 48573 23424 47566
do
	mkdir $SEED
	cd $SEED
	for K in 1 2 3 4 5 6 7 8 9
	do
		$ADMIXPATH -s $SEED -j10 $PLINKPATH $K | tee log${K}.out 
	done
	cd ..
done

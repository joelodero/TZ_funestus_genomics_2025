PFILE=/home/dennist/lstm_scratch/funestus/plink/2RL.50000.0.05.0.0.TZ.bed
ADMIX=~/software/admixture_linux-1.3.0/admixture

for SEED in 1234 4857 5854 2943 7234
do
	mkdir $SEED.TZ
	cd $SEED.TZ
		for K in 1 2 3 4 5 6 7 8 9
		do
			$ADMIX -j10 --cv=5 $PFILE $K | tee log${K}.out
		done
		cd ..
	done	

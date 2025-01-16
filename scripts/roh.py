#create conda env for malariagen
#install malariagen data, dask, pandas, multiprocessing


import malariagen_data
from multiprocessing import Process
from dask.distributed import Client
import pandas as pd


af1 = malariagen_data.Af1(pre=True)

sample_sets_tz='1236-VO-TZ-OKUMU-OKFR-TZ-2008', '1236-VO-TZ-OKUMU-VMF00248','1236-VO-TZ-OKUMU-VMF00252','1236-VO-TZ-OKUMU-VMF00261','1236-VO-TZ-OKUMU-VMF00090'
tz_sample_df=af1.sample_metadata(sample_sets_tz)

def chr2():
	print('converting chr2')
	if __name__ == "__main__":

		client = Client(n_workers=1, threads_per_worker=1)

		dflist = []
		for sample in tz_sample_df['sample_id']:
			roh_df = af1.roh_hmm(sample = sample, region = "2RL")
			roh_list.append(roh_df)

		roh_bigdf = pd.concat(roh_list)
		roh_bigdf.to_csv('/home/dennist/lstm_data/funestus/ROH_2RL.csv') #modify path here



def chr3():
	print('converting chr3')
	if __name__ == "__main__":

		client = Client(n_workers=1, threads_per_worker=1)

		dflist = []
		for sample in tz_sample_df['sample_id']:
			roh_df = af1.roh_hmm(sample = sample, region = "3RL")
			roh_list.append(roh_df)

		roh_bigdf = pd.concat(roh_list)
		roh_bigdf.to_csv('/home/dennist/lstm_data/funestus/ROH_3RL.csv') #modify path here



def chrX():
	print('converting chrX')
	if __name__ == "__main__":

		client = Client(n_workers=1, threads_per_worker=1)

		dflist = []
		for sample in tz_sample_df['sample_id']:
			roh_df = af1.roh_hmm(sample = sample, region = "X")
			roh_list.append(roh_df)

		roh_bigdf = pd.concat(roh_list)
		roh_bigdf.to_csv('/home/dennist/lstm_data/funestus/ROH_X.csv') #modify path here


#do in parallel
if __name__ == "__main__":
	p1 = Process(target=chr2)
	p1.start()
	p2 = Process(target=chr3)
	p2.start()
	p3 = Process(target=chrX)
	p3.start()

	p1.join()
	p2.join()
	p3.join()

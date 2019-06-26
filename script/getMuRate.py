import numpy as np
import glob
eps = 1e-8
files = glob.glob('../p17b/data_rmu/*.root')

n_range = 3

r_mu = np.zeros((3,4,n_range,3))

n_mu = np.zeros((3,4,n_range,3))
livetime = np.zeros((3,4))

for idx, f in enumerate(files):
	with open(f) as _f:
		lines = _f.readlines()
		site = int(lines[0].split()[0])
		if site == 4:
			site = 2
		else:
			site -= 1
		ads = int(lines[0].split()[1])
		for ad in range(ads):
			livetime[site][ad] += float(lines[1 + 2 * ad]) / 1e6
			n_tmp = [int(_tmp) for _tmp in lines[2 + 2 * ad].split()]
			for r in range(n_range):
				for i in range(3):
					n_mu[site][ad][r][i] += n_tmp[3 * r + i]				
	
for site in range(3):
	for ad in range(4):
		for r in range(n_range):
			for t in range(3):
				r_mu[site][ad][r][t] = n_mu[site][ad][r][t] / (livetime[site][ad] + eps)
		
for site in range(3):
	for r in range(n_range):
		for t in range(3):
			_sum = 0
			_sum_lt = 0
			for ad in range(4):
				_sum_lt += livetime[site][ad]
				_sum += n_mu[site][ad][r][t]
			print _sum / _sum_lt

			
for site in range(3):
	_sum_lt = 0.
	for ad in range(4):
		_sum_lt += livetime[site][ad]
	print _sum_lt / 1e3	/ 86400

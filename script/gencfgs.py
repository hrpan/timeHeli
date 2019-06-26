import itertools

fitMin = [ i for i in range(5, 101, 10) ]
fitMax = [5000]

fix_B12 = [0,]
fix_He8 = [0,]
bound_eps = [0, 1]
fix_lifetime = [1,]
fix_rmu = [0,]
use_eps_pull = [0, 1]

arg_lists = [
	fitMin,
	fitMax,
	fix_B12,
	fix_He8,
	bound_eps,
	fix_lifetime,
	fix_rmu,
	use_eps_pull
]


for _min, _max, _b12, _he8, _eps, _lt, _rmu, _pull in itertools.product(*arg_lists):
	if _pull == 1:
		_eps = 0
	filename = './cfit_cfgs/cfg_%d_%d_%d_%d_%d_%d_%d_%d' % (_min, _max, _b12, _he8, _eps, _lt, _rmu, _pull)

	#if ( _min >= 50 and _b12 == 0 ) or ( _min < 50 and _b12 == 1 ):
	#	continue
	with open(filename, 'w') as f:
		f.write(str(_min) + '\n')
		f.write(str(_max) + '\n')
		f.write(str(_b12) + '\n')
		f.write(str(_he8) + '\n')
		f.write(str(_eps) + '\n')
		f.write(str(_lt) + '\n')
		f.write(str(_rmu) + '\n')
		f.write(str(_pull) + '\n')
		f.close()

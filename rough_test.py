import csv
import numpy as np
from radar_rough_funcs import detrend
from radar_rough_funcs import calc_roughness_profile
from radar_rough_funcs import calc_hurst_rmsslope
from radar_rough_funcs import calc_fractalrough_DTM
from datetime import datetime
from datetime import timedelta

# Test script to make sure things are working correctly:

def test1():
	test_filename = './sample_prof.csv'
	test_file = open(test_filename, 'r')
	X = []
	Z = []
	dx_array = 15+15*np.array(range(20))
	with test_file:
		csv_reader = csv.reader(test_file)
		for row in csv_reader:
			X.append(row[0])
			Z.append(row[1])

	X = np.array(X, dtype='float')
	Z = np.array(Z, dtype='float')
	lam = 15 

	h_array, nu_array, slope_array = calc_roughness_profile(X, Z, dx_array)
	H, s_l, r2 = calc_hurst_rmsslope( dx_array, slope_array, lam, 1)

	for x in range(len(slope_array)):
		print('Delta X = {0} m: {1} Degrees, {2} m'.format(dx_array[x], slope_array[x], h_array[x]))

	print('H = {}'.format(H))
	print('s_l = {}'.format(s_l))
	print('R^2 = {}'.format(r2))

def main():
	dx_array = 15+50*np.array(range(30))
	lam = 15
	t1 = datetime.now()
	DTM_filename = '../Traces_Experiment/RawDTM/DTM066.tif'
	sl, sl_std, H, r2 =  calc_fractalrough_DTM( DTM_filename, dx_array)
	print('RMS Slope = {0} +- {1}; H = {2}, R2 = {3}'.format(sl, sl_std, H, r2))
	t2 = datetime.now()
	dt = t2-t1
	dt = dt.total_seconds()
	print('Time elapsed: {} s'.format(dt))

if __name__ == "__main__":
	main()

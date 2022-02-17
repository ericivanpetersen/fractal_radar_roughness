import csv
import numpy as np
from radar_rough_funcs import detrend
from radar_rough_funcs import calc_roughness_profile
from radar_rough_funcs import calc_hurst_rmsslope
from radar_rough_funcs import calc_fractalrough_DTM
from datetime import datetime
from datetime import timedelta

def fresnel_samples(object):
	"""Object containing information on DTM samples
	representative of individual SHARAD observations"""

	def __init__(self, filepath):
		"""Initiate by reading in summary file on
		SHARAD observations

		:param str filepath: file containing summary
			of SHARAD obs
		"""
	
		# Create fields to populate:
		self.sub_pow = []
		self.surf_pow = []
		self.depth = []
		self.point = []
		self.site = []
		self.daytime = []
		self.rough_param = []
		self.slope = []
		self.lon = []
		self.lat = []

		# Read in csv file:
		f = open(file_path, 'r')
		with f:
			content = f.readlines()
			line_num = 0
			for line in content:
				if line_num == 0: #Header
					line_num = line_num + 1
				else:
		# UNFINISHED! #
	

def main():
	dx_array = 15+15*np.array(range(20))
	lam = 15
	t1 = datetime.now()
	DTM_filename = '../Traces_Experiment/RawDTM/DTM001.tif'
	sl, sl_std, H, r2 =  calc_fractalrough_DTM( DTM_filename, dx_array)
	print('RMS Slope = {0} +- {1}; H = {2}, R2 = {3}'.format(sl, sl_std, H, r2))
	t2 = datetime.now()
	dt = t2-t1
	dt = dt.total_seconds()
	print('Time elapsed: {} s'.format(dt))

if __name__ == "__main__":
	main()

import numpy as np
import osr
import gdal
import matplotlib.pyplot as plt
from scipy import optimize

def detrend(X,Z):
	"""This function calculates the least-squares
	linear regression for the input profile and returns the 
	detrended elevations Z_detrend
	
	:param float X: array of distance along profile
	:param float Z: array of elevation (or other parameter)
		along profile
	:return float Z_detrend: array of detrended elevation (or other
		parameter) along profile
	:return float p: linear fit parameters corresponding to detrend
	"""
	
	p = np.polyfit(X,Z,1)
	Zp = np.polyval(p,X)
	Z_detrend = Z - Zp
	return Z_detrend, p

def movstd(X, dx):
	"""Calculates the standard deviation in a moving window 
	of length dx along the array X
	Returns array X_std which is of length X - dx, corresponding
		to elements dx/2 to len(X)-dx/2 of array X, thus 
		allowing all calculations of standard deviation 
		to be made with the same number of datapoints
	
	:param float X: values in an array
	:param int dx: size of moving window to calculate std in;
		units of array elements
	:return float X_std: array of std calculated in moving 
		window on array X
	"""
	
	XX = len(X)
	X_std = np.zeros(XX)
	for x in range(XX-dx):
		index = np.array(range(dx),dtype='int') + x
		X_std[x] = np.nanstd(X[index])
	return X_std

def calc_roughness_profile( X, Z, dx_array ):
	"""This function calculates roughness parameters 
	as a function of horizontal distance dx for 
	an input topographic profile 

	:param float X: distance along profile
	:param float Z: elevation of profile
	:param float dx_array: array of horizontal scales across which
			to calculate roughness parameters
	OUTPUTS (all are same size as dx_array):
	:return float nu_array: Allan deviation
	:return float slope_array: root-mean-square surface slope
	"""

	Z_dt, p_dt = detrend(X, Z) # Detrend topography
	LL = len(Z_dt) # Length of Profile
	XX = len(dx_array)
	DX = max(dx_array)
	
	min_length = 300 # Minimum profile length
	
	# Initialize arrays to fill:
	nu_array = np.zeros(XX)
	slope_array = np.zeros(XX)
	#h_array = np.zeros(XX)

	# Return NaNs if there's something fishy about this
	# 	profile, i.e.:
	# Range in detrended elevation > range in run
	#if ( np.ptp(Z_dt) > np.ptp(X) ):
	#	print("Profile has suspicious Z variation; skipping")
	#	exit()

	# Find "Spurs"
	spurs = np.absolute(Z_dt[1:] - Z_dt[:-1])
	spur_cut = 5 # "Spur" cut-off in meters
	spurs_ind = np.argwhere( np.absolute(spurs) > spur_cut )
	spurs_mask = np.where( np.absolute(spurs) > spur_cut )
	if len(spurs_ind) > 0: # If there are any spurs:
		rhs = LL-max(spurs_ind) #Length of spur-free profile on RHS
		lhs = min(spurs_ind) #Length of spur-free profile on LHS
		success = 0
		# If they are towards the left hand side of profile:
		if (rhs > lhs):
			# Determine if the clipped profile is long enough
			if ( (LL-max(spurs_ind)) > min_length):
				# Clip the profile and redo the detrend:
				X = X[max(spurs_ind)+1:]
				Z = Z[max(spurs_ind)+1:]
				Z_dt, p_dt = detrend(X, Z)
				success = 1
		# Or if they are towards the right hand side of profile:
		if (lhs > rhs):
			# Determine if the clipped profile is long enough
			if ( min(spurs_ind) > min_length):
				# Clip the profile and redo the detrend:
				X = X[:min(spurs)-1]
				Z = Z[:min(spurs)-1]
				Z_dt, p_dt = detrend(X, Z)
				success = 1
		# Or if the clipped profiles are not long enough
		if (success==0):
			print("Clipped Profile Not Long Enough; Continuing")
			return


	# Calculate roughness values for each value of dx in dx_array:
	for xx in range(XX):
		dx = dx_array[xx]
		dev = Z_dt[:-dx] - Z_dt[dx:] # Raw deviations
		dev = dev[ ( np.absolute(dev) < dx ) ] # Remove dev > dx
		element_slope = dev/dx # element-wise slope
		nu_array[xx] = np.std(dev)
		slope_array[xx] = np.std(element_slope)
		#prof_h_rms = movstd(Z_dt,dx)
		#prof_h_rms = prof_h_rms[ (prof_h_rms < dx ) ]
		#h_array[xx] = np.mean(prof_h_rms)

	return nu_array, slope_array

def fractal_eqn_rmsslope( dx_array, H, s_l, lam=15):
	"""This function defines the fractal equation for 
	rms slope as a function of horizontal scale dx.
	INPUTS:
		:param float H: Hurst Exponent
		:param float s_l: rms slope at reference scale lam
		:param int lam: reference scale/radar wavelength; 
			default = 15 (m). translates to array elements
		:param flaot dx_array: array of dx to calculate rms slope for
	OUTPUT: 
		:output float slope_array: rms slopes as function of 
			dx_array
	Written to be used as a fit function for fractal
		methods
	"""
	
	slope_array = s_l * (dx_array / lam)**(H-1)
	return slope_array

def calc_hurst_rmsslope( dx_array, slope_array, lam=15, plotnum='None' ):
	"""Function provides a power-law fit to rms slope
	as a function of delta x to calculate the hurst 
	exponent and rms slope at radar wavelenght lam
	INPUTS:
		:param float dx_array: horizontal scale delta x
		:param float slope_array: rms slope as a function of delta x
		:param float lam: radar wavelength lambda
		:param int plotnum: can use to plot results of fit
			set plotnum = 1 to plot
	OUTPUTS:
		:output float hurst: hurst exponent H
		:output float s_l: rms slope at radar wavelength lam
		:output float r2: R-squared value of power law fit"""

	# Calculate natural log of dx and rms slope:
	log_dx = np.log(dx_array)
	log_slope = np.log(slope_array)
	# Calculate least-squares fit:
	p0 = [0.5, slope_array[0] ] # Initial parameters to start the fitting
	params, params_cov = optimize.curve_fit(fractal_eqn_rmsslope, dx_array, slope_array, p0)
	slope_mod = fractal_eqn_rmsslope(dx_array, params[0], params[1]) #Slope resultant from model
	# Calculate R^2 residual of the fit:
	SSres = np.sum( (slope_array - slope_mod)**2 )
	SStot = np.sum( (slope_array - np.mean(slope_array))**2 )
	r2 = 1 - SSres / SStot * ( len(slope_array) - 1 ) / ( len(slope_array) - 2)

	# Plot fit if desired: 
	if (plotnum==1):
		plt.scatter(dx_array, slope_array)
		plt.plot(dx_array, slope_mod)
		plt.show()
	return params[0], params[1], r2

def calc_fractalrough_DTM( DTM_filename, dx_array, lam=15 ):
	"""Function reads a DTM in tiff format, 
	calculates roughness values for each profile that
	can be defined along X & Y axes, for each value dx
	in the dx_array, and finds the best fit hurst 
	exponent to the mean of those values as well as the
	reference roughness value at the radar wavelength
	lam
	
	:param str DTM_filename: DTM file to read in
	:param float dx_array: array of dx to assess roughness
	:param float lam: radar wavelength; default = 15m
	
	:output float sl_lam: mean rms slope at dx=lam
	:output float sl_lam_std: std of rms slope at dx=lam
	:output float H: Hurst exponent
	:output float r2: R^2 of fractal fit
	"""
	
	# Read in DTM file:
	df = gdal.Open(DTM_filename)
	band = df.GetRasterBand(1)
	DTM = band.ReadAsArray()
	DTM = np.array(DTM)
	# Define NaNs and Remove Empty Rows/Columns:
	nan_const = DTM[0,0]
	NaN_Mask = (DTM == nan_const) # Mask of NaN values
	DTM[ NaN_Mask] = np.nan # Replace values with NaN
	IsaN_Mask = np.logical_not(NaN_Mask) # Mask of elements with data
	cols = np.argwhere(np.sum(IsaN_Mask,0) > 0) # Columns with 1 or more data elements in them
	rows = np.argwhere(np.sum(IsaN_Mask,1) > 0) # Rows with 1 or more data elements in them
	rows = rows[:, np.newaxis]
	# Reduce DTM size to remove rows/cols of nan:
	DTM = DTM[ rows, cols ]
	# Recalculate nan and isan masks:
	NaN_Mask = np.isnan(DTM)
	IsaN_Mask = np.logical_not(NaN_Mask)

	# Calculate length of profiles and find those which are at least
	# 	2x the length of the max dx value:
	min_length = 300 # Minimum profile length
	prof_length_long = np.sum(IsaN_Mask,0)
	pind_long = np.argwhere(prof_length_long > min_length)
	prof_length_trans = np.sum(IsaN_Mask,1)
	pind_trans = np.argwhere(prof_length_trans > min_length)


	YY = len(pind_long)
	XX = len(pind_trans)
	DX = len(dx_array)
	# Initialize arrays:
	h_rms_x = np.zeros((XX,DX))
	h_rms_y = np.zeros((YY,DX))
	h_rms_arr = np.zeros((XX+YY,DX))
	sl_x = np.zeros((XX,DX))
	sl_y = np.zeros((YY,DX))
	sl_arr = np.zeros((XX+YY,DX))
	
	for xx in range(XX):
		X0 = range(YY)
		line = DTM[pind_trans[xx,0],X0]
		prof = line[ np.isfinite(line) ]
		X = range(len(prof))
		nu, sl = calc_roughness_profile(X, prof, dx_array) 
		#h_rms_x[xx,range(DX)] = h_rms
		#h_rms_arr[xx, range(DX)] = h_rms
		sl_x[xx, range(DX)] = sl
		sl_arr[xx, range(DX)] = sl
	for yy in range(YY):
		X0 = range(XX)
		line = DTM[X0,pind_long[yy,0]]
		prof = line[ np.isfinite(line) ]
		X = range(len(prof))
		nu, sl = calc_roughness_profile(X, prof, dx_array)
		#h_rms_y[yy, range(DX)] = h_rms
		#h_rms_arr[XX+yy, range(DX)] = h_rms
		sl_y[yy, range(DX)] = sl
		sl_arr[XX+yy, range(DX)] = sl

	mean_h_rms = np.nanmean(h_rms_arr,0)
	std_h_rms = np.nanstd(h_rms_arr,0)
	mean_sl = np.nanmean(sl_arr,0)
	std_sl = np.nanstd(sl_arr,0)

	xind = np.argwhere(dx_array == lam)
	h_rms_lam = mean_h_rms[xind]
	h_rms_lam_std = std_h_rms[xind]
	sl_lam = mean_sl[xind]
	sl_lam_std = std_sl[xind]

	[H, sl_lam_mod, r2] = calc_hurst_rmsslope(dx_array,mean_sl,lam)

	return sl_lam, sl_lam_std, H, r2

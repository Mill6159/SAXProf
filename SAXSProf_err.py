####################################
# import numpy as np
from SAXS8 import *
from PlotSAXS import *
from SAXProf_ErrCalcs import *
from matplotlib import pyplot as plt
# import SASM
#from scipy.stats.distributions import chi2
#from scipy.stats import chisqprob
# from scipy import stats
# from subprocess import check_output, CalledProcessError
# import colorsys
# import sys
import warnings
####################################

#############################################################################

# chisqprob is depricated. Replacement will be stats.distributions.chi2.sf

# derived from Fig3A_molmovdb_time
#
# from Fig3A_stats_energy_scan_detect_comp.py

### Version of code used to test self.sigma_profile_fit
# This code uses vac, window, and buffer subtracted profiles from
# the Nov 7 run to reproduce lysozyme data for the first time 
# in a way that you can change sample and window thickness correctly.

## This version has been modified for the atlastin models

#############################################################################
#############################################################################
# Do not change these parameters here. They need to be specific values that
# correspond to the experimental buffer, window, and vacuum profiles

###### Suppress RuntimeWarning from calculating PDDR ########

warnings.filterwarnings("ignore", category=RuntimeWarning)

#############################################################
# Initial User Inputs
#
#############################################################
energy = 9.962 # energy (keV)
P  = 8.4e11     # flux (photons/s)
t  = 1.0       # exposure time for individual snapshots (s)
snaps =  1      # number of snapshots averaged
a  = 150.47    # sample-to-detector distance (cm)
d  = 0.15       # sample cell path length (cm)

window_type = "mica"  # type of window material on sample cell

sensor_thickness = 0.032  # Pilatus sensor thickness

# These files below are only used when generating the "standard profiles."

# buffer with empty cell subtracted
buffer_model = "/S_A_nov07_028_lysbufnorm_flat.dat"
buffer_exposure = t  # exposure used for buffer model
buffer_flux = P      # flux used for buffer model

vac_model = "/A_nov07_072_lysbufnorm.dat"
vac_exposure = t
vac_flux = P

win_model = "/S_A_nov07_026_lysbufnorm.dat"
win_exposure = t
win_flux = P

#############################################################################
#############################################################################

sample_model_1 = "6lyz.pdb.dat"

c = 4.0 # concentration (mg.ml)
t = 10 # Exposure time (seconds)
mw1 = 14.3 # Molecular Weight (kDa)

saxs1 = SAXS(mw=mw1,a=a,d=d,energy=energy,P=P,t=t,total_t=t*snaps,c=c,shape='FoxS', detector='100K')
saxs1.set_window(window_type)
saxs1.sensor_thickness = sensor_thickness
saxs1.det_eff = saxs1.QE(saxs1.lam, saxs1.sensor_thickness)

#  determines the q-space range based on detector and mask. Default_q starts at q = 0, mask_q starts at q_min

saxs1.create_Mask(98, 3, 45, 14, wedge=360.0, type="rectangle") 

# buffer, vacuum, window, and model profiles read in and interpolated onto default_q
# saxs1.buf_q are the q points of default_q that fall within buf's q range
#
#saxs1.load_buf(buffer_model, t=buffer_exposure, P=buffer_flux, interpolate=True, q_array = saxs1.default_q)
#saxs1.load_vac(vac_model, t=vac_exposure, P=vac_flux, interpolate=True, q_array = saxs1.default_q)
#saxs1.load_win(win_model, t=win_exposure, P=win_flux, interpolate=True, q_array = saxs1.default_q)
saxs1.load_I(sample_model_1,interpolate=True,q_array = saxs1.default_q)

#############################################################################
# Here we load in
#############################################################################

saxs1.readStandardProfiles("M_Nov07_")

#############################################################################
#############################################################################

#saxs1.writeStandardProfiles("Nov07_")


#############################################################################
#############################################################################

saxs1.v = 0.72 # protein specific volume
saxs1.pp = 322e21/saxs1.v  # (from dry mass/specific volume)
saxs1.ps = 334.6e21  # electrons/cm3 for buffer (water)
saxs1.p = (saxs1.pp-saxs1.ps)*2.818e-13 # contrast

saxs1.d = 0.15     # microfluidic mixing chip
#saxs1.P = 1.6e11   # 8 um x 13 um CRL beam
saxs1.P = 3.8e12   # CHESS-U Flux, Ph/s

#saxs1.P = 3.8e11   # Oct 2017 30 um x 50 um CRL beam

# don't forget to set the beam dimensions (cm)

#h_beam = 0.005
#v_beam = 0.003

#h_beam = 0.0008
#v_beam = 0.0013

#############################################################################
#############################################################################

saxs1.set_energy(14.14) # this energy is energy for simulated data

# re-calculate visible q-range and "pixels per bin" array (NofPixels)
# RM! 04.26.20, mask parameters must be updated for Eiger4M detector.
saxs1.create_Mask(98, 3, 45, 14, wedge=360.0, type="rectangle")

# need to re-load model so we can re-interpolate onto new default_q
saxs1.load_I(sample_model_1,interpolate=True,q_array = saxs1.default_q)


# RM & REG! Apr2020
##### Rg Error Modeling #####

err_data = ERRORPROP(saxs1 = saxs1)
conc, rgError, log_sig = err_data.calc_errRg_conc()

err_data.plot_S1(conc, [x * 100 for x in rgError],
                 plotlabel = 'Simulated Error - Analytical model', savelabel = 'Simulated_Error_Func_Conc',
                 xlabel = 'Conc. ($\\frac{mg}{ml}$)', ylabel = '($\\frac{\sigma_{R_{g}}}{R_{g}}$) $\cdot 100$')


# Quick calculate model from initial points (slope) of the simulated data
inv_err = 1/np.array(rgError)
final_slope = (inv_err[1]-inv_err[0])/(conc[1]-conc[0])

# Technically this final_slope term should be empirically model as it may not be known apriori

err_data.plot_S1(conc, 1.0/(final_slope*np.array(conc)),
                 plotlabel= '($\\frac{1}{conc}$) Model', savelabel = 'Inv_c_Model',
                 xlabel = 'Conc. ($\\frac{mg}{ml}$)', ylabel = 'Model')

err_data.plot_S2(conc, rgError, 1.0/(final_slope*np.array(conc)),
                 plotlabel1 = 'Simulated Error - Analytical model', plotlabel2 = '($\\frac{1}{final \ slope \cdot conc}$)',
                 savelabel = 'Analytical_and_Inv_c_Rg_ErrorModel',
                 xlabel = 'Conc. ($\\frac{mg}{ml}$)', ylabel = '($\\frac{\sigma_{R_{g}}}{R_{g}}$)')

# RM! 04.28.2020
test = [30376724222.528214, 29959417761.80821, 28367436004.208202, 26530533976.208206, 23327730440.208195, 20586507413.8082, 19343066041.0082]
I = saxs1.I_of_q_variable_contrast(saxs1.c, saxs1.mw, saxs1.buf_model_q, test)
print ('testing function:', I[0])
######################################################

## RM! 04.15.2020 Edits to call PlotSAXS class
#saxsp = PlotSAXS()
#saxsp.plot_Kratky(saxs1, I_no_noise, I_w_noise)
#saxsp.plot_I_w_noise(saxs1, I_no_noise, I_w_noise)
#saxsp.plot_guinier_pure(saxs1, I_no_noise, I_w_noise)
## Does not work with FoxS input ##
# saxsp.plot_shape(saxs1)
#saxsp.plot_PDDR(saxs1)
#saxsp.print_output(saxs1)

######################################################


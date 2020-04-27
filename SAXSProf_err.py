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

err_data = ERRORPROP(saxs1 = saxs1)
conc, rgError, log_sig = err_data.calc_errRg()

err_data.plot_S1(conc, [x * 100 for x in rgError], 'ERROR', 'SAVELABEL')

# plt.plot(err_data.saxs1.buf_model_q, err_data.I_w_noise)
# plt.savefig("TEST.png", format = 'png')
# plt.show()
# print (rgError)
# err_data.plot_S1(conc, rgError, plotlabel='LABEL', savelabel='SAVELABEL')

# # generate buffer profile. simulate_buf uses trimmed mask_q for q-values
# try:
#     (synth_buf, MT_cell, Vac_cell, win_cell, buf_cell) = saxs1.simulate_buf(subtracted=True)
# except ValueError as e: # this essentially says if a ValueError occurs, call it e, extracts the args component
#     print (e.args) # of the instance, and has it printed, and then forces the ValueError to occur.
#     raise e
#
#
# conc = []
# errRg =[]
#
# for c in np.arange(0.05,5.0,0.1):
#     saxs1.c = c
#
#     # calculate synthetic curve on buf_model profile generated by simulate_buf
#     I_no_noise = saxs1.I_of_q(saxs1.c,saxs1.mw,saxs1.buf_model_q)
#
#     # calculate noisy profile on mask_q points that fall within both buf's q range and the specified default_q range (mask_q_short)
#     I_w_noise  = saxs1.t * saxs1.pixel_size ** 2 * saxs1.with_noise(saxs1.t,saxs1.buf_model_q,I_no_noise)
#
#     # calculated smooth I on default_q (goes all the way to q = 0)
#     I_no_noise = saxs1.t*saxs1.pixel_size**2*I_no_noise
#
#     # calculated smooth I and sigma on same q points as I_w_noise
#     I_no_noise_calc = saxs1.Icalc*saxs1.t*saxs1.pixel_size**2  # in register with buffer
#
#     imin = 0
#     imax = 90
#
#     log_sigma = np.abs(saxs1.sigma/I_w_noise)   # converts sigma of the variable to sigma for log of the variable
#
#     (inter,slope,sig2_inter,sig2_slope) = saxs1.lsqfity_with_sigmas(saxs1.buf_model_q[0:imax]**2, np.log(I_w_noise[0:imax]), log_sigma[0:imax])
#
#     Rg = np.sqrt(-3*slope)
#
#     sig_Rg = np.abs(Rg*np.sqrt(sig2_slope)/(2*slope))    # converts the sigma^2 for the slope to the sigma of Rg
#
#     #err_Rg = np.sqrt(-3*err_slope)
#
#     # print c,Rg,sig_Rg/Rg,Rg*saxs1.buf_model_q[imax]
#
#     conc.append(c)
#     errRg.append(sig_Rg/Rg)
#
# #    errRg.append()
#
#     #ax1.set_xlim(5.0,30.0)
#     #ax1.set_ylim(0.002,0.1)
#
# #    plt.plot(saxs1.buf_model_q[0:imax]**2,np.log(I_no_noise[0:imax]),lw=2,
# #                 label = 'No Noise')
#
# #    plt.plot(saxs1.buf_model_q[0:imax]**2,np.log(I_w_noise[0:imax]),'.',
# #                 label = 'Noise')
#
# #    plt.semilogy(saxs1.buf_model_q[0:imax],I_no_noise[0:imax],lw=2,
# #                 label = 'No Noise')
# #    plt.semilogy(saxs1.buf_model_q[0:imax],I_w_noise[0:imax],'.',
# #                 label = 'Noise')
#
# #ax1.set_yscale("log")
# #ax1.legend(loc="lower left", fontsize=14)
#
# inv_err = 1.0/np.array(errRg)
# plt.rc("axes",linewidth=2)
# #plt.rc("xticks",linewidth=3)
# plt.rc("lines", markeredgewidth=2)
# plt.rc('font', **{"sans-serif":["Arial"]})
# fig = plt.figure(figsize = (8,8))
# ax1 = fig.add_subplot(1,1,1)
# #plt.plot(conc, inv_err,label="SAXSProf")
# plt.plot(conc,errRg,'-o',label='SAXSProf')
# plt.savefig('Inv_Error_vs_Concentration2.png', format='png')
# plt.show()
#
# # fit data to 1/c by calculating scale factor from first two points
#
# final_slope = (inv_err[1]-inv_err[0])/(conc[1]-conc[0])
# plt.plot(conc,1.0/(final_slope*np.array(conc)),label='1/c fit')
# #plt.plot(conc,final_slope*np.array(conc),label='1/c fit')
#
#
# plt.ylabel("1/sig_Rg/Rg",size=22)
# plt.xlabel('$concentration$',size=22)
# # for tick in ax1.xaxis.get_major_ticks():
# #     tick.label1.set_fontsize(20)
# # #    tick.label1.set_fontname('Helvetica')
# # for tick in ax1.yaxis.get_major_ticks():
# #     tick.label1.set_fontsize(20)
# # plt.legend()
# plt.savefig('Inv_Error_vs_Concentration.png', format='png')
#
# plt.show()


# test = errProp()
# error.calc_errRg(imin = 1, imax = 125)



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


####################################
# import numpy as np
from SAXS8 import *
from PlotSAXS import *
from SAXSProf_ErrCalcs import *
from matplotlib import pyplot as plt
# import SASM
#from scipy.stats.distributions import chi2
#from scipy.stats import chisqprob
# from scipy import stats
# from subprocess import check_output, CalledProcessError
# import colorsys
# import sys
import warnings
#############################################################################
#############################################################################

# Do not change these parameters here. They need to be specific values that
# correspond to the experimental buffer, window, and vacuum profiles

###### Suppress RuntimeWarning from calculating PDDR ######## RM!

warnings.filterwarnings("ignore", category=RuntimeWarning)

#############################################################
# Initial User Inputs
#############################################################
energy = 9.962 # energy (keV)
P  = 8.4e11     # flux (photons/s)
t  = 1.0       # exposure time for individual snapshots (s)
snaps =  1      # number of snapshots averaged
a  = 150.47    # sample-to-detector distance (cm)
d  = 0.15       # sample cell path length (cm)

window_type = "mica"  # type of window material on sample cell

sensor_thickness = 0.032  # Pilatus sensor thickness

sample_model_1 = "6lyz.pdb.dat"

c = 5.0 # concentration (mg.ml)
t = 0.4450 # Exposure time (seconds)
mw1 = 14.3 # Molecular Weight (kDa)

saxs1 = SAXS(mw=mw1,a=a,d=d,energy=energy,P=P,t=t,total_t=t*snaps,c=c,shape='FoxS', detector='100K')
saxs1.set_window(window_type)
saxs1.sensor_thickness = sensor_thickness
saxs1.det_eff = saxs1.QE(saxs1.lam, saxs1.sensor_thickness)

#  determines the q-space range based on detector and mask. Default_q starts at q = 0, mask_q starts at q_min
saxs1.create_Mask(98, 3, 45, 14, wedge=360.0, type="rectangle")

#############################################################################
# Here we load in
#############################################################################

saxs1.readStandardProfiles("M_Nov07_")

#############################################################################

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

# set_energy() alone does NOT reset q-range or mask (see below).
saxs1.set_energy(14.14) # this energy is energy for simulated data

# re-calculate visible q-range and "pixels per bin" array (NofPixels)
# RM! Note the mask is energy dependent, because the observable angles are energy dependent.
# Therefore, when the energy is reset, the mask must be recalculated
# This will recalculate the q-range as well.
# RM! 04.26.20, mask parameters must be updated for Eiger4M detector.
saxs1.create_Mask(98, 3, 45, 14, wedge=360.0, type="rectangle")

# need to re-load model so we can re-interpolate onto new default_q
saxs1.load_I(sample_model_1,interpolate=True,q_array = saxs1.default_q)


# RM & REG! Apr2020
##### Rg Error Modeling #####

err_data = err_Calcs(saxs1 = saxs1)
conc, rgError, log_sig = err_data.calc_errRg_conc()

err_data.plot_S1(conc, [x * 100 for x in rgError],
                 plotlabel = 'Simulated Error - Analytical model',
                 savelabel = 'Simulated_Error_Func_Conc',
                 xlabel = 'Conc. ($\\frac{mg}{ml}$)',
                 ylabel = '($\\frac{\sigma_{R_{g}}}{R_{g}}$) $\cdot 100$')

# Quick calculate model from initial points (slope) of the simulated data
inv_err = 1/np.array(rgError)
final_slope = (inv_err[1]-inv_err[0])/(conc[1]-conc[0])

# Technically this final_slope term should be empirically model as it may not be known apriori

err_data.plot_S1(conc, 1.0/(final_slope*np.array(conc)),
                 plotlabel= '($\\frac{1}{conc}$) Model',
                 savelabel = 'Inv_c_Model',
                 xlabel = 'Conc. ($\\frac{mg}{ml}$)',
                 ylabel = 'Model')

err_data.plot_S2(conc, rgError, 1.0/(final_slope*np.array(conc)),
                 plotlabel1 = 'Simulated Error - Analytical model',
                 plotlabel2 = '($\\frac{1}{final \ slope \cdot conc}$)',
                 savelabel = 'Analytical_and_Inv_c_Rg_ErrorModel',
                 xlabel = 'Conc. ($\\frac{mg}{ml}$)',
                 ylabel = '($\\frac{\sigma_{R_{g}}}{R_{g}}$)')

# RM! 04.28.2020
# Contrast values taken from RM script.
density = [0.99707, 1.0015, 1.0184, 1.0379, 1.0719, 1.1010, 1.1142]
pressure = [0, 10, 50, 100, 200, 300, 350]
ps = []
rho = []
for i in range(len(density)):
    ps.append(density[i] * 3.3428E+23)  # conversion factor from g/ml to electrons per cm^3

for i in range(len(ps)):
    rho.append((saxs1.pp - ps[i]) * (2.818 * 10 ** -13))


I = saxs1.I_of_q_variable_contrast(saxs1.c, saxs1.mw, saxs1.buf_model_q, rho)
err_data.plot_S2(saxs1.buf_model_q, I[0], I[6],
                 plotlabel1= '0 MPa',
                 plotlabel2='350 MPa',
                 savelabel='Scattering_Curve_AtMultiplePressures',
                 xlabel='q($\\AA^{-1}$)',
                 ylabel='I(q)')


rho, Rg_error_contrast, sig2_Rg_out = err_data.calc_errRg_contrast()

err_data.plot_S1(rho, [x * 100 for x in Rg_error_contrast],
                 plotlabel= 'Simulated Error',
                 savelabel = 'Sim_err_Rg_func_of_Contrast',
                 xlabel = '$\Delta \\rho (cm^{-2})$',
                 ylabel = '($\\frac{\sigma_{R_{g}}}{R_{g}}$) $\cdot 100$')


model, c, Rg_error_contrast = err_data.rgErr_contrast_model()

err_data.plot_S2(rho, [x * 100 for x in Rg_error_contrast], [x * 100 for x in model],
                 plotlabel1 = 'Simulated Error - Analytical model',
                 plotlabel2 = '$\\frac{%s}{\\rho}$' % "{:.2e}".format(c[0]),
                 savelabel = 'Sim_err_Rg_func_of_Contrast_w_InvRho',
                 xlabel = '$\Delta \\rho (cm^{-2})$',
                 ylabel = '($\\frac{\sigma_{R_{g}}}{R_{g}}$) $\cdot 100$')

CytC_data = np.loadtxt("rgErr_Pressure_CytC.txt",
                       dtype={'names': ('Pressure', 'Rg', 'RgErr', 'RgErrRel', 'RgErrRelPercent'),
           'formats': (np.float, np.float, np.float, np.float, np.float)}, skiprows=2)
rho.pop(2)
popRho = rho
Rg_error_contrast.pop(2)
popRgErr = Rg_error_contrast
exptData = np.ndarray.tolist(CytC_data['RgErrRelPercent'])
exptData.pop(2)

err_data.plot_S2(popRho, [x * 100 for x in popRgErr], exptData,
                 plotlabel1 = 'Simulated Error - Analytical model',
                 plotlabel2 = 'Experimental Cytochrome C Data',
                 savelabel = 'SimandExpt_err_Rg_func_of_Contrast',
                 xlabel = '$\Delta \\rho (cm^{-2})$',
                 ylabel = '($\\frac{\sigma_{R_{g}}}{R_{g}}$) $\cdot 100$')

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


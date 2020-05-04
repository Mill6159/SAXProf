# chisqprob is depricated. Replacement will be stats.distributions.chi2.sf

# derived from Fig3A_molmovdb_time
#
# from Fig3A_stats_energy_scan_detect_comp.py

### Version of code used to test self.sigma_profile_fit
# This code uses vac, window, and buffer subtracted profiles from
# the Nov 7 run to reproduce lysozyme data for the first time
# in a way that you can change sample and window thickness correctly.

## This version has been modified for the atlastin models

# These files below are only used when generating the "standard profiles."

# buffer with empty cell subtracted
# buffer_model = "/S_A_nov07_028_lysbufnorm_flat.dat"
# buffer_exposure = t  # exposure used for buffer model
# buffer_flux = P      # flux used for buffer model
#
# vac_model = "/A_nov07_072_lysbufnorm.dat"
# vac_exposure = t
# vac_flux = P
#
# win_model = "/S_A_nov07_026_lysbufnorm.dat"
# win_exposure = t
# win_flux = P

#############################################################################
#############################################################################
# buffer, vacuum, window, and model profiles read in and interpolated onto default_q
# saxs1.buf_q are the q points of default_q that fall within buf's q range
#
# saxs1.load_buf(buffer_model, t=buffer_exposure, P=buffer_flux, interpolate=True, q_array = saxs1.default_q)
# saxs1.load_vac(vac_model, t=vac_exposure, P=vac_flux, interpolate=True, q_array = saxs1.default_q)
# saxs1.load_win(win_model, t=win_exposure, P=win_flux, interpolate=True, q_array = saxs1.default_q)
# saxs1.load_I(sample_model_1,interpolate=True,q_array = saxs1.default_q)
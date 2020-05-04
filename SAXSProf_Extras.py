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
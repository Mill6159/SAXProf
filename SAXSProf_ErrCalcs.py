import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

class err_Calcs:

    def __init__(self, saxs1 = [], I_no_noise = [], I_w_noise = [], I_no_noise_calc = [], c = [], synth_buf = [],
                 MT_cell = [], Vac_cell = [], win_cell = [], buf_cell = []):
        """
        Description: This initializes the class. It imports several important parameters, particularly from the SAXS() class.
        This class also assumes an object has already been assigned to the SAXS() class. In this case, it is saxs1.

        """
        print('err_Calcs class was called.')
        print "--------------"
        self.saxs1 = saxs1
        # self.c = c
        self.c = saxs1.c
        self.synth_buf = synth_buf
        self.MT_cell = MT_cell

        # generate buffer profile. simulate_buf uses trimmed mask_q for q-values
        try:
            (synth_buf, MT_cell, Vac_cell, win_cell, buf_cell) = saxs1.simulate_buf(subtracted=True)
        except ValueError as e:  # this essentially says if a ValueError occurs, call it e, extracts the args component
            print(e.args)  # of the instance, and has it printed, and then forces the ValueError to occur.
            raise e

        self.synth_buf = synth_buf
        self.MT_cell = MT_cell
        self.Vac_cell = Vac_cell
        self.win_cell = win_cell
        self.buf_cell = buf_cell

        # calculate synthetic curve on buf_model profile generated by simulate_buf
        I_no_noise = saxs1.I_of_q(self.c, saxs1.mw, saxs1.buf_model_q)

        # calculate noisy profile on mask_q points that fall within both buf's q range and the specified default_q range (mask_q_short)
        I_w_noise = saxs1.t * saxs1.pixel_size ** 2 * saxs1.with_noise(saxs1.t, saxs1.buf_model_q, I_no_noise)

        # calculated smooth I on default_q (goes all the way to q = 0)
        I_no_noise = saxs1.t * saxs1.pixel_size ** 2 * I_no_noise

        # calculated smooth I and sigma on same q points as I_w_noise
        I_no_noise_calc = saxs1.Icalc * saxs1.t * saxs1.pixel_size ** 2  # in register with buffer

        self.I_no_noise = I_no_noise
        self.I_w_noise = I_w_noise
        self.I_no_noise_calc = I_no_noise_calc

    def calc_errRg_conc(self, imin = 0, imax = 90):
        saxs1 = self.saxs1
        I_w_noise = self.I_w_noise
        conc = []
        errRg = []
        sig2_Rg_out = []
        for c in np.arange(0.05,5.0,0.1):
            saxs1.c = c
            # calculate synthetic curve on buf_model profile generated by simulate_buf
            I_no_noise = saxs1.I_of_q(saxs1.c, saxs1.mw, saxs1.buf_model_q)
            # calculate noisy profile on mask_q points that fall within both buf's q range and the specified default_q range (mask_q_short)
            I_w_noise = saxs1.t * saxs1.pixel_size ** 2 * saxs1.with_noise(saxs1.t, saxs1.buf_model_q, I_no_noise)

            log_sigma = np.empty(len(saxs1.sigma))
            for i in range(len(saxs1.sigma) - 1):
                log_sigma[i] = np.abs(saxs1.sigma[i] / I_w_noise[i])  # converts sigma of the variable to sigma for log of the variable
                # return log_sigma

            (inter, slope, sig2_inter, sig2_slope) = saxs1.lsqfity_with_sigmas(saxs1.buf_model_q[imin:imax] ** 2, np.log(I_w_noise[imin:imax]),
                                                                               log_sigma[imin:imax])

            Rg = np.sqrt(-3 * slope)
            sig_Rg = np.abs(
                Rg * np.sqrt(sig2_slope) / (2 * slope))  # converts the sigma^2 for the slope to the sigma of Rg

            # err_Rg = np.sqrt(-3*err_slope)
            # print c,Rg,sig_Rg/Rg,Rg*saxs1.buf_model_q[imax]

            ## lists to append to each iteration of the loop.
            sig2_Rg_out.append(sig_Rg)
            conc.append(saxs1.c)
            errRg.append(sig_Rg / Rg)
        return conc, errRg, sig2_Rg_out
        print('Default imin/imax values used: 0 and 90')



    def calc_errRg_contrast(self, imin = [], imax = [], density = [],
                            pressure = []):
        """
        Density must be in g/ml
        pressure must be in MPa (default 0-350 MPa). Not actually a necessary input.
        Vectors must be of equal length.
        """
        if imin == [] and imax == [] and density == [] and pressure == []:
            # Input general parameters #
            saxs1 = self.saxs1
            c = self.saxs1.c # use default conc.
            pp = self.saxs1.pp
            mw = self.saxs1.mw # default molecular weight
            t = self.saxs1.t # default exposure time
            pixel_size = self.saxs1.pixel_size
            q = self.saxs1.buf_model_q
            density = [0.99707, 1.0015, 1.0184, 1.0379, 1.0719, 1.1010, 1.1142]
            pressure = [0, 10, 50, 100, 200, 300, 350]
            # Calculate contrast term based on supplied densities #
            ps = []
            rho = []
            for i in range(len(density)):
                ps.append(density[i] * 3.3428E+23)  # conversion factor from g/ml to electrons per cm^3

            for i in range(len(ps)):
                rho.append((pp - ps[i]) * (2.818 * 10 ** -13))

            ## Now iterate over the error calculation for the SAXS curve at each contrast value
            # This is the 'crucial' component of the function.
            errRg = []
            sig2_Rg_out = []
            if len(ps) == len(rho): # To make sure calculations went well.
                for i in range(len(ps)):
                    imin = 0
                    imax = 90
                # calculate synthetic curve on buf_model profile generated by simulate_buf
                    I_no_noise = []
                    contrast = rho[i]
                    I_no_noise.append(saxs1.I_of_q_input_rho(c, mw, q, contrast))
                # calculate noisy profile on mask_q points that fall within both buf's q
                # calculate noisy profile on mask_q points that fall within both buf's q range and the specified default_q range (mask_q_short)
                    I_w_noise = []
                    I = I_no_noise[0]
                    I_w_noise.append(t * pixel_size ** 2 * saxs1.with_noise(t, q, I))
                    I_w_noise = I_w_noise[0]

                    # with_noise() should recalculate the sigma values
                    # therefore, must call sigma for each iteration of the loop.
                    sigma = saxs1.sigma
                    log_sigma = []
                    for i in range(len(sigma)):
                        log_sigma.append((np.abs(sigma[i]/I_w_noise[i])))

                    (inter, slope, sig2_inter, sig2_slope) = saxs1.lsqfity_with_sigmas(saxs1.buf_model_q[imin:imax] ** 2, np.log(I_w_noise[imin:imax]),
                                                                               np.array(log_sigma[imin:imax]))

                    Rg = np.sqrt(-3 * slope)
                    sig_Rg = np.abs(
                        Rg * np.sqrt(sig2_slope) / (2 * slope))  # converts the sigma^2 for the slope to the sigma of Rg

                # err_Rg = np.sqrt(-3*err_slope)
                # print c,Rg,sig_Rg/Rg,Rg*saxs1.buf_model_q[imax]

                ## lists to append to each iteration of the loop.
                    errRg.append(sig_Rg / Rg)
                    sig2_Rg_out.append(sig_Rg)
                print('Rg relative error as calculated as a function of contrast (i.e pressure) via the calc_errRg_contrast() function')
                return rho, errRg, sig2_Rg_out
            else:
                print('Something went wrong calculating the contrast term.')
        print('Default imin/imax values used: 0 and 90')

    def rgErr_contrast_model(self):
        """
        RM! Currently uses default settings for calc_errRg_contrast() function. I will need to add functionality to allow the user
        to define those inputs as they will severely impact the model
        """
        model = []
        rho, Rg_error_contrast, sig2_Rg_out = self.calc_errRg_contrast()

        def cont_model(rho, constant):
            return (constant / (rho))
        g = [4 * 10 **7]
        # x = [x ** 2 for x in rho]
        c, cov = curve_fit(cont_model, rho, Rg_error_contrast, g)

        model = cont_model(rho, c[0])
        return model, c, Rg_error_contrast


    def plot_S1(self, X, Y, plotlabel = '', savelabel = '', xlabel = '', ylabel = ''):
        if plotlabel == '' and savelabel == '':
            plotlabel = 'No label provided'
            savelabel = 'No_Label_Provided'
            plt.rc("axes", linewidth=2)
            plt.rc("lines", markeredgewidth=2)
            plt.rc('font', **{"sans-serif": ["Helvetica"]})
            fig = plt.figure(figsize=(8, 8))
            ax1 = fig.add_subplot(1, 1, 1)
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
                tick.label1.set_fontname('Helvetica')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
            plt.ylabel('No Y label provided',size=22)
            plt.xlabel('No X label provided',size=22)
            plt.plot(X, Y, '-o', label=plotlabel)
            plt.legend(numpoints=1, fontsize=18, loc="best")
            plt.savefig(savelabel + ".png", format='png',
                        bbox_inches = 'tight')
            plt.show()
        else:
            plt.rc("axes", linewidth=2)
            plt.rc("lines", markeredgewidth=2)
            plt.rc('font', **{"sans-serif": ["Helvetica"]})
            fig = plt.figure(figsize=(8, 8))
            ax1 = fig.add_subplot(1, 1, 1)
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
                tick.label1.set_fontname('Helvetica')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
            plt.ylabel(ylabel,size=22)
            plt.xlabel(xlabel,size=22)
            plt.plot(X, Y, '-o', label=plotlabel)
            plt.legend(numpoints=1, fontsize=18, loc="best")
            plt.savefig(savelabel + ".png", format='png',
                        bbox_inches = 'tight')
            plt.show()

    def plot_S2(self, X, Y1, Y2, plotlabel1 = '', plotlabel2 = '', savelabel = '', xlabel = '', ylabel = ''):
        if plotlabel1 == '' and plotlabel2 == '' and savelabel == '':
            plotlabel1 = 'No label provided'
            plotlabel2 = 'No label provided'
            savelabel = 'No_Label_Provided'
            plt.rc("axes", linewidth=2)
            plt.rc("lines", markeredgewidth=2)
            plt.rc('font', **{"sans-serif": ["Helvetica"]})
            fig = plt.figure(figsize=(8, 8))
            ax1 = fig.add_subplot(1, 1, 1)
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
                tick.label1.set_fontname('Helvetica')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
            plt.ylabel('No Y label provided',size=22)
            plt.xlabel('No X label provided',size=22)
            plt.plot(X, Y1, '-o', label=plotlabel1)
            plt.plot(X, Y2, '-o',
                     color = 'k',
                     label=plotlabel2)
            plt.legend(numpoints=1, fontsize=18, loc="best")
            plt.savefig(savelabel + ".png", format='png',
                        bbox_inches = 'tight')
            plt.show()
        else:
            plt.rc("axes", linewidth=2)
            plt.rc("lines", markeredgewidth=2)
            plt.rc('font', **{"sans-serif": ["Helvetica"]})
            fig = plt.figure(figsize=(8, 8))
            ax1 = fig.add_subplot(1, 1, 1)
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
                tick.label1.set_fontname('Helvetica')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
            plt.ylabel(ylabel,size=22)
            plt.xlabel(xlabel,size=22)
            plt.plot(X, Y1, '-o', label=plotlabel1)
            plt.plot(X, Y2, '--',
                     linewidth = 4,
                     color = 'k',
                     label=plotlabel2)
            plt.legend(numpoints=1, fontsize=18, loc="best")
            plt.savefig(savelabel + ".png", format='png',
                        bbox_inches = 'tight')
            plt.show()

    def tester(self):
        print(self.saxs1.sigma)



#!/usr/bin/env python
import os, StringIO, tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
from matplotlib import cm
import sys
import matplotlib
matplotlib.use('agg') 
from mpl_toolkits.mplot3d import axes3d
import time
from matplotlib.backends.backend_agg import FigureCanvasAgg 
from matplotlib.figure import Figure 
from matplotlib import pyplot as plt
#from matplotlib import rc
from scipy import special
from scipy import integrate
from scipy import polyfit
from scipy import stats
from scipy.interpolate import interpolate
from math import *
import numpy as np
from array import array
import re

FIGX = 6
FIGY = 3

class PlotSAXS:
    def __init__ (self,results_path="",name=""):
        """ Create canvas and other set ups for drawing the graph """
        plt.ioff()     #turn off interactive mode so the graph does not refresh automatically
        self.chi_sqr = -1
        # self.results_path=results_path
        self.results_path = os.getcwd() # RM!
        self.results_url=""
        self.name = name 
        
	#determine current dictory, as well as directory to store plots, relatively
	"""
#       self.app_path = os.path.abspath(os.path.join(__file__, '..'))+"/"
        self.app_path = os.path.dirname(__file__)
	if self.app_path not in sys.path:
	    sys.path.append(self.app_path)
	self.app_url = re.sub(r'.*/wsgi/', '/', self.app_path)
#	self.results_path = re.sub(r'wsgi', 'html', self.app_path)
#	self.results_path+= "results/"
	self.results_path = self.app_path+"/../../../html/macchess/Web_SAXS/results/"
	self.results_url = self.app_url+"/results/"
	os.environ['MPLCONFIGDIR'] = self.app_path+".matplotlib/"
	"""
        
    

    def plot_mask(self,saxs):
        """Make a plot of real and simulated masks
        Important: This function assumes that saxs object has 
        loaded the real mask, 
        created a simulated mask. 
        """
    
        fig = plt.figure(figsize = (FIGX,FIGY))
        # create a q array for the actual mask from experiment
        # retrieve saxs.realMask and plot it
        q = np.linspace(0.0,0.2813,len(saxs.realMask))
        plt.plot(q,saxs.realMask,label = 'Actual mask')
        
        # plot the simulated mask
        plt.plot(saxs.default_q, saxs.NofPixels,'.-',label=saxs.mask_model)
        #print len(saxs.default_q), len(saxs.NofPixels)
            
        plt.gcf().subplots_adjust(bottom=0.15)
        #labels and annotations        
        plt.xlabel('q')
        plt.ylabel('pixels')
        ymin, ymax = plt.ylim()
        plt.axis([0.0,0.3,ymin,ymax],emit=True)
        #plt.axis([saxs.alt_q_min,0.3,ymin,ymax],emit=True)      # emit ?
        plt.legend()
        plt.show()

    def plot_count(self,saxs,q, N_noise,Sig_noise):
        """ Plot the simulated photon counts per pixel on the detector.
        This function doesn't have capability to plot 2 saxs objects."""
       	fig=plt.figure(figsize=(FIGX,FIGY))
        plt.semilogy(q ,N_noise,'r.',markersize = 3)
        
        #=raise the plot, such that label of the x-axis is not cut off
        plt.gcf().subplots_adjust(bottom=0.15)
        #labels and annotations        
        plt.xlabel('q')
        plt.ylabel('Intensity(avg. photon counts/pixel)')
        ymin, ymax = plt.ylim()
        plt.axis([0,0.3,ymin,ymax],emit=True)      # emit ?
        plt.axvline(x=(np.pi/saxs.D_max))
        plt.text(x=(np.pi/saxs.D_max),y=(ymin+ymax)/2,s="Shannon limit")
	#add the second y-axis for the absolute scale
        plt.twinx()
        ymin,ymax = ymin/saxs.det_eff/saxs.I_I0_ratio/saxs.P*(saxs.a**2)/saxs.d, ymax/saxs.t/(0.0172**2)/saxs.det_eff/saxs.I_I0_ratio/saxs.P*(saxs.a**2)/saxs.d
        plt.semilogy([0],[1])
        plt.ylabel('absolute scale cm^-1')
        plt.ylim(ymin,ymax)
        plt.show()

    def plot_bufs(self,saxs,saxs2):
        (Ibuff,Iwind,Ivac) = saxs.simulate_buf()
        Ibuff = Ibuff*(saxs.pixel_size)**2
        fig = plt.figure(figsize = (8,4))
        plt.semilogy(saxs.buf_q, saxs.buf_I, label="actual buffer")
        plt.semilogy(saxs2.exp_q,saxs2.exp_I)
        plt.semilogy(saxs.default_q,Ibuff,label="sim buffer")
        plt.semilogy(saxs.default_q,Iwind,label="windows")
        plt.legend()
        plt.show()

    def plot_I_sig(self,saxs1,I1,I2,saxs2,title="test"):
        """  compare sig values for two profiles """
        fig = plt.figure(figsize = (16,8))
        f1 = fig.add_subplot(2,1,1)

        plt.semilogy(saxs1.short_q,I1,label = "sim")
        plt.semilogy(saxs1.exp_q,I2,label = "no noise")
        plt.semilogy(saxs2.exp_q,saxs2.exp_I, label = "exp")
        plt.ylabel("I")
        plt.legend(loc=3)
        plt.title(title)

        f2 = fig.add_subplot(2,1,2)
        plt.semilogy(saxs1.short_q,I1/saxs1.sigma, label = "sim")
        plt.semilogy(saxs2.exp_q,saxs2.exp_I/saxs2.exp_Sig, label = "exp")
        plt.ylabel("I/sig")
        plt.legend()
        plt.show()

    def plot_I(self,saxs1,I1,I1_noise = None, saxs2 = None, I2 = None, I2_noise = None, noise = False):
        """Plot I, with options to show noise. Default is without noise. 
        Have the capability to plot 2 saxs objects. 
        User must provide saxs2, I2 if only want smooth curve;
        must provide saxs2, I2, I2_noise if want simulated noise"""

        fig = plt.figure(figsize = (FIGX,FIGY))
        plt.semilogy(saxs1.exp_q,I1, label = "Intensity 1")
        if not (saxs2 == None or I2 == None):
            plt.semilogy(saxs2.exp_q,I2, label = "Intensity 2")

        if noise == True:

            if not I1_noise==None:
                plt.semilogy(saxs1.buf_q[:len(I1_noise)],I1_noise,'g.', label="Sim. Noise of I1")
                # Ensure that first saxs object noise is plotted 
                # even in the case where I2_noise isn't provided.
                # Also ensure that second saxs object noise is plotted when provided.
            if not I2_noise == None:
                plt.semilogy(saxs2.buf_q[:len(I2_noise)],I2_noise,'r.',label = "Sim. Noise of I2")
        
        elif noise == False:
            print "plot_I(): You didn't plot noise."
        else:
            print "plot_I(): Check your flags and provided args."
        plt.xlabel('q')
        plt.ylabel('Intensity [photons/cm^2*s]')
        plt.legend()
        plt.show()
        
    def plot_I_w_noise(self,saxs1,I_no_noise,I_w_noise):
        """Plot I, with noise and the pure theoretical smooth curve. """

        fig = plt.figure(figsize = (FIGX,FIGY))
        plt.semilogy(saxs1.buf_model_q,I_no_noise, color = 'r', label = "no noise")
        plt.semilogy(saxs1.mask_q_short,I_w_noise, color = 'b', label = "with noise")

        #plt.axvline(x=(np.pi/saxs1.D_max),color='r',label ="Shannon limit" )

        plt.xlabel('q')
        plt.ylabel('Intensity [photons/cm^2*s]')
        plt.legend(loc=3,prop={'size':8})
        #plt.show()
        plt.savefig("%s/%s_I_w_noise.png"%(self.results_path,'SAXProf'), format='png', dpi=500, bbox_inches='tight')
       
    def plot_Kratky (self,saxs1, I_no_noise, I_w_noise):
        """ Plot a Kratky plot, using the curves with simulated noise.
        """
        fig=plt.figure(figsize=(FIGX,FIGY))
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        ax = fig.add_subplot(1, 1, 1)
        plt.rc("axes", linewidth=2)
        # plt.rc("lines", markeredgewidth=2)
        plt.rc('font', **{"sans-serif": ["Arial"]})
        plt.plot(saxs1.buf_model_q, I_no_noise*saxs1.buf_model_q**2, color='r')
        plt.plot(saxs1.mask_q_short, I_w_noise*saxs1.mask_q_short**2, color='g')
        
        #=raise the plot, such that label of the x-axis is not cut off
        plt.gcf().subplots_adjust(bottom=0.15)
        #labels and annotations        
        plt.xlabel('q ($\\AA^{-1}$)')
        plt.ylabel('Intensity$^{2}$')
        ymin, ymax = plt.ylim()
        plt.axis([0,0.3,ymin,ymax],emit=True) # fix limits RM!
        ax.xaxis.set_tick_params(which='both', width=2)
        ax.yaxis.set_tick_params(which='both', width=2)

        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(12)
            tick.label1.set_fontname('Arial')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(12)
            tick.label1.set_fontname('Arial')
        plt.legend(numpoints=1,fontsize=14, loc="best")
        fig.tight_layout()
        plt.savefig("%s/%s_Kratky.png" % (self.results_path, 'SAXProf'), format='png', dpi=500, bbox_inches='tight')
        plt.show()  # RM !
        # plt.savefig("%s/%s_Kratky.png"%(self.results_path,self.name))
   

    def plot_guinier_pure(self,saxs1,  I_no_noise, I_w_noise):
        """ Plot guinier region of the theoretical intensity curve """
        #rc('text', usetex=True)
        fig = plt.figure(figsize=(FIGX,FIGY))
        ax = fig.add_subplot(2,1,1)
        plt.plot(saxs1.buf_model_q**2, np.log(I_no_noise),color='r')
        plt.plot(saxs1.mask_q_short**2, np.log(I_w_noise),'.',color='g')

        # determine the Rg from the first few points on the noise-free (positive) curve

        k_noise,r_noise, ry, smy, sby = saxs1.lsqfity(saxs1.buf_model_q[0:10]**2,np.log(I_no_noise[0:10]))

        Rg = sqrt(-k_noise*3)
        saxs1.Rg_pure = Rg   # store for output

        #self.Rg_noise = sqrt(-k_noise*3)
        #self.I_0 = np.exp(r_noise)
        #self.sigI_0 = np.exp(sby)

        qmax = 2.0/Rg
        #Imin = np.log(I_no_noise[10])
        Imax = 1.1*np.log(I_no_noise[0])

        # find I(qmax == 2/Rg)
        indx = (np.abs(saxs1.buf_model_q-qmax)).argmin()
        Imin = 0.9*np.log(I_no_noise[indx])
        
        # calculate Guinier values of both datasets on improved interval

        k_noise,r_noise, ry, smy, sby = saxs1.lsqfity(saxs1.buf_model_q[0:indx]**2,np.log(I_no_noise[0:indx]))

        Rg = sqrt(-k_noise*3)

        indx_n = (np.abs(saxs1.mask_q_short-qmax)).argmin()

        k_noise,r_noise, ry, smy, sby = saxs1.lsqfity(saxs1.mask_q_short[0:indx_n]**2,np.log(I_w_noise[0:indx_n]))

        Rg_noise = sqrt(-k_noise*3)
        saxs1.Rg_noise = Rg_noise   # store for output

        # construct line from noise fit

        x_fit = saxs1.mask_q_short**2
        y_fit = x_fit*k_noise+r_noise

        plt.plot(x_fit, y_fit,'--',color='k')

        plt.gcf().subplots_adjust(bottom=0.15)
        # Adjust the axes
        ymin, ymax = plt.ylim()
        plt.xlim(0, qmax*qmax)
        plt.ylim(Imin,Imax)
        #plt.xlabel('q**2')
        ax.set_xticklabels([])
        plt.ylabel('ln(I)')
        plt.legend()

        ax1 = fig.add_subplot(2,1,2)
        plt.axhline(color='k')
        plt.plot(saxs1.mask_q_short**2, I_w_noise-np.exp(y_fit),'.')
        plt.xlim(0, qmax*qmax)
        plt.ylabel("$\Delta$")
        plt.xlabel('$q^{2}$')

        plt.gcf().subplots_adjust(bottom=0.15)
        plt.gcf().subplots_adjust(left=0.15)
        plt.savefig("%s/%s_guinier_noise.png" % (self.results_path, 'SAXProf'), format='png', dpi=500, bbox_inches='tight')
        plt.show()

    def plot_guinier_noise(self, saxs1, q_noise1, guinier_noise1, saxs2=None, q_noise2=None, guinier_noise2=None, fit=False, x1 = None, y1 = None, x2 = None, y2= None):
        """ Plot 2 guinier region, and capable of 2 fits of the simulated noise curves
        Have capability to plot 2 saxs objects. 
        User  must provide saxs2, q_noise2, guinier_noise2 if want the second graph.
        If want to use fit, user should provide x1, y1 arrays, or x2, y2 arrays, or both sets."""
        fig = plt.figure(figsize=(FIGX,FIGY))
        
        plt.plot(q_noise1**2, np.log(guinier_noise1),'r.',markersize=3)                 #identify the limits of coordinates, for the use of placing texts(labels)
        plt.axvline(x=q_noise1[-1]**2, color = 'g',label = "q*Rg=%.1f"%saxs1.qRg)
        plt.axvline(x=(np.pi/saxs1.D_max)**2,color = 'b', label = "Shannon limit ")  
        xmax = max(1.1*(np.pi/saxs1.D_max)**2, 1.1*q_noise1[-1]**2)
        
        if not (saxs2==None or q_noise2==None or guinier_noise2==None): 
            plt.plot(q_noise2**2, np.log(guinier_noise2),'g.',markersize=3) 
            xmax = max(xmax, 1.1*(np.pi/saxs2.D_max)**2, 1.1*q_noise2[-1]**2)
            plt.axvline(x=q_noise2[-1]**2, color = 'b',label = "Curve 2: q*Rg=%.1f"%saxs2.qRg)

            plt.axvline(x=(np.pi/saxs2.D_max)**2,color = 'm', label = "Shannon limit of curve 2")
        elif saxs2==None and q_noise2==None and guinier_noise2==None:
            print "plot_guinier_noise(): Plotting only one guinier plot."
        else:
            print "plot_guinier_noise(): Check your flags and args."



        plt.gcf().subplots_adjust(bottom=0.15)
        ymin, ymax = plt.ylim()
        plt.xlim(0, xmax)
        xmin, xmax = plt.xlim()        
        plt.xlabel('$q^{2} (\AA^{-2})$')
        plt.ylabel('ln(I)')


        if fit ==True:
            
            if not(x1==None or y1==None):
                plt.plot(x1, y1,color = 'k', label = "Guinier Fit")
 
            if not (x2==None or y2==None):
                plt.plot(x2, y2,color = 'g', label = "Fit 2")

        elif fit == False:
            print "plot_guinier_noise(): You didn't plot the fit."
        else: 
            print "plot_guinier_noise(): Check your flags and args.\nIf fit = True, one must provide (x1 and y1) or (x2 and y2), or both sets of args."

        plt.savefig("%s/%s_guinier_noise.png" % (self.results_path, 'SAXProf'), format='png', dpi=500, bbox_inches='tight')
        plt.legend()
        #plt.show()





    def plot_guinier_rel_err(self,q_noise1,guinier_noise1,guinier_fit1,q_noise2=None,guinier_noise2=None,guinier_fit2=None):
        """ Plot 2 guinier relative error graphs, (data point - fit)/ fit
        Have capability to plot 2 saxs objects. 
        User  must provide q_noise2, guinier_noise2, and guinier_fit2 if want the second graph.
        """
        fig = plt.figure(figsize=(FIGX,FIGY))
        plt.plot(q_noise1,(guinier_noise1-guinier_fit1)/guinier_fit1*100, color = 'r', label= "Rel. err. of Guinier fit")
        # Only plot the second graph when all data sets required are in place
        if not (q_noise2==None or guinier_noise2==None or guinier_fit2==None):
            plt.plot(q_noise2,(guinier_noise2-guinier_fit2)/guinier_fit2*100, color = 'g', label = " Rel. err. of curve 2")
        elif q_noise2==None and guinier_noise2==None and guinier_fit2==None:
            print "plot_guinier_rel_err(): Plotting only one relative error curve."
        plt.axhline(y=0)
        plt.xlabel('q')
        plt.ylabel('(I_noise - I_Guinier)/I_Guinier(%)')
        plt.legend()
        #plt.show()
        plt.savefig("%s/%s_guinier_rel_err.png"%(self.results_path,self.name))

    def plot_shape(self,saxs):
        """plot a model of the molecule in indicate its shape,
        Cannot plot two saxs objects"""
        #reference:http://msemac.redwoods.edu/~darnold/html/cylinders.html
        if saxs.shape=="Cylin":
            fig=plt.figure(figsize=(3,3))
            ax=axes3d.Axes3D(fig,azim=30,elev=30)
            saxsR=saxs.R
            saxsH=saxs.H
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = saxsR * np.outer(np.cos(u), np.sin(v))
            y = saxsR * np.outer(np.sin(u), np.sin(v))
            z = saxsH/2 * np.outer(np.ones(100), np.ones(100))
            ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')
            x=np.linspace(-saxsR,saxsR,100)
            z=np.linspace(-saxsH/2,saxsH/2,100)
            X, Z=np.meshgrid(x,z)
            Y=np.sqrt(saxsR**2-X**2)
            ax.plot_surface(X,Y,Z, rstride=5, cstride=5, color='b')
            ax.plot_surface(X,-Y,Z, rstride=5, cstride=5, color='b')
            ax.set_xlabel("$\AA$")
            ax.set_ylabel("$\AA$")
            ax.set_zlabel("$\AA$")
            for tl in ax.w_xaxis.get_ticklabels(): 
                tl.set_size(6)
            for tl in ax.w_yaxis.get_ticklabels(): 
                tl.set_size(6) 
            for tl in ax.w_zaxis.get_ticklabels(): 
                tl.set_size(6)   
                #reference:http://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
                # Create cubic bounding box to simulate equal aspect ratio
            max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
            Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
            Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
            Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
            # Comment or uncomment following both lines to test the fake bounding box:
            for xb, yb, zb in zip(Xb, Yb, Zb):
                ax.plot([xb], [yb], [zb], 'w')
            plt.show()
        elif saxs.shape=="Sphere":
            fig = plt.figure(figsize=(3,3))
            ax=axes3d.Axes3D(fig,azim=30,elev=30)
            saxsR=saxs.R
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)

            x = saxsR * np.outer(np.cos(u), np.sin(v))
            y = saxsR * np.outer(np.sin(u), np.sin(v))
            z = saxsR * np.outer(np.ones(np.size(u)), np.cos(v))
            ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')
            ax.set_xlabel("$\AA$")
            ax.set_ylabel("$\AA$")
            ax.set_zlabel("$\AA$")
            for tl in ax.w_xaxis.get_ticklabels(): 
                tl.set_size(6)
            for tl in ax.w_yaxis.get_ticklabels(): 
                tl.set_size(6) 
            for tl in ax.w_zaxis.get_ticklabels(): 
                tl.set_size(6)   
                plt.savefig("%s/MoleculeShape.png"%self.results_path)
            plt.show()
        elif saxs.shape == "FoxS":
            print "plot_shape(): Unable to plot shape. Input is a FoxS theoretical curve.\nPlease use other program to determine shape."

        else:
            print "plot_shape(): Unable to plot shape."

       
    def plot_PDDR(self, saxs1,saxs2=None):
        """Plot Pair distance distribution function (capable of plotting 2 saxs objects)
        Must provide saxs2 if want the second graph."""
        fig=plt.figure(figsize=(FIGX,FIGY))
        r1 = np.arange(0,saxs1.D_max*1.4,saxs1.D_max/50.0)
        PDDF1 = saxs1.PDDF(r1)
        #print PDDF1
        PDDF1[0] = 0
        #print PDDF1
        plt.plot(r1,PDDF1, color ='#FB900F', label = 'PDDF')
        if saxs2 != None:
            r2 = np.arange(0,saxs2.D_max*1.2,saxs2.D_max/50.0)
            PDDF2 = saxs2.PDDF(r2)
            #print PDDF2
        #set P(r=0) value to zero manually, as the algorithm often returns nan for this value
            PDDF2[0] = 0
            #print PDDF2
            plt.plot(r2,PDDF2, color ='g', label = 'PDDF of Curve 2')
        else:
            print "plot_PDDR(): Plotting only one PDDF plot."
        
        plt.legend(numpoints=1,fontsize=14, loc="best")
        plt.savefig("%s/%s_PDDR.png" % (self.results_path, 'SAXProf'), format='png', dpi=500, bbox_inches='tight')
        #plt.show()



    def get_chi_sqr(self, theo_saxs, exp_saxs):
        """Given 2 curves, calculate the chi_sqr
        Two saxs object must share the same exp_q array. It's used with loaded data,
        Not with theoretical curve for sphere or cylindrical shapes
        
        """
        N1, Sig1 = theo_saxs.get_count_noise(theo_saxs.exp_q)
        N2, Sig2 = exp_saxs.get_count_noise(exp_saxs.exp_q)
 
        totalN = len(N2)

        self.chi_sqr = np.sum((N1-N2)**2/(Sig1**2+Sig2**2))/totalN

        #print N1[100],N2[100],Sig1[100],Sig2[100]

        return self.chi_sqr

    def plot_I_sig_ratio(self, saxs1, saxs2):
        """Given 2 curves, calculate the signal to noise ratio. 
        Two saxs object must share the same exp_q and buf_q array. It's used with loaded data,
        Not with theoretical curve for sphere or cylindrical shapes"""
        fig=plt.figure(figsize=(FIGX,FIGY))
        N1, Sig1 = saxs1.get_count_noise(saxs1.exp_q)
        N2, Sig2 = saxs2.get_count_noise(saxs2.exp_q)
        ItoSig = (N1-N2)/np.sqrt(N1+N2)
        
        plt.plot(saxs1.buf_q, ItoSig, 'g.', label = "Signal to Noise Ratio")
        plt.xlabel('q')
        plt.ylabel('I/sigma')
        plt.axhline(y=2.0)
        plt.legend()
        plt.show()
        
    
        
    
    def print_output(self,saxs1,saxs2 = None):
        """ Print the text output of the simulation.
        Good for stand alone mode, without HTML formatting
        """
        out = "\n"
        out+= "Simulation parameters and results:\n"
        out+= "Curve 1:\n"
        if saxs1.shape=='Sphere':
            out+="Radius of molecule = %.2f A\n>" % (saxs1.R)
        elif saxs1.shape == 'Cylin':
            out+="Height H = %.2f A\n" % saxs1.H
            out+="Radius R = %.2f A\n" % saxs1.R
            
        elif saxs1.shape == 'FoxS':
            out += "Shape is encoded in FoxS profile.\n"
        out+="Concentrantion = %.2f mg/ml\n" % saxs1.c
        out+="Molecular weight = %.2f kDa\n" % saxs1.mw
        if saxs1.Rg_noise == -1:
            out+= "print_output(): Sim. Rg and I_0 not calculated. \nRun SAXS.cal_guinier first.\n"
        else:
            out+="Sim. Rg (with noise)= %.2f +- %.4f A\n" % (saxs1.Rg_noise,saxs1.rms_error_noise)
            out+="Sim. Rg (no noise)=%.2f +- %.4f A\n" % (saxs1.Rg_pure,saxs1.rms_error_pure)
            # out+="Ideal Radius of gyration (directly from protein size&shape)=%.2f A\n" % (self.qRg/self.q_fit(self.R))
            #out+="q_min*Rg=%.2f\n"% (self.q_min*self.qRg/self.q_fit(self.R))
            out+="I(0)=%.2f photons/pixel\n" % (saxs1.I_0)
        
        if self.chi_sqr != -1:
            out+= "Chi square = %.9f\n" %(self.chi_sqr)
        else:
            out+= "Chi square not calculated."
            
        if saxs2!=None:
            out += "Curve 2: \n"
            if saxs2.shape=='Sphere':
                out+="Radius of molecule = %.2f A\n>" % (saxs2.R)
            elif saxs2.shape == 'Cylin':
                out+="Height H = %.2f A\n" % saxs2.H
                out+="Radius R = %.2f A\n" % saxs2.R
                
            elif saxs2.shape == 'FoxS':
                out += "Shape is encoded in FoxS profile.\n"
            out+="Concentrantion = %.2f mg/ml\n" % saxs2.c
            out+="Molecular weight = %.2f kDa\n" % saxs2.mw
            if saxs2.Rg_noise == -1:
                out+= "print_output(): Sim. Rg and I_0 not calculated. \nRun SAXS.cal_guinier first.\n"
            else:
                out+="Sim. Rg (with noise)= %.2f +- %.4f A\n" % (saxs2.Rg_noise,saxs2.rms_error_noise)
                out+="Sim. Rg (no noise)=%.2f +- %.4f A\n" % (saxs2.Rg_pure,saxs2.rms_error_pure)
                # out+="Ideal Radius of gyration (directly from protein size&shape)=%.2f A\n" % (self.qRg/self.q_fit(self.R))
            #out+="q_min*Rg=%.2f\n"% (self.q_min*self.qRg/self.q_fit(self.R))
                out+="I(0)=%.2f photons/pixel\n" % (saxs2.I_0)
        
            if self.chi_sqr != -1:
                out+= "Chi square = %.9f\n" %(self.chi_sqr)
            else:
                out+= "Chi square not calculated."
        print out 

    def get_output(self,saxs1, saxs2=None):
        #Need to be modified#
        """Get the text output of the simulation, in HTML formatting."""
        if saxs2==None:
            out="Simulation parameters and results:<br>"
            out+= r"<table border='1' width='600'><tr><th>Parameter</th><th>Value</th></tr>"
            if saxs1.shape=='Sphere':
                out+="<tr><th>Size</th><th>Radius = %.2f &#197</th></tr>" % (saxs1.R)
            elif saxs1.shape == 'Cylin':
                out+="<tr><th>Size</th><th>Height H = %.2f &#197<br>" % saxs1.H
                out+="Radius R = %.2f &#197</th></tr>" % saxs1.R
            ############   HOW DO WE FORMAT THIS???   #############
            elif saxs1.shape == 'FoxS':
                out += "Shape is encoded in FoxS profile.\n"
        

            out+="<tr><th>Concentrantion</th><th>%.2f mg/ml</th></tr>" % saxs1.c
            out+="<tr><th>Molecular weight</th><th>%.2f kDa</th></tr>" % saxs1.mw
            #out+="<tr><th>Sim. Rg (with noise)</th><th>%.2f &#177 %.4f &#197</th></tr>" % (saxs1.Rg_noise,saxs1.rms_error_noise)
            #out+="<tr><th>Sim. Rg (no noise)</th><th>%.2f &#177 %.4f &#197</th></tr>" % (saxs1.Rg_pure,saxs1.rms_error_pure)

            out+="<tr><th>Sim. Rg (with noise)</th><th>%.2f &#197</th></tr>" % (saxs1.Rg_noise)
            #out+="<tr><th>Dose</th><th>%.0f Gy</tr>" % (saxs1.dose)


            #out+="<tr><th>Sim. Rg (no noise)</th><th>%.2f &#197</th></tr>" % (saxs1.Rg_pure)

            #out+="<tr><th>Ideal Rg</th><th>%.2f &#197</th></tr>" % (saxs1.qRg/saxs1.q_fit(saxs1.R))
            #out+="<tr><th>q_min*Rg</th><th>%.2f</th></tr>"% (saxs1.q_min*saxs1.qRg/saxs1.q_fit(saxs1.R))
            #out+="<tr><th>I(0)</th><th>%.2f photons/pixel</th></tr>" % (saxs1.I_0)
            out+="</table>"

        else:
            out="Simulation parameters and results:<br>"
            out+= r"<table border='1'><tr><th>Parameter</th><th>Value1</th><th>Value2</th</tr>"
            out+= "<tr><th>Shape</th><th>%s </th><th>%s </th></tr>" % (saxs1.shape, saxs2.shape)
            if saxs1.shape == 'FoxS':
                out += "Shape is encoded in FoxS profile.\n"
            else:

                out+="<tr><th>Size</th><th>Radius1 = %.2f &#197</th><th>Radius2 = %.2f &#197</th></tr>" % (saxs1.R, saxs2.R)
                if saxs1.shape == 'Cylin' and saxs2.shape == 'Cylin':
                    out+="<tr><th>Size</th><th>Height1 = %.2f &#197</th><th>Height2 = %.2f &#197</th></tr>" % (saxs1.H, saxs1.H)
                elif saxs1.shape == "Cylin" and saxs2.shape != "Cylin":
                    out+="<tr><th>Size</th><th>Height1 = %.2f &#197</th><th>Height2 = none</th></tr>" % (saxs1.H)
                elif saxs1.shape != "Cylin" and saxs2.shape == "Cylin":
                    out+="<tr><th>Size</th><th>Height1 = none</th><th>Height2 = %.2f &#197</th></tr>" % (saxs2.H)
 
        

            out+="<tr><th>Concentrantion</th><th>%.2f mg/ml</th><th>%.2f mg/ml</th></tr>" % (saxs1.c, saxs2.c)
            out+="<tr><th>Molecular weight</th><th>%.2f kDa</th><th>%.2f kDa</th></tr>" % (saxs1.mw,saxs2.mw)
            out+="<tr><th>Sim. Rg (with noise)</th><th>%.2f &#177 %.4f &#197</th><th>%.2f &#177 %.4f &#197</th></tr>" % (saxs1.Rg_noise,saxs1.rms_error_noise, saxs2.Rg_noise, saxs2.rms_error_noise)
            out+="<tr><th>Sim. Rg (no noise)</th><th>%.2f &#177 %.4f &#197</th><th>%.2f &#177 %.4f &#197</th></tr>" % (saxs1.Rg_pure,saxs1.rms_error_pure, saxs2.Rg_pure, saxs2.rms_error_pure)
            out+="<tr><th>Ideal Rg</th><th>%.2f &#197</th><th>%.2f &#197</th></tr>" % (saxs1.qRg/saxs1.q_fit(saxs1.R), saxs2.qRg/saxs2.q_fit(saxs2.R))
            out+="<tr><th>q_min*Rg</th><th>%.2f</th><th>%.2f</th></tr>"% (saxs1.q_min*saxs1.qRg/saxs1.q_fit(saxs1.R),saxs2.q_min*saxs2.qRg/saxs2.q_fit(saxs2.R))
            out+="<tr><th>I(0)</th><th>%.2f photons/pixel</th><th>%.2f photons/pixel</th></tr>" % (saxs1.I_0, saxs2.I_0)
            out+="<tr><th>Chi Square</th><th>%.10f</th><th></th></tr>" %(self.chi_sqr)
            out+="</table>"
    
        return out 

   

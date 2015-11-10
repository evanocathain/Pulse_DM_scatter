#!/usr/bin/python

## Import various important packages
from sympy import *
import numpy as np
import scipy as sp
from scipy.stats import skew
#import matplotlib.pylab as plt
import argparse
import sys

## Parse command line arguments & set default values
parser = argparse.ArgumentParser()
parser.add_argument('-ftop', type=float, dest='ftop', help='set the top frequency in GHz (default: 1.500)', default=1.500)
parser.add_argument('-df', type=float, dest='df', help='set the channel bandwidth in MHz (default: 1.000)', default=1.000)
parser.add_argument('-nchans', type=int, dest='nchans', help='set the number of frequency channels (default: 100)', default=100)
parser.add_argument('-tmax', type=float, dest='tmax', help='set the time integral upper limit in ms (default: 100.0)', default=100.0)
parser.add_argument('-tsamp', type=float, dest='tsamp', help='set the output signal time sampling in us (default: 10.0)', default=10.0)
parser.add_argument('-nterms', type=int, dest='nterms', help='set the number of terms in intergal sum approximation (default: 100)', default=100)
parser.add_argument('-sigma', type=float, dest='sigma', help='set the width of the Gaussian profile in ms (default: 1.0)', default=1.0)
parser.add_argument('-dm', type=float, dest='dm', help='set the DM (default: 10.0)', default=10.0)
parser.add_argument('-tau0', type=float, dest='tau0', help='set the scattering width at fref in ms (default: 1.0)', default=1.0)
parser.add_argument('-fref', type=float, dest='fref', help='set the scattering fref in GHz (default: 1.500)', default=1.500)
parser.add_argument('-scat_idx', type=float, dest='scat_idx', help='set the scattering power-law index (default: -4.0)', default=-4.0)
parser.add_argument('-output', dest='output', help='set the output to file OR screen (default: file)', default='file')
parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')
args = parser.parse_args()
## Set some values - system
ftop=args.ftop              # in GHz
df = args.df*0.001          # in GHz
Nchans=args.nchans
fbot=ftop-Nchans*df         # in GHz
## Set some values - integration
tmax=args.tmax              # in ms
tsamp=args.tsamp*0.001      # in ms
nsamps=tmax/tsamp
nterms=args.nterms          # Need 100 terms (slow) for tmax of 100.0, for 10 terms can't see exp tail.
## Set some values - pulse & ISM
sigma=args.sigma            # in ms
dm = args.dm                # in cm^-3*pc
offset=-4.150*dm*ftop**-2   # in ms
tau0=args.tau0              # in ms
fref=args.fref              # in GHz
scat_idx=args.scat_idx
## Set some values - output
output=args.output
if output == 'screen':
    print "Outputting to screen using matplotlib"
    import matplotlib.pylab as plt
if output == 'file':
    print "Outputting to file"
    outf = file("outf", "w")
    np.set_printoptions(threshold=np.nan)

## Define some sympy symbols, t is in ms, f is in GHz
t, f, tt, ff = symbols('t f tt ff')

## Define some functions - a Gaussian profile and an exponential scattering function
gauss_prof = Lambda((t,f), exp(-(t-offset-(4.15*dm*f**-2))**2/(2*sigma**2)))
exp_scat = Lambda((t,f), Heaviside(t)*exp(-(t/tau0)*(fref/f)**scat_idx))
#pretty_print(gauss_prof)
#pretty_print(exp_scat)
## Create the signal which arrives at the telescope
convolve = Integral(gauss_prof(tt-t,f)*exp_scat(t,f), (t,0,tmax)).as_sum(nterms,method="midpoint")

# Create the channelised signal
for i in range(0,Nchans):
    chan = Integral(convolve, (f,ftop-(i+1)*df,ftop-i*df)).as_sum(nterms,method="midpoint")    
    chan_np = lambdify(tt,chan,"numpy")
    a = sp.linspace(-tmax,tmax,nsamps+1)
    ar = chan_np(a)
    if (i==0):
        area = ar.sum()
    ar = (ar/ar.sum())*area  # Make sure that the area under the curve is conserved
#    print ftop-(i+0.5)*df, ar.sum(), tau0*((ftop-(i+0.5)*df)/fref)**(-scat_idx), skew(ar)
    if (output == "screen"):
        plt.plot(ar,'black')
        plt.xlabel('Time')
        plt.ylabel('Amplitude')
    gauss_prof_assume = Lambda(t, exp(-(t-offset-4.15*dm*(ftop-(i+0.5)*df)**-2)**2/(2*sigma**2)))
    exp_scat_assume = Lambda(t, Heaviside(t)*exp(-(t/tau0)*(fref/(ftop-(i+0.5)*df))**scat_idx))
    chan_assume = Integral(gauss_prof_assume(tt-t)*exp_scat_assume(t), (t,0,tmax)).as_sum(nterms,method="midpoint")
    chan_assume_np = lambdify(tt,chan_assume,"numpy")
    br = chan_assume_np(a)
    br = (br/br.sum())*area
    if (output == "screen"):
        plt.plot(br,'blue')
        plt.plot(br-ar,'green')
    if (output == "file"):
        out_arr = np.append(np.transpose(ar[np.newaxis]),np.transpose(br[np.newaxis]),axis=1) # appended array with real function, and assumed model
        outf.write("%s %s\n" % (ftop-(i+0.5)*df, tsamp*1000.0))
        outf.write("%s" % (out_arr))
if (output == "screen"):
    plt.show()

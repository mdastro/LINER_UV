# Libraries
library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)

# Libraries and other fitting parameters
emiles   <- data("EMILES")
dale     <- data("Dale_NormTot")
wl_pivot <- data("pivwave")

# Selecting the filters that I will use
filters = c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 
            'z_SDSS', 
            #'Z_VISTA', 
            'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 
            'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE', 'P100_Herschel', 
            'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')

filtout=foreach(i = filters)%do%{approxfun(getfilt(i))}
print(filtout)

temppiv=pivwave[pivwave$filter %in% filters,]

# Setting fit parameters
agemax=13.3e9-cosdistTravelTime(z=redshift, H0 = 67.8, OmegaM = 0.308)*1e9

print(agemax/10e9)

inpar=c(mSFR = 0,            #log-space
        mpeak = 0.7,         #log-space
        mperiod = 0.3,       #log-space
        mskew = 0.3,
        tau_birth = 0,       #log-space
        tau_screen = -0.5,   #log-space
        alpha_SF_birth = 1,
        alpha_SF_screen = 3
)

# Let's plot the SFH with the following Prospect plot method
plotSFH=function(par, agemax=13.3, add=FALSE, col='black', ylim=NULL,...){
  magcurve(massfunc_snorm_trunc(age=x, 
                                mSFR=10^par[1], 
                                mpeak=10^par[2],
                                mperiod=10^par[3],
                                mskew=par[4], 
                                magemax=agemax),
           0, 13.8e9, add=add, col=col, ylim=ylim, xlab='Age (Yr)', 
           ylab='SFR (Msol / Yr)',...)
}

plotSFH(inpar)

genSED=ProSpectSED(massfunc=massfunc_snorm_trunc,
                   mSFR=10^inpar[1],
                   mpeak=10^inpar[2],
                   mperiod=10^inpar[3],
                   mskew=inpar[4],
                   tau_birth=10^inpar[5], 
                   tau_screen=10^inpar[6], 
                   alpha_SF_birth=inpar[7], 
                   alpha_SF_screen=inpar[8],
                   z=0.1,
                   Z=Zfunc_massmap_lin,
                   filtout=filtout,
                   Dale=Dale_NormTot,
                   speclib=EMILES,
                   agemax=agemax
)

one_galaxy = 
  read.csv('/home/mlldantas/Projects/LINER_UV/Data/SEDfitTest/one_galaxy.csv')

print(one_galaxy)

flux_input=data.frame(filter=temppiv$filter, pivwave=temppiv$pivwave, 
                      flux=one_galaxy$fluxes, fluxerr=one_galaxy$errors)
print(flux_input)

LumDist_Mpc = cosdistLumDist(z=0.1, H0 = 67.8, OmegaM = 0.308)


Data=list(flux=flux_input,
          arglist=list(z=0.1, massfunc=massfunc_snorm_trunc, agemax=agemax, 
                       Z=Zfunc_massmap_lin, LumDist_Mpc=LumDist_Mpc),
          speclib=EMILES, 
          Dale=Dale_NormTot, 
          filtout=filtout, 
          SFH=SFHfunc, 
          # the preferred functional form of the SFH (eg either SFHfunc, 
          # SFHburst)
          parm.names=c('mSFR','mpeak','mperiod','mskew','tau_birth','tau_screen',
                       'alpha_SF_birth','alpha_SF_screen'), 
          # which parameters to fit for
          logged=c(T,T,T,F,T,T,F,F), # fit parameters in logged or linear space
          intervals=list(lo=c(-4,-2,-1,-0.5,-2.5,-2.5,0,0), 
                         hi=c(3,1,1,1,1.5,1,4,4)), # fitting range for parameters
          fit = 'LD', 
          # specifies the way in which the SED should be fitted ('LD', 'optim', 
          #'CMA', or 'check')
          mon.names=c('LP','masstot','SFRburst',
                      paste('flux.',flux_input$filter,sep='')),
          N=length(filters), # number of observed filters
          like='norm',
          verbose=FALSE
)

set.seed(1)
LDout = LaplacesDemon(Model=ProSpectSEDlike, Data=Data,  Initial.Values=inpar,
                      control=list(abstol=0.1), Iterations=1e4, 
                      Algorithm='CHARM', Thinning=1)

set.seed(1)
Data$fit = 'optim'
optout = optim(par=inpar, fn=ProSpectSEDlike, Data=Data, method = 'BFGS')
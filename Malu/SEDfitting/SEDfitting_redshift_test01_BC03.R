# Libraries---------------------------------------------------------------------
library(ProSpect)
library(LaplacesDemon)
library(foreach)
library(celestial)
library(magicaxis)
library(cmaeshpc)
require(reshape2)

# Libraries and other fitting parameters----------------------------------------
bc03     <- data("BC03lr")
dale     <- data("Dale_NormTot")
wl_pivot <- data("pivwave")

# Selecting the filters that I will use-----------------------------------------
filters = c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 
            'z_SDSS', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 
            'W1_WISE' , 'W2_WISE')

filtout=foreach(i = filters)%do%{approxfun(getfilt(i))}
print(filtout)

temppiv=pivwave[pivwave$filter %in% filters,]

# Setting fit parameters--------------------------------------------------------
redshift <- 0.1
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

# Let's plot the SFH with the following Prospect plot method--------------------
## Creating a function for this plot- - - - - - - - - - - - - - - - - - - - - - 
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

## Actually plotting the SFH- - - - - - - - - - - - - - - - - - - - - - - - - - 
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
                   speclib=BC03lr,
                   agemax=agemax
)

# Setting the data--------------------------------------------------------------
one_galaxy = 
  read.csv('/home/mlldantas/Projects/LINER_UV/Data/SEDfitTest/one_galaxy.csv')

new_galaxy = one_galaxy[c(1:13),]

print(new_galaxy)

flux_input=data.frame(filter=temppiv$filter, pivwave=temppiv$pivwave, 
                      flux=new_galaxy$fluxes, fluxerr=new_galaxy$errors)
for(i in temppiv$filter){
  print(i)
}

print(temppiv$filter)
print(flux_input)

LumDist_Mpc = cosdistLumDist(z=0.1, H0 = 67.8, OmegaM = 0.308)


Data_z=list(flux=flux_input,
            arglist=list(massfunc=massfunc_snorm_trunc, agemax=agemax, Z=Zfunc_massmap_lin),
            speclib=BC03lr, 
            Dale=Dale_NormTot, 
            filtout=filtout, 
            SFH=SFHfunc, # the preferred functional form of the SFH 
                         #(eg either SFHfunc, SFHburst, SFHp4, SFHp5)
            parm.names=c('z', 'mSFR','mpeak','mperiod','mskew','tau_birth',
                         'tau_screen', 'alpha_SF_birth','alpha_SF_screen'), 
                                                  # which parameters to fit for
            logged=c(T,T,T,T,F,T,T,F,F), 
                                      # fit parameters in logged or linear space
            intervals=list(lo=c(-3, -4,-2,-1,-0.5,-2.5,-2.5,0,0), 
                           hi=c(0, 3,1,1,1,1.5,1,4,4)), 
                           # fitting range for parameters
            fit = 'CMA', # specifies the way in which the SED should be fitted 
                         #('LD', 'optim', 'CMA', or 'check')
            mon.names=c('LP','masstot','SFRburst',
                        paste('flux.',flux_input$filter,sep='')),
            N=length(filters), # number of observed filters
            like='norm',
            verbose=FALSE
)

badpar_z = (Data_z$intervals$lo + Data_z$intervals$hi) / 2 
#CMA is pretty tolerant of terrible initial guesses, unlike optim and LD.

CMAout_z = cmaeshpc(par=badpar, fn=ProSpectSEDlike, Data=Data_z, 
                    lower=Data_z$intervals$lo,
                    upper=Data_z$intervals$hi, 
                    control=list(trace=TRUE, maxwalltime=2))

# Running the fit!
set.seed(1)
LDout = LaplacesDemon(Model=ProSpectSEDlike, Data=Data,  Initial.Values=inpar,
                      control=list(abstol=0.1), Iterations=1e4, 
                      Algorithm='CHARM', Thinning=1)

set.seed(1)
Data$fit = 'optim'
optout = optim(par=inpar, fn=ProSpectSEDlike, Data=Data, method = 'BFGS')

set.seed(1)
Data$fit = 'CMA'
#CMA is pretty tolerant of terrible initial guesses, unlike optim and LD.
badpar = (Data$intervals$lo + Data$intervals$hi) / 2 

print(CMAout$par)

# Exploring the stellar mass "formed" during the SFH----------------------------
maghist(LDout$Monitor[,"masstot"], verbose = FALSE, xlab='Stellar Mass / Msol', 
        ylab='PDF')
abline(v=genSED$Stars$masstot, col='red')

# Plotting the fitting results--------------------------------------------------
## Fluxes generated by MCMC sampler---------------------------------------------
magplot(flux_input$pivwave, LDout$Monitor[1,4:16], type='l', log='xy', grid=TRUE,
        xlab="Wavelength (Ang)", ylab='Flux Density / Jy' 
        #xlim=c(1e3,1e5), ylim=c(1e-7,1e-3)
) #LDout$Monitor gives the chosen bands results (if I understood correctly)

for(i in 2:1e4){
  lines(flux_input$pivwave, LDout$Monitor[i,4:16], col=hsv(alpha=0.01))
} # plotting all the 10^4 results!!

points(flux_input[,c("pivwave","flux")])
magerr(flux_input$pivwave, flux_input$flux, ylo=flux_input$fluxerr)

print(LDout$Summary1)
print(LDout$Summary2)
print(LDout$Summary2['mSFR'])
print(LDout$Monitor[1,]['LP'])
print(LDout$Monitor)


maghist(LDout$Monitor[,"masstot"], verbose = FALSE, xlab='Stellar Mass / Msol', ylab='PDF')
abline(v=genSED$Stars$masstot, col='red')



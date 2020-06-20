# We will use examine the effects that Horndeski models of 
# modified gravity/dark energy have on the growth rate of 
# structure fsigma8(z) from redshift z=6 to z=0.

# Parameterization of background: 'lcdm'

# Parameterization of gravity: 'hill_scale'

import numpy as np
import matplotlib.pyplot as plt
#To use LaTeX and select Helvetica as the default font:
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
from classy import Class


def D(z):
    """ 
    Linear growth function
        D(z) = ( P(k_ls, z) / P(k_ls, 0) )^(1 / 2)
    where k_ls is a large-scale mode.
    """
    k=0.01 #Our choice of large-scale mode
    mPk=cosmo.pk(k,z)
    mPk_norm=cosmo.pk(k,0) #Normalize at z=0
    D=np.sqrt(mPk/mPk_norm)
    return D

def f(z):
    """
    Linear growth rate f(z) where
        f(a) = d log D / d log a
    where a = 1 / (1 + z) is the scale factor.
    """
    a=1./(1.+z)
    #da=0.01
    da=0.01*a
    #da=1e-7
    gp,g,gm=[D(1./ia-1.) for ia in [a+da,a,a-da]]
    f=a*(gp-gm)/(2*g*da)
    #dz=0.01
    #gp,g,gm=[D(zi) for zi in [z+dz,z,z-dz]]
    #f=(z)*(gp-gm)/(2.*g*dz)
    return f

def sigma8(z):
    """
    sigma8(z) = sigma8*D(z)
    """
    s8=cosmo.sigma8()*D(z)
    return s8

def fsigma8(z):
    """
    Growth rate of structure
        fsigma8(z) = f(z) * sigma8(z)
    """
    fs8=f(z)*sigma8(z)
    return fs8


def reldev(y1, y2):
    """
    Calculate the relative percent deviation between LCDM 
    and non-LCDM outputs
    """
    return 100.*(y2-y1)/y1



### Our redshift range ###

# We want to see the behavior of fsigma8 from early
# to late times i.e. redshifts 0 to 1000.
# particularly, we want to show that our MG models of
# fsigma8 reduce to LCDM at early times.
zlist = np.linspace(0.011, 1000., 5000.)
#zstr=','.join(map(str,zlist+zlist[-1]+2))
zstr=str(zlist[-1]+2.)

### The LCDM model ###

lcdmpars = {'output': 'mPk',
            'P_k_max_h/Mpc':20,
            'z_pk': zstr,
            'background_verbose': 1, #Info
            'tau_reio': 0.07,
            'omega_cdm': 0.11987,     
            'A_s': 2.204e-9, 
            'h': 0.6715918, #computed from 100*theta_s=1.042143 
            'N_ur': 3.046-1.,  
            'N_ncdm': 1.,           
            'm_ncdm': 0.06,       
            'omega_b': 0.022252,     
            'n_s': 0.96475, 
           }
cosmo = Class() #create universe
cosmo.set(lcdmpars) #feed params to cosmos
cosmo.compute() 

lcdm_fsigma8=np.array([fsigma8(z) for z in zlist])

### Modified Gravity models ###

## Modifying gravity: alpha_T = 0, ^alpha_K = 0.1, ^alpha_M = 0.1, varying ^alpha_B = 0 ##

pars = {'output': 'mPk',
        'P_k_max_h/Mpc':20,
        'z_pk': zstr,
        'background_verbose': 1, #Info
        'tau_reio': 0.07,
        'omega_cdm': 0.11987,     
        'A_s': 2.204e-9,     
        'h': 0.6715918, #computed from 100*theta_s=1.042143 
        'N_ur': 3.046-1.,  
        'N_ncdm': 1.,           
        'm_ncdm': 0.06,       
        'omega_b': 0.022252,     
        'n_s': 0.96475, 
        'Omega_Lambda': 0, #no cosmological constant
        'Omega_fld': 0, #no perfect fluid DE
        'Omega_smg': -1, #Omega_DE -- use closure relation to 
                         #find the density of modified gravity
        'gravity_model': 'hill_scale',
        'expansion_model': 'lcdm' #model for rho_DE(tau)
       }

# Fiducial Horndeski coupling values and M_*^2,ini
M=1.
a_K = 0.001#0.1#1.e-6 #Seems to actually have an effect. See: 0.1, and 1.0 difference
a_B = 0.#0.005
a_M = 0.1 #Note: 1 and greater is insane.
a_T = 0.#1e-30

# Transition scale factor and rapidity
a_trans = 1. / (1. + 7.)
a_rapidity = 4.#3./2.#0.01 

# Varying a_B
a_B_mg_fsigma8={}

a_B_mg_list=[0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16] #a_M = 0.1, zt=7, a_r=4.

for a_B in a_B_mg_list:#, 0.4, 0.6, 0.8]:
    pars['parameters_smg'] = "{}, {}, {}, {}, {}, {}, {}".format(
        M, a_K, a_B, a_M, a_T, a_trans, a_rapidity)
    cosmo=Class()
    cosmo.set(pars)
    cosmo.compute()
    a_B_mg_fsigma8[a_B] = np.array([fsigma8(z) for z in zlist])


## Reproducing LCDM:  ^alpha_T = 0, ^alpha_K = 0.1, ^alpha_M = ^alpha_B = 0 ##

pars = {'output': 'mPk',
        'P_k_max_h/Mpc':20,
        'z_pk': zstr,
        'background_verbose': 1, #Info
        'tau_reio': 0.07,
        'omega_cdm': 0.11987,     
        'A_s': 2.204e-9,     
        'h': 0.6715918, #computed from 100*theta_s=1.042143 
        'N_ur': 3.046-1.,  
        'N_ncdm': 1.,           
        'm_ncdm': 0.06,       
        'omega_b': 0.022252,     
        'n_s': 0.96475, 
        'Omega_Lambda': 0, #no cosmological constant
        'Omega_fld': 0, #no perfect fluid DE
        'Omega_smg': -1, #Omega_DE -- use closure relation to 
                         #find the density of modified gravity
        'gravity_model': 'hill_scale',
        'expansion_model': 'lcdm' #model for rho_DE(tau)
       }

# Fiducial Horndeski coupling values and M_*^2,ini
_M=1.
_a_K = 0.001#1.e-6
_a_B = 0.
_a_M = 0.
_a_T = 0.

pars['parameters_smg'] = "{}, {}, {}, {}, {}, {}, {}".format(
                         _M, _a_K, _a_B, _a_M, _a_T, a_trans, a_rapidity)
cosmo=Class()
cosmo.set(pars)
cosmo.compute()
mg_lcdm_fsigma8 = np.array([fsigma8(z) for z in zlist])


### Plot fsigma8(z) for LCDM and Horndeski gravity ###

f, ax = plt.subplots(2, 1, figsize=(15/1.5, 10/1.5), sharex='col')

# First plot LCDM results

ax[0].semilogx(zlist, lcdm_fsigma8, 
               dashes=[6,4], color='black', 
               zorder=2)
ax[0].set_ylabel(r'$f\sigma_8$')

ax[1].semilogx(zlist, reldev(lcdm_fsigma8, lcdm_fsigma8), 
               dashes=[6,4], color='black', 
               label=r'$\Lambda\mathrm{CDM}$', 
               zorder=2)
ax[1].set_ylabel('$\mathrm{rel. dev. [\%]$}')
ax[1].set_xlabel(r'$z$')

ax[1].set_ylim(-0.25, 0.05) ###

# Now plot the MG/DE results

# Reproducing LCDM

ax[0].plot(zlist, mg_lcdm_fsigma8, 
           zorder=1)
ax[1].plot(zlist, reldev(lcdm_fsigma8, mg_lcdm_fsigma8), 
           label=r'$\hat{{\alpha}}_\mathrm{{B}}={}, \hat{{\alpha}}_\mathrm{{M}}={}$'.format(_a_B, _a_M), 
           zorder=1)
            
# MG
for a_B in a_B_mg_list:
    ax[0].plot(zlist, a_B_mg_fsigma8[a_B], 
               zorder=1) 
        
    ax[1].plot(zlist, reldev(lcdm_fsigma8, a_B_mg_fsigma8[a_B]), 
               label=r'$\hat{{\alpha}}_\mathrm{{B}}={}, \hat{{\alpha}}_\mathrm{{M}}={}$'.format(a_B, a_M), 
               zorder=1)

plt.tight_layout()

# Remove whitespace between upper and lower plots
plt.subplots_adjust(hspace=0, wspace=0.3) 
# Tick marks on all sides of each plot
for j in range(2):
    #for i in range(2):
    axes=ax[j]
    axes.tick_params(axis='both', which='both', direction='in',
                         top=True, right=True)
    legend=axes.legend(framealpha=0)
        
#plt.savefig('fsigma8_Horndeski_hill_scale_braiding_planck_mass_run_rate_z_0_to_1000.pdf', format='pdf')
plt.savefig('fsigma8_Horndeski_hill_scale_braiding_planck_mass_run_rate_z_0_to_1000_close_up.pdf', format='pdf')
plt.show()

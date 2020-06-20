import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

from CMB_lensing_potential import CMB_lensing_potential

CMBlp = CMB_lensing_potential()

llist = np.array(range(2, 2501))

# Make plot
CMBlp.make_plot()

# Plot LCDM Case
CMBlp.plot_Cl_kappa_kappa_GR(llist)

#-------------------------------------------------------------------------- 
# Plot MG models
#--------------------------------------------------------------------------

# Reproduce LCDM
CMBlp.plot_Cl_kappa_kappa_MG(llist, 
                             a_B = 0., a_M = 0., a_trans = 1./(1.+7.), r = 4.)

# Some simple models that have desired effects at z = 2 to 6:
CMBlp.plot_Cl_kappa_kappa_MG(llist, 
                             a_B = 0.1, 
                             a_M = 0.1, 
                             a_trans = 1./(1.+7.), 
                             r = 4.)
CMBlp.plot_Cl_kappa_kappa_MG(llist, 
                             a_B = 0.13, 
                             a_M = 0.1, 
                             a_trans = 1./(1.+7.), 
                             r = 4.)
CMBlp.plot_Cl_kappa_kappa_MG(llist,
                             a_B = 0.15, 
                             a_M = 0.1, 
                             a_trans = 1./(1.+7.), 
                             r = 4.)
CMBlp.plot_Cl_kappa_kappa_MG(llist, 
                             a_B = 0.16, 
                             a_M = 0.1, 
                             a_trans = 1./(1.+7.), 
                             r = 4.)

# Best fit varying alpha_b only. 
# We set alpha_m = 0.1, a_trans = 1./(1.+7.), and r = 4.
# Finds a_B = 0.00130177 which has chisqr = 3.45158072065
# This works -- it has a smaller chisqr than 0.00130176 and 0.00130178
CMBlp.plot_Cl_kappa_kappa_MG(llist, 
                             a_B = 0.00130177, 
                             a_M = 0.1, 
                             a_trans = 1./(1.+7.), 
                             r = 4.)

# Best fit varying alpha_B and a_trans
# We set alpha_m = 0.1 and r = 4.
# Finds a_B = 0.0001971 and a_trans = 0.12738292 
#  which has chisqr = 3.44268933614
# This works -- it has a smaller chisqr than
#  a_B = 0.0001971 and a_trans = 0.12738291
#a_B = 0.0001971 and a_trans = 0.12738293
#a_B = 0.0001970 and a_trans = 0.12738292   
#a_B = 0.0001972 and a_trans = 0.12738292
CMBlp.plot_Cl_kappa_kappa_MG(llist, 
                             a_B = 0.0001971, 
                             a_M = 0.1, 
                             a_trans = 0.12738292, 
                             r = 4.)

"""
# Best fit varying alpha_B, a_trans, and r
# We set alpha_m = 0.1
# Finds a_B = 1.85595399e-03, a_trans = 1.29077787e-01, and r = 3.98879999e+00
#  which has chisqr = 3.43401073595
fsigma8.plot_fsigma8_MG(zlist, fsigma8_datafile, 
                        a_B = 1.85595399e-03, 
                        a_M = 0.1, 
                        a_trans = 1.29077787e-01, 
                        r = 3.98879999)
"""

# Best fit varying alpha_B, a_trans, and r
# We set alpha_m = 0.1
# Finds a_B = 1.70409951e-03, a_trans = 1.28222985e-01, and r = 3.99463473e+00
#  which has chisqr = 
CMBlp.plot_Cl_kappa_kappa_MG(llist, 
                             a_B = 1.70409951e-03, 
                             a_M = 0.1, 
                             a_trans = 1.28222985e-01, 
                             r = 3.99463473)


#--------------------------------------------------------------------------

# Show and save plot
CMBlp.show_plot()

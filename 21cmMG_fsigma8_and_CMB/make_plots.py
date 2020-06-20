import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

from fsigma8_plotter import Fsigma8
from fs8s8ClTTClkk_plotter import Fs8s8ClTTClkk

# Our relevant redshift range is from 0 to 6
zlist = np.linspace(0.011, 6., 100.)

# Let's make our ell list from 2 to 2500:
llist = np.array(range(2, 2501))
#llist = np.linspace(2., 2500., 5000.)

"""
#----------------------------------------------------------------
# fsigma8 plot
#----------------------------------------------------------------
fsigma8 = Fsigma8()

# Datafile
fsigma8_datafile = 'fsigma8_measurements.dat'

# Make plot
fsigma8.make_plot()

# Plot measurement data
fsigma8.plot_fsigma8_data(fsigma8_datafile)

# Plot LCDM Case
fsigma8.plot_fsigma8_GR(zlist, fsigma8_datafile)

# MG model reproducing LCDM
fsigma8.plot_fsigma8_MG(zlist, fsigma8_datafile, 
                         a_B = 0., a_M = 0., a_trans = 1./(1.+7.), r = 4.)

# MG model with desired effects at z = 2 to 6:
fsigma8.plot_fsigma8_MG(zlist, fsigma8_datafile, 
                        a_B = 0.15, 
                        a_M = 0.1, 
                        a_trans = 1./(1.+7.), 
                        r = 4.)

# MG best fit varying alpha_b only. 
# We set alpha_m = 0.1, a_trans = 1./(1.+7.), and r = 4.
# Finds a_B = 0.00130177 which has chisqr = 3.45158072065
# This works -- it has a smaller chisqr than 0.00130176 and 0.00130178
fsigma8.plot_fsigma8_MG(zlist, fsigma8_datafile, 
                        a_B = 0.00130177, 
                        a_M = 0.1, 
                        a_trans = 1./(1.+7.), 
                        r = 4.)

# Show and save plot
fsigma8.show_plot()
"""
#----------------------------------------------------------------
# fsigma8 + sigma8 + C_ell^TT + C_ell^kappakappa plot
#----------------------------------------------------------------
fs8s8ClTTClkk = Fs8s8ClTTClkk()

fs8s8ClTTClkk.make_plot()

# Plot LCDM Case
fs8s8ClTTClkk.plot_fsigma8_GR(zlist)
fs8s8ClTTClkk.plot_sigma8_GR(zlist)
fs8s8ClTTClkk.plot_Cl_TT_GR(zlist, llist)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_GR(zlist, llist)

# MG model reproducing LCDM

fs8s8ClTTClkk.plot_fsigma8_MG(zlist,  
                              a_B = 0., 
                              a_M = 0., 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,  
                              a_B = 0., 
                              a_M = 0., 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                              a_B = 0., 
                              a_M = 0., 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist, 
                                    a_B = 0., 
                                    a_M = 0., 
                                    a_trans = 1./(1.+7.), 
                                    r = 4.)



# MG model with desired effects at z = 2 to 6:
fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = 0.15, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = 0.15, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                            a_B = 0.15, 
                            a_M = 0.1, 
                            a_trans = 1./(1.+7.), 
                            r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                                    a_B = 0.15, 
                                    a_M = 0.1, 
                                    a_trans = 1./(1.+7.), 
                                    r = 4.)


"""
# Additional model
fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = 0.14, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = 0.14, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                            a_B = 0.14, 
                            a_M = 0.1, 
                            a_trans = 1./(1.+7.), 
                            r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                                    a_B = 0.14, 
                                    a_M = 0.1, 
                                    a_trans = 1./(1.+7.), 
                                    r = 4.)

fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = 0.13, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = 0.13, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                            a_B = 0.13, 
                            a_M = 0.1, 
                            a_trans = 1./(1.+7.), 
                            r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                                    a_B = 0.13, 
                                    a_M = 0.1, 
                                    a_trans = 1./(1.+7.), 
                                    r = 4.)
fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = 0.12, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = 0.12, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                            a_B = 0.12, 
                            a_M = 0.1, 
                            a_trans = 1./(1.+7.), 
                            r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                                    a_B = 0.12, 
                                    a_M = 0.1, 
                                    a_trans = 1./(1.+7.), 
                                    r = 4.)
fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = 0.11, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = 0.11, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                            a_B = 0.11, 
                            a_M = 0.1, 
                            a_trans = 1./(1.+7.), 
                            r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                                    a_B = 0.11, 
                                    a_M = 0.1, 
                                    a_trans = 1./(1.+7.), 
                                    r = 4.)

fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = 0.1, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = 0.1, 
                              a_M = 0.1, 
                              a_trans = 1./(1.+7.), 
                              r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                            a_B = 0.1, 
                            a_M = 0.1, 
                            a_trans = 1./(1.+7.), 
                            r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                                    a_B = 0.1, 
                                    a_M = 0.1, 
                                    a_trans = 1./(1.+7.), 
                                    r = 4.)
"""


# No slip models a_B = -2 * a_M
a_M = 0.1
a_B = -2. * a_M
fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = a_B, 
                              a_M = a_M, 
                              a_trans = 1./(1.+7.), 
                              r = 1.5)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = a_B, 
                              a_M = a_M, 
                              a_trans = 1./(1.+7.), 
                              r = 1.5)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                              a_B = a_B, 
                              a_M = a_M, 
                              a_trans = 1./(1.+7.), 
                              r = 1.5)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                              a_B = a_B, 
                              a_M = a_M,  
                              a_trans = 1./(1.+7.), 
                            r = 1.5)

"""
a_M = -0.1
a_B = -2. * a_M
fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                              a_B = a_B, 
                              a_M = a_M, 
                              a_trans = 1./(1.+7.), 
                              r = 1.5)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                              a_B = a_B, 
                              a_M = a_M, 
                              a_trans = 1./(1.+7.), 
                              r = 1.5)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                              a_B = a_B, 
                              a_M = a_M, 
                              a_trans = 1./(1.+7.), 
                              r = 1.5)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                              a_B = a_B, 
                              a_M = a_M,  
                              a_trans = 1./(1.+7.), 
                            r = 1.5)
"""

# MG best fit varying alpha_b only. 
# We set alpha_m = 0.1, a_trans = 1./(1.+7.), and r = 4.
# Finds a_B = 0.00130177 which has chisqr = 3.45158072065
# This works -- it has a smaller chisqr than 0.00130176 and 0.00130178


fs8s8ClTTClkk.plot_fsigma8_MG(zlist,
                        a_B = 0.00130177, 
                        a_M = 0.1, 
                        a_trans = 1./(1.+7.), 
                        r = 4.)
fs8s8ClTTClkk.plot_sigma8_MG(zlist,
                        a_B = 0.00130177, 
                        a_M = 0.1, 
                        a_trans = 1./(1.+7.), 
                        r = 4.)
fs8s8ClTTClkk.plot_Cl_TT_MG(zlist, llist,
                        a_B = 0.00130177, 
                        a_M = 0.1, 
                        a_trans = 1./(1.+7.), 
                        r = 4.)
fs8s8ClTTClkk.plot_Cl_kappa_kappa_MG(zlist, llist,
                        a_B = 0.00130177, 
                        a_M = 0.1, 
                        a_trans = 1./(1.+7.), 
                        r = 4.)


fs8s8ClTTClkk.show_plot()

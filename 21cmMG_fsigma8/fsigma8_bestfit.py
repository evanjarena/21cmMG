import numpy as np
from scipy.optimize import curve_fit

from peturbations import Peturbations

peturb = Peturbations()

"""
Current observational measurements on fsigma8 

We read in current measurements between redshifts 0 and 2
from fsigma8_measurements.dat
Column format: Index_number, z, fs8, fs8_error
"""
data = np.loadtxt('fsigma8_measurements.dat')
meas_fsigma8 = data[:,2]
meas_fsigma8_error = data[:,3]
meas_zlist = data[:,1]
"""
popt, pcov = curve_fit(peturb.fsigma8_MG_fit, meas_zlist, meas_fsigma8, 
                       #p0=[0.1, 1./(1.+7.), 4.],
                       sigma=1./(meas_fsigma8_error*meas_fsigma8_error), 
                       bounds=([0., 1./(1.+10.), 1.5], [0.15, 1./(1.+0.), 4.]))
#perr = np.sqrt(np.diag(pcov))
"""

"""
# Varying alpha_b only. 
# We set alpha_m = 0.1, a_trans = 1./(1.+7.), and r = 4.
# Finds a_B = 0.00130177 which has chisqr = 3.45158072065
# This works -- it has a smaller chisqr than 0.00130176 and 0.00130178
popt, pcov = curve_fit(peturb.fsigma8_MG_fit, meas_zlist, meas_fsigma8, 
                       p0=[0.],
                       sigma=meas_fsigma8_error, 
                       bounds=([-0.01], [0.15]))
"""
"""
# Varying alpha_b and a_trans
# We set alpha_m = 0.1 and r = 4.
# Finds a_b = 0.0001971 and a_trans = 0.12738292 
#  which has chisqr = 3.44268933614
# 
popt, pcov = curve_fit(peturb.fsigma8_MG_fit, meas_zlist, meas_fsigma8, 
                       p0=[0., 1./(1.+7.)],
                       sigma=meas_fsigma8_error, 
                       bounds=([-0.01, 1./(1.+10.)], [0.15, 1./(1.+0.)]))
"""

# Varying alpha_b, a_trans, and r
# We set alpha_m = 0.1.
# 
popt, pcov = curve_fit(peturb.fsigma8_MG_fit, meas_zlist, meas_fsigma8, 
                       p0=[0., 1./(1.+7.), 4.],
                       sigma=meas_fsigma8_error, 
                       bounds=([-0.01, 1./(1.+10.), 1.5], [0.15, 1./(1.+0.), 5.]))

print popt

popt=np.array((popt))

np.savetxt('fsigma8_bestfit.dat', popt, delimiter='    ')

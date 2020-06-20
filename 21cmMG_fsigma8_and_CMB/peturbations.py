import numpy as np
from classy import Class

class Peturbations:

    def __init__(self):
        self.cosmo = Class()

    def cosmo_GR(self, zlist):
        """
        Compute cosmology for General Relativity 
        with Planck15 parameters.
        """
        zstr=','.join(map(str,zlist+zlist[-1]+2))
        lcdmpars = {'output': 'mPk, tCl, pCl, lCl',
                    'lensing': 'yes',
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

        self.cosmo.set(lcdmpars)
        self.cosmo.compute()

    def cosmo_MG(self, zlist, alpha_b, alpha_m, a_trans, r):
        """
        Compute cosmology for Modified gravity with a 
        LCDM background expansion (with Planck15
        parameters) and 'hill_scale' parameterization 
        for the property functions.
        """
        zstr=','.join(map(str,zlist+zlist[-1]+2))
        mgpars =   {'output': 'mPk, tCl, pCl, lCl',
                    'lensing': 'yes',
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

                    'Omega_Lambda': 0,
                    'Omega_fld': 0,
                    'Omega_smg': -1,
                    'expansion_model': 'lcdm',
                    'gravity_model': 'hill_scale',
                    }
        
        # Planck Mass
        Mstarsq = 1.
        # Property Functions
        alpha_t = 0.
        alpha_k = 0.001

        mgpars['parameters_smg'] = "{}, {}, {}, {}, {}, {}, {}".format(
                                   Mstarsq, 
                                   alpha_k, alpha_b, alpha_m, alpha_t,
                                   a_trans, r)
        self.cosmo.set(mgpars)
        self.cosmo.compute()

    def D(self, z):
        """ 
        Linear growth function
            D(z) = ( P(k_ls, z) / P(k_ls, 0) )^(1 / 2)
        where k_ls is a large-scale mode.
        """
        k=0.01 #Our choice of large-scale mode
        mPk=self.cosmo.pk(k,z)
        mPk_norm=self.cosmo.pk(k,0) #Normalize at z=0
        D=np.sqrt(mPk/mPk_norm)
        return D

    def f(self, z):
        """
        Linear growth rate f(z) where
            f(a) = d log D / d log a
        where a = 1 / (1 + z) is the scale factor.
        """
        a=1./(1.+z)
        da=0.01#*a
        gp,g,gm=[self.D(1./ia-1.) for ia in [a+da,a,a-da]]
        f=a*(gp-gm)/(2*g*da)
        return f

    def sigma8(self, z):
        """
        sigma8(z) = sigma8*D(z)
        """
        s8=self.cosmo.sigma8()*self.D(z)
        return s8

    def fsigma8(self, z):
        """
        Growth rate of structure
            fsigma8(z) = f(z) * sigma8(z)
        """
        fs8=self.f(z)*self.sigma8(z)
        return fs8

    def fsigma8_GR(self, zlist):
        """
        Get fsigma8 for General Relativity
        """
        self.cosmo_GR(zlist) 
        fs8_GR = np.array([self.fsigma8(z) for z in zlist])
        return fs8_GR

    def fsigma8_MG(self, zlist, alpha_b, alpha_m, a_trans, r):
        """
        Get fsigma8 for Modified Gravity with a LCDM background 
        expansion and hill_scale parameterization for the property
        functions.
        """
        self.cosmo_MG(zlist, alpha_b, alpha_m, a_trans, r)
        fs8_MG = np.array([self.fsigma8(z) for z in zlist])
        return fs8_MG

    def fsigma8_MG_fit(self, zlist, alpha_b, a_trans, r):
        """
        Get fsigma8 for Modified Gravity with a LCDM background 
        expansion and hill_scale parameterization for the property
        functions.
        """
        alpha_m = 0.1
        #a_trans = 1./(1.+7.)
        #r = 4.
        self.cosmo_MG(zlist, alpha_b, alpha_m, a_trans, r)
        fs8_MG = np.array([self.fsigma8(z) for z in zlist])
        return fs8_MG

    def sigma8_GR(self, zlist):
        """
        Get fsigma8 for General Relativity
        """
        self.cosmo_GR(zlist) 
        s8_GR = np.array([self.sigma8(z) for z in zlist])
        return s8_GR

    def sigma8_MG(self, zlist, alpha_b, alpha_m, a_trans, r):
        """
        Get sigma8 for Modified Gravity with a LCDM background 
        expansion and hill_scale parameterization for the property
        functions.
        """
        self.cosmo_MG(zlist, alpha_b, alpha_m, a_trans, r)
        s8_MG = np.array([self.sigma8(z) for z in zlist])
        return s8_MG

    def Cl_TT_GR(self, zlist, llist):
        """
        Get C_ell^TT for General Relativity
        """
        self.cosmo_GR(zlist) 
        lmax = llist[-1]
        Cl_TT_GR = self.cosmo.lensed_cl(lmax)['tt'][2:]
        Cl_TT_GR = (1./(2.*np.pi)) * (llist * (llist + 1.)) * Cl_TT_GR
        return Cl_TT_GR

    def Cl_TT_MG(self, zlist, llist, alpha_b, alpha_m, a_trans, r):
        """
        Get C_ell^TT for Modified Gravity with a LCDM background 
        expansion and hill_scale parameterization for the property
        functions.
        """
        self.cosmo_MG(zlist, alpha_b, alpha_m, a_trans, r)
        lmax = llist[-1]
        Cl_TT_MG = self.cosmo.lensed_cl(lmax)['tt'][2:]
        Cl_TT_MG = (1./(2.*np.pi)) * (llist * (llist + 1.)) * Cl_TT_MG
        return Cl_TT_MG

    def Cl_kappa_kappa_GR(self, zlist, llist):
        """
        Get C_ell^kappakappa for General Relativity
        """
        self.cosmo_GR(zlist) 
        lmax = llist[-1]
        Cl_phi_phi_GR = self.cosmo.lensed_cl(lmax)['pp'][2:]
        Cl_kap_kap_GR = (1./4.) * ((llist * (llist + 1.))**2.) * Cl_phi_phi_GR
        return Cl_kap_kap_GR

    def Cl_kappa_kappa_MG(self, zlist, llist, alpha_b, alpha_m, a_trans, r):
        """
        Get C_ell^kappakappa for Modified Gravity with a LCDM background 
        expansion and hill_scale parameterization for the property
        functions.
        """
        self.cosmo_MG(zlist, alpha_b, alpha_m, a_trans, r)
        lmax = llist[-1]
        Cl_phi_phi_MG = self.cosmo.lensed_cl(lmax)['pp'][2:]
        Cl_kap_kap_MG = (1./4.) * ((llist * (llist + 1.))**2.) * Cl_phi_phi_MG
        return Cl_kap_kap_MG


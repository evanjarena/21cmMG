/** @file background.c Documented background module
 *
 * Julien Lesgourgues, 17.04.2011
 * routines related to ncdm written by T. Tram in 2011
 *
 * Deals with the cosmological background evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the background, i.e. to integrate
 *    the background equations, and store all background quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondance between redhsift and conformal time.
 *
 *
 * The overall logic in this module is the following:
 *
 * 1. most background parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of a few variables that we will call {B} (in simplest
 * models, of the scale factor 'a'; in extended cosmologies, of 'a'
 * plus e.g. (phi, phidot) for quintessence, or some temperature for
 * exotic particles, etc...).
 *
 * 2. in turn, quantitites {B} can be found as a function of conformal
 * time by integrating the background equations.
 *
 * 3. some other quantitites that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be infered analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background
 *    quantitites {A} as a function of quantitites {B}.
 *
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for conformal time 'tau'; one or
 *    more for quantitites {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most procise way is to call directly
 * background_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'tau' and want
 * any other qunatity, we can call background_at_tau(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'tau' for a given redshift 'z': this can be done with
 * background_tau_of_z(). So if we are somewhere in the code, knowing
 * z and willing to get background quantitites, we should call first
 * background_tau_of_z() and then background_at_tau().
 *
 *
 * In order to save time, background_at_tau() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantitites, or rarely useful
 * quantities). Each line in the interpolation table is a vector which
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantitites like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantitites needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning
 * -# background_at_tau(), background_tau_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

#include "background.h"

/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table ant interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Ouput: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_tau(
                      struct background *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size comptible with return_format) */
                      ) {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that tau is in the pre-computed range */

  class_test(tau < pba->tau_table[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table[0]);

  class_test(tau > pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dtau2_table,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table,
                                                        pba->bt_size,
                                                        pba->background_table,
                                                        pba->d2background_dtau2_table,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;
}

/**
 * Conformal time at given redhsift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

int background_tau_of_z(
                        struct background *pba,
                        double z,
                        double * tau
                        ) {

  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index;

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < pba->z_table[pba->bt_size-1],
             pba->error_message,
             "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->z_table[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  class_call(array_interpolate_spline(
                                      pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      pba->d2tau_dz2_table,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Background quantities at given a.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of {B} parameters such as the scale factor 'a'
 * (see discussion at the begining of this file). In extended
 * comsological models, the pvecback_B vector contains other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background structure
 * @param a             Input: value of scale factor
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assmued to be already allocated)
 * @return the error status
 */

int background_functions(
                         struct background *pba,
                         double * pvecback_B, /* Vector containing all {B} quantities. */
                         short return_format,
                         double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size comptible with return_format) */
                         ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* total non-relativistic density */
  double rho_m;
  /* scale factor relative to scale factor today */
  double a_rel;
  /* background ncdm quantities */
  double rho_ncdm,p_ncdm,pseudo_p_ncdm;
  /* index for n_ncdm species */
  int n_ncdm;
  /* scale factor */
  double a;

  /** - initialize local variables */
  a = pvecback_B[pba->index_bi_a];
  rho_tot = 0.;
  p_tot = 0.;
  rho_r=0.;
  rho_m=0.;
  a_rel = a / pba->a_today;

  class_test(a_rel <= 0.,
             pba->error_message,
             "a = %e instead of strictly positive",a_rel);

  /** - pass value of a to output */
  pvecback[pba->index_bg_a] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a_rel,4);
  rho_tot += pvecback[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback[pba->index_bg_rho_g];
  rho_r += pvecback[pba->index_bg_rho_g];

  /* baryons */
  pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  rho_tot += pvecback[pba->index_bg_rho_b];
  p_tot += 0;
  rho_m += pvecback[pba->index_bg_rho_b];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }

  /* dcdm */
  if (pba->has_dcdm == _TRUE_) {
    /* Pass value of rho_dcdm to output */
    pvecback[pba->index_bg_rho_dcdm] = pvecback_B[pba->index_bi_rho_dcdm];
    rho_tot += pvecback[pba->index_bg_rho_dcdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_dcdm];
  }

  /* dr */
  if (pba->has_dr == _TRUE_) {
    /* Pass value of rho_dr to output */
    pvecback[pba->index_bg_rho_dr] = pvecback_B[pba->index_bi_rho_dr];
    rho_tot += pvecback[pba->index_bg_rho_dr];
    p_tot += (1./3.)*pvecback[pba->index_bg_rho_dr];
    rho_r += pvecback[pba->index_bg_rho_dr];
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {

    /* Loop over species: */
    for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){

      /* function returning background ncdm[n_ncdm] quantities (only
         those for which non-NULL pointers are passed) */
      class_call(background_ncdm_momenta(
                                         pba->q_ncdm_bg[n_ncdm],
                                         pba->w_ncdm_bg[n_ncdm],
                                         pba->q_size_ncdm_bg[n_ncdm],
                                         pba->M_ncdm[n_ncdm],
                                         pba->factor_ncdm[n_ncdm],
                                         1./a_rel-1.,
                                         NULL,
                                         &rho_ncdm,
                                         &p_ncdm,
                                         NULL,
                                         &pseudo_p_ncdm),
                 pba->error_message,
                 pba->error_message);

      pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[pba->index_bg_p_ncdm1+n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm] = pseudo_p_ncdm;

      /* (3 p_ncdm1) is the "relativistic" contrinution to rho_ncdm1 */
      rho_r += 3.* p_ncdm;

      /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution
         to rho_ncdm1 */
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_) {
    pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
    rho_tot += pvecback[pba->index_bg_rho_lambda];
    p_tot -= pvecback[pba->index_bg_rho_lambda];
  }

  /* fluid with w=w0+wa(1-a/a0) and constant cs2 */
  if (pba->has_fld == _TRUE_) {
    pvecback[pba->index_bg_rho_fld] = pba->Omega0_fld * pow(pba->H0,2)
      / pow(a_rel,3.*(1.+pba->w0_fld+pba->wa_fld))
      * exp(3.*pba->wa_fld*(a_rel-1.));
    rho_tot += pvecback[pba->index_bg_rho_fld];
    p_tot += (pba->w0_fld+pba->wa_fld*(1.-a_rel)) * pvecback[pba->index_bg_rho_fld];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_ur == _TRUE_) {
    pvecback[pba->index_bg_rho_ur] = pba->Omega0_ur * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback[pba->index_bg_rho_ur];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_ur];
    rho_r += pvecback[pba->index_bg_rho_ur];
  }

 /** - compute expansion rate H from Friedmann equation: this is the
      unique place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of [3c^2/8piG], ie
      rho_class = [8 pi G rho_physical / 3 c^2]
      NOTE: different computation if scalar field is present */  
  if (pba->has_smg == _FALSE_){
    
    pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);
    
    /** - compute derivative of H with respect to conformal time */
    pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;
  }
  /* Scalar field */
  /** if (has_smg) do all the mess to compute the Hubble rate, etc... */
  else{
    
//       /** - save rho_tot and p_tot without smg (it is important that smg is the last thing to be computed) */
      pvecback[pba->index_bg_rho_tot_wo_smg] = rho_tot;
      pvecback[pba->index_bg_p_tot_wo_smg] = p_tot;
      //NOTE: add the field energy and pressure after the debug has been added
     
    class_call(background_gravity_functions(pba,
					    pvecback_B,
					    return_format,
					    pvecback),
	       pba->error_message,
               pba->error_message);
    
    rho_tot += pvecback[pba->index_bg_rho_smg];
    p_tot += pvecback[pba->index_bg_p_smg];
    //divide relativistic & nonrelativistic (not very meaningful for oscillatory models)
    
    // TODO: where is this used? consider removing altogether
    rho_r += 3.*pvecback[pba->index_bg_p_smg]; //field pressure contributes radiation
    rho_m += pvecback[pba->index_bg_rho_smg] - 3.* pvecback[pba->index_bg_p_smg]; //the rest contributes matter     
  }

  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_tot;

  /** - compute other quantities in the exhaustive, redundent format: */
  if (return_format == pba->long_info) {

    /** - compute critical density */
    pvecback[pba->index_bg_rho_crit] = rho_tot-pba->K/a/a;
    class_test(pvecback[pba->index_bg_rho_crit] <= 0.,
               pba->error_message,
               "rho_crit = %e instead of strictly positive",pvecback[pba->index_bg_rho_crit]);

    /** - compute Omega_m */
    pvecback[pba->index_bg_Omega_m] = rho_m / rho_tot;

    /* one can put other variables here */
    /*  */
    /*  */

  }

  return _SUCCESS_;

}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input/Output : pointer to initialized background structure
 * @return the error status
 */

int background_init(
                    struct precision * ppr,
                    struct background * pba
                    ) {

  /** Summary: */

  /** - local variables : */
  int n_ncdm;
  double rho_ncdm_rel,rho_nu_rel;
  double Neff;
  int filenum=0;

  /** - in verbose mode, provide some information */
  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");

    /* below we want to inform the user about ncdm species*/
    if (pba->N_ncdm > 0) {

      Neff = pba->Omega0_ur/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;

      /* loop over ncdm species */
      for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

        /* inform if p-s-d read in files */
        if (pba->got_files[n_ncdm] == _TRUE_) {
          printf(" -> ncdm species i=%d read from file %s\n",n_ncdm+1,pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
          filenum++;
        }

        /* call this function to get rho_ncdm */
        background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                pba->w_ncdm_bg[n_ncdm],
                                pba->q_size_ncdm_bg[n_ncdm],
                                0.,
                                pba->factor_ncdm[n_ncdm],
                                0.,
                                NULL,
                                &rho_ncdm_rel,
                                NULL,
                                NULL,
                                NULL);

        /* inform user of the contribution of each species to
           radiation density (in relativistic limit): should be
           between 1.01 and 1.02 for each active neutrino species;
           evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
           density of one neutrino in the instantaneous decoupling
           limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
           from the definition of N_eff) */
        rho_nu_rel = 56.0/45.0*pow(_PI_,6)*pow(4.0/11.0,4.0/3.0)*_G_/pow(_h_P_,3)/pow(_c_,7)*
          pow(_Mpc_over_m_,2)*pow(pba->T_cmb*_k_B_,4);

        printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
               n_ncdm+1,
               pba->q_size_ncdm_bg[n_ncdm],
               pba->q_size_ncdm[n_ncdm],
               rho_ncdm_rel/rho_nu_rel);

        Neff += rho_ncdm_rel/rho_nu_rel;

      }

      printf(" -> total N_eff = %g (sumed over ultra-relativistic and ncdm species)\n",Neff);

    }
  }

  /** - if shooting failed during input, catch the error here: */
  class_test(pba->shooting_failed == _TRUE_,
             pba->error_message,
             "Shooting failed, try optimising input_get_guess(). Error message:\n\n%s",
             pba->shooting_error);

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  class_call(background_indices(pba),
             pba->error_message,
             pba->error_message);

  /** - control that cosmological parameter values make sense */

  /* H0 in Mpc^{-1} */
  class_test((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_),
             pba->error_message,
             "H0=%g out of bounds (%g<H0<%g) \n",pba->H0,_H0_SMALL_,_H0_BIG_);

  class_test(fabs(pba->h * 1.e5 / _c_  / pba->H0 -1.)>ppr->smallest_allowed_variation,
             pba->error_message,
             "inconsistency between Hubble and reduced Hubble parameters: you have H0=%f/Mpc=%fkm/s/Mpc, but h=%f",pba->H0,pba->H0/1.e5* _c_,pba->h);

  /* T_cmb in K */
  class_test((pba->T_cmb < _TCMB_SMALL_)||(pba->T_cmb > _TCMB_BIG_),
             pba->error_message,
             "T_cmb=%g out of bounds (%g<T_cmb<%g)",pba->T_cmb,_TCMB_SMALL_,_TCMB_BIG_);

  /* H0 in Mpc^{-1} */
  class_test((pba->Omega0_k < _OMEGAK_SMALL_)||(pba->Omega0_k > _OMEGAK_BIG_),
             pba->error_message,
             "Omegak = %g out of bounds (%g<Omegak<%g) \n",pba->Omega0_k,_OMEGAK_SMALL_,_OMEGAK_BIG_);

  /* fluid equation of state */
  if (pba->has_fld == _TRUE_) {
    class_test(pba->w0_fld+pba->wa_fld>=1./3.,
               pba->error_message,
               "Your choice for w0_fld+wa_fld=%g is suspicious, there would not be radiation domination at early times\n",
               pba->w0_fld+pba->wa_fld);
  }

  /* in verbose mode, inform the user about the value of the ncdm
     masses in eV and about the ratio [m/omega_ncdm] in eV (the usual
     93 point something)*/
  if ((pba->background_verbose > 0) && (pba->has_ncdm == _TRUE_)) {
    for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
      printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
             n_ncdm+1,
             pba->m_ncdm_in_eV[n_ncdm],
             pba->m_ncdm_in_eV[n_ncdm]*pba->deg_ncdm[n_ncdm]/pba->Omega0_ncdm[n_ncdm]/pba->h/pba->h);
    }
  }

  /* check other quantities which would lead to segmentation fault if zero */
  class_test(pba->a_today <= 0,
             pba->error_message,
             "input a_today = %e instead of strictly positive",pba->a_today);

  class_test(_Gyr_over_Mpc_ <= 0,
             pba->error_message,
             "_Gyr_over_Mpc = %e instead of strictly positive",_Gyr_over_Mpc_);

  /** - this function integrates the background over time, allocates
      and fills the background table */
  class_call(background_solve(ppr,pba),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 *
 *
 * @param pba Input : pointer to background structure (to be freed)
 * @return the error status
 */

int background_free(
                    struct background *pba
                    ) {
  int err;

  free(pba->tau_table);
  free(pba->z_table);
  free(pba->d2tau_dz2_table);
  free(pba->background_table);
  free(pba->d2background_dtau2_table);

  err = background_free_input(pba);

  return err;
}

/**
 * Free pointers inside background structure which were
 * allocated in input_read_parameters()
 *
 * @param pba Input : pointer to background structure
 * @return the error status
 */

int background_free_input(
                          struct background *pba
                          ) {

  int k;
  if (pba->Omega0_ncdm_tot != 0.){
    for(k=0; k<pba->N_ncdm; k++){
      free(pba->q_ncdm[k]);
      free(pba->w_ncdm[k]);
      free(pba->q_ncdm_bg[k]);
      free(pba->w_ncdm_bg[k]);
      free(pba->dlnf0_dlnq_ncdm[k]);
    }

    free(pba->q_ncdm);
    free(pba->w_ncdm);
    free(pba->q_ncdm_bg);
    free(pba->w_ncdm_bg);
    free(pba->dlnf0_dlnq_ncdm);
    free(pba->q_size_ncdm);
    free(pba->q_size_ncdm_bg);
    free(pba->M_ncdm);
    free(pba->T_ncdm);
    free(pba->ksi_ncdm);
    free(pba->deg_ncdm);
    free(pba->Omega0_ncdm);
    free(pba->m_ncdm_in_eV);
    free(pba->factor_ncdm);
    if(pba->got_files!=NULL)
      free(pba->got_files);
    if(pba->ncdm_psd_files!=NULL)
      free(pba->ncdm_psd_files);
    if(pba->ncdm_psd_parameters!=NULL)
      free(pba->ncdm_psd_parameters);
  }

  if (pba->Omega0_scf != 0.){
    if (pba->scf_parameters != NULL)
      free(pba->scf_parameters);
  }
  if (pba->Omega0_smg != 0.){
    if (pba->parameters_smg != NULL)
      free(pba->parameters_smg);
    //dealocate parameters_2_smg only for parameterizations
    if (pba->field_evolution_smg == _FALSE_ && pba->parameters_2_smg != NULL)
      free(pba->parameters_2_smg);
  }
  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input : pointer to background structure
 * @return the error status
 */

int background_indices(
                       struct background *pba
                       ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - initialize all flags: which species are present? */

  pba->has_cdm = _FALSE_;
  pba->has_ncdm = _FALSE_;
  pba->has_dcdm = _FALSE_;
  pba->has_dr = _FALSE_;
  pba->has_scf = _FALSE_;
  pba->has_smg = _FALSE_; /*Scalar field*/    
  pba->has_lambda = _FALSE_;
  pba->has_fld = _FALSE_;
  pba->has_ur = _FALSE_;
  pba->has_curvature = _FALSE_;

  if (pba->Omega0_cdm != 0.)
    pba->has_cdm = _TRUE_;

  if (pba->Omega0_ncdm_tot != 0.)
    pba->has_ncdm = _TRUE_;

  if (pba->Omega0_dcdmdr != 0.){
    pba->has_dcdm = _TRUE_;
    if (pba->Gamma_dcdm != 0.)
      pba->has_dr = _TRUE_;
  }
  
  if (pba->Omega0_smg != 0.)
    pba->has_smg = _TRUE_;    

  if (pba->Omega0_lambda != 0.)
    pba->has_lambda = _TRUE_;

  if (pba->Omega0_fld != 0.)
    pba->has_fld = _TRUE_;

  if (pba->Omega0_ur != 0.)
    pba->has_ur = _TRUE_;

  if (pba->sgnK != 0)
    pba->has_curvature = _TRUE_;

  /** - intialization of all indices */

  index_bg=0;

  /* index for scale factor */
  class_define_index(pba->index_bg_a,_TRUE_,index_bg,1);

  /* - indices for H and its conformal-time-derivative */
  class_define_index(pba->index_bg_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H_prime,_TRUE_,index_bg,1);

  /* - end of indices in the short vector of background values */
  pba->bg_size_short = index_bg;

  /* - index for rho_g (photon density) */
  class_define_index(pba->index_bg_rho_g,_TRUE_,index_bg,1);

  /* - index for rho_b (baryon density) */
  class_define_index(pba->index_bg_rho_b,_TRUE_,index_bg,1);

  /* - index for rho_cdm */
  class_define_index(pba->index_bg_rho_cdm,pba->has_cdm,index_bg,1);

  /* - indices for ncdm. We only define the indices for ncdm1
     (density, pressure, pseudo-pressure), the other ncdm indices
     are contiguous */
  class_define_index(pba->index_bg_rho_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_pseudo_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);

  /* - index for dcdm */
  class_define_index(pba->index_bg_rho_dcdm,pba->has_dcdm,index_bg,1);

  /* - index for dr */
  class_define_index(pba->index_bg_rho_dr,pba->has_dr,index_bg,1);
  
  /* - indices for scalar field (modified gravity) */
  class_define_index(pba->index_bg_M2_smg,pba->has_smg,index_bg,1);      

  class_define_index(pba->index_bg_rho_smg,pba->has_smg,index_bg,1);   
  class_define_index(pba->index_bg_p_smg,pba->has_smg,index_bg,1);  
  
  class_define_index(pba->index_bg_kineticity_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_braiding_smg,pba->has_smg,index_bg,1);  
  class_define_index(pba->index_bg_tensor_excess_smg,pba->has_smg,index_bg,1);   

  class_define_index(pba->index_bg_mpl_running_smg,pba->has_smg,index_bg,1);  
  class_define_index(pba->index_bg_kineticity_prime_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_braiding_prime_smg,pba->has_smg,index_bg,1);  
  class_define_index(pba->index_bg_mpl_running_prime_smg,pba->has_smg,index_bg,1);  
  class_define_index(pba->index_bg_tensor_excess_prime_smg,pba->has_smg,index_bg,1);    
  
  class_define_index(pba->index_bg_cs2_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_cs2num_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_kinetic_D_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_kinetic_D_prime_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_1_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_2_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_3_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_4_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_5_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_6_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_7_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_lambda_8_smg,pba->has_smg,index_bg,1);

  class_define_index(pba->index_bg_rho_tot_wo_smg,pba->has_smg,index_bg,1);  
  class_define_index(pba->index_bg_p_tot_wo_smg,pba->has_smg,index_bg,1);   
  
  class_define_index(pba->index_bg_H_prime_prime,pba->has_smg,index_bg,1);  
  class_define_index(pba->index_bg_p_tot_wo_prime_smg,pba->has_smg,index_bg,1);
  class_define_index(pba->index_bg_p_prime_smg,pba->has_smg,index_bg,1);
  

  /* - index for Lambda */
  class_define_index(pba->index_bg_rho_lambda,pba->has_lambda,index_bg,1);

  /* - index for fluid */
  class_define_index(pba->index_bg_rho_fld,pba->has_fld,index_bg,1);

  /* - index for ultra-relativistic neutrinos/species */
  class_define_index(pba->index_bg_rho_ur,pba->has_ur,index_bg,1);

  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_r,_TRUE_,index_bg,1);

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */
  /*    */

  /* - end of indices in the normal vector of background values */
  pba->bg_size_normal = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  class_define_index(pba->index_bg_rho_crit,_TRUE_,index_bg,1);

  /* - index for Omega_m (non-relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_m,_TRUE_,index_bg,1);

  /* -> conformal distance */
  class_define_index(pba->index_bg_conf_distance,_TRUE_,index_bg,1);

  /* -> angular diameter distance */
  class_define_index(pba->index_bg_ang_distance,_TRUE_,index_bg,1);

  /* -> luminosity distance */
  class_define_index(pba->index_bg_lum_distance,_TRUE_,index_bg,1);

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bg_time,_TRUE_,index_bg,1);

  /* -> conformal sound horizon */
  class_define_index(pba->index_bg_rs,_TRUE_,index_bg,1);

  /* -> density growth factor in dust universe */
  class_define_index(pba->index_bg_D,_TRUE_,index_bg,1);

  /* -> velocity growth factor in dust universe */
  class_define_index(pba->index_bg_f,_TRUE_,index_bg,1);

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  pba->bg_size = index_bg;

  /* - now, indices in vector of variables to integrate.
     First {B} variables, then {C} variables. */

  index_bi=0;

  /* -> scale factor */
  class_define_index(pba->index_bi_a,_TRUE_,index_bi,1);

  /* -> energy density in DCDM */
  class_define_index(pba->index_bi_rho_dcdm,pba->has_dcdm,index_bi,1);

  /* -> energy density in DR */
  class_define_index(pba->index_bi_rho_dr,pba->has_dr,index_bi,1);
  
  /* -> scalar field and its derivative wrt conformal time (only if needs to evolve the field)
   * plus other parameters that might be integrated in certain parameterizations
   */
  if (pba->has_smg == _TRUE_){    
    //if model needs to integrate M_pl from alpha_M, declare an index
    class_define_index(pba->index_bi_M_pl_smg,
		       pba->M_pl_evolution_smg,index_bi,1);
  }  

  /* End of {B} variables, now continue with {C} variables */
  pba->bi_B_size = index_bi;

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bi_time,_TRUE_,index_bi,1);

  /* -> sound horizon */
  class_define_index(pba->index_bi_rs,_TRUE_,index_bi,1);

  /* -> integral for growth factor */
  class_define_index(pba->index_bi_growth,_TRUE_,index_bi,1);

  /* -> index for conformal time in vector of variables to integrate */
  class_define_index(pba->index_bi_tau,_TRUE_,index_bi,1);

  /* -> end of indices in the vector of variables to integrate */
  pba->bi_size = index_bi;

  /* index_bi_tau must be the last index, because tau is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over tau itself) */
  class_test(pba->index_bi_tau != index_bi-1,
             pba->error_message,
             "background integration requires index_bi_tau to be the last of all index_bi's");

  /* flags for calling the interpolation routine */

  pba->short_info=0;
  pba->normal_info=1;
  pba->long_info=2;

  pba->inter_normal=0;
  pba->inter_closeby=1;

  return _SUCCESS_;

}

/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified (it is the only place to modify if you
 * need a partlar f0(q))
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_distribution(
                                 void * pbadist,
                                 double q,
                                 double * f0
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;
  int n_ncdm,lastidx;
  double ksi;
  double qlast,dqlast,f0last,df0last;
  double *param;
  /** Variables corresponing to entries in param: */
  //double square_s12,square_s23,square_s13;
  //double mixing_matrix[3][3];
  //int i;

  /** - extract from the input structure pbadist all the relevant information */
  pbadist_local = pbadist;          /* restore actual format of pbadist */
  pba = pbadist_local->pba;         /* extract the background structure from it */
  param = pba->ncdm_psd_parameters; /* extract the optional parameter list from it */
  n_ncdm = pbadist_local->n_ncdm;   /* extract index of ncdm species under consideration */
  ksi = pba->ksi_ncdm[n_ncdm];      /* extract chemical potential */

  /** - shall we interpolate in file, or shall we use analytical formula below? */

  /** -> deal first with the case of interpolating in files */
  if (pba->got_files[n_ncdm]==_TRUE_) {

    lastidx = pbadist_local->tablesize-1;
    if(q<pbadist_local->q[0]){
      //Handle q->0 case:
      *f0 = pbadist_local->f0[0];
    }
    else if(q>pbadist_local->q[lastidx]){
      //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast=qlast - pbadist_local->q[lastidx-1];
      df0last=f0last - pbadist_local->f0[lastidx-1];

      *f0 = f0last*exp(-(qlast-q)*df0last/f0last/dqlast);
    }
    else{
      //Do interpolation:
      class_call(array_interpolate_spline(
                                          pbadist_local->q,
                                          pbadist_local->tablesize,
                                          pbadist_local->f0,
                                          pbadist_local->d2f0,
                                          1,
                                          q,
                                          &pbadist_local->last_index,
                                          f0,
                                          1,
                                          pba->error_message),
                 pba->error_message,     pba->error_message);
    }
  }

  /** -> deal now with case of reading analytical function */
  else{
    /**
       Enter here your analytic expression(s) for the p.s.d.'s. If
       you need different p.s.d.'s for different species, put each
       p.s.d inside a condition, like for instance: if (n_ncdm==2) {
       *f0=...}.  Remember that n_ncdm = 0 refers to the first
       species.
    */

    /**************************************************/
    /*    FERMI-DIRAC INCLUDING CHEMICAL POTENTIALS   */
    /**************************************************/

    *f0 = 1.0/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));

    /**************************************************/

    /** This form is only appropriate for approximate studies, since in
        reality the chemical potential are associated with flavor
        eigenstates, not mass eigenstates. It is easy to take this into
        account by introducing the mixing angles. In the later part
        (not read by the code) we illustrate how to do this */

    if (_FALSE_) {

      /* We must use the list of extra parameters read in input, stored in the
         ncdm_psd_parameter list, extracted above from the structure
         and now called param[..] */

      /* check that this list has been read */
      class_test(param == NULL,
                 pba->error_message,
                 "Analytic expression wants to use 'ncdm_psd_parameters', but they have not been entered!");

      /* extract values from the list (in this example, mixing angles) */
      double square_s12=param[0];
      double square_s23=param[1];
      double square_s13=param[2];

      /* infer mixing matrix */
      double mixing_matrix[3][3];
      int i;

      mixing_matrix[0][0]=pow(fabs(sqrt((1-square_s12)*(1-square_s13))),2);
      mixing_matrix[0][1]=pow(fabs(sqrt(square_s12*(1-square_s13))),2);
      mixing_matrix[0][2]=fabs(square_s13);
      mixing_matrix[1][0]=pow(fabs(sqrt((1-square_s12)*square_s13*square_s23)+sqrt(square_s12*(1-square_s23))),2);
      mixing_matrix[1][1]=pow(fabs(sqrt(square_s12*square_s23*square_s13)-sqrt((1-square_s12)*(1-square_s23))),2);
      mixing_matrix[1][2]=pow(fabs(sqrt(square_s23*(1-square_s13))),2);
      mixing_matrix[2][0]=pow(fabs(sqrt(square_s12*square_s23)-sqrt((1-square_s12)*square_s13*(1-square_s23))),2);
      mixing_matrix[2][1]=pow(sqrt((1-square_s12)*square_s23)+sqrt(square_s12*square_s13*(1-square_s23)),2);
      mixing_matrix[2][2]=pow(fabs(sqrt((1-square_s13)*(1-square_s23))),2);

      /* loop over flavor eigenstates and compute psd of mass eigenstates */
      *f0=0.0;
      for(i=0;i<3;i++){

    	*f0 += mixing_matrix[i][n_ncdm]*1.0/pow(2*_PI_,3)*(1./(exp(q-pba->ksi_ncdm[i])+1.) +1./(exp(q+pba->ksi_ncdm[i])+1.));

      }
    } /* end of region not used, but shown as an example */
  }

  return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weigths. The logic is: if we can convolve accurately
 * f0(q) with this function, then we can convolve it accuractely with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_,4));
  double e = 2.0/(45.0*_zeta5_);

  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as 1/r or 1/r^2 for r->0 */
  *test = pow(2.0*_PI_,3)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q);

  return _SUCCESS_;
}

/**
 * This function finds optimal quadrature weights for each ncdm
 * species
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_ncdm_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  int index_q, k,tolexp,row,status,filenum;
  double f0m2,f0m1,f0,f0p1,f0p2,dq,q,df0dq,tmp1,tmp2;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  /* Allocate pointer arrays: */
  class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);

  /* Allocate pointers: */
  class_alloc(pba->q_size_ncdm,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_size_ncdm_bg,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->factor_ncdm,sizeof(double)*pba->N_ncdm,pba->error_message);

  for(k=0, filenum=0; k<pba->N_ncdm; k++){
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    /*Do we need to read in a file to interpolate the distribution function? */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      psdfile = fopen(pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(psdfile == NULL,pba->error_message,
                 "Could not open file %s!",pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
      // Find size of table:
      for (row=0,status=2; status==2; row++){
        status = fscanf(psdfile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row-1;

      /*Allocate room for interpolation table: */
      class_alloc(pbadist.q,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.d2f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      for (row=0; row<pbadist.tablesize; row++){
        status = fscanf(psdfile,"%lf %lf",
                        &pbadist.q[row],&pbadist.f0[row]);
        //		printf("(q,f0) = (%g,%g)\n",pbadist.q[row],pbadist.f0[row]);
      }
      fclose(psdfile);
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pbadist.q,
                                          pbadist.tablesize,
                                          pbadist.f0,
                                          1,
                                          pbadist.d2f0,
                                          _SPLINE_EST_DERIV_,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
      filenum++;
    }

    /* Handle perturbation qsampling: */
    class_alloc(pba->q_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_ncdm[k],
                             pba->w_ncdm[k],
                             &(pba->q_size_ncdm[k]),
                             _QUADRATURE_MAX_,
                             ppr->tol_ncdm,
                             pbadist.q,
                             pbadist.tablesize,
                             background_ncdm_test_function,
                             background_ncdm_distribution,
                             &pbadist,
                             pba->error_message),
               pba->error_message,
               pba->error_message);
    pba->q_ncdm[k]=realloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));
    pba->w_ncdm[k]=realloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));


    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
             k+1,
             pba->q_size_ncdm[k]);

    /* Handle background q_sampling: */
    class_alloc(pba->q_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_ncdm_bg[k],
                             pba->w_ncdm_bg[k],
                             &(pba->q_size_ncdm_bg[k]),
                             _QUADRATURE_MAX_BG_,
                             ppr->tol_ncdm_bg,
                             pbadist.q,
                             pbadist.tablesize,
                             background_ncdm_test_function,
                             background_ncdm_distribution,
                             &pbadist,
                             pba->error_message),
               pba->error_message,
               pba->error_message);


    pba->q_ncdm_bg[k]=realloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));
    pba->w_ncdm_bg[k]=realloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));

    /** - in verbose mode, inform user of number of sampled momenta
        for background quantities */
    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
             k+1,
             pba->q_size_ncdm_bg[k]);

    class_alloc(pba->dlnf0_dlnq_ncdm[k],
                pba->q_size_ncdm[k]*sizeof(double),
                pba->error_message);


    for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
      q = pba->q_ncdm[k][index_q];
      class_call(background_ncdm_distribution(&pbadist,q,&f0),
                 pba->error_message,pba->error_message);

      //Loop to find appropriate dq:
      for(tolexp=_PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++){

        if (index_q == 0){
          dq = MIN((0.5-ppr->smallest_allowed_variation)*q,2*exp(tolexp)*(pba->q_ncdm[k][index_q+1]-q));
        }
        else if (index_q == pba->q_size_ncdm[k]-1){
          dq = exp(tolexp)*2.0*(pba->q_ncdm[k][index_q]-pba->q_ncdm[k][index_q-1]);
        }
        else{
          dq = exp(tolexp)*(pba->q_ncdm[k][index_q+1]-pba->q_ncdm[k][index_q-1]);
        }

        class_call(background_ncdm_distribution(&pbadist,q-2*dq,&f0m2),
                   pba->error_message,pba->error_message);
        class_call(background_ncdm_distribution(&pbadist,q+2*dq,&f0p2),
                   pba->error_message,pba->error_message);

        if (fabs((f0p2-f0m2)/f0)>sqrt(ppr->smallest_allowed_variation)) break;
      }

      class_call(background_ncdm_distribution(&pbadist,q-dq,&f0m1),
                 pba->error_message,pba->error_message);
      class_call(background_ncdm_distribution(&pbadist,q+dq,&f0p1),
                 pba->error_message,pba->error_message);
      //5 point estimate of the derivative:
      df0dq = (+f0m2-8*f0m1+8*f0p1-f0p2)/12.0/dq;
      //printf("df0dq[%g] = %g. dlf=%g ?= %g. f0 =%g.\n",q,df0dq,q/f0*df0dq,
      //Avoid underflow in extreme tail:
      if (fabs(f0)==0.)
        pba->dlnf0_dlnq_ncdm[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
      else
        pba->dlnf0_dlnq_ncdm[k][index_q] = q/f0*df0dq;
    }

    pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow(pba->T_cmb*pba->T_ncdm[k]*_k_B_,4)*8*_PI_*_G_
      /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_;

    /* If allocated, deallocate interpolation table:  */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      free(pbadist.q);
      free(pbadist.f0);
      free(pbadist.d2f0);
    }
  }


  return _SUCCESS_;
}

/**
 * For a given ncdm sepcies: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 *
 * @param qvec     Input: smapled momenta
 * @param wvec     Input: quadrature weigths
 * @param qsize    Input: number of momenta/weigths
 * @param M        Input: mass
 * @param factor   Input: normalization factor for the p.s.d.
 * @param z        Input: redhsift
 * @param n        Output: number density
 * @param rho      Output: energy density
 * @param p        Output: pressure
 * @param drho_dM  Output: derivative used in next function
 * @param pseudo_p Ouput: pseudo-pressure used in perturbation module for fluid approx
 *
 */

int background_ncdm_momenta(
                            /* Only calculate for non-NULL pointers: */
                            double * qvec,
                            double * wvec,
                            int qsize,
                            double M,
                            double factor,
                            double z,
                            double * n,
                            double * rho, // density
                            double * p,   // pressure
                            double * drho_dM,  // d rho / d M used in next function
                            double * pseudo_p  // pseudo-p used in ncdm fluid approx
                            ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;

  /** - rescale normalization at given redshift */
  factor2 = factor*pow(1+z,4);

  /** - initialize quantities */
  if (n!=NULL) *n = 0.;
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dM!=NULL) *drho_dM = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;

  /** - loop over momenta */
  for (index_q=0; index_q<qsize; index_q++) {

    /* squared momentum */
    q2 = qvec[index_q]*qvec[index_q];

    /* energy */
    epsilon = sqrt(q2+M*M/(1.+z)/(1.+z));

    /* integrand of the various quantities */
    if (n!=NULL) *n += q2*wvec[index_q];
    if (rho!=NULL) *rho += q2*epsilon*wvec[index_q];
    if (p!=NULL) *p += q2*q2/3./epsilon*wvec[index_q];
    if (drho_dM!=NULL) *drho_dM += q2*M/(1.+z)/(1.+z)/epsilon*wvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += pow(q2/epsilon,3)/3.0*wvec[index_q];
  }

  /** - ajust normalization */
  if (n!=NULL) *n *= factor2*(1.+z);
  if (rho!=NULL) *rho *= factor2;
  if (p!=NULL) *p *= factor2;
  if (drho_dM!=NULL) *drho_dM *= factor2;
  if (pseudo_p!=NULL) *pseudo_p *=factor2;

  return _SUCCESS_;
}

/**
 * When the user passed in input the density fraction Omeha_ncdm or
 * omega_ncdm but not the mass, infer the mass with Newton iteration method.
 *
 * @param ppr    Input: precision structure
 * @param pba    Input/Output: background structure
 * @param n_ncdm Input: index of ncdm species
 */

int background_ncdm_M_from_Omega(
                                 struct precision *ppr,
                                 struct background *pba,
                                 int n_ncdm
                                 ) {
  double rho0,rho,n,M,deltaM,drhodM;
  int iter,maxiter=50;

  rho0 = pba->H0*pba->H0*pba->Omega0_ncdm[n_ncdm]; /*Remember that rho is defined such that H^2=sum(rho_i) */
  M = 0.0;

  background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                          pba->w_ncdm_bg[n_ncdm],
                          pba->q_size_ncdm_bg[n_ncdm],
                          M,
                          pba->factor_ncdm[n_ncdm],
                          0.,
                          &n,
                          &rho,
                          NULL,
                          NULL,
                          NULL);

  /* Is the value of Omega less than a massless species?*/
  class_test(rho0<rho,pba->error_message,
             "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
             n_ncdm,pba->Omega0_ncdm[n_ncdm],pba->Omega0_ncdm[n_ncdm]*rho/rho0);

  /* In the strict NR limit we have rho = n*(M) today, giving a zero'th order guess: */
  M = rho0/n; /* This is our guess for M. */
  for (iter=1; iter<=maxiter; iter++){

    /* Newton iteration. First get relevant quantities at M: */
    background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                            pba->w_ncdm_bg[n_ncdm],
                            pba->q_size_ncdm_bg[n_ncdm],
                            M,
                            pba->factor_ncdm[n_ncdm],
                            0.,
                            NULL,
                            &rho,
                            NULL,
                            &drhodM,
                            NULL);

    deltaM = (rho0-rho)/drhodM; /* By definition of the derivative */
    if ((M+deltaM)<0.0) deltaM = -M/2.0; /* Avoid overshooting to negative M value. */
    M += deltaM; /* Update value of M.. */
    if (fabs(deltaM/M)<ppr->tol_M_ncdm){
      /* Accuracy reached.. */
      pba->M_ncdm[n_ncdm] = M;
      break;
    }
  }
  class_test(iter>=maxiter,pba->error_message,
             "Newton iteration could not converge on a mass for some reason.");
  return _SUCCESS_;
}

/**
 *  This function integrates the background over time, allocates and
 *  fills the background table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_solve(
                     struct precision *ppr,
                     struct background *pba
                     ) {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  struct background_parameters_and_workspace bpaw;
  /* a growing table (since the number of time steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* initial conformal time */
  double tau_start;
  /* final conformal time */
  double tau_end;
  /* an index running over bi indices */
  int i;
  /* vector of quantities to be integrated */
  double * pvecback_integration;
  /* vector of all background quantities */
  double * pvecback;
  /* necessary for calling array_interpolate(), but never used */
  int last_index=0;
  /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
  double comoving_radius=0.;

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;

  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */

  /* Size of vector to integrate is (pba->bi_size-1) rather than
   * (pba->bi_size), since tau is not integrated.
   */
  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
             gi.error_message,
             pba->error_message);

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);

  /* here tau_end is in fact the initial time (in the next loop
     tau_start = tau_end) */
  tau_end=pvecback_integration[pba->index_bi_tau];

  /** - create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             pba->error_message);

  /* initialize the counter for the number of steps */
  pba->bt_size=0;

  /** - loop over integration steps : call background_functions(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of tau */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    tau_start = tau_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      tau_end = tau_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }
    else {
      tau_end = tau_start + (pba->a_today/pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    /* -> store value of tau */
    pvecback_integration[pba->index_bi_tau]=tau_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
             gTable.error_message,
             pba->error_message);
  pba->bt_size++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pba->error_message);

  /** - interpolate to get quantities precisely today with array_interpolate() */
  class_call(array_interpolate(
                               pData,
                               pba->bi_size,
                               pba->bt_size,
                               pba->index_bi_a,
                               pba->a_today,
                               &last_index,
                               pvecback_integration,
                               pba->bi_size,
                               pba->error_message),
             pba->error_message,
             pba->error_message);

  /* substitute last line with quantities today */
  for (i=0; i<pba->bi_size; i++)
    pData[(pba->bt_size-1)*pba->bi_size+i]=pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_tau];
  /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
  if (pba->has_dcdm == _TRUE_){
    pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm]/pba->H0/pba->H0;
  }
  if (pba->has_dr == _TRUE_){
    pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr]/pba->H0/pba->H0;
  }


  /** - allocate background tables */
  class_alloc(pba->tau_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->z_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2tau_dz2_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->background_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_dtau2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  /** - In a loop over lines, fill background table using the result of the integration plus background_functions() */
  for (i=0; i < pba->bt_size; i++) {

    /* -> establish correspondance between the integrated variable and the bg variables */

    pba->tau_table[i] = pData[i*pba->bi_size+pba->index_bi_tau];

    class_test(pData[i*pba->bi_size+pba->index_bi_a] <= 0.,
               pba->error_message,
               "a = %e instead of strictly positiv",pData[i*pba->bi_size+pba->index_bi_a]);

    pba->z_table[i] = pba->a_today/pData[i*pba->bi_size+pba->index_bi_a]-1.;

    pvecback[pba->index_bg_time] = pData[i*pba->bi_size+pba->index_bi_time];
    pvecback[pba->index_bg_conf_distance] = pba->conformal_age - pData[i*pba->bi_size+pba->index_bi_tau];

    if (pba->sgnK == 0) comoving_radius = pvecback[pba->index_bg_conf_distance];
    else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(pba->K);
    else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(-pba->K);

    pvecback[pba->index_bg_ang_distance] = pba->a_today*comoving_radius/(1.+pba->z_table[i]);
    pvecback[pba->index_bg_lum_distance] = pba->a_today*comoving_radius*(1.+pba->z_table[i]);
    pvecback[pba->index_bg_rs] = pData[i*pba->bi_size+pba->index_bi_rs];

    /* -> compute all other quantities depending only on {B} variables.
       The value of {B} variables in pData are also copied to pvecback.*/
    class_call(background_functions(pba,pData+i*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

    /* -> compute growth functions (valid in dust universe) */

    /* D = H \int [da/(aH)^3] = H \int [dtau/(aH^2)] = H * growth */
    pvecback[pba->index_bg_D] = pvecback[pba->index_bg_H]*pData[i*pba->bi_size+pba->index_bi_growth];

    /* f = [dlnD]/[dln a] = 1/(aH) [dlnD]/[dtau] = H'/(aH^2) + 1/(a^2 H^3 growth) */
    pvecback[pba->index_bg_f] = pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_a]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_H] + 1./(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pData[i*pba->bi_size+pba->index_bi_growth]);

    /* -> write in the table */
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_table + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_table");
  }

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);

  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      1,
                                      pba->d2tau_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table,
                                      pba->bt_size,
                                      pba->background_table,
                                      pba->bg_size,
                                      pba->d2background_dtau2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);
  
/** - second loop over lines, overwrite derivatives that can't be analytically computed from background_functions 
   * Fill the derivatives of the Bellini-Sawicki functions in pvecback
   * This is done just by overwriting the pvecback entries corresponding to the relevant indice
   */ 
  
  double * pvecback_derivs;
  class_alloc(pvecback_derivs,pba->bg_size*sizeof(double),pba->error_message);  
  
  for (i=0; i < pba->bt_size; i++) {

    // write the derivatives in the structure
    class_call(array_derivate_spline(pba->tau_table, // x_array 
				     pba->bt_size, // int n_lines
				     pba->background_table, // array
				     pba->d2background_dtau2_table, // double * array_splined
				     pba->bg_size, // n_columns
				     pba->tau_table[i], // double x -> tau
				     &last_index, // int* last_index // this is something for the interpolation to talk to each other when using a loop
				     pvecback_derivs, // double * result 
				     pba->bg_size, //result_size, from 1 to n_columns
				     pba->error_message),
		pba->error_message,
		pba->error_message);  
            
      /* -> write in the table (overwrite the alpha time derivatives, which were set to nan in background_functions)
       * direction of copy: add the corresponding indices to the coordinates
       * thing to be copied: the address (&) to the element of pvecback corresponding to the index we want
       * size: just a single double number
       * -> repeat for all necessary quantities
       */
    if (pba->has_smg == _TRUE_){
      
      //Need to update pvecback
      class_call(background_at_tau(pba,
				   pba->tau_table[i],
				   pba->long_info,
				   pba->inter_normal,
				   &last_index, //should be no problem to use the same one as for the derivatives
				   pvecback),
		 pba->error_message,
		 pba->error_message);
      
      /*NOTE: here we compute the derivatives of quantities coputed in background_gravity_functions during the integration.
       * for quantities that depend on these derivatives (e.g. the gamma_i functions determining the effective mass)
       * there is an additional loop at the end of background_solve
       */
      
      // Kineticity'
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_kineticity_prime_smg,
			      &pvecback_derivs[pba->index_bg_kineticity_smg],
			      1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_kineticity_prime_smg,
               pba->error_message,
               "cannot copy data back to pba->background_table");
      
      //Braiding'
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_braiding_prime_smg,
			      &pvecback_derivs[pba->index_bg_braiding_smg],
			      1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_braiding_prime_smg,
               pba->error_message,
               "cannot copy data back to pba->background_table");
      
      //Planck mass run rate'
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_prime_smg,
			      &pvecback_derivs[pba->index_bg_mpl_running_smg],
			      1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_prime_smg,
               pba->error_message,
               "cannot copy data back to pba->background_table");
      
      //Tensor excess'
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_tensor_excess_prime_smg,
			      &pvecback_derivs[pba->index_bg_tensor_excess_smg],
			      1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_tensor_excess_prime_smg,
               pba->error_message,
               "cannot copy data back to pba->background_table");

      //H''
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_H_prime_prime,
         		&pvecback_derivs[pba->index_bg_H_prime],
         		1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_H_prime_prime,
              pba->error_message,
              "cannot copy data back to pba->background_table");
      
      // p_tot_wo_smg'
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_p_tot_wo_prime_smg,
            &pvecback_derivs[pba->index_bg_p_tot_wo_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_p_tot_wo_prime_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");
      
      // p_smg'
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_p_prime_smg,
            &pvecback_derivs[pba->index_bg_p_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_p_prime_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");

      //D'
      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_prime_smg,
			      &pvecback_derivs[pba->index_bg_kinetic_D_smg],
			      1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_prime_smg,
               pba->error_message,
               "cannot copy data back to pba->background_table");

      // Planck's mass running
      // Only need to compute it if neither self consistent field evolution nor evolving M_pl in terms of alpha_M
      // check equation 3.3 of Bellini & Sawicki 2014 
      
      if (pba->field_evolution_smg == _FALSE_ && pba->M_pl_evolution_smg == _FALSE_){
	
	double alpha_M = pvecback_derivs[pba->index_bg_M2_smg]/pvecback[pba->index_bg_M2_smg]/pvecback[pba->index_bg_a]/pvecback[pba->index_bg_H];

	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_smg,
				&alpha_M, //write using the address
				1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_smg,
		   pba->error_message,
		   "cannot copy data back to pba->background_table");
      }      
      
      double a = pvecback[pba->index_bg_a];
      double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
      double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
      double p_smg = pvecback[pba->index_bg_p_smg];
      double rho_smg = pvecback[pba->index_bg_rho_smg];
      
      if(pba->background_verbose > 15 && fabs(1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H])>1e-8)
	printf("a = %g, (delta H')/H' = %g \n", a, 1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H]);
      
      
      //TODO: clean up this!!
      // speed of sound, wanted for output
      
      double H = pvecback[pba->index_bg_H];
      double H_p = pvecback[pba->index_bg_H_prime];

      double M2 = pvecback[pba->index_bg_M2_smg];
      double kin = pvecback[pba->index_bg_kineticity_smg];
      double bra = pvecback[pba->index_bg_braiding_smg];
      double run = pvecback[pba->index_bg_mpl_running_smg];
      double ten = pvecback[pba->index_bg_tensor_excess_smg];
      
      //need to update the time derivatives of the interesting functions 
      
      double kin_p = pvecback_derivs[pba->index_bg_kineticity_smg];
      double bra_p = pvecback_derivs[pba->index_bg_braiding_smg];
      double run_p = pvecback_derivs[pba->index_bg_mpl_running_smg];
      double ten_p = pvecback_derivs[pba->index_bg_tensor_excess_smg];
      double p_tot_p = pvecback_derivs[pba->index_bg_p_tot_wo_smg];
      double p_smg_p = pvecback_derivs[pba->index_bg_p_smg];
      double D = pvecback[pba->index_bg_kinetic_D_smg];
      
      pvecback[pba->index_bg_lambda_1_smg] = (run + (-1.)*ten)*(-3.)*bra + (1. + ten)*kin;
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_1_smg,
            &pvecback[pba->index_bg_lambda_1_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_1_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      

      
      pvecback[pba->index_bg_lambda_2_smg] = (2.*(1. + (-1.)*M2) + bra*M2)*(rho_tot + p_tot)*(-3.)/2.*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/2.*pow(H,-2) + pow(H,-1)*bra_p*pow(a,-1);
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_2_smg,
            &pvecback[pba->index_bg_lambda_2_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_2_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      

      
      pvecback[pba->index_bg_lambda_3_smg] = (2. + run)*(-1.)/2.*D + (-3.)/4.*bra*pvecback[pba->index_bg_lambda_2_smg];
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_3_smg,
            &pvecback[pba->index_bg_lambda_3_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_3_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      

      
      pvecback[pba->index_bg_lambda_4_smg] = kin*pvecback[pba->index_bg_lambda_2_smg] + (2.*kin*bra_p + (-1.)*bra*kin_p)*(-1.)*pow(H,-1)*pow(a,-1);
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_4_smg,
            &pvecback[pba->index_bg_lambda_4_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_4_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      

      
      pvecback[pba->index_bg_lambda_5_smg] = (bra + 2.*run + (-2.)*ten + bra*ten)*3./2.*bra + (run + (-1.)*ten)*D + 3./2.*bra*pvecback[pba->index_bg_lambda_2_smg];
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_5_smg,
            &pvecback[pba->index_bg_lambda_5_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_5_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");


      pvecback[pba->index_bg_lambda_6_smg] = 3./2.*(((9./2.*bra + kin)*(M2 - 1.)*pow(M2,-1) + (-9.)/4.*pow(bra,2) - bra*kin/2. + D*run)*pow(rho_tot,2) + ((9.*bra + kin)*(M2 - 1.)*pow(M2,-1) + (-9.)/2.*pow(bra,2) - bra*kin/2. + D*run)*rho_tot*p_tot + 9./2.*bra*(M2 - 1. - M2*bra/2.)*pow(M2,-1)*pow(p_tot,2) + (kin*(M2 - 1.)*pow(M2,-1) - bra*kin/2. + D*run)*(rho_tot + p_tot)*rho_smg + ((kin - bra*kin/2. + D*run)*rho_smg + ((9.*bra + kin)*(2. - bra)/2. + D*run - 9./2.*bra*pow(M2,-1))*rho_tot + 9.*bra*(1. - bra/2. - pow(M2,-1)/2.)*p_tot)*(rho_smg + p_smg) + 9./2.*bra*(1. - bra/2.)*pow(rho_smg + p_smg,2))*pow(H,-4) + (((9.*bra*(rho_tot + p_tot) - 2.*kin*(rho_tot + rho_smg)) + (rho_smg + p_smg)*9.*bra)*bra_p/2. + (rho_tot + rho_smg)*bra*kin_p + (-2.*(1. - M2)*kin + 3.*pow(bra,2)*M2)*3./2.*pow(M2,-1)*p_tot_p + 3.*D*p_smg_p)*pow(H,-3)*pow(a,-1)/2.;
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_6_smg,
            &pvecback[pba->index_bg_lambda_6_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_6_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      

      
      pvecback[pba->index_bg_lambda_7_smg] = ((-2.) + bra)*(4. + run)*(-1.)/8.*D + ((-2.)*(1. + M2) + bra*M2)*(rho_tot + p_tot)*3./16.*pow(H,-2)*D*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*3./16.*pow(H,-2)*D + (D*bra_p + ((-2.) + bra)*((-3.)*bra*bra_p + (-1.)*kin_p))*1./8.*pow(H,-1)*pow(a,-1);
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_7_smg,
            &pvecback[pba->index_bg_lambda_7_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_7_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      

      
      pvecback[pba->index_bg_lambda_8_smg] = ((-2.) + bra)*(4. + run)*1./8.*D + 3./8.*(rho_tot + p_tot)*(((-9.)*bra + (-2.)*D*(1. + 2.*M2 - bra*M2))*(-1.)/2. + (rho_tot*(1. - M2) - (p_smg + rho_smg*M2))*9.*pow(H,-2)*pow(M2,-1))*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/8.*pow(H,-2)*D + (2.*(1. - M2) + bra*M2)*(rho_tot + p_tot)*(p_tot + p_smg)*27./16.*pow(H,-4)*pow(M2,-2) + ((-9.)*(rho_tot + p_tot) + (-6.)*bra*pow(H,2)*M2 + 3.*pow(bra,2)*pow(H,2)*M2 + (-1.)*pow(H,2)*D*M2)*1./8.*pow(H,-3)*pow(M2,-1)*bra_p*pow(a,-1) + ((-2.) + bra)*1./8.*pow(H,-1)*kin_p*pow(a,-1) + ((-2.) + bra)*9./16.*bra*pow(H,-3)*pow(M2,-1)*p_tot_p*pow(a,-1);
      
            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_8_smg,
            &pvecback[pba->index_bg_lambda_8_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_8_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      

      
      pvecback[pba->index_bg_cs2num_smg] = ((-2.) + bra)*((-1.)*bra + (-2.)*run + 2.*ten + (-1.)*bra*ten)*1./2. + pvecback[pba->index_bg_lambda_2_smg];

            memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_cs2num_smg,
            &pvecback[pba->index_bg_cs2num_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_cs2num_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");       


      pvecback[pba->index_bg_cs2_smg] = pvecback[pba->index_bg_cs2num_smg]/D;

      memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_cs2_smg,
            &pvecback[pba->index_bg_cs2_smg],
            1*sizeof(double));
      class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_cs2_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");      


      /* Here we update the minimum values of the stability quantities
       * test will be performed based on the lowest values
       */
      if (a > 0. && pba->parameters_tuned_smg == _TRUE_){
	if (D < pba->min_D_smg){
	  pba->min_D_smg = D;
	}
	
	if (pvecback[pba->index_bg_cs2_smg] <= pba->min_cs2_smg){
	  pba->min_cs2_smg = pvecback[pba->index_bg_cs2_smg];
	}
      
	if (pvecback[pba->index_bg_M2_smg] < pba->min_M2_smg){
	  pba->min_M2_smg = pvecback[pba->index_bg_M2_smg];
	}
      
	if (pvecback[pba->index_bg_tensor_excess_smg] + 1. < pba->min_ct2_smg){
	  pba->min_ct2_smg = 1. + pvecback[pba->index_bg_tensor_excess_smg];
	}
      }
     
    }//end of has_smg   

  }

  /* Horndeski stability tests
  * only if not overriden
  * and model is tuned!
  * TODO: commented out till properly checked!
  */
  if ((pba->has_smg == _TRUE_) && 
      (pba->parameters_tuned_smg == _TRUE_) &&
      (pba->skip_stability_tests_smg == _FALSE_)){
    
    class_test(pba->min_D_smg <= -fabs(pba->D_safe_smg),
        pba->error_message,
        "Ghost instability for scalar field perturbations with minimum D=%g \n",pba->min_D_smg);
    class_test(pba->min_cs2_smg < -fabs(pba->cs2_safe_smg),
        pba->error_message,
        "Gradient instability for scalar field perturbations with minimum c_s^2=%g \n",pba->min_cs2_smg);
    class_test(pba->min_M2_smg < -fabs(pba->M2_safe_smg),
        pba->error_message,
        "Ghost instability for metric tensor perturbations with minimum M*^2=%g \n",pba->min_M2_smg);
    class_test(pba->min_ct2_smg < -fabs(pba->ct2_safe_smg),
        pba->error_message,
        "Gradient instability for metric tensor perturbations with minimum c_t^2=%g \n",pba->min_ct2_smg);
     
   } 

   /* Yet another (third!) loop to make sure the background table makes sense
   */
  for (i=0; i < pba->bt_size; i++) {
      
     class_call(background_at_tau(pba,
				   pba->tau_table[i],
				   pba->long_info,
				   pba->inter_normal,
				   &last_index, //should be no problem to use the same one as for the derivatives
				   pvecback),
		 pba->error_message,
		 pba->error_message);
     
     // check if any of the values becomes nan
     int j = 0;
        while (j < pba->bg_size){
	  class_test( isnan(pvecback[j]) && (pba->parameters_tuned_smg == _TRUE_),
               pba->error_message,
               "pvecback[%i] = %e at a = %e in background!",j,pvecback[j],pvecback[pba->index_bg_a]);
	 j++; 
	}
  
  }
  
  free(pvecback_derivs);  //free the structure  
  

  /** - compute remaining "related parameters" */

  /** -> so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */
  pba->Neff = (pba->background_table[pba->index_bg_Omega_r]
               *pba->background_table[pba->index_bg_rho_crit]
               -pba->background_table[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pba->background_table[pba->index_bg_rho_g]);

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);
  }

  if (pba->background_verbose > 2) {
    if ((pba->has_dcdm == _TRUE_)&&(pba->has_dr == _TRUE_)){
      printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
      printf("     -> Omega0_dcdm = %f\n",pba->Omega0_dcdm);
      printf("     -> Omega0_dr = %f\n",pba->Omega0_dr);
      printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n",
             pba->Omega0_dr+pba->Omega0_dcdm,pba->Omega0_dcdmdr);
      printf("     -> Omega_ini_dcdm/Omega_b = %f\n",pba->Omega_ini_dcdm/pba->Omega0_b);
    }
    if (pba->has_smg == _TRUE_){
      printf(" -> Omega_smg = %f, whished %f ",pvecback[pba->index_bg_rho_smg]/pvecback[pba->index_bg_rho_crit], pba->Omega0_smg); 
      if(pba->has_lambda == _TRUE_)
	printf(", Omega_Lambda = %f", pba->Omega0_lambda);
      printf("\n");
      if (pba->background_verbose > 3) {
	printf("Minimal stability values: cs2 = %g, ct2 = %g, D = %g, M2 = %g \n",pba->min_cs2_smg,pba->min_ct2_smg,pba->min_D_smg,pba->min_M2_smg);
      }
      background_gravity_parameters(pba);
	      
    }    
  }

  free(pvecback);
  free(pvecback_integration);
  
  return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 *
 * @param ppr                  Input : pointer to precision structure
 * @param pba                  Input : pointer to background structure
 * @param pvecback             Input : vector of background quantitites used as workspace
 * @param pvecback_integration Output : vector of background quantitites to be integrated, returned with proper initial values
 * @return the error status
 */

int background_initial_conditions(
                                  struct precision *ppr,
                                  struct background *pba,
                                  double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
                                  double * pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
                                  ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  double rho_ncdm, p_ncdm, rho_ncdm_rel_tot=0.;
  double f,Omega_rad, rho_rad;
  int counter,is_early_enough,n_ncdm;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today;

  /** If we have ncdm species, perhaps we need to start earlier
      than the standard value for the species to be relativistic.
      This could happen for some WDM models.
  */

  if (pba->has_ncdm == _TRUE_) {

    for (counter=0; counter < _MAX_IT_; counter++) {

      is_early_enough = _TRUE_;
      rho_ncdm_rel_tot = 0.;

      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

	class_call(background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
					   pba->w_ncdm_bg[n_ncdm],
					   pba->q_size_ncdm_bg[n_ncdm],
					   pba->M_ncdm[n_ncdm],
					   pba->factor_ncdm[n_ncdm],
					   pba->a_today/a-1.0,
					   NULL,
					   &rho_ncdm,
					   &p_ncdm,
					   NULL,
					   NULL),
                   pba->error_message,
                   pba->error_message);
	rho_ncdm_rel_tot += 3.*p_ncdm;
	if (fabs(p_ncdm/rho_ncdm-1./3.)>ppr->tol_ncdm_initial_w)
	  is_early_enough = _FALSE_;
      }
      if (is_early_enough == _TRUE_)
	break;
      else
	a *= _SCALE_BACK_;
    }
    class_test(counter == _MAX_IT_,
	       pba->error_message,
	       "Search for initial scale factor a such that all ncdm species are relativistic failed.");
  }

  pvecback_integration[pba->index_bi_a] = a;

  /* Set initial values of {B} variables: */
  Omega_rad = pba->Omega0_g;
  if (pba->has_ur == _TRUE_)
    Omega_rad += pba->Omega0_ur;
  rho_rad = Omega_rad*pow(pba->H0,2)/pow(a/pba->a_today,4);
  if (pba->has_ncdm == _TRUE_){
    /** We must add the relativistic contribution from NCDM species: */
    rho_rad += rho_ncdm_rel_tot;
  }
  if (pba->has_dcdm == _TRUE_){
    /* Remember that the critical density today in CLASS conventions is H0^2 */
    pvecback_integration[pba->index_bi_rho_dcdm] =
      pba->Omega_ini_dcdm*pba->H0*pba->H0*pow(pba->a_today/a,3);
    if (pba->background_verbose > 3)
      printf("Density is %g. a_today=%g. Omega_ini=%g\n",pvecback_integration[pba->index_bi_rho_dcdm],pba->a_today,pba->Omega_ini_dcdm);
  }

  if (pba->has_dr == _TRUE_){
    if (pba->has_dcdm == _TRUE_){
      /** f is the critical density fraction of DR. The exact solution is
	  f = -Omega_rad+pow(pow(Omega_rad,3./2.)+0.5*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3),2./3.);
	  but it is not numerically stable for very small f which is always the case.
	  Instead we use the Taylor expansion of this equation, which is equivalent to
	  ignoring f(a) in the Hubble rate.
      */
      f = 1./3.*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3)/sqrt(Omega_rad);
      pvecback_integration[pba->index_bi_rho_dr] = f*pba->H0*pba->H0/pow(a/pba->a_today,4);
    }
    else{
      /** This is reserved for a future case where dr is not sourced by dcdm */
      pvecback_integration[pba->index_bi_rho_dr] = 0.0;
    }
  }

 
  
/** - fix initial value of modified gravity
   * run over all possible model cases
   */
  if(pba->has_smg == _TRUE_){
    
    pba->initial_conditions_set_smg = _FALSE_;
    
    //default value, can override latter
    if (pba->M_pl_evolution_smg ==_TRUE_){
      pvecback_integration[pba->index_bi_M_pl_smg] = 1.;
    }
       
    switch (pba->gravity_model_smg) {
	
      case propto_omega:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[4];
	break;	
      
      case propto_scale:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[4];
	break;	

      case eft_alphas_power_law:
	pvecback_integration[pba->index_bi_M_pl_smg] = 1. + pba->parameters_2_smg[0]*pow(a, pba->parameters_2_smg[4]);
	break;

      case eft_gammas_power_law:
	pvecback_integration[pba->index_bi_M_pl_smg] = 1. + pba->parameters_2_smg[0]*pow(a,pba->parameters_2_smg[4]) + pba->parameters_2_smg[3]*pow(a,pba->parameters_2_smg[7]);
	break;
      case eft_gammas_exponential:
	pvecback_integration[pba->index_bi_M_pl_smg] = exp(pba->parameters_2_smg[0]*pow(a,pba->parameters_2_smg[4])) + exp(pba->parameters_2_smg[3]*pow(a,pba->parameters_2_smg[7])) -1.;
	break;

      case constant:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[4];
	break;

      case tanh_scale:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[0];
	break;

      case hill_scale:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[0];
	break;

      case oscillation_scale:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[0];
	break;

      case gaussian_scale:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[0];
	break;

      case scale_early_late:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[0];
	break;

      case arctan:
	pvecback_integration[pba->index_bi_M_pl_smg] = pba->parameters_2_smg[0];
	break;

    }
    
      if (pba->M_pl_evolution_smg == _TRUE_)
	if (pba->background_verbose>3) 
	  printf(" -> Initial conditions: M_pl = %e \n",pvecback_integration[pba->index_bi_M_pl_smg]);    
  }  

  /* Infer pvecback from pvecback_integration */
  class_call(background_functions(pba, pvecback_integration, pba->long_info, pvecback),
	     pba->error_message,
	     pba->error_message);

  /* Just checking that our initial time indeed is deep enough in the radiation
     dominated regime */
  class_test(fabs(pvecback[pba->index_bg_Omega_r]-1.) > ppr->tol_initial_Omega_r,
	     pba->error_message,
	     "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to start from radiation domination.",
	     pvecback[pba->index_bg_Omega_r]);

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[pba->index_bg_H] <= 0.,
             pba->error_message,
             "H = %e instead of strictly positive",pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_time] = 1./(2.* pvecback[pba->index_bg_H]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[pba->index_bi_tau] = 1./(a * pvecback[pba->index_bg_H]);

  /** - compute initial sound horizon, assuming c_s=1/sqrt(3) initially */
  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_tau]/sqrt(3.);

  /** - compute initial value of the integral over dtau/(aH^2),
      assumed to be proportional to a^4 during RD, but with arbitrary
      normalization */
  pvecback_integration[pba->index_bi_growth] = 1./(4.*a*a*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]);

  return _SUCCESS_;

}

/**
 * Subroutine for formatting background output
 */

int background_output_titles(struct background * pba,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){

  /** Length of the columntitle should be less than _OUTPUTPRECISION_+6
      to be indented correctly, but it can be as long as . */
  int n;
  char tmp[20];

  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"proper time [Gyr]",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"H [1/Mpc]",_TRUE_);
  class_store_columntitle(titles,"comov. dist.",_TRUE_);
  class_store_columntitle(titles,"ang.diam.dist.",_TRUE_);
  class_store_columntitle(titles,"lum. dist.",_TRUE_);
  class_store_columntitle(titles,"comov.snd.hrz.",_TRUE_);
  class_store_columntitle(titles,"(.)rho_g",_TRUE_);
  class_store_columntitle(titles,"(.)rho_b",_TRUE_);
  class_store_columntitle(titles,"(.)rho_cdm",pba->has_cdm);
  if (pba->has_ncdm == _TRUE_){
    for (n=0; n<pba->N_ncdm; n++){
      sprintf(tmp,"(.)rho_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
    }
  }
  class_store_columntitle(titles,"(.)rho_lambda",pba->has_lambda);
  class_store_columntitle(titles,"(.)rho_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)rho_ur",pba->has_ur);
  class_store_columntitle(titles,"(.)rho_crit",_TRUE_);
  class_store_columntitle(titles,"(.)rho_dcdm",pba->has_dcdm);
  class_store_columntitle(titles,"(.)rho_dr",pba->has_dr);

  class_store_columntitle(titles,"gr.fac. D",_TRUE_);
  class_store_columntitle(titles,"gr.fac. f",_TRUE_);
  
  class_store_columntitle(titles,"(.)rho_smg",pba->has_smg);   
  class_store_columntitle(titles,"(.)p_smg",pba->has_smg);   
  class_store_columntitle(titles,"M*^2_smg",pba->has_smg);       
  class_store_columntitle(titles,"kineticity_smg",pba->has_smg);   
  class_store_columntitle(titles,"braiding_smg",pba->has_smg);
  class_store_columntitle(titles,"tensor_excess_smg",pba->has_smg);   
  class_store_columntitle(titles,"Mpl_running_smg",pba->has_smg);       
  class_store_columntitle(titles,"c_s^2",pba->has_smg);        
  class_store_columntitle(titles,"kin (D)",pba->has_smg);        

  return _SUCCESS_;
}

int background_output_data(
                           struct background *pba,
                           int number_of_titles,
                           double *data){
  int index_tau, storeidx, n;
  double *dataptr, *pvecback;

  /** Store quantities: */
  for (index_tau=0; index_tau<pba->bt_size; index_tau++){
    dataptr = data + index_tau*number_of_titles;
    pvecback = pba->background_table + index_tau*pba->bg_size;
    storeidx = 0;

    class_store_double(dataptr,pba->a_today/pvecback[pba->index_bg_a]-1.,_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_time]/_Gyr_over_Mpc_,_TRUE_,storeidx);
    class_store_double(dataptr,pba->conformal_age-pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_H],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ang_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_lum_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rs],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_g],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_b],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_cdm],pba->has_cdm,storeidx);
    if (pba->has_ncdm == _TRUE_){
      for (n=0; n<pba->N_ncdm; n++)
        class_store_double(dataptr,pvecback[pba->index_bg_rho_ncdm1+n],_TRUE_,storeidx);
    }
    class_store_double(dataptr,pvecback[pba->index_bg_rho_lambda],pba->has_lambda,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_ur],pba->has_ur,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_crit],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dcdm],pba->has_dcdm,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dr],pba->has_dr,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_D],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_f],_TRUE_,storeidx);
    
    class_store_double(dataptr,pvecback[pba->index_bg_rho_smg],pba->has_smg,storeidx);  
    class_store_double(dataptr,pvecback[pba->index_bg_p_smg],pba->has_smg,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_M2_smg],pba->has_smg,storeidx); 
    class_store_double(dataptr,pvecback[pba->index_bg_kineticity_smg],pba->has_smg,storeidx);   
    class_store_double(dataptr,pvecback[pba->index_bg_braiding_smg],pba->has_smg,storeidx);   
    class_store_double(dataptr,pvecback[pba->index_bg_tensor_excess_smg],pba->has_smg,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_mpl_running_smg],pba->has_smg,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_cs2_smg],pba->has_smg,storeidx);     
    class_store_double(dataptr,pvecback[pba->index_bg_kinetic_D_smg],pba->has_smg,storeidx);
    
  }

  return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which are passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and wokspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input : conformal time
 * @param y                        Input : vector of variable
 * @param dy                       Output : its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output : error message
 */
int background_derivs(
                      double tau,
                      double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
                      double* dy, /* vector with argument dy[index_bi]
                                     (must be already allocated with
                                     size pba->bi_size) */
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      ) {

  /** Summary: */

  /** - define local variables */

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - Calculates functions of /f$ a /f$ with background_functions() */
  class_call(background_functions(pba, y, pba->normal_info, pvecback),
             pba->error_message,
             error_message);

  /** - calculate /f$ a'=a^2 H /f$ */
  dy[pba->index_bi_a] = y[pba->index_bi_a] * y[pba->index_bi_a] * pvecback[pba->index_bg_H];

  /** - calculate /f$ t' = a /f$ */
  dy[pba->index_bi_time] = y[pba->index_bi_a];

  class_test(pvecback[pba->index_bg_rho_g] <= 0.,
             error_message,
             "rho_g = %e instead of strictly positive",pvecback[pba->index_bg_rho_g]);

  /** - calculate rs' = c_s */
  dy[pba->index_bi_rs] = 1./sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]))*sqrt(1.-pba->K*y[pba->index_bi_rs]*y[pba->index_bi_rs]); // TBC: curvature correction

  /** calculate growth' = 1/(aH^2) */
  dy[pba->index_bi_growth] = 1./(y[pba->index_bi_a] * pvecback[pba->index_bg_H] * pvecback[pba->index_bg_H]);

  if (pba->has_dcdm == _TRUE_){
    /** compute dcdm density rho' = -3aH rho - a Gamma rho*/
    dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dcdm]-
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if (pba->has_dcdm == _TRUE_){
    /** compute dcdm density rho' = -3aH rho - a Gamma rho*/
    dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dcdm]-
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if ((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_)){
    /** Compute dr density rho' = -4aH rho - a Gamma rho*/
    dy[pba->index_bi_rho_dr] = -4.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dr]+
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }
  
  /** - Scalar field equation: \f$ \phi'' + 2 a H \phi' + a^2 dV = 0 \f$  (note H is wrt cosmic time)**/
  if (pba->has_smg == _TRUE_){
    /** - Planck mass equation (if parameterization in terms of alpha_m **/
    if (pba->M_pl_evolution_smg == _TRUE_)
      dy[pba->index_bi_M_pl_smg] = y[pba->index_bi_a]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_mpl_running_smg]*y[pba->index_bi_M_pl_smg];   //in this case the running has to be integrated (eq 3.3 of 1404.3713 yields M2' = aH\alpha_M)
  }  


  return _SUCCESS_;

}


/** This function fills the modified gravity part of background_functions.
 * First all the Horndeski functions G_i(X,phi) are computed. 
 * A loop is used to allow different implementations without erasing previous ones
 * Note that in CLASS units a canonical field has G2 = X ~ [Mpc]^-2
 * This implies that phi_prime ~ [Mpc]^-1 
 * and phi is dimensionless (we can think of phi as given in units of the Planck mass
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added. 
 * Numerical derivatives may further serve as a consistency check.
 * TODO: Add a background_write_alpha_primes
 */

int background_gravity_functions(
				 struct background *pba,
				 double * pvecback_B,
				 short return_format,
				 double * pvecback
				 ){
  
      // scalar field + curvature not yet implemented
  class_test(pba->K !=0 ,
	     pba->error_message,
	     "has_smg with curvature K = %e not yet implemented",pba->K);  

  if (pba->field_evolution_smg == _TRUE_) {
  
    class_test(pba->field_evolution_smg == _TRUE_,
	       pba->error_message,
	       "self consistent evolution not implemented in this version. You should choose a parameterization\n");
      
  }// end of if pba->field_evolution_smg
  else{
    
    //TODO: this structure will need revision if we have more complicated parameterizations
    //TODO: separate parameterization from background evolution
    double a, M_pl;
    double rho_tot, p_tot;
    double Omega_smg;
    
    a = pvecback_B[pba->index_bi_a];
    M_pl = pvecback_B[pba->index_bi_M_pl_smg];    
    
    rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
    p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
    
    //initialize the values to the defaults
    pvecback[pba->index_bg_kineticity_smg] = 0;
    pvecback[pba->index_bg_braiding_smg] = 0.;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_M2_smg] = 1.;
    pvecback[pba->index_bg_mpl_running_smg] = 0.;
    
     if (pba->expansion_model_smg == lcdm){
      
      double Omega_const_smg = pba->parameters_smg[0];
      
      pvecback[pba->index_bg_rho_smg] = Omega_const_smg*pow(pba->H0,2);
      pvecback[pba->index_bg_p_smg] = -Omega_const_smg*pow(pba->H0,2);
    }
    
    if (pba->expansion_model_smg == wowa){
      
      double Omega_const_smg = pba->parameters_smg[0];
      double w0 = pba->parameters_smg[1];
      double wa = pba->parameters_smg[2];
      
      pvecback[pba->index_bg_rho_smg] = Omega_const_smg * pow(pba->H0,2)/pow(a,3.*(1. + w0 + wa)) * exp(3.*wa*(a-1.));
      pvecback[pba->index_bg_p_smg] = (w0+(1-a)*wa) * Omega_const_smg * pow(pba->H0,2)/pow(a,3.*(1.+w0+wa)) * exp(3.*wa*(a-1.));
    }

    rho_tot += pvecback[pba->index_bg_rho_smg];
    p_tot += pvecback[pba->index_bg_p_smg];    

    
    Omega_smg = pvecback[pba->index_bg_rho_smg]/rho_tot; //used for some parameterizations
    

    if (pba->gravity_model_smg == propto_omega) {	
      
      double c_k = pba->parameters_2_smg[0];
      double c_b = pba->parameters_2_smg[1];
      double c_m = pba->parameters_2_smg[2];
      double c_t = pba->parameters_2_smg[3];      
      
      pvecback[pba->index_bg_kineticity_smg] = c_k*Omega_smg;
      pvecback[pba->index_bg_braiding_smg] = c_b*Omega_smg;	
      pvecback[pba->index_bg_tensor_excess_smg] = c_t*Omega_smg;
      pvecback[pba->index_bg_mpl_running_smg] = c_m*Omega_smg;
      pvecback[pba->index_bg_M2_smg] = M_pl;
    }
    else if (pba->gravity_model_smg == propto_scale) {	
      
      double c_k = pba->parameters_2_smg[0];
      double c_b = pba->parameters_2_smg[1];
      double c_m = pba->parameters_2_smg[2];
      double c_t = pba->parameters_2_smg[3];      
      
      pvecback[pba->index_bg_kineticity_smg] = c_k*a;
      pvecback[pba->index_bg_braiding_smg] = c_b*a;	
      pvecback[pba->index_bg_tensor_excess_smg] = c_t*a;
      pvecback[pba->index_bg_mpl_running_smg] = c_m*a;
      pvecback[pba->index_bg_M2_smg] = M_pl;
    }
    else if (pba->gravity_model_smg == eft_alphas_power_law) {	
      
      double M_0 = pba->parameters_2_smg[0];
      double c_k = pba->parameters_2_smg[1];
      double c_b = pba->parameters_2_smg[2];
      double c_t = pba->parameters_2_smg[3];      
      double M_0_exp = pba->parameters_2_smg[4];
      double c_k_exp = pba->parameters_2_smg[5];
      double c_b_exp = pba->parameters_2_smg[6];
      double c_t_exp = pba->parameters_2_smg[7];      
      
      pvecback[pba->index_bg_kineticity_smg] = c_k*pow(a, c_k_exp);
      pvecback[pba->index_bg_braiding_smg] = c_b*pow(a, c_b_exp);	
      pvecback[pba->index_bg_tensor_excess_smg] = c_t*pow(a, c_t_exp);
      pvecback[pba->index_bg_mpl_running_smg] = M_0*M_0_exp*pow(a, M_0_exp)/(1. + M_0*pow(a, M_0_exp));
      pvecback[pba->index_bg_M2_smg] = M_pl;
    }
    else if ((pba->gravity_model_smg == eft_gammas_power_law) || (pba->gravity_model_smg == eft_gammas_exponential)) {
      
      double Omega=0., g1=0., g2=0., g3=0., Omega_p=0., Omega_pp=0., g3_p=0.;
      double Omega_0=0., g1_0=0., g2_0=0., g3_0=0.;
      double Omega_exp=0., g1_exp=0., g2_exp=0., g3_exp=0.;
      
      Omega_0 = pba->parameters_2_smg[0];
      g1_0 = pba->parameters_2_smg[1];
      g2_0 = pba->parameters_2_smg[2];
      g3_0 = pba->parameters_2_smg[3];
      Omega_exp = pba->parameters_2_smg[4];
      g1_exp = pba->parameters_2_smg[5];
      g2_exp = pba->parameters_2_smg[6];
      g3_exp = pba->parameters_2_smg[7];
      
      if (pba->gravity_model_smg == eft_gammas_power_law) {
        Omega = Omega_0*pow(a,Omega_exp);
        Omega_p = Omega_0*Omega_exp*pow(a,Omega_exp-1.); // Derivative w.r.t. the scale factor
        Omega_pp = Omega_0*Omega_exp*(Omega_exp-1.)*pow(a,Omega_exp-2.); // Derivative w.r.t. the scale factor
        g1 = g1_0*pow(a,g1_exp);
        g2 = g2_0*pow(a,g2_exp);
        g3 = g3_0*pow(a,g3_exp);
        g3_p = g3_0*g3_exp*pow(a,g3_exp-1.); // Derivative w.r.t. the scale factor
      }
      else { //(pba->gravity_model_smg == eft_gammas_exponential)
        Omega = exp(Omega_0*pow(a,Omega_exp))-1.;
        Omega_p = Omega_0*Omega_exp*pow(a,Omega_exp-1.)*exp(Omega_0*pow(a,Omega_exp)); // Derivative w.r.t. the scale factor
        Omega_pp = Omega_0*Omega_exp*pow(a,Omega_exp-2.)*exp(Omega_0*pow(a,Omega_exp))*(Omega_exp-1.+Omega_0*Omega_exp*pow(a,Omega_exp)); // Derivative w.r.t. the scale factor
        g1 = exp(g1_0*pow(a,g1_exp))-1.;
        g2 = exp(g2_0*pow(a,g2_exp))-1.;
        g3 = exp(g3_0*pow(a,g3_exp))-1.;
        g3_p = g3_0*g3_exp*pow(a,g3_exp-1.)*exp(g3_0*pow(a,g3_exp)); // Derivative w.r.t. the scale factor
      }
      
      
      double c_over_H2 = (-pow(a,2.)*Omega_pp + 3.*(Omega + a*Omega_p/2.)*(rho_tot + p_tot)/rho_tot + 3.*(pvecback[pba->index_bg_rho_smg]+pvecback[pba->index_bg_p_smg])/rho_tot)/2.;
      
      pvecback[pba->index_bg_kineticity_smg] = 2.*(2.*g1*pow(pba->H0,2.)/rho_tot + c_over_H2)/(1. + Omega + g3);
      pvecback[pba->index_bg_braiding_smg] = -(g2*pba->H0/sqrt(rho_tot) + a*Omega_p)/(1. + Omega + g3);
      pvecback[pba->index_bg_tensor_excess_smg] = -g3/(1. + Omega + g3);
      pvecback[pba->index_bg_mpl_running_smg] = a*(Omega_p + g3_p)/(1. + Omega + g3);
      pvecback[pba->index_bg_M2_smg] = M_pl;

    }

    else if (pba->gravity_model_smg == constant) {	
      
      double c_k = pba->parameters_2_smg[0];
      double c_b = pba->parameters_2_smg[1];
      double c_m = pba->parameters_2_smg[2];
      double c_t = pba->parameters_2_smg[3];      
      
      pvecback[pba->index_bg_kineticity_smg] = c_k;
      pvecback[pba->index_bg_braiding_smg] = c_b;	
      pvecback[pba->index_bg_tensor_excess_smg] = c_t;
      pvecback[pba->index_bg_mpl_running_smg] = c_m;
      pvecback[pba->index_bg_M2_smg] = M_pl;

    }

    else if (pba->gravity_model_smg == tanh_scale) {	
      
      double c_k_early = pba->parameters_2_smg[1];
      double c_k_late = pba->parameters_2_smg[2];
      double c_b_early = pba->parameters_2_smg[3];
      double c_b_late = pba->parameters_2_smg[4];
      double c_m_early = pba->parameters_2_smg[5];
      double c_m_late = pba->parameters_2_smg[6];
      double c_t_early = pba->parameters_2_smg[7];    
      double c_t_late = pba->parameters_2_smg[8];
      double a_trans = pba->parameters_2_smg[9];
      double da = pba->parameters_2_smg[10];
      
      pvecback[pba->index_bg_kineticity_smg] = c_k_late + (1./2.)*(c_k_early - c_k_late)  * (1 - tanh((a - a_trans)/(2.* da))) ;
      pvecback[pba->index_bg_braiding_smg] = c_b_late + (1./2.)*(c_b_early - c_b_late)   * (1 - tanh((a - a_trans)/(2.* da))) ;    
      pvecback[pba->index_bg_tensor_excess_smg] = c_t_late + (1./2.)*(c_t_early - c_t_late)  * (1 - tanh((a - a_trans)/(2.* da))) ;
      pvecback[pba->index_bg_mpl_running_smg] = c_m_late + (1./2.)*(c_m_early - c_m_late)  * (1 - tanh((a - a_trans)/(2.* da))) ;
      pvecback[pba->index_bg_M2_smg] = M_pl;

      /*pvecback[pba->index_bg_kineticity_smg] = c_k_late + (1./2.)*(c_k_early - c_k_late)  * (1 + tanh((a - a_trans)/(2.* da))) ;
      pvecback[pba->index_bg_braiding_smg] = c_b_late + (1./2.)*(c_b_early - c_b_late)   * (1 + tanh((a - a_trans)/(2.* da))) ;    
      pvecback[pba->index_bg_tensor_excess_smg] = c_t_late + (1./2.)*(c_t_early - c_t_late)  * (1 + tanh((a - a_trans)/(2.* da))) ;
      pvecback[pba->index_bg_mpl_running_smg] = c_m_late + (1./2.)*(c_m_early - c_m_late)  * (1 + tanh((a - a_trans)/(2.* da))) ;
      pvecback[pba->index_bg_M2_smg] = M_pl;*/

      /*pvecback[pba->index_bg_kineticity_smg] = ( c_k_late + (1./2.)*(c_k_early - c_k_late)  * (1 - tanh((a - a_trans)/(2.* da))) )*(1.-Omega_smg) ;
      pvecback[pba->index_bg_braiding_smg] = ( c_b_late + (1./2.)*(c_b_early - c_b_late)   * (1 - tanh((a - a_trans)/(2.* da))) )*(1.-Omega_smg);    
      pvecback[pba->index_bg_tensor_excess_smg] = ( c_t_late + (1./2.)*(c_t_early - c_t_late)  * (1 - tanh((a - a_trans)/(2.* da))) )*(1.-Omega_smg);
      pvecback[pba->index_bg_mpl_running_smg] = ( c_m_late + (1./2.)*(c_m_early - c_m_late)  * (1 - tanh((a - a_trans)/(2.* da))) )*(1.-Omega_smg);
      pvecback[pba->index_bg_M2_smg] = M_pl;*/

      /* pvecback[pba->index_bg_kineticity_smg] = c_k_late + (1./2.)*(c_k_early - c_k_late)  * (1 - exp(2.*(a-a_trans)/(2.*da))) / (1 - exp(2.*(a-a_trans)/(2.*da)));
      pvecback[pba->index_bg_braiding_smg] = c_b_late + (1./2.)*(c_b_early - c_b_late)   * (1 - exp(2.*(a-a_trans)/(2.*da))) / (1 - exp(2.*(a-a_trans)/(2.*da)));    
      pvecback[pba->index_bg_tensor_excess_smg] = c_t_late + (1./2.)*(c_t_early - c_t_late)  * (1 - exp(2.*(a-a_trans)/(2.*da))) / (1 - exp(2.*(a-a_trans)/(2.*da)));
      pvecback[pba->index_bg_mpl_running_smg] = c_m_late + (1./2.)*(c_m_early - c_m_late)  * (1 - exp(2.*(a-a_trans)/(2.*da))) / (1 - exp(2.*(a-a_trans)/(2.*da)));
      pvecback[pba->index_bg_M2_smg] = M_pl; */
     
    }

    else if (pba->gravity_model_smg == hill_scale) {	
      
      double c_k = pba->parameters_2_smg[1];
      double c_b = pba->parameters_2_smg[2];
      double c_m = pba->parameters_2_smg[3];
      double c_t = pba->parameters_2_smg[4];    
      double a_trans = pba->parameters_2_smg[5];
      double a_rapidity = pba->parameters_2_smg[6];
      
      pvecback[pba->index_bg_kineticity_smg] = c_k * (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.);
      pvecback[pba->index_bg_braiding_smg] =  c_b * (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.);    
      pvecback[pba->index_bg_tensor_excess_smg] =  c_t * (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.);
      pvecback[pba->index_bg_mpl_running_smg] =  c_m * (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.);
      pvecback[pba->index_bg_M2_smg] = M_pl;
     
    }

    else if (pba->gravity_model_smg == oscillation_scale) {	
      
      double c_k = pba->parameters_2_smg[1];
      double c_b = pba->parameters_2_smg[2];
      double c_m = pba->parameters_2_smg[3];
      double c_t = pba->parameters_2_smg[4];    
      double a_trans = pba->parameters_2_smg[5];
      double da = pba->parameters_2_smg[6];
      double a_rapidity = 3./2.;
      
      pvecback[pba->index_bg_kineticity_smg] = c_k * ( 0.5*(1./0.4) * sin((a - a_trans) / da) * exp(-pow( ((a - a_trans) / da), 2)) + (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.));
      pvecback[pba->index_bg_braiding_smg] =  c_b * ( 0.5*(1./0.4) * sin((a - a_trans) / da) * exp(-pow( ((a - a_trans) / da), 2)) + (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.));    
      pvecback[pba->index_bg_tensor_excess_smg] =  c_t * ( 0.5*(1./0.4) * sin((a - a_trans) / da) * exp(-pow( ((a - a_trans) / da), 2)) + (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.));
      pvecback[pba->index_bg_mpl_running_smg] =  c_m * ( 0.5*(1./0.4) * sin((a - a_trans) / da) * exp(-pow( ((a - a_trans) / da), 2)) + (4. * pow((a / a_trans), a_rapidity)) / pow((pow((a / a_trans), a_rapidity) + 1.), 2.));
      pvecback[pba->index_bg_M2_smg] = M_pl;
     
    }

    else if (pba->gravity_model_smg == gaussian_scale) {	
      
      double c_k = pba->parameters_2_smg[1];
      double c_b = pba->parameters_2_smg[2];
      double c_m = pba->parameters_2_smg[3];
      double c_t = pba->parameters_2_smg[4];    
      double a_trans = pba->parameters_2_smg[5];
      double da = pba->parameters_2_smg[6];
      
      pvecback[pba->index_bg_kineticity_smg] = c_k * (1./sqrt(2.*3.14159265*da)) * exp(-(pow((a-a_trans),2.))/(2.*da));
      pvecback[pba->index_bg_braiding_smg] =  c_b * (1./sqrt(2.*3.14159265*da)) * exp(-(pow((a-a_trans),2.))/(2.*da));    
      pvecback[pba->index_bg_tensor_excess_smg] =  c_t * (1./sqrt(2.*3.14159265*da)) * exp(-(pow((a-a_trans),2.))/(2.*da));
      pvecback[pba->index_bg_mpl_running_smg] =  c_m * (1./sqrt(2.*3.14159265*da)) * exp(-(pow((a-a_trans),2.))/(2.*da));
      pvecback[pba->index_bg_M2_smg] = M_pl;
     
    }

    else if (pba->gravity_model_smg == scale_early_late) {	
      
      double c_k_early = pba->parameters_2_smg[1];
      double c_k_late = pba->parameters_2_smg[2];
      double c_b_early = pba->parameters_2_smg[3];
      double c_b_late = pba->parameters_2_smg[4];
      double c_m_early = pba->parameters_2_smg[5];
      double c_m_late = pba->parameters_2_smg[6];
      double c_t_early = pba->parameters_2_smg[7];    
      double c_t_late = pba->parameters_2_smg[8];
      double a_trans = pba->parameters_2_smg[9];
      double da = pba->parameters_2_smg[10];
      
      pvecback[pba->index_bg_kineticity_smg] = (1./2.) * (c_k_early + c_k_late) + (c_k_late - c_k_early) * (1./3.14159265) * ((a - a_trans) / da);
      pvecback[pba->index_bg_braiding_smg] = (1./2.) * (c_b_early + c_b_late) + (c_b_late - c_b_early) * (1./3.14159265) * ((a - a_trans) / da);	
      pvecback[pba->index_bg_tensor_excess_smg] = (1./2.) * (c_t_early + c_t_late) + (c_t_late - c_t_early) * (1./3.14159265) * ((a - a_trans) / da);
      pvecback[pba->index_bg_mpl_running_smg] = (1./2.) * (c_m_early + c_m_late) + (c_m_late - c_m_early) * (1./3.14159265) * ((a - a_trans) / da);
      pvecback[pba->index_bg_M2_smg] = M_pl;
     
    }

    else if (pba->gravity_model_smg == arctan) {	
      
      double c_k_early = pba->parameters_2_smg[1];
      double c_k_late = pba->parameters_2_smg[2];
      double c_b_early = pba->parameters_2_smg[3];
      double c_b_late = pba->parameters_2_smg[4];
      double c_m_early = pba->parameters_2_smg[5];
      double c_m_late = pba->parameters_2_smg[6];
      double c_t_early = pba->parameters_2_smg[7];    
      double c_t_late = pba->parameters_2_smg[8];
      double a_trans = pba->parameters_2_smg[9];
      double da = pba->parameters_2_smg[10];
      
      pvecback[pba->index_bg_kineticity_smg] = (1./2.) * (c_k_early + c_k_late) + (c_k_late - c_k_early) * (1./3.14159265) * atan((a - a_trans) / da);
      pvecback[pba->index_bg_braiding_smg] = (1./2.) * (c_b_early + c_b_late) + (c_b_late - c_b_early) * (1./3.14159265) * atan((a - a_trans) / da);	
      pvecback[pba->index_bg_tensor_excess_smg] = (1./2.) * (c_t_early + c_t_late) + (c_t_late - c_t_early) * (1./3.14159265) * atan((a - a_trans) / da);
      pvecback[pba->index_bg_mpl_running_smg] = (1./2.) * (c_m_early + c_m_late) + (c_m_late - c_m_early) * (1./3.14159265) * atan((a - a_trans) / da);
      pvecback[pba->index_bg_M2_smg] = M_pl;
     
    }
    
    
    pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);
    /** - compute derivative of H with respect to conformal time */
    pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

  }//end of parameterized mode  
    
  // add a value to the kineticity to avoid problems with perturbations in certain models.
  // NOTE: this needs to be done here to avoid interfering with the equations
  pvecback[pba->index_bg_kineticity_smg] += pba->kineticity_safe_smg;
  
  // kinetic term
  pvecback[pba->index_bg_kinetic_D_smg] = 3./2.*pow(pvecback[pba->index_bg_braiding_smg],2) + pvecback[pba->index_bg_kineticity_smg];
  
  //Derivatives of the BS functions and others. Set to zero here and computed numerically once the background is integrated (needed so that debuggers don't complain).
      
  pvecback[pba->index_bg_kineticity_prime_smg] = 0.;
  pvecback[pba->index_bg_braiding_prime_smg] = 0.;
  pvecback[pba->index_bg_mpl_running_prime_smg] = 0.;
  pvecback[pba->index_bg_tensor_excess_prime_smg] = 0.;
  pvecback[pba->index_bg_H_prime_prime] = 0.;
  pvecback[pba->index_bg_p_tot_wo_prime_smg] = 0.;
  pvecback[pba->index_bg_p_prime_smg] = 0.;
  pvecback[pba->index_bg_cs2_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_prime_smg] = 0.;
  pvecback[pba->index_bg_cs2num_smg] = 0.;
  pvecback[pba->index_bg_lambda_1_smg] = 0.;
  pvecback[pba->index_bg_lambda_2_smg] = 0.;
  pvecback[pba->index_bg_lambda_3_smg] = 0.;
  pvecback[pba->index_bg_lambda_4_smg] = 0.;
  pvecback[pba->index_bg_lambda_5_smg] = 0.;
  pvecback[pba->index_bg_lambda_6_smg] = 0.;
  pvecback[pba->index_bg_lambda_7_smg] = 0.;
  pvecback[pba->index_bg_lambda_8_smg] = 0.;
  if (pba->field_evolution_smg == _FALSE_ && pba->M_pl_evolution_smg == _FALSE_){
    pvecback[pba->index_bg_mpl_running_smg] = 0.;
  }

  return _SUCCESS_;
  
}

int background_gravity_parameters(
				  struct background *pba
				  ){
 
  switch (pba->gravity_model_smg) {
     
   case propto_omega:
     printf("Modified gravity: propto_omega with parameters: \n");
     printf("-> c_K = %g, c_B = %g, c_M = %g, c_T = %g, M_*^2_init = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	    pba->parameters_2_smg[4]);
     break;
     
   case propto_scale:
     printf("Modified gravity: propto_scale with parameters: \n");
     printf("-> c_K = %g, c_B = %g, c_M = %g, c_T = %g, M_*^2_init = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	    pba->parameters_2_smg[4]);
     break;

   case eft_alphas_power_law:
     printf("Modified gravity: eft_alphas_power_law with parameters: \n");
     printf("-> M_*^2_0 = %g, c_K = %g, c_B = %g, c_T = %g, M_*^2_exp = %g, c_K_exp = %g, c_B_exp = %g, c_T_exp = %g\n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	    pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7]);
     break;

   case eft_gammas_power_law:
     printf("Modified gravity: eft_gammas_power_law with parameters: \n");
     printf("-> Omega_0 = %g, gamma_1 = %g, gamma_2 = %g, gamma_3 = %g, Omega_0_exp = %g, gamma_1_exp = %g, gamma_2_exp = %g, gamma_3_exp = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7]);
     break;

   case eft_gammas_exponential:
     printf("Modified gravity: eft_gammas_exponential with parameters: \n");
     printf("-> Omega_0 = %g, gamma_1 = %g, gamma_2 = %g, gamma_3 = %g, Omega_0_exp = %g, gamma_1_exp = %g, gamma_2_exp = %g, gamma_3_exp = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7]);
     break;

   case constant:
     printf("Modified gravity: constant with parameters: \n");
     printf("-> c_K = %g, c_B = %g, c_M = %g, c_T = %g, M_*^2_init = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	    pba->parameters_2_smg[4]);
     break;

   case tanh_scale:
     printf("Modified gravity: tanh with parameters: \n");
     printf("-> M_*^2_init, c_K_early = %g, c_K_late = %g, c_B_early = %g, c_B_late = %g, c_M_early = %g, c_M_late = %g, c_T_early = %g, c_T_early = %g, a_trans = %g, da = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7],pba->parameters_2_smg[8],pba->parameters_2_smg[9],pba->parameters_2_smg[10]);
     break;

   case hill_scale:
     printf("Modified gravity: tanh with parameters: \n");
     printf("-> M_*^2_init, c_K = %g, c_B = %g, c_M = %g, c_T = %g, a_trans = %g, a_rapidity = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6]);
     break;

   case oscillation_scale:
     printf("Modified gravity: tanh with parameters: \n");
     printf("-> M_*^2_init, c_K = %g, c_B = %g, c_M = %g, c_T = %g, a_trans = %g, da = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6]);
     break;

   case gaussian_scale:
     printf("Modified gravity: tanh with parameters: \n");
     printf("-> M_*^2_init, c_K = %g, c_B = %g, c_M = %g, c_T = %g, a_trans = %g, da = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6]);
     break;

   case scale_early_late:
     printf("Modified gravity: scale_early_late with parameters: \n");
     printf("-> M_*^2_init, c_K_early = %g, c_K_late = %g, c_B_early = %g, c_B_late = %g, c_M_early = %g, c_M_late = %g, c_T_early = %g, c_T_early = %g, a_trans = %g, da = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7],pba->parameters_2_smg[8],pba->parameters_2_smg[9],pba->parameters_2_smg[10]);
     break;

   case arctan:
     printf("Modified gravity: arctan with parameters: \n");
     printf("-> M_*^2_init, c_K_early = %g, c_K_late = %g, c_B_early = %g, c_B_late = %g, c_M_early = %g, c_M_late = %g, c_T_early = %g, c_T_early = %g, a_trans = %g, da = %g \n",
	    pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7],pba->parameters_2_smg[8],pba->parameters_2_smg[9],pba->parameters_2_smg[10]);
     break;

    
  }
  
  if(pba->field_evolution_smg==_FALSE_) {
    switch (pba->expansion_model_smg) {
      
    case lcdm:
      printf("Parameterized model with LCDM expansion \n");
      printf("-> Omega_smg = %f \n",pba->parameters_smg[0]);
      break;

      case wowa:
      printf("Parameterized model with CPL expansion \n");
      printf("-> Omega_smg = %f, w0 = %f, wa = %e \n",
	     pba->parameters_smg[0],pba->parameters_smg[1],pba->parameters_smg[2]);
      break;
    }
    
  }
  
  return _SUCCESS_;
  
}

/**
 * Scalar field potential and its derivatives with respect to the field _scf
 * For Albrecht & Skordis model: 9908085
 * \f$ V = V_p_scf*V_e_scf \f$
 * \f$ V_e =  \exp(-\lambda \phi) (exponential) \f$
 * \f$ V_p = (\phi - B)^\alpha + A (polynomial bump) \f$
 * TODO: -Add some functionality to include different models/potentials (tuning would be difficult, though)
 * - Generalize to Kessence/Horndeski/PPF and/or couplings
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added.
 * Numerical derivatives may further serve as a consistency check.
 */

/** The units of phi, tau in the derivatives and the potential V are the following:
    --> phi is given in units of the reduced Planck mass m_pl = (8 pi G)^(-1/2)
    --> tau in the derivative is given in units of Mpc.
    --> the potential V(phi) is given in units of m_pl^2/Mpc^2.
    With this convention, we have
    rho^{class} = (8 pi G)/3 rho^{physical} = 1/(3 m_pl^2) rho^{physical}
                = 1/3 * [ 1/(2a^2) (phi')^2 + V(phi) ]
    and rho^{class} has the proper dimension Mpc^-2.
 */

double V_e_scf(struct background *pba,
               double phi
               ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return  exp(-scf_lambda*phi);
}

double dV_e_scf(struct background *pba,
                double phi
                ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return -scf_lambda*V_scf(pba,phi);
}

double ddV_e_scf(struct background *pba,
                 double phi
                 ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return pow(-scf_lambda,2)*V_scf(pba,phi);
}


/** parameters and functions for the polynomial coefficient
 * \f$ V_p = (\phi - B)^\alpha + A (polynomial bump) \f$
 * double scf_alpha = 2;
 * double scf_B = 34.8;
 * double scf_A = 0.01; (values for their Figure 2)
 */

double V_p_scf(
               struct background *pba,
               double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  pow(phi - scf_B,  scf_alpha) +  scf_A;
}

double dV_p_scf(
                struct background *pba,
                double phi) {

  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return   scf_alpha*pow(phi -  scf_B,  scf_alpha - 1);
}

double ddV_p_scf(
                 struct background *pba,
                 double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  scf_alpha*(scf_alpha - 1.)*pow(phi -  scf_B,  scf_alpha - 2);
}

/** now the overall potential \f$ V = V_p*V_e \f$
 */

double V_scf(
             struct background *pba,
             double phi) {
  return  V_e_scf(pba,phi)*V_p_scf(pba,phi);
}

double dV_scf(
              struct background *pba,
	      double phi) {
  return dV_e_scf(pba,phi)*V_p_scf(pba,phi) + V_e_scf(pba,phi)*dV_p_scf(pba,phi);
}

double ddV_scf(
               struct background *pba,
               double phi) {
  return ddV_e_scf(pba,phi)*V_p_scf(pba,phi) + 2*dV_e_scf(pba,phi)*dV_p_scf(pba,phi) + V_e_scf(pba,phi)*ddV_p_scf(pba,phi);
}

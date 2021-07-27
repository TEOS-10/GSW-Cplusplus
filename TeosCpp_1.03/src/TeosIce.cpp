#include "TeosIce.h"

/*****************************************
  "TeosIce.cpp" Version 1.03
  by Randall Kent Whited
  rkwhited@gmail.com
  --------------------------------
  This is one modification of the TEOS-10
  file "gsw_oceanographic_toolbox.c"
  as adapted from the TEOS-10
  C version 3.05 http://www.teos-10.org
  ---------------------------------------
  All copyrights and all license issues
  are thes same as for the C version
******************************************/

/*****************************************
  The C++ methods herein are simply
  implementations of C functions in the
  TEOS-10 "gsw_oceanographic_toolbox.c"
  source code file.
*****************************************/

/** descendant class constructor : base class constructor */
TeosIce::TeosIce(): TeosBase()
{
}

/** descendant class destructor */
TeosIce::~TeosIce()
{
}

/**
==========================================================================
  method: gsw_alpha_wrt_t_ice (t, p)
==========================================================================

   Calculates the thermal expansion coefficient of ice with respect to
   in-situ temperature.

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   alpha_wrt_t_ice  =  thermal expansion coefficient of ice with respect
             to in-situ temperature             [ 1/K ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_alpha_wrt_t_ice(double t, double p)
{
    return (gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(0,1,t,p));
}



/**
==========================================================================
  method: gsw_chem_potential_water_ice (t, p)
==========================================================================

   Calculates the chemical potential of water in ice from in-situ
   temperature and pressure.

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   chem_potential_water_ice  =  chemical potential of ice     [ J/kg ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_chem_potential_water_ice(double t, double p)
{
    return (gsw_gibbs_ice(0,0,t,p));
}

/**
==========================================================================
  method: gsw_deltasa_from_sp(sp,p,lon,lat)
==========================================================================

  Calculates Absolute Salinity Anomaly, deltaSA, from Practical Salinity, SP.

  sp     : Practical Salinity               [unitless]
  p      : sea pressure                [dbar]
  lon    : longitude                   [deg E]
  lat    : latitude               [deg N]

  gsw_deltasa_from_sp : Absolute Salinty Anomaly      [g/kg]
*/
double TeosIce::gsw_deltasa_from_sp(double sp, double p, double lon, double lat)
{
    double  res;

    res = gsw_sa_from_sp(sp,p,lon,lat) - gsw_sr_from_sp(sp);

    if (res > cppGSW_ERROR_LIMIT)
    {
      res = cppGSW_INVALID_VALUE;
    }
    return (res);
}



/**
==========================================================================
  method: gsw_enthalpy_diff (sa, ct, p_shallow, p_deep)
==========================================================================

   Calculates the difference of the specific enthalpy of seawater between
   two different pressures, p_deep (the deeper pressure) and p_shallow
   (the shallower pressure), at the same values of SA and CT.  This
   method uses the computationally-efficient expression for specific
   volume in terms of SA, CT and p (Roquet et al., 2014).  The output
   (enthalpy_diff_CT) is the specific enthalpy evaluated at (SA,CT,p_deep)
   minus the specific enthalpy at (SA,CT,p_shallow).

   SA    =  Absolute Salinity             [ g/kg ]
   CT    =  Conservative Temperature (ITS-90)      [ deg C ]
   p_shallow  =  upper sea pressure            [ dbar ]
       ( i.e. shallower absolute pressure - 10.1325 dbar )
   p_deep     =  lower sea pressure            [ dbar ]
       ( i.e. deeper absolute pressure - 10.1325 dbar )

   enthalpy_diff_CT  =  difference of specific enthalpy       [ J/kg ]
              (deep minus shallow)
--------------------------------------------------------------------------
*/
double TeosIce::gsw_enthalpy_diff(double sa, double ct, double p_shallow, double p_deep)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

    double  dynamic_enthalpy_shallow, dynamic_enthalpy_deep,
       part_1, part_2, part_3, part_4, part_5, xs, ys,
       z_deep, z_shallow;

    xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
    ys = ct*0.025;
    z_shallow = p_shallow*1e-4;
    z_deep = p_deep*1e-4;

    part_1 = gsvco.h001
       + xs*(gsvco.h101 + xs*(gsvco.h201 + xs*(gsvco.h301 + xs*(gsvco.h401
       + xs*(gsvco.h501 + gsvco.h601*xs))))) + ys*(gsvco.h011 + xs*(gsvco.h111
       + xs*(gsvco.h211+ xs*(gsvco.h311
       + xs*(gsvco.h411 + gsvco.h511*xs)))) + ys*(gsvco.h021 + xs*(gsvco.h121
       + xs*(gsvco.h221 + xs*(gsvco.h321
       + gsvco.h421*xs))) + ys*(gsvco.h031 + xs*(gsvco.h131 + xs*(gsvco.h231
       + gsvco.h331*xs)) + ys*(gsvco.h041
       + xs*(gsvco.h141 + gsvco.h241*xs) + ys*(gsvco.h051 + gsvco.h151*xs
       + gsvco.h061*ys)))));

    part_2 = gsvco.h002
       + xs*( gsvco.h102 + xs*( gsvco.h202 + xs*( gsvco.h302
       + xs*( gsvco.h402 + gsvco.h502*xs))))
       + ys*( gsvco.h012 + xs*( gsvco.h112 + xs*( gsvco.h212
       + xs*( gsvco.h312 + gsvco.h412*xs))) + ys*( gsvco.h022
       + xs*( gsvco.h122 + xs*( gsvco.h222 + gsvco.h322*xs))
       + ys*( gsvco.h032 + xs*( gsvco.h132 + gsvco.h232*xs)
       + ys*( gsvco.h042 + gsvco.h142*xs + gsvco.h052*ys))));

    part_3 = gsvco.h003
       + xs*(gsvco.h103 + xs*(gsvco.h203 + xs*(gsvco.h303
       + gsvco.h403*xs))) + ys*(gsvco.h013
       + xs*(gsvco.h113 + xs*(gsvco.h213 + gsvco.h313*xs))
       + ys*(gsvco.h023 + xs*(gsvco.h123 + gsvco.h223*xs)
       + ys*(gsvco.h033 + gsvco.h133*xs + gsvco.h043*ys)));

    part_4 = gsvco.h004 + xs*(gsvco.h104 + gsvco.h204*xs)
       + ys*(gsvco.h014 + gsvco.h114*xs + gsvco.h024*ys);

    part_5 = gsvco.h005 + gsvco.h105*xs + gsvco.h015*ys;

    dynamic_enthalpy_shallow =  z_shallow*(part_1 + z_shallow*(part_2
       + z_shallow*(part_3 + z_shallow*(part_4 + z_shallow*(part_5
       + z_shallow*(gsvco.h006 + gsvco.h007*z_shallow))))));

    dynamic_enthalpy_deep = z_deep*(part_1 + z_deep*(part_2 + z_deep*(part_3
       + z_deep*(part_4+z_deep*(part_5 + z_deep*(gsvco.h006 + gsvco.h007*z_deep))))));

    return ((dynamic_enthalpy_deep - dynamic_enthalpy_shallow)*gtc.db2pa*1e4);
}



/**
==========================================================================
  method: gsw_geo_strf_dyn_height_pc (sa, ct, delta_p, &
                   geo_strf_dyn_height_pc, p_mid)
==========================================================================

   Calculates dynamic height anomaly as the integral of specific volume
   anomaly from the the sea surface pressure (0 Pa) to the pressure p.
   This method, gsw_geo_strf_dyn_height_pc, is to used when the
   Absolute Salinity and Conservative Temperature are piecewise constant in
   the vertical over sucessive pressure intervals of delta_p (such as in
   a forward "z-coordinate" ocean model).  "geo_strf_dyn_height_pc" is
   the dynamic height anomaly with respect to the sea surface.  That is,
   "geo_strf_dyn_height_pc" is the geostrophic streammethod for the
   difference between the horizontal velocity at the pressure concerned, p,
   and the horizontal velocity at the sea surface.  Dynamic height anomaly
   is the geostrophic streammethod in an isobaric surface.  The reference
   values used for the specific volume anomaly are SA = SSO = 35.16504 g/kg
   and CT = 0 deg C.  The output values of geo_strf_dyn_height_pc are
   given at the mid-point pressures, p_mid, of each layer in which SA and
   CT are vertically piecewice constant (pc).  This method calculates
   enthalpy using the computationally-efficient 75-term expression for
   specific volume of Roquet et al., (2015).

   SA       =  Absolute Salinity               [ g/kg ]
   CT       =  Conservative Temperature (ITS-90)        [ deg C ]
   delta_p  =  difference in sea pressure between the deep and     [ dbar ]
          shallow extents of each layer in which SA and CT
          are vertically constant. delta_p must be positive.

   Note. sea pressure is absolute pressure minus 10.1325 dbar.

   geo_strf_dyn_height_pc =  dynamic height anomaly        [ m^2/s^2 ]
   p_mid        =  mid-point pressure in each layer      [ dbar ]
---------------------------------------------------------------------------
*/
double *TeosIce::gsw_geo_strf_dyn_height_pc(double *sa, double *ct, double *delta_p,
                   int n_levels, double *geo_strf_dyn_height_pc,
                   double *p_mid)
{
    int i, np;
    double *delta_h=NULL, delta_h_half, dyn_height_deep=0.0,
      *p_deep=NULL, *p_shallow=NULL;

    for (i=0; i<n_levels; i++)
    {
      if (delta_p[i] < 0.0)
      {
         return (NULL);
      }
    }

    np = n_levels;

    delta_h = new double[3*np];
    p_deep = delta_h+np;
    p_shallow = p_deep+np;

    for (i=0; i<np; i++)
    {
      p_deep[i] = (i==0)? delta_p[0] : p_deep[i-1] + delta_p[i];
      p_shallow[i] = p_deep[i] - delta_p[i];
      delta_h[i] = gsw_enthalpy_diff(sa[i],ct[i],p_shallow[i],p_deep[i]);
    }

    for (i=0; i<np; i++)
    {
      dyn_height_deep = dyn_height_deep - delta_h[i];
      /** This is Phi minus Phi_0 of Eqn. (3.32.2) of IOC et al. (2010).*/
      p_mid[i] = 0.5*(p_shallow[i]  + p_deep[i]);
      delta_h_half = gsw_enthalpy_diff(sa[i],ct[i],p_mid[i],p_deep[i]);

      geo_strf_dyn_height_pc[i] = gsw_enthalpy_sso_0(p_mid[i]) +
               dyn_height_deep + delta_h_half;
    }

    if (delta_h) delete []delta_h;

    return (geo_strf_dyn_height_pc);
}

/**
==========================================================================
  method: gsw_helmholtz_energy_ice (t, p)
==========================================================================

   Calculates the Helmholtz energy of ice.

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   Helmholtz_energy_ice  =  Helmholtz energy of ice      [ J/kg ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_helmholtz_energy_ice(double t, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    return (gsw_gibbs_ice(0,0,t,p) -
      (gtc.db2pa*p + gtc.gsw_p0)*
          gsw_gibbs_ice(0,1,t,p));
}
/**
==========================================================================
elemental subroutine gsw_ice_fraction_to_freeze_seawater (sa, ct, p, &
                 t_ih, sa_freeze, ct_freeze, w_ih)
==========================================================================

   Calculates the mass fraction of ice (mass of ice divided by mass of ice
   plus seawater), which, when melted into seawater having (SA,CT,p) causes
   the final dilute seawater to be at the freezing temperature.  The other
   outputs are the Absolute Salinity and Conservative Temperature of the
   final diluted seawater.

   SA   =  Absolute Salinity of seawater            [ g/kg ]
   CT   =  Conservative Temperature of seawater (ITS-90)     [ deg C ]
   p    =  sea pressure                   [ dbar ]
        ( i.e. absolute pressure - 10.1325d0 dbar )
   t_Ih =  in-situ temperature of the ice at pressure p (ITS-90)  [ deg C ]

   SA_freeze = Absolute Salinity of seawater after the mass fraction of
          ice, ice_fraction, at temperature t_Ih has melted into the
          original seawater, and the final mixture is at the freezing
          temperature of seawater.             [ g/kg ]

   CT_freeze = Conservative Temperature of seawater after the mass
          fraction, w_Ih, of ice at temperature t_Ih has melted into
          the original seawater, and the final mixture is at the
          freezing temperature of seawater.        [ deg C ]

   w_Ih      = mass fraction of ice, having in-situ temperature t_Ih,
          which, when melted into seawater at (SA,CT,p) leads to the
          final diluted seawater being at the freezing temperature.
          This output must be between 0 and 1.         [unitless]
---------------------------------------------------------------------------
*/
void TeosIce::gsw_ice_fraction_to_freeze_seawater(double sa, double ct, double p,
                      double t_ih, double *sa_freeze,
                      double *ct_freeze, double *w_ih)
{
    int no_iter;
    double ctf, ctf_mean, ctf_old, ctf_plus1, ctf_zero,
      dfunc_dsaf, func, func_plus1, func_zero, h, h_ih,
      saf, saf_mean, saf_old, tf, h_hat_sa, h_hat_ct, ctf_sa;
    double  sa0 = 0.0, saturation_fraction = 0.0;

    ctf = gsw_ct_freezing(sa,p,saturation_fraction);
    if (ct < ctf)
    {
      /**The seawater ct input is below the freezing temp*/
      *sa_freeze = cppGSW_INVALID_VALUE;
      *ct_freeze = *sa_freeze;
      *w_ih = *sa_freeze;

      return;
    }

    tf = gsw_t_freezing(sa0,p,saturation_fraction);
    if (t_ih > tf)
    {
      /**The input, t_Ih, exceeds the freezing temperature at sa = 0*/
      *sa_freeze = cppGSW_INVALID_VALUE;
      *ct_freeze = *sa_freeze;
      *w_ih = *sa_freeze;

      return;
    }

    h = gsw_enthalpy_ct_exact(sa,ct,p);
    h_ih = gsw_enthalpy_ice(t_ih,p);

    ctf_zero = gsw_ct_freezing(sa0,p,saturation_fraction);
    func_zero = sa*(gsw_enthalpy_ct_exact(sa0,ctf_zero,p) - h_ih);

    ctf_plus1 = gsw_ct_freezing(sa+1.0,p,saturation_fraction);
    func_plus1 = sa*(gsw_enthalpy_ct_exact(sa+1.0,ctf_plus1,p) - h) - (h - h_ih);

    saf = -(sa+1.0)*func_zero/(func_plus1 - func_zero);   /**initial guess*/
    ctf = gsw_ct_freezing(saf,p,saturation_fraction);
    gsw_enthalpy_first_derivatives_ct_exact(saf,ctf,p,&h_hat_sa,&h_hat_ct);
    gsw_ct_freezing_first_derivatives(saf,p,1.0,&ctf_sa,NULL);

    dfunc_dsaf = sa*(h_hat_sa + h_hat_ct*ctf_sa) - (h - h_ih);

    for (no_iter = 1; no_iter <= 2; no_iter++)
    {
      saf_old = saf;
      ctf_old = ctf;
      func = sa*(gsw_enthalpy_ct_exact(saf_old,ctf_old,p) - h)
         - (saf_old - sa)*(h - h_ih);
      saf = saf_old - func/dfunc_dsaf;
      saf_mean = 0.5*(saf + saf_old);
      ctf_mean = gsw_ct_freezing(saf_mean,p,saturation_fraction);
      gsw_enthalpy_first_derivatives_ct_exact(saf_mean,ctf_mean,p,
         &h_hat_sa, &h_hat_ct);
      gsw_ct_freezing_first_derivatives(saf_mean,p,saturation_fraction,
         &ctf_sa, NULL);
      dfunc_dsaf = sa*(h_hat_sa + h_hat_ct*ctf_sa) - (h - h_ih);
      saf = saf_old - func/dfunc_dsaf;
      ctf = gsw_ct_freezing(saf,p,saturation_fraction);
    }
    /**
        After these 2 iterations of this modified Newton-Raphson method, the
        error in SA_freeze is less than 1.3d0x10^-13 g/kg, in CT_freeze is
        less than   4x10^-13 deg C and in w_Ih is less than 3.8d0x10^-15
        which represent machine precision for these calculations.
    */

    *sa_freeze = saf;
    *ct_freeze = ctf;
    *w_ih = (h - gsw_enthalpy_ct_exact(*sa_freeze,*ct_freeze,p))/(h - h_ih);
}

/**
==========================================================================
  method: gsw_internal_energy_ice (t, p)
==========================================================================

   Calculates the specific internal energy of ice.

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   internal_energy_ice  =  specific internal energy (u)         [J/kg]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_internal_energy_ice(double t, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    return (gsw_gibbs_ice(0,0,t,p)
           - (gtc.gsw_t0 + t)*gsw_gibbs_ice(1,0,t,p)
           - (gtc.db2pa*p + gtc.gsw_p0)*gsw_gibbs_ice(0,1,t,p));
}
/**
==========================================================================
subroutine gsw_ipv_vs_fnsquared_ratio(sa,ct,p,pref,nz,ipv_vs_fnsquared_ratio,p_mid)
==========================================================================
   Calculates the ratio of the vertical gradient of potential density to
   the vertical gradient of locally-referenced potential density.  This
   ratio is also the ratio of the planetary Isopycnal Potential Vorticity
   (IPV) to f times N^2, hence the name for this variable,
   IPV_vs_fNsquared_ratio (see Eqn. (3.20.5) of IOC et al. (2010)).
   The reference sea pressure, p_ref, of the potential density surface must
   have a constant value.

   IPV_vs_fNsquared_ratio is evaluated at the mid pressure between the
   individual data points in the vertical.

  sa      : Absolute Salinity    (a profile (length nz))     [g/kg]
  ct      : Conservative Temperature  (a profile (length nz))     [deg C]
  p       : sea pressure         (a profile (length nz))     [dbar]
  p_ref   : reference sea pressure of the potential density surface
    ( i.e. absolute reference pressure - 10.1325 dbar )      [dbar]
  nz      : number of bottles
  IPV_vs_fNsquared_ratio
     : The ratio of the vertical gradient of potential density
       referenced to p_ref, to the vertical gradient of locally-
       referenced potential density.  It is ouput on the same
       vertical (M-1)xN grid as p_mid.
       IPV_vs_fNsquared_ratio is dimensionless.     [ unitless ]
  p_mid   : Mid pressure between p grid  (length nz-1)      [dbar]
*/
void TeosIce::gsw_ipv_vs_fnsquared_ratio(double *sa, double *ct, double *p,
                   double p_ref, int nz,
                   double *ipv_vs_fnsquared_ratio,
                   double *p_mid)
{
    int     k;
    double  dsa, sa_mid, dct, ct_mid;
    double  alpha_mid, beta_mid;
    double  alpha_pref, beta_pref, numerator, denominator;

    if (nz < 2)
    {
      *p_mid = *ipv_vs_fnsquared_ratio = cppGSW_INVALID_VALUE;

      return;
    }

    for (k = 0; k < nz-1; k++)
    {
      dsa = (sa[k] - sa[k+1]);
      dct = (ct[k] - ct[k+1]);
      sa_mid = 0.5*(sa[k] + sa[k+1]);
      ct_mid = 0.5*(ct[k] + ct[k+1]);
      p_mid[k] = 0.5*(p[k] + p[k+1]);

      alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
      beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);
      alpha_pref = gsw_alpha(sa_mid,ct_mid,p_ref);
      beta_pref = gsw_beta(sa_mid,ct_mid,p_ref);

      numerator = dct*alpha_pref - dsa*beta_pref;
      denominator = dct*alpha_mid - dsa*beta_mid;

      if (denominator == 0.0)
      {
         ipv_vs_fnsquared_ratio[k] = cppGSW_INVALID_VALUE;
      }
      else
         ipv_vs_fnsquared_ratio[k] = numerator/denominator;
   }
}

/**
==========================================================================
  method: gsw_kappa_const_t_ice (t, p)
==========================================================================

   Calculates isothermal compressibility of ice.
   Note. This is the compressibility of ice AT CONSTANT IN-SITU
     TEMPERATURE

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   kappa_const_t_ice  =  isothermal compressibility      [ 1/Pa ]
    Note. The output units are 1/Pa not 1/dbar.
--------------------------------------------------------------------------
*/
double TeosIce::gsw_kappa_const_t_ice(double t, double p)
{
   return (-gsw_gibbs_ice(0,2,t,p)/gsw_gibbs_ice(0,1,t,p));
}
/**
==========================================================================
  method: gsw_kappa_ice (t, p)
==========================================================================

   Calculates the isentropic compressibility of ice.

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   kappa_ice  =  isentropic compressibility         [ 1/Pa ]
    Note. The output units are 1/Pa not 1/dbar.
--------------------------------------------------------------------------
*/
double TeosIce::gsw_kappa_ice(double t, double p)
{
    double gi_tp, gi_tt;

    gi_tt = gsw_gibbs_ice(2,0,t,p);
    gi_tp = gsw_gibbs_ice(1,1,t,p);

    return ((gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))/
        (gsw_gibbs_ice(0,1,t,p)*gi_tt));
}
/**
==========================================================================
  method: gsw_melting_ice_equilibrium_sa_ct_ratio (sa, p)
==========================================================================

   Calculates the ratio of SA to CT changes when ice melts into seawater
   with both the seawater and the seaice temperatures being almost equal to
   the equilibrium freezing temperature.  It is assumed that a small mass
   of ice melts into an infinite mass of seawater.  If indeed the
   temperature of the seawater and the ice were both equal to the freezing
   temperature, then no melting or freezing would occur an imbalance
   between these three temperatures is needed for freezing or melting to
   occur (the three temperatures being (1) the seawater temperature,
   (2) the ice temperature, and (3) the freezing temperature.

   The output, melting_ice_equilibrium_SA_CT_ratio, is dSA/dCT rather than
   dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is zero
   whereas dCT/dSA would be infinite.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   p  :  sea pressure at which the melting occurs       [ dbar ]
     ( i.e. absolute pressure - 10.1325d0 dbar )

   melting_ice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
             changes when ice melts into seawater, with
             the seawater and seaice being close to the
             freezing temperature.    [ g/(kg K) ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_ice_equilibrium_sa_ct_ratio(double sa, double p)
{
	double ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
	double saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	t_seaice = gsw_t_freezing(sa,p,saturation_fraction);

	h = gsw_enthalpy_ct_exact(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	/**note that h_hat_ct is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa*h_hat_ct/(h - h_ih - sa*h_hat_sa));
}
/**
==========================================================================
  method: gsw_melting_ice_equilibrium_sa_ct_ratio_poly (sa, p)
==========================================================================
!
   Calculates the ratio of SA to CT changes when ice melts into seawater
   with both the seawater and the seaice temperatures being almost equal to
   the equilibrium freezing temperature.  It is assumed that a small mass
   of ice melts into an infinite mass of seawater.  If indeed the
   temperature of the seawater and the ice were both equal to the freezing
   temperature, then no melting or freezing would occur an imbalance
   between these three temperatures is needed for freezing or melting to
   occur (the three temperatures being (1) the seawater temperature,
   (2) the ice temperature, and (3) the freezing temperature.

   The output, melting_ice_equilibrium_SA_CT_ratio, is dSA/dCT rather than
   dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is zero
   whereas dCT/dSA would be infinite.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   p  :  sea pressure at which the melting occurs       [ dbar ]
     ( i.e. absolute pressure - 10.1325d0 dbar )

   melting_ice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
             changes when ice melts into seawater, with
             the seawater and seaice being close to the
             freezing temperature.    [ g/(kg K) ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_ice_equilibrium_sa_ct_ratio_poly(double sa, double p)
{
	double  ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	t_seaice = gsw_t_freezing_poly(sa,p,saturation_fraction);

	h = gsw_enthalpy(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);
	gsw_enthalpy_first_derivatives(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	/**note that h_hat_ct is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa*h_hat_ct / (h - h_ih - sa*h_hat_sa));
}
/**
==========================================================================
elemental subroutine gsw_melting_ice_into_seawater (sa, ct, p, w_ih, t_ih,&
                   sa_final, ct_final, w_ih_final)
==========================================================================

   Calculates the final Absolute Salinity, final Conservative Temperature
   and final ice mass fraction that results when a given mass fraction of
   ice melts and is mixed into seawater whose properties are (SA,CT,p).
   This code takes the seawater to contain no dissolved air.

   When the mass fraction w_Ih_final is calculated as being a positive
   value, the seawater-ice mixture is at thermodynamic equlibrium.

   This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
   is sufficiently large (i.e. sufficiently "warm") so that there is no ice
   present in the final state.  In this case the final state consists of
   only seawater rather than being an equlibrium mixture of seawater and
   ice which occurs when w_Ih_final is positive.  Note that when
   w_Ih_final = 0, the final seawater is not at the freezing temperature.

   SA   =  Absolute Salinity of seawater            [ g/kg ]
   CT   =  Conservative Temperature of seawater (ITS-90)     [ deg C ]
   p    =  sea pressure at which the melting occurs      [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   w_Ih =  mass fraction of ice, that is the mass of ice divided by the
      sum of the masses of ice and seawater.  That is, the mass of
      ice divided by the mass of the final mixed fluid.
      w_Ih must be between 0 and 1.             [ unitless ]
   t_Ih =  the in-situ temperature of the ice (ITS-90)       [ deg C ]

   SA_final    =  Absolute Salinity of the seawater in the final state,
        whether or not any ice is present.          [ g/kg ]
   CT_final    =  Conservative Temperature of the seawater in the the final
        state, whether or not any ice is present.       [ deg C ]
   w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
        If this ice mass fraction is positive, the system is at
        thermodynamic equilibrium.  If this ice mass fraction is
        zero there is no ice in the final state which consists
        only of seawater which is warmer than the freezing
        temperature.               [unitless]
---------------------------------------------------------------------------
*/
void TeosIce::gsw_melting_ice_into_seawater(double sa, double ct, double p,
                   double w_ih, double t_ih, double *sa_final,
                   double *ct_final, double *w_ih_final)
{
    double  ctf, h_bulk, sa_bulk, tf_ih;
    double  saturation_fraction = 0.0;

    ctf = gsw_ct_freezing(sa,p,saturation_fraction);
    if (ct < ctf)
    {
      /**The seawater ct input is below the freezing temp*/
      *sa_final = cppGSW_INVALID_VALUE;
      *ct_final = *sa_final;
      *w_ih_final = *sa_final;

      return;
    }

    tf_ih = gsw_t_freezing(0.0,p,saturation_fraction) - 1e-6;
    if (t_ih > tf_ih)
    {
      /**
         t_ih input exceeds the freezing temp.
         The 1e-6 C buffer in the allowable
         t_Ih is to ensure that there is some ice Ih in the sea ice.
      */
      *sa_final = cppGSW_INVALID_VALUE;
      *ct_final = *sa_final;
      *w_ih_final = *sa_final;

      return;
    }

    sa_bulk = (1.0 - w_ih)*sa;
    h_bulk = (1.0 - w_ih)*gsw_enthalpy_ct_exact(sa,ct,p)
           + w_ih*gsw_enthalpy_ice(t_ih,p);

    gsw_frazil_properties(sa_bulk,h_bulk,p,sa_final,ct_final,w_ih_final);

    if (*sa_final > cppGSW_ERROR_LIMIT)
    {
      *sa_final = cppGSW_INVALID_VALUE;
      *ct_final = *sa_final;
      *w_ih_final = *sa_final;

      return;
    }
}
/**
==========================================================================
  method: gsw_melting_ice_sa_ct_ratio (sa, ct, p, t_ih)
==========================================================================

   Calculates the ratio of SA to CT changes when ice melts into seawater.
   It is assumed that a small mass of ice melts into an infinite mass of
   seawater.  Because of the infinite mass of seawater, the ice will always
   melt.

   The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
   This is done so that when SA = 0, the output, dSA/dCT is zero whereas
   dCT/dSA would be infinite.

   SA   =  Absolute Salinity of seawater            [ g/kg ]
   CT   =  Conservative Temperature of seawater (ITS-90)     [ deg C ]
   p    =  sea pressure at which the melting occurs      [ dbar ]
     ( i.e. absolute pressure - 10.1325d0 dbar )
   t_Ih =  the in-situ temperature of the ice (ITS-90)       [ deg C ]

   melting_ice_SA_CT_ratio = the ratio of SA to CT changes when ice melts
              into a large mass of seawater
                        [ g kg^-1 K^-1 ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_ice_sa_ct_ratio(double sa, double ct, double p, double t_ih)
{
	double  ctf, h, h_ih, tf, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);

	if (ct < ctf)
	{
		/**the seawater ct input is below the freezing temperature*/
		return (cppGSW_INVALID_VALUE);
	}

	tf = gsw_t_freezing(0.0,p,saturation_fraction);
	if (t_ih > tf)
	{
		/**t_ih exceeds the freezing temperature at sa = 0*/
		return (cppGSW_INVALID_VALUE);
	}

	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_ih,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,&h_hat_sa,&h_hat_ct);
	/**Note that h_hat_CT is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa*h_hat_ct/(h - h_ih - sa*h_hat_sa));
}
/**
==========================================================================
  method: gsw_melting_ice_sa_ct_ratio_poly (sa, ct, p, t_ih)
==========================================================================

   Calculates the ratio of SA to CT changes when ice melts into seawater.
   It is assumed that a small mass of ice melts into an infinite mass of
   seawater.  Because of the infinite mass of seawater, the ice will always
   melt.

   The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
   This is done so that when SA = 0, the output, dSA/dCT is zero whereas
   dCT/dSA would be infinite.

   SA   =  Absolute Salinity of seawater            [ g/kg ]
   CT   =  Conservative Temperature of seawater (ITS-90)     [ deg C ]
   p    =  sea pressure at which the melting occurs      [ dbar ]
     ( i.e. absolute pressure - 10.1325d0 dbar )
   t_Ih =  the in-situ temperature of the ice (ITS-90)       [ deg C ]

   melting_ice_SA_CT_ratio = the ratio of SA to CT changes when ice melts
              into a large mass of seawater
                        [ g kg^-1 K^-1 ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_ice_sa_ct_ratio_poly(double sa, double ct, double p, double t_ih)
{
	double  ctf, h, h_ih, tf, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**the seawater ct input is below the freezing temperature*/
		return (cppGSW_INVALID_VALUE);
	}

	tf = gsw_t_freezing_poly(0.0,p,saturation_fraction);
	if (t_ih > tf)
	{
		/**t_ih exceeds the freezing temperature at sa = 0*/
		return (cppGSW_INVALID_VALUE);
	}

	h = gsw_enthalpy(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_ih,p);
	gsw_enthalpy_first_derivatives(sa,ct,p,&h_hat_sa,&h_hat_ct);
	/**Note that h_hat_CT is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa*h_hat_ct/(h - h_ih - sa*h_hat_sa));
}
/**
==========================================================================
  method: gsw_melting_seaice_equilibrium_sa_ct_ratio (sa, p)
==========================================================================

   Calculates the ratio of SA to CT changes when sea ice melts into
   seawater with both the seawater and the sea ice temperatures being
   almost equal to the equilibrium freezing temperature.  It is assumed
   that a small mass of seaice melts into an infinite mass of seawater.  If
   indeed the temperature of the seawater and the sea ice were both equal
   to the freezing temperature, then no melting or freezing would occur; an
   imbalance between these three temperatures is needed for freezing or
   melting to occur (the three temperatures being (1) the seawater
   temperature, (2) the sea ice temperature, and (3) the freezing
   temperature.

   Note that the output of this method, dSA/dCT is independent of the
   sea ice salinity, SA_seaice.  That is, the output applies equally to
   pure ice Ih and to sea ice with seaice salinity, SA_seaice.  This result
   is proven in the manuscript, McDougall et al. (2013).

   The output, melting_seaice_equilibrium_SA_CT_ratio, is dSA/dCT rather
   than dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is
   zero whereas dCT/dSA would be infinite.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   p  :  sea pressure at which the melting occurs       [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   melting_seaice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
              changes when sea ice melts into seawater, with
              the seawater and sea ice being close to the
              freezing temperature.        [ g/(kg K) ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_seaice_equilibrium_sa_ct_ratio(double sa, double p)
{
   double  ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
   double  saturation_fraction = 0.0;

   ctf = gsw_ct_freezing(sa,p,saturation_fraction);
   t_seaice = gsw_t_freezing(sa,p,saturation_fraction);

   h = gsw_enthalpy_ct_exact(sa,ctf,p);
   h_ih = gsw_enthalpy_ice(t_seaice,p);
   gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,&h_hat_sa,&h_hat_ct);

   return (sa*h_hat_ct / (h - h_ih - sa*h_hat_sa));
}
/**
==========================================================================
  method: gsw_melting_seaice_equilibrium_sa_ct_ratio_poly (sa, p)
==========================================================================

   Calculates the ratio of SA to CT changes when sea ice melts into
   seawater with both the seawater and the sea ice temperatures being
   almost equal to the equilibrium freezing temperature.  It is assumed
   that a small mass of seaice melts into an infinite mass of seawater.  If
   indeed the temperature of the seawater and the sea ice were both equal
   to the freezing temperature, then no melting or freezing would occur; an
   imbalance between these three temperatures is needed for freezing or
   melting to occur (the three temperatures being (1) the seawater
   temperature, (2) the sea ice temperature, and (3) the freezing
   temperature.

   Note that the output of this method, dSA/dCT is independent of the
   sea ice salinity, SA_seaice.  That is, the output applies equally to
   pure ice Ih and to sea ice with seaice salinity, SA_seaice.  This result
   is proven in the manuscript, McDougall et al. (2013).

   The output, melting_seaice_equilibrium_SA_CT_ratio, is dSA/dCT rather
   than dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT is
   zero whereas dCT/dSA would be infinite.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   p  :  sea pressure at which the melting occurs       [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )

   melting_seaice_equilibrium_SA_CT_ratio = the ratio dSA/dCT of SA to CT
              changes when sea ice melts into seawater, with
              the seawater and sea ice being close to the
              freezing temperature.        [ g/(kg K) ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(double sa, double p)
{
    double  ctf, h, h_ih, t_seaice, h_hat_sa, h_hat_ct;
    double  saturation_fraction = 0.0;

    ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
    t_seaice = gsw_t_freezing_poly(sa,p,saturation_fraction);

    h = gsw_enthalpy(sa,ctf,p);
    h_ih = gsw_enthalpy_ice(t_seaice,p);
    gsw_enthalpy_first_derivatives(sa,ctf,p,&h_hat_sa,&h_hat_ct);

    return (sa*h_hat_ct / (h - h_ih - sa*h_hat_sa));
}
/**
==========================================================================
elemental subroutine gsw_melting_seaice_into_seawater (sa, ct, p, &
          w_seaice, sa_seaice, t_seaice, sa_final, ct_final)
==========================================================================

   Calculates the Absolute Salinity and Conservative Temperature that
   results when a given mass of sea ice (or ice) melts and is mixed into a
   known mass of seawater (whose properties are (SA,CT,p)).

   If the ice contains no salt (e.g. if it is of glacial origin), then the
   input 'SA_seaice' should be set to zero.

   Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
   and 12 g/kg of salt (defined as the mass of salt divided by the mass of
   ice Ih plus brine) and this programme returns NaN's if the input
   SA_seaice is greater than 15 g/kg.  If the SA_seaice input is not zero,
   usually this would imply that the pressure p should be zero, as sea ice
   only occurs near the sea surface.  The code does not impose that p = 0
   if SA_seaice is non-zero.  Rather, this is left to the user.

   The Absolute Salinity, SA_brine, of the brine trapped in little pockets
   in the sea ice, is in thermodynamic equilibrium with the ice Ih that
   surrounds these pockets.  As the sea ice temperature, t_seaice, may be
   less than the freezing temperature, SA_brine is usually greater than the
   Absolute Salinity of the seawater at the time and place when and where
   the sea ice was formed.  So usually SA_brine will be larger than SA.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   CT    :  Conservative Temperature of seawater (ITS-90)      [ deg C ]
   p  :  sea pressure at which the melting occurs       [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   w_seaice  =  mass fraction of sea ice, that is the mass of sea ice
           divided by the sum of the masses of sea ice and seawater.
           That is, the mass of sea ice divided by the mass of the
           final mixed fluid.  w_seaice must be between 0 and 1.
                            [ unitless ]
   SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction of
           salt in sea ice, expressed in g of salt per kg of sea ice.
                           [ g/kg ]
   t_seaice  =  the in-situ temperature of the sea ice (or ice) (ITS-90)
                          [ deg C ]

   SA_final  =  Absolute Salinity of the mixture of the melted sea ice
           (or ice) and the orignal seawater        [ g/kg ]
   CT_final  =  Conservative Temperature of the mixture of the melted
           sea ice (or ice) and the orignal seawater    [ deg C ]
--------------------------------------------------------------------------
*/
void TeosIce::gsw_melting_seaice_into_seawater(double sa, double ct, double p,
                  double w_seaice, double sa_seaice,
                  double t_seaice, double *sa_final,
                  double *ct_final)
{
	double  ctf, h, h_brine, h_final, h_ih, sa_brine, tf_sa_seaice;
	double  saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**The seawater ct input is below the freezing temp*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;

		return;
	}

	tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction) - 1e-6;

	if (t_seaice > tf_sa_seaice)
	{
		/**
		  The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
		  some ice Ih in the sea ice. Without this buffer, that is if t_seaice
		  is allowed to be exactly equal to tf_sa_seaice, the seaice is
		  actually 100% brine at Absolute Salinity of SA_seaice.
		*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;

		return;
	}

	sa_brine = gsw_sa_freezing_from_t(t_seaice,p,saturation_fraction);

	if (sa_brine >= cppGSW_ERROR_LIMIT)
	{
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;

		return;
	}

	h_brine = gsw_enthalpy_t_exact(sa_brine,t_seaice,p);

	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);

	h_final = h - w_seaice*(h - h_ih - (h_brine - h_ih)*sa_seaice/sa_brine);

	*sa_final = sa - w_seaice*(sa - sa_seaice);
	/**
	ctf = gsw_ct_freezing(sa_final,p,saturation_fraction)

	if (h_final .lt. gsw_enthalpy_ct_exact(sa_final,ctf,p)) then
	       Melting this much seaice is not possible as it would result in
	       frozen seawater
	     sa_final = gsw_error_code(4,func_name)
	     ct_final = sa_final
	     return
	end if
	*/
	*ct_final = gsw_ct_from_enthalpy_exact(*sa_final,h_final,p);

	if (*ct_final > cppGSW_ERROR_LIMIT)
	{
		*sa_final = *ct_final;
	}
}
/**
==========================================================================
  method: gsw_melting_seaice_sa_ct_ratio (sa, ct, p, sa_seaice, &
                     t_seaice)
==========================================================================

  Calculates the ratio of SA to CT changes when sea ice melts into seawater.
  It is assumed that a small mass of sea ice melts into an infinite mass of
  seawater.  Because of the infinite mass of seawater, the sea ice will
  always melt.

  Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
  and 12 g/kg of salt (defined as the mass of salt divided by the mass of
  ice Ih plus brine) and this programme returns NaN's if the input
  SA_seaice is greater than 15 g/kg.  If the SA_seaice input is not zero,
  usually this would imply that the pressure p should be zero, as sea ice
  only occurs near the sea surface.  The code does not impose that p = 0 if
  SA_seaice is non-zero.  Rather, this is left to the user.

  The Absolute Salinity, SA_brine, of the brine trapped in little pockets
  in the sea ice, is in thermodynamic equilibrium with the ice Ih that
  surrounds these pockets.  As the seaice temperature, t_seaice, may be
  less than the freezing temperature, SA_brine is usually greater than the
  Absolute Salinity of the seawater at the time and place when and where
  the sea ice was formed.  So usually SA_brine will be larger than SA.

  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
  This is done so that when (SA - seaice_SA) = 0, the output, dSA/dCT is
  zero whereas dCT/dSA would be infinite.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   CT    :  Conservative Temperature of seawater (ITS-90)      [ deg C ]
   p  :  sea pressure at which the melting occurs       [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   SA_seaice  :  Absolute Salinity of sea ice, that is, the mass fraction
       of salt in sea ice expressed in g of salt per kg of
       sea ice                  [ g/kg ]
   t_seaice : the in-situ temperature of the sea ice (ITS-90)     [ deg C ]

   melting_seaice_SA_CT_ratio : the ratio dSA/dCT of SA to CT changes when
       sea ice melts into a large mass of seawater   [ g/(kg K) ]
---------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_seaice_sa_ct_ratio(double sa, double ct, double p,
   double sa_seaice, double t_seaice)
{
   double  ctf, delsa, h, h_brine, h_ih, sa_brine,
      tf_sa_seaice, h_hat_sa, h_hat_ct;
   double  saturation_fraction = 0.0;

   if (sa_seaice < 0.0 || sa_seaice > 15.0) {
      return (cppGSW_INVALID_VALUE);
   }

   ctf = gsw_ct_freezing(sa,p,saturation_fraction);
   if (ct < ctf) {    /*the seawater ct input is below the freezing temp*/
      return (cppGSW_INVALID_VALUE);
   }

   tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction) - 1e-6;
   if (t_seaice > tf_sa_seaice) {   /*t_seaice exceeds the freezing sa*/
      return (cppGSW_INVALID_VALUE);
   }
   /**
   ------------------------------------------------------------------------
   The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
   some ice Ih in the sea ice.  Without this buffer, that is if t_seaice
   is allowed to be exactly equal to tf_sa_seaice, the sea ice is actually
   100% brine at Absolute Salinity of SA_seaice.
   ------------------------------------------------------------------------
    */
   h = gsw_enthalpy_ct_exact(sa,ct,p);
   h_ih = gsw_enthalpy_ice(t_seaice,p);
   gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,&h_hat_sa,&h_hat_ct);

   sa_brine = gsw_sa_freezing_from_t(t_seaice,p,saturation_fraction);
   if (sa_brine > cppGSW_ERROR_LIMIT) {
      return (cppGSW_INVALID_VALUE);
   }
   h_brine = gsw_enthalpy_t_exact(sa_brine,t_seaice,p);
   delsa = sa - sa_seaice;

   return (h_hat_ct*delsa /
      (h - h_ih - delsa*h_hat_sa - sa_seaice*(h_brine - h_ih)/sa_brine));
}
/**
==========================================================================
  method: gsw_melting_seaice_sa_ct_ratio_poly (sa, ct, p, &
                         sa_seaice, t_seaice)
==========================================================================

  Calculates the ratio of SA to CT changes when sea ice melts into seawater.
  It is assumed that a small mass of sea ice melts into an infinite mass of
  seawater.  Because of the infinite mass of seawater, the sea ice will
  always melt.

  Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
  and 12 g/kg of salt (defined as the mass of salt divided by the mass of
  ice Ih plus brine) and this programme returns NaN's if the input
  SA_seaice is greater than 15 g/kg.  If the SA_seaice input is not zero,
  usually this would imply that the pressure p should be zero, as sea ice
  only occurs near the sea surface.  The code does not impose that p = 0 if
  SA_seaice is non-zero.  Rather, this is left to the user.

  The Absolute Salinity, SA_brine, of the brine trapped in little pockets
  in the sea ice, is in thermodynamic equilibrium with the ice Ih that
  surrounds these pockets.  As the seaice temperature, t_seaice, may be
  less than the freezing temperature, SA_brine is usually greater than the
  Absolute Salinity of the seawater at the time and place when and where
  the sea ice was formed.  So usually SA_brine will be larger than SA.

  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA.
  This is done so that when (SA - seaice_SA) = 0, the output, dSA/dCT is
  zero whereas dCT/dSA would be infinite.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   CT    :  Conservative Temperature of seawater (ITS-90)      [ deg C ]
   p  :  sea pressure at which the melting occurs       [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   SA_seaice  :  Absolute Salinity of sea ice, that is, the mass fraction
       of salt in sea ice expressed in g of salt per kg of
       sea ice                  [ g/kg ]
   t_seaice : the in-situ temperature of the sea ice (ITS-90)     [ deg C ]

   melting_seaice_SA_CT_ratio : the ratio dSA/dCT of SA to CT changes when
       sea ice melts into a large mass of seawater   [ g/(kg K) ]
---------------------------------------------------------------------------
*/
double TeosIce::gsw_melting_seaice_sa_ct_ratio_poly(double sa, double ct, double p,
   double sa_seaice, double t_seaice)
{
   double  ctf, delsa, h, h_brine, h_ih, sa_brine, ct_brine,
      tf_sa_seaice, h_hat_sa, h_hat_ct;
   double  saturation_fraction = 0.0;

   if (sa_seaice < 0.0 || sa_seaice > 15.0) {
      return (cppGSW_INVALID_VALUE);
   }

   ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
   if (ct < ctf)
   { /**the seawater ct input is below the freezing temp*/
      return (cppGSW_INVALID_VALUE);
   }

   tf_sa_seaice = gsw_t_freezing_poly(sa_seaice,p,saturation_fraction) - 1e-6;
   if (t_seaice > tf_sa_seaice)
   {   /**t_seaice exceeds the freezing sa*/
      return (cppGSW_INVALID_VALUE);
   }
   /**
   -----------------------------------------------------------------------
   The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
   some ice Ih in the sea ice.  Without this buffer, that is if t_seaice
   is allowed to be exactly equal to tf_sa_seaice, the sea ice is actually
   100% brine at Absolute Salinity of SA_seaice.
   -----------------------------------------------------------------------
   */
   h = gsw_enthalpy(sa,ct,p);
   h_ih = gsw_enthalpy_ice(t_seaice,p);
   gsw_enthalpy_first_derivatives(sa,ct,p,&h_hat_sa,&h_hat_ct);

   sa_brine = gsw_sa_freezing_from_t_poly(t_seaice,p,saturation_fraction);
   if (sa_brine > cppGSW_ERROR_LIMIT)
   {
      return (cppGSW_INVALID_VALUE);
   }
   ct_brine = gsw_ct_from_t(sa_brine,t_seaice,p);
   h_brine = gsw_enthalpy(sa_brine,ct_brine,p);
   delsa = sa - sa_seaice;

   return (h_hat_ct*delsa /
        (h - h_ih - delsa*h_hat_sa - sa_seaice*(h_brine - h_ih)/sa_brine));
}
/**
==========================================================================
method: gsw_nsquared(*sa,*ct,*p,*lat,nz,*n2,*p_mid)
==========================================================================
--------------------------------------------------------------------------
  (water column properties, based on the 48-term expression for density)
--------------------------------------------------------------------------
   Calculates the buoyancy frequency squared (N^2)(i.e. the Brunt-Vaisala
   frequency squared) at the mid pressure from the equation,


       2      2        beta.d(SA) - alpha.d(CT)
     N   =  g  .rho_local. -------------------------
                  dP

   The pressure increment, dP, in the above formula is in Pa, so that it is
   10^4 times the pressure increment dp in dbar.

  sa     : Absolute Salinity    (a profile (length nz))     [g/kg]
  ct     : Conservative Temperature  (a profile (length nz))     [deg C]
  p      : sea pressure         (a profile (length nz))     [dbar]
  lat    : latitude        (a profile (length nz))     [deg N]
  nz     : number of levels in the profile
  n2     : Brunt-Vaisala Frequency squared  (length nz-1)   [s^-2]
  p_mid  : Mid pressure between p grid      (length nz-1)   [dbar]
*/
void TeosIce::gsw_nsquared(double *sa, double *ct, double *p, double *lat, int nz,
   double *n2, double *p_mid)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */

   int     k;
   double  p_grav, n_grav, grav_local, dsa, sa_mid, dct, ct_mid,
      dp, rho_mid, alpha_mid, beta_mid;

   if (nz < 2)
       return;

   p_grav  = gsw_grav(lat[0],p[0]);

   for (k = 0; k < nz-1; k++)
   {
       n_grav = gsw_grav(lat[k+1],p[k+1]);
       grav_local = 0.5*(p_grav + n_grav);
       dsa = (sa[k+1] - sa[k]);
       sa_mid = 0.5*(sa[k] + sa[k+1]);
       dct = (ct[k+1] - ct[k]);
       ct_mid      = 0.5*(ct[k] + ct[k+1]);
       dp = (p[k+1] - p[k]);
       p_mid[k] = 0.5*(p[k] + p[k+1]);
       rho_mid = gsw_rho(sa_mid,ct_mid,p_mid[k]);
       alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
       beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);

       n2[k] = (grav_local*grav_local)*(rho_mid/(gtc.db2pa*dp))*
           (beta_mid*dsa - alpha_mid*dct);

       p_grav = n_grav;
   }
}


/**
==========================================================================
  method: gsw_pot_enthalpy_from_pt_ice_poly (pt0_ice)
==========================================================================

   Calculates the potential enthalpy of ice from potential temperature of
   ice (whose reference sea pressure is zero dbar).  This is a
   compuationally efficient polynomial fit to the potential enthalpy of
   ice.

   pt0_ice  :  potential temperature of ice (ITS-90)         [ deg C ]
   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_pot_enthalpy_from_pt_ice_poly(double pt0_ice)
{
   int iteration;
   double df_dt, f, pot_enthalpy_ice,
      pot_enthalpy_ice_mid, pot_enthalpy_ice_old,
      p0 = -3.333601570157700e5,
      p1 =  2.096693916810367e3,
      p2 =  3.687110754043292,
      p3 =  4.559401565980682e-4,
      p4 = -2.516011957758120e-6,
      p5 = -1.040364574632784e-8,
      p6 = -1.701786588412454e-10,
      p7 = -7.667191301635057e-13;

   /**initial estimate of the potential enthalpy.*/
   pot_enthalpy_ice = p0 + pt0_ice*(p1 + pt0_ice*(p2 + pt0_ice*(p3
            + pt0_ice*(p4 + pt0_ice*(p5 + pt0_ice*(p6
            + pt0_ice*p7))))));

   df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice);

   for (iteration = 1; iteration <= 5; iteration++) {
       pot_enthalpy_ice_old = pot_enthalpy_ice;
       f = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old)
         - pt0_ice;
       pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt;
       pot_enthalpy_ice_mid = 0.5*(pot_enthalpy_ice+pot_enthalpy_ice_old);
       df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice_mid);
       pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt;
   }
   /**
   The error of this fit ranges between -6e-3 and 6e-3 J/kg over the
   potential temperature range of -100 to 2 deg C, or the potential
   enthalpy range of -5.7 x 10^5 to -3.3 x 10^5 J/kg.
   */
   return (pot_enthalpy_ice);
}



/**
==========================================================================
  method: gsw_pressure_coefficient_ice (t, p)
==========================================================================

   Calculates pressure coefficient of ice.

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  :  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   pressure_coefficient_ice  :  pressure coefficient of ice     [Pa/K]
    Note. The output units are Pa/K NOT dbar/K.
--------------------------------------------------------------------------
*/
double TeosIce::gsw_pressure_coefficient_ice(double t, double p)
{
   return (-gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(0,2,t,p));
}
/**
==========================================================================
  method: gsw_pressure_freezing_ct (sa, ct, saturation_fraction)
==========================================================================

   Calculates the pressure (in dbar) of seawater at the freezing
   temperature.  That is, the output is the pressure at which seawater,
   with Absolute Salinity SA, Conservative Temperature CT, and with
   saturation_fraction of dissolved air, freezes.  If the input values are
   such that there is no value of pressure in the range between 0 dbar and
   10,000 dbar for which seawater is at the freezing temperature, the
   output, pressure_freezing, is put equal to NaN.

   SA    :  Absolute Salinity of seawater             [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)             [ deg C ]
   saturation_fraction : the saturation fraction of dissolved air in
          seawater
   pressure_freezing : sea pressure at which the seawater freezes  [ dbar ]
    ( i.e. absolute pressure - 10.1325 dbar )
---------------------------------------------------------------------------
*/
double TeosIce::gsw_pressure_freezing_ct(double sa, double ct, double saturation_fraction)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */

   int     i_iter, number_of_iterations = 3;
   double  ct_freezing_p0, ct_freezing_p10000, dctf_dp, f,
      pf, pfm, pf_old, ctfreezing_p;

   /**
     rec_Pa2dbar is to have dCTf_dp in units of K/dbar rather than K/Pa
     Find CT > CT_freezing_p0.  If this is the case, the input CT value
     represent seawater that will not be frozen at any positive p.
   */
   ct_freezing_p0 = gsw_ct_freezing(sa,0.0,saturation_fraction);
   if (ct > ct_freezing_p0) {
       return (cppGSW_INVALID_VALUE);
   }
   /**
     Find CT < CT_freezing_p10000.  If this is the case, the input CT value
     represent seawater that is frozen even at p = 10,000 dbar.
   */
   ct_freezing_p10000 = gsw_ct_freezing(sa,1e4,saturation_fraction);
   if (ct < ct_freezing_p10000)
   {
       return (cppGSW_INVALID_VALUE);
   }
   /**
     This is the initial (linear) guess of the freezing pressure, in dbar.
   */
   pf = gtc.rec_pa2db*(ct_freezing_p0 - ct)/
         (ct_freezing_p0 - ct_freezing_p10000);

   gsw_ct_freezing_first_derivatives(sa,pf,saturation_fraction,
                      NULL,&ctfreezing_p);
   dctf_dp = gtc.rec_pa2db*ctfreezing_p;
       /**
          this dctf_dp is the initial value of the partial derivative of
          ct_freezing with respect to pressure (in dbar) at fixed sa,
          assuming that the saturation_fraction is zero.
       */
   for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
   {
       pf_old = pf;
       f = gsw_ct_freezing(sa,pf_old,saturation_fraction) - ct;
       pf = pf_old - f/dctf_dp;
       pfm = 0.5*(pf + pf_old);
       gsw_ct_freezing_first_derivatives(sa,pfm,saturation_fraction,
                     NULL,&ctfreezing_p);
       dctf_dp = gtc.rec_pa2db*ctfreezing_p;
       pf = pf_old - f/dctf_dp;
   }

   if (gsw_sa_p_inrange(sa,pf))
       return (pf);

   return (cppGSW_INVALID_VALUE);
}
/**
==========================================================================
  method: gsw_pt0_cold_ice_poly (pot_enthalpy_ice)
==========================================================================

   Calculates an initial estimate of pt0_ice when it is less than about
   -100 deg C.

   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]

   pt0_cold_ice_poly  :  initial estimate of potential temperatur
          of very cold ice in dgress C (not K)     [ deg C ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_pt0_cold_ice_poly(double pot_enthalpy_ice)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */
   double  log_abs_theta0, log_h_diff,
      /**h00 = gsw_enthalpy_ice(-gtc.gsw_t0,0)*/
      h00 = -6.320202333358860e5,

      s0 =  1.493103204647916,
      s1 =  2.372788609320607e-1,
      s2 = -2.014996002119374e-3,
      s3 =  2.640600197732682e-6,
      s4 =  3.134706016844293e-5,
      s5 =  2.733592344937913e-6,
      s6 =  4.726828010223258e-8,
      s7 = -2.735193883189589e-9,
      s8 = -8.547714991377670e-11;

   log_h_diff = log(pot_enthalpy_ice - h00);

   log_abs_theta0 = s0 + log_h_diff*(s1 + log_h_diff*(s2 + log_h_diff*(s3
         + log_h_diff*(s4 + log_h_diff*(s5 + log_h_diff*(s6
         + log_h_diff*(s7 + log_h_diff*s8)))))));

   return (exp(log_abs_theta0) - gtc.gsw_t0);
}



/**
==========================================================================
  method: gsw_pt_from_pot_enthalpy_ice (pot_enthalpy_ice)
==========================================================================

   Calculates the potential temperature of ice from the potential enthalpy
   of ice.  The reference sea pressure of both the potential temperature
   and the potential enthalpy is zero dbar.

   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]
   pt0_ice  :  potential temperature of ice (ITS-90)         [ deg C ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_pt_from_pot_enthalpy_ice(double pot_enthalpy_ice)
{
	int iteration;
	double df_dt, f, mod_pot_enthalpy_ice, pt0_cold_ice, recip_df_dt,
	pt0_cold_ice_old, pt0_ice, pt0_ice_old, ptm_cold_ice, ptm_ice;
	double  h00 = -6.320202333358860e5, /*gsw_enthalpy_ice(-gtc.gsw_t0,0)*/
	p0 = 0.0;

	mod_pot_enthalpy_ice = gsw_max(pot_enthalpy_ice,h00);

	if (mod_pot_enthalpy_ice >= -5.1e5)
	{
		/**
		   For input potential enthalpies greater than -5.1e-5, the above part of
         the code gives the output potential temperature of ice accurate to
         1e-13 degrees C.
		*/
		pt0_ice = gsw_pt_from_pot_enthalpy_ice_poly(mod_pot_enthalpy_ice);
		/**
		The variable "df_dt" below is the derivative of the above polynomial
		with respect to pot_enthalpy_ice.  This is the initial value of the
		derivative of the method f.
		*/
		recip_df_dt =
		gsw_pt_from_pot_enthalpy_ice_poly_dh(mod_pot_enthalpy_ice);

		pt0_ice_old = pt0_ice;
		f = gsw_pot_enthalpy_from_pt_ice(pt0_ice_old)
		- mod_pot_enthalpy_ice;
		pt0_ice = pt0_ice_old - f*recip_df_dt;
		ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
		recip_df_dt = 1.0/gsw_cp_ice(ptm_ice,p0);
		pt0_ice = pt0_ice_old - f*recip_df_dt;

	}
	else
	{
		/**
		   For  pot_enthalpy_ice < -5.1e5 (or pt0_ice less than about -100 deg c)
         these temperatures are less than those found in nature on planet earth
		*/
		pt0_cold_ice = gsw_pt0_cold_ice_poly(mod_pot_enthalpy_ice);

		df_dt = gsw_cp_ice(pt0_cold_ice+0.02,p0);
		/**         the heat capacity, cp, is
		  evaluated at 0.02 c greater than usual in order to avoid stability
		  issues and to ensure convergence near zero absolute temperature.
		*/
		for (iteration = 1; iteration <= 6; iteration++)
		{
			pt0_cold_ice_old = pt0_cold_ice;
			f = gsw_pot_enthalpy_from_pt_ice(pt0_cold_ice_old)
			- mod_pot_enthalpy_ice;
			pt0_cold_ice = pt0_cold_ice_old - f/df_dt;
			ptm_cold_ice = 0.5*(pt0_cold_ice + pt0_cold_ice_old);
			df_dt = gsw_cp_ice(ptm_cold_ice+0.02,p0);
			/**note the extra 0.02 c here as well*/
			pt0_cold_ice = pt0_cold_ice_old - f/df_dt;
		}
		pt0_ice = pt0_cold_ice;
	}
	/**
	The potential temerature has a maximum error as listed in the table below.

	   potential temerature error (deg C)  |  @ potential temerature (deg C)
	--------------------------------------|---------------------------------
	       0.012       |     -273.15 to -273.12
	          4 x 10^-4          |     -232.12 to -273.0
	         2.5 x 10^-6         |     -273
	          7 x 10^-9          |     -272
	        3.7 x 10^-10         |     -270
	          6 x 10^-11         |     -268
	         2.5 x 10^11         |     -266
	         3 x 10^-12          |     -260
	         7 x 10^-13          |     -250
	        2.2 x 10^-13         |     -220
	        1.7 x 10^-13         |    >= -160

	  Note.  The above errors in each temperature range are machine precisions
	  for this calculation.
	*/

	return (pt0_ice);
}
/**
==========================================================================
  method: gsw_pt_from_pot_enthalpy_ice_poly (pot_enthalpy_ice)
==========================================================================

   Calculates the potential temperature of ice (whose reference sea
   pressure is zero dbar) from the potential enthalpy of ice.  This is a
   compuationally efficient polynomial fit to the potential enthalpy of
   ice.

   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]
   pt0_ice  :  potential temperature of ice (ITS-90)         [ deg C ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_pt_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice)
{
   double q0 = 2.533588268773218e2,
      q1 = 2.594351081876611e-3,
      q2 = 1.765077810213815e-8,
      q3 = 7.768070564290540e-14,
      q4 = 2.034842254277530e-19,
      q5 = 3.220014531712841e-25,
      q6 = 2.845172809636068e-31,
      q7 = 1.094005878892950e-37;

   /**
      The error of this fit ranges between -5e-5 and 2e-4 deg C over the potential
      temperature range of -100 to 2 deg C, or the potential enthalpy range of
      -5.7 x 10^5 to -3.3 x 10^5 J/kg.
   */

   return (q0
    + pot_enthalpy_ice*(q1 + pot_enthalpy_ice*(q2 + pot_enthalpy_ice*(q3
    + pot_enthalpy_ice*(q4 + pot_enthalpy_ice*(q5 + pot_enthalpy_ice*(q6
    + pot_enthalpy_ice*q7)))))));
}
/**
==========================================================================
  method: gsw_pt_from_pot_enthalpy_ice_poly_dh (pot_enthalpy_ice)
==========================================================================

   Calculates the derivative of potential temperature of ice with respect
   to potential enthalpy.  This is based on the compuationally-efficient
   polynomial fit to the potential enthalpy of ice.

   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]
   dpt0_ice_dh  :  derivative of potential temperature of ice
         with respect to potential enthalpy        [ deg C ]
-------------------------------------------------------------------------
*/
double TeosIce::gsw_pt_from_pot_enthalpy_ice_poly_dh(double pot_enthalpy_ice)
{
   double q1 = 2.594351081876611e-3,
      p2 = 3.530155620427630e-8,
      p3 = 2.330421169287162e-13,
      p4 = 8.139369017110120e-19,
      p5 = 1.610007265856420e-24,
      p6 = 1.707103685781641e-30,
      p7 = 7.658041152250651e-37;

   return (q1 + pot_enthalpy_ice*(p2 + pot_enthalpy_ice*(p3
       + pot_enthalpy_ice*(p4 + pot_enthalpy_ice*(p5 + pot_enthalpy_ice*(p6
       + pot_enthalpy_ice*p7))))));
}
/**
=========================================================================
  method: gsw_pt_from_t_ice (t, p, p_ref)
=========================================================================

   Calculates potential temperature of ice Ih with the general reference
   pressure, p_ref, from in-situ temperature, t.

   A faster gsw routine exists if p_ref is indeed zero dbar.  This routine
   is "gsw_pt0_from_t_ice(t,p)".

   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  :  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   p_ref  :  reference pressure                [ dbar ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_pt_from_t_ice(double t, double p, double p_ref)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */
   int number_of_iterations;
   double dentropy, dentropy_dt, dp,
      pt_ice, pt_ice_old, ptm_ice, true_entropy,

      p1 = -2.259745637898635e-4,
      p2 =  1.486236778150360e-9,
      p3 =  6.257869607978536e-12,
      p4 = -5.253795281359302e-7,
      p5 =  6.752596995671330e-9,
      p6 =  2.082992190070936e-11,

      q1 = -5.849191185294459e-15,
      q2 =  9.330347971181604e-11,
      q3 =  3.415888886921213e-13,
      q4 =  1.064901553161811e-12,
      q5 = -1.454060359158787e-10,
      q6 = -5.323461372791532e-13;
      /**This is the starting polynomial for pt of ice Ih.*/

   dp = p - p_ref;

   pt_ice = t + dp*(p1 + (p + p_ref)*(p2 + p3*t) + t*(p4 + t*(p5 + p6*t)));

   if (pt_ice < -gtc.gsw_t0) pt_ice = -gtc.gsw_t0;

   if (pt_ice < -273.0) pt_ice = pt_ice + 0.05;
   /**
     we add 0.05 to the initial estimate of pt_ice at temps less than
     -273 to ensure that it is never less than -273.15.
   */
   dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice,p_ref);

   true_entropy = -gsw_gibbs_ice_part_t(t,p);

   for (number_of_iterations = 1; number_of_iterations <= 3;
       number_of_iterations++)
   {
       pt_ice_old = pt_ice;
       dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref) - true_entropy;
       pt_ice = pt_ice_old - dentropy/dentropy_dt;
       ptm_ice = 0.5*(pt_ice + pt_ice_old);
       dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref);
       pt_ice = pt_ice_old - dentropy/dentropy_dt;
   }

   if (pt_ice < -273.0)
   {
      pt_ice = t + (p - p_ref)*(q1 + (p + p_ref)*(q2 + q3*t)
         + t*(q4 + t*(q5 + q6*t)));

      dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice+0.01,p_ref);
      /**
         we add 0.01 to the initial estimate of pt_ice used in the derivative
         to ensure that it is never less than -273.15 because the derivative
         approaches zero at absolute zero.
      */
      for (number_of_iterations = 1; number_of_iterations <= 3;
      number_of_iterations++)
      {
         pt_ice_old = pt_ice;
         dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref)
            - true_entropy;
         pt_ice = pt_ice_old - dentropy/dentropy_dt;
         ptm_ice = 0.5*(pt_ice + pt_ice_old);
         ptm_ice = ptm_ice + 0.01;
         dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref);
         pt_ice = pt_ice_old - dentropy/dentropy_dt;
      }

   }
   /**
     For temperatures less than -273.1 degsC the maximum error is less than
     2x10^-7 degsC. For temperatures between -273.1 and 273 the maximum error
     is less than 8x10^-8 degsC, and for temperatures greater than -273 degsC the
     maximum error is 1.5x10^-12 degsC.  These errors are over the whole
     ocean depths with both p and pref varying independently between 0 and
     10,000 dbar, while the in-situ temperature varied independently between
     -273.15 and +2 degsC.
   */
   return (pt_ice);
}


/**
==========================================================================
  method: gsw_rho_ice (t, p)
==========================================================================

   Calculates in-situ density of ice from in-situ temperature and pressure.
   Note that the output, rho_ice, is density, not density anomaly;  that
   is, 1000 kg/m^3 is not subracted from it.

   t   :  in-situ temperature (ITS-90)             [ deg C ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   rho_ice  :  in-situ density of ice (not density anomaly)      [ kg/m^3 ]
---------------------------------------------------------------------------
*/
double TeosIce::gsw_rho_ice(double t, double p)
{
   return (1.0/gsw_gibbs_ice(0,1,t,p));
}

/**
==========================================================================
  method: gsw_sa_freezing_estimate (p, saturation_fraction,*ct,*t)
==========================================================================

  From an estimate of SA from a polynomial in CT and p

--------------------------------------------------------------------------
*/
double TeosIce::gsw_sa_freezing_estimate(double p, double saturation_fraction, double *ct,
   double *t)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */
   double ctx, ctsat, sa,
      /**note that aa = 0.502500117621d0/35.16504*/
      aa = 0.014289763856964,

      bb = 0.057000649899720,

      p0  =  2.570124672768757e-1,
      p1  = -1.917742353032266e1,
      p2  = -1.413382858617969e-2,
      p3  = -5.427484830917552e-1,
      p4  = -4.126621135193472e-4,
      p5  = -4.176407833276121e-7,
      p6  =  4.688217641883641e-5,
      p7  = -3.039808885885726e-8,
      p8  = -4.990118091261456e-11,
      p9  = -9.733920711119464e-9,
      p10 = -7.723324202726337e-12,
      p11 =  7.121854166249257e-16,
      p12 =  1.256474634100811e-12,
      p13 =  2.105103897918125e-15,
      p14 =  8.663811778227171e-19;

   /** A very rough estimate of sa to get the saturated ct */
   if (ct != NULL)
   {
      sa = gsw_max(-(*ct + 9e-4*p)/0.06, 0.0);
      ctx = *ct;
   }
   else if (t != NULL)
   {
      sa = gsw_max(-(*t + 9e-4*p)/0.06, 0.0);
      ctx = gsw_ct_from_t(sa,*t,p);
   }
   else
   {
      return (0.0);
   }
   /**
     CTsat is the estimated value of CT if the seawater were saturated with
     dissolved air, recognizing that it actually has the air fraction
     saturation_fraction; see McDougall, Barker and Feistel, 2014).
   */
   ctsat = ctx - (1.0-saturation_fraction)*
      (1e-3)*(2.4-aa*sa)*(1.0+bb*(1.0-sa/gtc.gsw_sso));

   return (p0 + p*(p2 + p4*ctsat + p*(p5 + ctsat*(p7 + p9*ctsat)
       + p*(p8  + ctsat*(p10 + p12*ctsat) + p*(p11 + p13*ctsat + p14*p))))
       + ctsat*(p1 + ctsat*(p3 + p6*p)));
}
/**
==========================================================================
  method: gsw_sa_freezing_from_ct (ct, p, saturation_fraction)
==========================================================================

   Calculates the Absolute Salinity of seawater at the freezing temperature.
   That is, the output is the Absolute Salinity of seawater, with
   Conservative Temperature CT, pressure p and the fraction
   saturation_fraction of dissolved air, that is in equilibrium
   with ice at the same in situ temperature and pressure.  If the input
   values are such that there is no positive value of Absolute Salinity for
   which seawater is frozen, the output is made a NaN.

   CT    :  Conservative Temperature of seawater (ITS-90)      [ deg C ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction  =  the saturation fraction of dissolved air in
            seawater

   sa_freezing_from_ct  =  Absolute Salinity of seawater when it freezes,
        for given input values of its Conservative Temperature,
        pressure and air saturation fraction.       [ g/kg ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_sa_freezing_from_ct(double ct, double p, double saturation_fraction)
{
   int i_iter, number_of_iterations = 3;
   double ct_freezing_zero_sa, f, ctfreezing_sa, sa, sa_mean, sa_old;
   /**
     This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
     differently in calculating the initial values of both SA and dCT_dSA.
   */
   double sa_cut_off = 2.5;
   /**
     Find CT > CT_freezing_zero_SA.  If this is the case, the input values
     represent seawater that is not frozen (for any positive SA).
   */
   ct_freezing_zero_sa = gsw_ct_freezing(0.0,p,saturation_fraction);
   if (ct > ct_freezing_zero_sa)
       return (cppGSW_INVALID_VALUE);

   /**Form the first estimate of SA from a polynomial in CT and p*/
   sa = gsw_sa_freezing_estimate(p,saturation_fraction,&ct,NULL);
   if (sa < -sa_cut_off)
       return (cppGSW_INVALID_VALUE);
   /**
     Form the first estimate of CTfreezing_SA,
     the derivative of CT_freezing with respect to SA at fixed p.
   */
   sa = gsw_max(sa,0.0);
   gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,
                      &ctfreezing_sa, NULL);
   /**
     For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
     with one based on (CT_freezing_zero_SA - CT).
   */
   if (fabs(sa) < sa_cut_off)
       sa = (ct - ct_freezing_zero_sa)/ctfreezing_sa;
   /**
   ------------------------------------------------------
     Begin the modified Newton-Raphson method to solve
     f = (CT_freezing - CT) = 0 for SA.
   ------------------------------------------------------
   */
   for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
   {
       sa_old = sa;
       f = gsw_ct_freezing(sa,p,saturation_fraction) - ct;
       sa = sa_old - f/ctfreezing_sa;
       sa_mean = 0.5*(sa + sa_old);
       gsw_ct_freezing_first_derivatives(sa_mean,p,saturation_fraction,
               &ctfreezing_sa, NULL);

       sa = sa_old - f/ctfreezing_sa;
   }

   if (gsw_sa_p_inrange(sa,p))
       return (sa);

   return (cppGSW_INVALID_VALUE);
}
/**
==========================================================================
  method: gsw_sa_freezing_from_ct_poly (ct, p, saturation_fraction)
==========================================================================

   Calculates the Absolute Salinity of seawater at the freezing temperature.
   That is, the output is the Absolute Salinity of seawater, with the
   fraction saturation_fraction of dissolved air, that is in equilibrium
   with ice at Conservative Temperature CT and pressure p.  If the input
   values are such that there is no positive value of Absolute Salinity for
   which seawater is frozen, the output is put equal to Nan.

   CT    :  Conservative Temperature (ITS-90)             [ deg C ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction  :  the saturation fraction of dissolved air in
            seawater
   sa_freezing_from_ct  :  Absolute Salinity of seawater when it freezes,
        for given input values of Conservative Temperature
        pressure and air saturation fraction.       [ g/kg ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_sa_freezing_from_ct_poly(double ct, double p, double saturation_fraction)
{
   int i_iter, number_of_iterations = 2;
   double ct_freezing, ct_freezing_zero_sa, dct_dsa, sa, sa_old, sa_mean;
   /**
     This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
     differently in calculating the initial values of both SA and dCT_dSA.
   */
   double sa_cut_off = 2.5;
   /**
     Find CT > CT_freezing_zero_SA.  If this is the case, the input values
     represent seawater that is not frozen (at any positive SA).
   */
   ct_freezing_zero_sa = gsw_ct_freezing_poly(0.0,p,saturation_fraction);
   if (ct > ct_freezing_zero_sa)
       return (cppGSW_INVALID_VALUE);

   /**Form the first estimate of SA from a polynomial in CT and p */
   sa = gsw_sa_freezing_estimate(p,saturation_fraction,&ct,NULL);
   if (sa < -sa_cut_off)
       return (cppGSW_INVALID_VALUE);
   /**
     Form the first estimate of dCT_dSA, the derivative of CT with respect
     to SA at fixed p.
   */
   sa = gsw_max(sa,0.0);
   gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,
                      &dct_dsa, NULL);
   /**
     For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
     with one based on (CT_freezing_zero_SA - CT).
   */
   if (fabs(sa) < sa_cut_off)
       sa = (ct - ct_freezing_zero_sa)/dct_dsa;
   /**
   -----------------------------------------------------------------------
     Begin the modified Newton-Raphson method to solve the root of
     CT_freezing = CT for SA.
   -----------------------------------------------------------------------
   */
   for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
   {
      sa_old = sa;
      ct_freezing = gsw_ct_freezing_poly(sa_old,p,saturation_fraction);
      sa = sa_old - (ct_freezing - ct)/dct_dsa;
      sa_mean = 0.5*(sa + sa_old);
      gsw_ct_freezing_first_derivatives_poly(sa_mean,p,
            saturation_fraction, &dct_dsa, NULL);

      sa = sa_old - (ct_freezing - ct)/dct_dsa;
   }

   if (gsw_sa_p_inrange(sa,p))
       return (sa);

   return (cppGSW_INVALID_VALUE);
}
/**
==========================================================================
  method: gsw_sa_freezing_from_t (t, p, saturation_fraction)
==========================================================================

   Calculates the Absolute Salinity of seawater at the freezing temperature.
   That is, the output is the Absolute Salinity of seawater, with the
   fraction saturation_fraction of dissolved air, that is in equilibrium
   with ice at in-situ temperature t and pressure p.  If the input values
   are such that there is no positive value of Absolute Salinity for which
   seawater is frozen, the output is set to NaN.

   t  :  in-situ Temperature (ITS-90)              [ deg C ]
   p  :  sea pressure                     [ dbar ]
    ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
          seawater
   (i.e., saturation_fraction must be between 0 and 1, and the default
     is 1, completely saturated)
   sa_freezing_from_t  :  Absolute Salinity of seawater when it freezes, for
       given input values of in situ temperature, pressure and
       air saturation fraction.           [ g/kg ]
-----------------------------------------------------------------------------
*/
double TeosIce::gsw_sa_freezing_from_t(double t, double p, double saturation_fraction)
{
   int i_iter, number_of_iterations = 2;
   double f, sa, sa_mean, sa_old, t_freezing_zero_sa, tfreezing_sa;
   /**
     This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
     differently in calculating the initial values of both SA and dCT_dSA.
   */
   double  sa_cut_off = 2.5;
   /**
     Find t > t_freezing_zero_SA.  If this is the case, the input values
     represent seawater that is not frozen (at any positive SA).
   */
   t_freezing_zero_sa = gsw_t_freezing(0.0,p,saturation_fraction);
   if (t > t_freezing_zero_sa)
       return (cppGSW_INVALID_VALUE);
   /**
     This is the inital guess of SA using a purpose-built
     polynomial in CT and p
   */
   sa = gsw_sa_freezing_estimate(p,saturation_fraction,NULL,&t);
   if (sa < -sa_cut_off)
       return (cppGSW_INVALID_VALUE);
   /**
     Form the first estimate of tfreezing_SA, the derivative of
     CT_freezing with respect to SA at fixed p.
   */
   sa = gsw_max(sa,0.0);
   gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
                     &tfreezing_sa,NULL);
   /**     For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
     with one based on (t_freezing_zero_SA - t).
   */
   if (fabs(sa) < sa_cut_off)
       sa = (t - t_freezing_zero_sa)/tfreezing_sa;
   /**
   -----------------------------------------------------------------------
     Begin the modified Newton-Raphson method to find the root of
     f = (t_freezing - t) = 0 for SA.
   -----------------------------------------------------------------------
   */
   for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
   {
      sa_old = sa;
      f = gsw_t_freezing(sa_old,p,saturation_fraction) - t;
      sa = sa_old - f/tfreezing_sa;
      sa_mean = 0.5*(sa + sa_old);
      gsw_t_freezing_first_derivatives(sa_mean,p,saturation_fraction,
                    &tfreezing_sa, NULL);

      sa = sa_old - f/tfreezing_sa;
   }

   if (gsw_sa_p_inrange(sa,p))
       return (sa);

   return (cppGSW_INVALID_VALUE);
}
/**
==========================================================================
  method: gsw_sa_freezing_from_t_poly (t, p, saturation_fraction)
==========================================================================

   Calculates the Absolute Salinity of seawater at the freezing temperature.
   That is, the output is the Absolute Salinity of seawater, with the
   fraction saturation_fraction of dissolved air, that is in equilibrium
   with ice at in-situ temperature t and pressure p.  If the input values
   are such that there is no positive value of Absolute Salinity for which
   seawater is frozen, the output is put equal to Nan.

   t  :  in-situ Temperature (ITS-90)              [ deg C ]
   p  :  sea pressure                     [ dbar ]
    ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
          seawater
   sa_freezing_from_t_poly  :  Absolute Salinity of seawater when it freezes,
       for given input values of in situ temperature, pressure and
       air saturation fraction.           [ g/kg ]
-----------------------------------------------------------------------------
*/
double TeosIce::gsw_sa_freezing_from_t_poly(double t, double p, double saturation_fraction)
{
   int i_iter, number_of_iterations = 5;
   double dt_dsa, sa, sa_old, sa_mean, t_freezing, t_freezing_zero_sa;
   /**
     This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
     differently in calculating the initial values of both SA and dCT_dSA.
   */
   double  sa_cut_off = 2.5;
   /**
     Find t > t_freezing_zero_SA.  If this is the case, the input values
     represent seawater that is not frozen (at any positive SA).
   */
   t_freezing_zero_sa = gsw_t_freezing_poly(0.0,p,saturation_fraction);
   if (t > t_freezing_zero_sa)
       return (cppGSW_INVALID_VALUE);
   /**
     This is the inital guess of SA using a purpose-built
     polynomial in CT and p
   */
   sa = gsw_sa_freezing_estimate(p,saturation_fraction,NULL,&t);
   if (sa < -sa_cut_off)
       return (cppGSW_INVALID_VALUE);
   /**
      Form the first estimate of dt_dSA, the derivative of t with respect
      to SA at fixed p.
   */
   sa = gsw_max(sa,0.0);
   gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction,
                     &dt_dsa, NULL);
   /**
     For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
     with one based on (t_freezing_zero_SA - t).
   */
   if (fabs(sa) < sa_cut_off)
       sa = (t - t_freezing_zero_sa)/dt_dsa;
   /**
   -----------------------------------------------------------------------
     Begin the modified Newton-Raphson method to find the root of
     t_freezing = t for SA.
   -----------------------------------------------------------------------
   */
   for (i_iter = 1; i_iter <= number_of_iterations; i_iter++)
   {
      sa_old = sa;
      t_freezing = gsw_t_freezing_poly(sa_old,p,saturation_fraction);
      sa = sa_old - (t_freezing - t)/dt_dsa;
      sa_mean = 0.5*(sa + sa_old);
      gsw_t_freezing_first_derivatives_poly(sa_mean,p,saturation_fraction,
                         &dt_dsa, NULL);

       sa = sa_old - (t_freezing - t)/dt_dsa;
   }

   if (gsw_sa_p_inrange(sa,p))
       return (sa);

   return (cppGSW_INVALID_VALUE);
}

/**
==========================================================================
  method: gsw_sa_p_inrange (sa, p)
==========================================================================

   Check for any values that are out of the TEOS-10 range ...

   SA    :  Absolute Salinity               [ g/kg ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
---------------------------------------------------------------------------
*/
int TeosIce::gsw_sa_p_inrange(double sa, double p)
{
   if (p > 10000.0 || sa > 120.0 ||
      (p + sa*71.428571428571402) > 13571.42857142857)
   {
       return (0);
   }

   return (1);
}
/**
==========================================================================
method: gsw_seaice_fraction_to_freeze_seawater(sa,ct,p,sa_seaice,t_seaice,
*sa_freeze, *ct_freeze,*w_seaice)
==========================================================================

   Calculates the mass fraction of sea ice (mass of sea ice divided by mass
   of sea ice plus seawater), which, when melted into seawater having the
   properties (SA,CT,p) causes the final seawater to be at the freezing
   temperature.  The other outputs are the Absolute Salinity and
   Conservative Temperature of the final seawater.

   SA   :  Absolute Salinity of seawater            [ g/kg ]
   CT   :  Conservative Temperature of seawater (ITS-90)     [ deg C ]
   p    :  sea pressure                   [ dbar ]
        ( i.e. absolute pressure - 10.1325 dbar )
   SA_seaice :  Absolute Salinity of sea ice, that is, the mass fraction of
           salt in sea ice, expressed in g of salt per kg of sea ice.[ g/kg ]
   t_seaice  :  in-situ temperature of the sea ice at pressure p (ITS-90)[ deg C ]
   SA_freeze  :  Absolute Salinity of seawater after the mass fraction of
       sea ice, w_seaice, at temperature t_seaice has melted into
       the original seawater, and the final mixture is at the
       freezing temperature of seawater.       [ g/kg ]
   CT_freeze  :  Conservative Temperature of seawater after the mass
       fraction, w_seaice, of sea ice at temperature t_seaice has
       melted into the original seawater, and the final mixture
       is at the freezing temperature of seawater.      [ deg C ]
   w_seaice   :  mass fraction of sea ice, at SA_seaice and t_seaice,
       which, when melted into seawater at (SA,CT,p) leads to the
       final mixed seawater being at the freezing temperature.
       This output is between 0 and 1.       [unitless]
--------------------------------------------------------------------------
*/
void TeosIce::gsw_seaice_fraction_to_freeze_seawater(double sa, double ct, double p,
   double sa_seaice, double t_seaice, double *sa_freeze, double *ct_freeze,
   double *w_seaice)
{
	int number_of_iterations;
	double ctf, ctf_mean, ctf_old, ctf_plus1, ctf_zero,
	   dfunc_dsaf, func, func_plus1, func_zero, h, h_brine,
	   h_ih, sa_freezing, saf, saf_mean, saf_old,
	   salt_ratio, tf_sa_seaice, h_hat_sa, h_hat_ct, ctf_sa,
	   sa0 = 0.0, saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**The seawater ct input is below the freezing temp*/
		*sa_freeze = *ct_freeze = *w_seaice = cppGSW_INVALID_VALUE;
		return;
	}

	tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction) - 1e-6;
	if (t_seaice > tf_sa_seaice)
	{
		/**
		The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
		some ice Ih in the sea ice.   Without this buffer, that is if t_seaice
		is allowed to be exactly equal to tf_sa_seaice, the sea ice is
		actually 100% brine at Absolute Salinity of SA_seaice.
		*/
		*sa_freeze = *ct_freeze = *w_seaice = cppGSW_INVALID_VALUE;
		return;
	}

	sa_freezing = gsw_sa_freezing_from_t(t_seaice,p,saturation_fraction);
	if (sa_freezing > cppGSW_ERROR_LIMIT)
	{
		*sa_freeze = *ct_freeze = *w_seaice = cppGSW_INVALID_VALUE;
		return;
	}
	h_brine = gsw_enthalpy_t_exact(sa_freezing,t_seaice,p);
	salt_ratio = sa_seaice/sa_freezing;

	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);

	ctf_plus1 = gsw_ct_freezing(sa+1.0,p,saturation_fraction);
	func_plus1 = (sa - sa_seaice)
	   *(gsw_enthalpy_ct_exact(sa+1.0,ctf_plus1,p)
      - h) - (h - h_ih) + salt_ratio*(h_brine - h_ih);

	ctf_zero = gsw_ct_freezing(sa0,p,saturation_fraction);
	func_zero = (sa - sa_seaice)
	   *(gsw_enthalpy_ct_exact(sa0,ctf_zero,p) - h)
	   + sa*((h - h_ih) - salt_ratio*(h_brine - h_ih));

	saf = -(sa+1.0)*func_zero/(func_plus1 - func_zero);
	/**initial guess of saf*/
	ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	gsw_enthalpy_first_derivatives_ct_exact(saf,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_ct_freezing_first_derivatives(saf,p,saturation_fraction,
												 &ctf_sa, NULL);

	dfunc_dsaf = (sa - sa_seaice)*(h_hat_sa + h_hat_ct*ctf_sa)
	- (h - h_ih) + salt_ratio*(h_brine - h_ih);

	for (number_of_iterations = 1; number_of_iterations <= 4;
			number_of_iterations++)
	{
		saf_old = saf;
		ctf_old = ctf;
		func = (sa - sa_seaice)
		   *(gsw_enthalpy_ct_exact(saf_old,ctf_old,p) - h)
		   - (saf_old - sa)*((h - h_ih) - salt_ratio*(h_brine - h_ih));

		saf = saf_old - func/dfunc_dsaf;
		saf_mean = 0.5*(saf + saf_old);
		ctf_mean = gsw_ct_freezing(saf_mean,p,saturation_fraction);
		gsw_enthalpy_first_derivatives_ct_exact(saf_mean,ctf_mean,p,
															 &h_hat_sa,&h_hat_ct);

		gsw_ct_freezing_first_derivatives(saf_mean,p,saturation_fraction,
													 &ctf_sa, NULL);

		dfunc_dsaf = (sa - sa_seaice)*(h_hat_sa + h_hat_ct*ctf_sa)
		- (h - h_ih) + salt_ratio*(h_brine - h_ih);
		saf = saf_old - func/dfunc_dsaf;
		ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	}
	/**
	 After these 4 iterations of this modified Newton-Raphson method, the
	  errors in SA_freeze is less than 1.5x10^-12 g/kg, in CT_freeze is less than
	  2x10^-13 deg C and in w_seaice is less than 2.8x10^-13 which represent machine
	  precision for these calculations.
	*/
	*sa_freeze = saf;
	*ct_freeze = ctf;

	*w_seaice = (h - gsw_enthalpy_ct_exact(*sa_freeze,*ct_freeze,p)) /
	(h - h_ih - salt_ratio*(h_brine - h_ih));

	return;
}

/**
==========================================================================
  method: gsw_sound_speed_ice (t, p)
==========================================================================

   Calculates the compression speed of sound in ice.

   t   :  in-situ temperature (ITS-90)             [ deg C ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   sound_speed_ice  :  compression speed of sound in ice       [ m/s ]
--------------------------------------------------------------------------
*/
double TeosIce::gsw_sound_speed_ice(double t, double p)
{
   double gi_tp, gi_tt;

   gi_tt = gsw_gibbs_ice(2,0,t,p);
   gi_tp = gsw_gibbs_ice(1,1,t,p);

   return (gsw_gibbs_ice(0,1,t,p) *
          sqrt(gi_tt/(gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))));
}
/**
=========================================================================
  method: gsw_t_from_pt0_ice (pt0_ice, p)
=========================================================================

   Calculates in-situ temperature from the potential temperature of ice Ih
   with reference pressure, p_ref, of 0 dbar (the surface), and the
   in-situ pressure.

   pt0_ice  :  potential temperature of ice Ih with reference pressure of
          zero dbar (ITS-90)             [ deg C ]
   p  :  sea pressure               [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
--------------------------------------------------------------------------
*/
double TeosIce::gsw_t_from_pt0_ice(double pt0_ice, double p)
{
   double  p0 = 0.0;

   return (gsw_pt_from_t_ice(pt0_ice,p0,p));
}
/**
==========================================================================
method: gsw_t_freezing_exact (sa,p,saturation_fraction)
==========================================================================

   Calculates the in-situ temperature at which seawater freezes. The
   in-situ temperature freezing point is calculated from the exact
   in-situ freezing temperature which is found by a modified Newton-Raphson
   iteration (McDougall and Wotherspoon, 2013) of the equality of the
   chemical potentials of water in seawater and in ice.

   An alternative GSW function, gsw_t_freezing_poly, it is based on a
   computationally-efficient polynomial, and is accurate to within -5e-4 K
   and 6e-4 K, when compared with this function.

   SA    :  Absolute Salinity               [ g/kg ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
          seawater
   (i.e., saturation_fraction must be between 0 and 1, and the default
     is 1, completely saturated)
   t_freezing : in-situ temperature at which seawater freezes.    [ deg C ]
---------------------------------------------------------------------------
*/
double TeosIce::gsw_t_freezing_exact (double sa, double p, double saturation_fraction)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */

	double df_dt, tf, tfm, tf_old, f, return_value;
	int polynomial=1;

	/** The initial value of t_freezing_exact (for air-free seawater) */
	tf = gsw_t_freezing_poly(sa,p,saturation_fraction,polynomial);

	df_dt = 1e3*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p) -
	gsw_gibbs_ice(1,0,tf,p);
	/**
	  df_dt here is the initial value of the derivative of the function f whose
	  zero (f = 0) we are finding (see Eqn. (3.33.2) of IOC et al (2010)).
	*/

	tf_old = tf;
	f = 1e3*gsw_chem_potential_water_t_exact(sa,tf_old,p) -
	gsw_gibbs_ice(0,0,tf_old,p);

	tf = tf_old - f/df_dt;
	tfm = 0.5*(tf + tf_old);

	df_dt = 1e3*gsw_t_deriv_chem_potential_water_t_exact(sa,tfm,p) -
	gsw_gibbs_ice(1,0,tfm,p);

	tf = tf_old - f/df_dt;

	tf_old = tf;

	f = 1e3*gsw_chem_potential_water_t_exact(sa,tf_old,p) -
	gsw_gibbs_ice(0,0,tf_old,p);

	tf = tf_old - f/df_dt;

	/** Adjust for the effects of dissolved air */
	return_value = tf -
	saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gtc.gsw_sso));

	return (return_value);
}
/**
==========================================================================
method: gsw_ct_freezing_exact(sa,p,saturation_fraction)
==========================================================================

   Calculates the Conservative Temperature at which seawater freezes.  The
   Conservative Temperature freezing point is calculated from the exact
   in-situ freezing temperature which is found by a modified Newton-Raphson
   iteration (McDougall and Wotherspoon, 2013) of the equality of the
   chemical potentials of water in seawater and in ice.

   An alternative GSW function, gsw_CT_freezing_poly, it is based on a
   computationally-efficient polynomial, and is accurate to within -5e-4 K
   and 6e-4 K, when compared with this function.

   SA    :  Absolute Salinity               [ g/kg ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
          seawater
   CT_freezing : Conservative Temperature at freezing of seawater [ deg C ]
---------------------------------------------------------------------------
*/
double TeosIce::gsw_ct_freezing_exact(double sa, double p, double saturation_fraction)
{
	double t_freezing;

	t_freezing = gsw_t_freezing_exact(sa,p,saturation_fraction);

	return (gsw_ct_from_t(sa,t_freezing,p));
}

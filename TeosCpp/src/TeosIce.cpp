#include "TeosIce.h"

/******************************************
  TeosIce.h      Version 1.05
  by Randall Kent Whited
  rkwhited@gmail.com
  ---------------------------------------
  All copyrights and all license issues
  are the same as previous versions
  ---------------------------------------
  OR converted from .m Matlab files in
  the Teos-10 Matlab version at
  http://www.teos-10.org
*******************************************/

/** class constructor : base class constructor */
TeosIce::TeosIce(): TeosBase()
{
   try
   {
      statusOk = true;
   }
   catch(const unsigned e)
   {
      statusOk = !(e == badBaseAlloc);
   }
   catch (...)
   {
      statusOk = false;
   }
}

/** descendant class destructor */
TeosIce::~TeosIce()
{
}

/**************************************************************************
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
*************************************************************************/

double TeosIce::gsw_ct_freezing_exact(double sa,double p,double saturation_fraction)
{
	double t_freezing;
	t_freezing = gsw_t_freezing_exact(sa,p,saturation_fraction);
	return (gsw_ct_from_t(sa,t_freezing,p));
}


/**************************************************************************
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
*************************************************************************/

double *TeosIce::gsw_geo_strf_dyn_height_pc(double *sa,double *ct,double *delta_p,int n_levels,double *geo_strf_dyn_height_pc,double *p_mid)
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

/**************************************************************************
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
*************************************************************************/

void TeosIce::gsw_ipv_vs_fnsquared_ratio(double *sa,double *ct,double *p,double p_ref,int nz,double *ipv_vs_fnsquared_ratio,double *p_mid)
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

/**************************************************************************
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
*************************************************************************/

double TeosIce::gsw_melting_ice_equilibrium_sa_ct_ratio_poly(double sa,double p)
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

/**************************************************************************
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
*************************************************************************/

void TeosIce::gsw_melting_ice_into_seawater(double sa,double ct,double p,double w_ih,double t_ih,double *sa_final,double *ct_final,double *w_ih_final)
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


/**************************************************************************
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
*************************************************************************/

void TeosIce::gsw_nsquared(double *sa,double *ct,double *p,double *lat,int nz,double *n2,double *p_mid)
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




/**************************************************************************
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
*************************************************************************/

double TeosIce::gsw_t_freezing_exact(double sa,double p,double saturation_fraction)
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
/** moved here from .m Matlab conversions to C++ */





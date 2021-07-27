#ifndef TEOSICE_H
#define TEOSICE_H


/******************************************
  "TeosIce.h" Version 1.05
  by Randall Kent Whited
  rkwhited@gmail.com
  --------------------------------
  This is a modification of the
  TEOS-10 file "gswteos-10.h"
  as adapted from the TEOS-10
  C version 3.05 http://www.teos-10.org
  ---------------------------------------
  OR converted from .m Matlab files in
  the Teos-10 Matlab version at
  http://www.teos-10.org
  ---------------------------------------
  All copyrights and all license issues
  are the same as for the C version
*******************************************/

/** Base-class object header */
#include "TeosBase.h"

using namespace std;

/****************************
  "TeosIce" C++ Class
  -------------------
  The methods of this class
  focus more on ice related
  issues (the "TeosSea"
  class focuses more on
  non-frozen seawater).
  ---------------------
  TEOS descriptions of the
  various methods are in
  the relevant ".cpp" file.
****************************/
class TeosIce:  public TeosBase
{
public:
    TeosIce();
    virtual ~TeosIce();

		double gsw_ct_freezing_exact(double sa,double p,double saturation_fraction);
		double *gsw_geo_strf_dyn_height_pc(double *sa,double *ct,double *delta_p,int n_levels,double *geo_strf_dyn_height_pc,double *p_mid);
		void gsw_ipv_vs_fnsquared_ratio(double *sa,double *ct,double *p,double p_ref,int nz,double *ipv_vs_fnsquared_ratio,double *p_mid);
		double gsw_melting_ice_equilibrium_sa_ct_ratio_poly(double sa,double p);
		void gsw_melting_ice_into_seawater(double sa,double ct,double p,double w_ih,double t_ih,double *sa_final,double *ct_final,double *w_ih_final);
		void gsw_nsquared(double *sa,double *ct,double *p,double *lat,int nz,double *n2,double *p_mid);
		double gsw_t_freezing_exact(double sa,double p,double saturation_fraction);

		/** .m Matlab version conversions to C++ */



};
#endif

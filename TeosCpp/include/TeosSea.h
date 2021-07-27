#ifndef TEOSSEA_H
#define TEOSSEA_H


/******************************************
  "TeosSea.h" Version 1.05
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

/***************************
  TEOS C++ CLASS
  ------------------------
  TEOS descriptions of the
  methods are in the .cpp
  file.
***************************/
class TeosSea:  public TeosBase
{
	public:
		TeosSea();
		virtual ~TeosSea();

		void gsw_ct_second_derivatives(double sa,double pt,double *ct_sa_sa,double *ct_sa_pt,double *ct_pt_pt);
		double *gsw_geo_strf_dyn_height(double *sa,double *ct,double *p,double p_ref,int n_levels,double *dyn_height);
		void gsw_p_sequence(double p1,double p2,double max_dp_i,double *pseq,int *nps);
		void gsw_pt_second_derivatives(double sa,double ct,double *pt_sa_sa,double *pt_sa_ct,double *pt_ct_ct);
		int gsw_refine_grid_for_dh(double *p,double p_ref,int nz,double dp,double *p_i,int ni_max,int *p_indices,int *p_ref_ind_ptr);
		void gsw_rho_first_derivatives(double sa,double ct,double p,double *drho_dsa,double *drho_dct,double *drho_dp);
		void gsw_rho_second_derivatives_wrt_enthalpy(double sa,double ct,double p,double *rho_sa_sa,
				double *rho_sa_h,double *rho_h_h);

};
#endif

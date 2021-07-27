#ifndef TEOSBASE_H
#define TEOSBASE_H

/*****************************************
  "TeosBase.h" Version 1.0
  by Randall Kent Whited
  rkwhited@gmail.com
  --------------------------------
  This is a modification of the
  TEOS-10 file "gswteos-10.h"
  as adapted from the TEOS-10
  C version 3.05 http://www.teos-10.org
  ---------------------------------------
  All copyrights and all license issues
  are the same as for the C version
******************************************/

/*************************************
  "TeosBase" is the base or ancestor
  class for the "TeosSea" and the
  "TeosIce" descendant classes.
**************************************/

/** std headers */
#include <string.h>
#include <math.h>
#include <complex>
#include <errno.h>

using namespace std;

/****************************
  The "gsw_structs.h" file
  has been converted to C++
  encapsulated format by
  placing data within
  structs.
*****************************/
#include "gsw_structs.h"

/************************
  Header file to the
  SAAR object for access
  to SAAR data arrays.
************************/
#include "SaarDataHandler.h"

/** error return values from various gsw methods */
const double cppGSW_INVALID_VALUE = 9e15;
const double cppGSW_ERROR_LIMIT = 1e10;

/********************************
  TEOS-10 data in the "TeosBase"
  class is handled by an object
  of the "SaarDataHandler" class
  and various C++ structs.
  ------------------------------
  Descendant classes have full
  access to the "gsw_structs"
  and "gsw_saar_data"
*********************************/

class TeosBase
{
	public:
		TeosBase();
		virtual ~TeosBase();

		/** these complex numbers are initialized in the constructor */
		complex<double> t1;
		complex<double> t2;
		complex<double> r1;
		complex<double> r20;
		complex<double> r21;
		complex<double> r22;

		/** Specialized data handler */
		SaarDataHandler *pSaarData;

		/********************************
		  struct data used by both
		  "TeosSea" & "TeosIce" objects
		  are initialized in the
		  "TeosBase" constructor.
		********************************/
		GSW_TEOS10_CONSTANTS gtc;
		GSW_SPECVOL_COEFFICIENTS gsvco;
		GSW_SP_COEFFICIENTS gspc;
		GSW_BALTIC_DATA gbalt;
		GSW_SAAR_DATA gsaar;
		GSW_FREEZING_POLY_COEFFICIENTS gfpc;
		GSW_GIBBS_ICE_COEFFICIENTS gic;

		/*************************************
		  Previous "C functions" have become
		  "C++ methods" that are used by both
		  the "TeosSea" & "TeosIce" objects
		  via the "TeosBase" base class or
		  their own methods.
		  ----------------------------------
		  Their C library names as well as
		  functionality, as described in
		  the "TeosSea.cpp" and "TeosIce.cpp"
		  files (descriptions by TEOS org),
		  have been preserved.
		**************************************/

		/** 'A' */
		void   gsw_add_barrier(double *input_data, double lon, double lat,
									  double long_grid, double lat_grid, double dlong_grid,
									  double dlat_grid, double *output_data);

		void   gsw_add_mean(double *data_in, double *data_out);
		double gsw_adiabatic_lapse_rate_ice(double t, double p);
		double gsw_alpha(double sa, double ct, double p);
		/** 'B' */
		double gsw_beta(double sa, double ct, double p);
		/** 'C' */
		double gsw_chem_potential_water_t_exact(double sa, double t, double p);
		double gsw_cp_ice(double t, double p);

		void gsw_ct_first_derivatives_wrt_t_exact(double sa, double t, double p,
				double *ct_sa_wrt_t, double *ct_t_wrt_t,
				double *ct_p_wrt_t);

		double gsw_ct_freezing(double sa, double p, double saturation_fraction);

		void gsw_ct_freezing_first_derivatives(double sa, double p,
															double saturation_fraction,
															double *ctfreezing_sa,
															double *ctfreezing_p);

		void gsw_ct_freezing_first_derivatives_poly(double sa, double p,
				double saturation_fraction,
				double *ctfreezing_sa,
				double *ctfreezing_p);

		double gsw_ct_freezing_poly(double sa, double p, double saturation_fraction);
		double gsw_ct_from_enthalpy(double sa, double h, double p);
		double gsw_ct_from_enthalpy_exact(double sa, double h, double p);
		double gsw_ct_from_pt(double sa, double pt);
		void gsw_ct_from_rho(double rho, double sa, double p, double *ct, double *ct_multiple);
		double gsw_ct_from_t(double sa, double t, double p);
		double gsw_ct_maxdensity(double sa, double p);
		/** 'D' */
		double gsw_deltasa_atlas(double p, double lon, double lat);
		double gsw_dilution_coefficient_t_exact(double sa, double t, double p);
		double gsw_dynamic_enthalpy(double sa, double ct, double p);
		double gsw_enthalpy(double sa, double ct, double p);
		/** 'E' */
		double gsw_enthalpy_ct_exact(double sa, double ct, double p);

		void gsw_enthalpy_first_derivatives(double sa, double ct, double p,
														double *h_sa, double *h_ct);

		void gsw_enthalpy_first_derivatives_ct_exact(double sa, double ct, double p,
				double *h_sa, double *h_ct);

		double gsw_enthalpy_ice(double t, double p);

		void gsw_enthalpy_second_derivatives_ct_exact(double sa, double ct, double p,
				double *h_sa_sa, double *h_sa_ct,
				double *h_ct_ct);

		double gsw_enthalpy_sso_0(double p);
		double gsw_enthalpy_t_exact(double sa, double t, double p);

		void gsw_entropy_first_derivatives(double sa, double ct,
													  double *eta_sa, double *eta_ct);

		double gsw_entropy_from_ct(double sa, double ct);
		double gsw_entropy_ice(double t, double p);
		double gsw_entropy_part(double sa, double t, double p);
		double gsw_entropy_part_zerop(double sa, double pt0);

		void gsw_entropy_second_derivatives(double sa, double ct, double *eta_sa_sa,
														double *eta_sa_ct, double *eta_ct_ct);

		/** 'F' */
		double gsw_fdelta(double p, double lon, double lat);

		void gsw_frazil_properties(double sa_bulk, double h_bulk, double p,
											double *sa_final, double *ct_final, double *w_ih_final);

		void gsw_frazil_properties_potential(double sa_bulk, double h_pot_bulk,
														 double p, double *sa_final, double *ct_final,
														 double *w_ih_final);

		void gsw_frazil_properties_potential_poly(double sa_bulk, double h_pot_bulk,
				double p, double *sa_final, double *ct_final,
				double *w_ih_final);

		void gsw_frazil_ratios_adiabatic(double sa, double p, double w_ih,
													double *dsa_dct_frazil, double *dsa_dp_frazil,
													double *dct_dp_frazil);

		void gsw_frazil_ratios_adiabatic_poly(double sa, double p, double w_ih,
														  double *dsa_dct_frazil, double *dsa_dp_frazil,
														  double *dct_dp_frazil);

		/** 'G' */
		double gsw_gibbs(int ns, int nt, int np, double sa, double t, double p);
		double gsw_gibbs_ice (int nt, int np, double t, double p);
		double gsw_gibbs_ice_part_t(double t, double p);
		double gsw_gibbs_ice_pt0(double pt0);
		double gsw_gibbs_ice_pt0_pt0(double pt0);
		double gsw_gibbs_pt0_pt0(double sa, double pt0);
		double gsw_grav(double lat, double p);
		/** 'H' */
		double gsw_hill_ratio_at_sp2(double t);
		/** 'I' */
		double gsw_internal_energy(double sa, double ct, double p);
		/** 'K' */
		double gsw_kappa(double sa, double ct, double p);
		double gsw_kappa_t_exact(double sa, double t, double p);
		/** 'L' */
		double gsw_latentheat_melting(double sa, double p);

		void gsw_linear_interp_sa_ct(double *sa, double *ct, double *p,
											  int np, double *p_i, int npi,
											  double *sa_i, double *ct_i);

		/** 'O' */
		double gsw_o2sol(double sa, double ct, double p, double lon, double lat);
		double gsw_o2sol_sp_pt(double sp, double pt);

		/** 'P' */
		double gsw_p_from_z(double z, double lat,
		                     double geo_strf_dyn_height=0.0, /** default value */
								   double sea_surface_geopotential=0.0); /** default value */

		int gsw_pchip_derivs(double *x, double *y, int n, double *d);
		double gsw_pchip_edge_case(double h0, double h1, double m0, double m1);
		double gsw_pot_enthalpy_from_pt_ice(double pt0_ice);
		double gsw_pot_enthalpy_ice_freezing(double sa, double p);

		void gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa, double p,
				double *pot_enthalpy_ice_freezing_sa,
				double *pot_enthalpy_ice_freezing_p);

		void gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(double sa, double p,
				double *pot_enthalpy_ice_freezing_sa,
				double *pot_enthalpy_ice_freezing_p);

		double gsw_pot_enthalpy_ice_freezing_poly(double sa, double p);
		double gsw_pt_from_ct(double sa, double ct);
		double gsw_pt_from_t(double sa, double t, double p, double p_ref);
		double gsw_pt0_from_t(double sa, double t, double p);
		double gsw_pt0_from_t_ice(double t, double p);
		/** 'R' */
		double gsw_rho(double sa, double ct, double p);

		void gsw_rho_alpha_beta (double sa, double ct, double p,
										 double *rho, double *alpha, double *beta);

		void gsw_rr68_interp_section(int sectnum, double *sa, double *ct, double *p,
											  int mp, int nsect, double *ip_sect, int *ip_isect,
											  double *p_i, double *sa_i, double *ct_i);

		/** 'S' */
		double gsw_sa_from_sp(double sp, double p, double lon, double lat);
		double gsw_sa_from_sp_baltic(double sp, double lon, double lat);
		double gsw_sa_from_sstar(double sstar, double p,double lon,double lat);
		double gsw_saar(double p, double lon, double lat);

		int gsw_sgn(double x);

		double gsw_sp_from_sa(double sa, double p, double lon, double lat);
		double gsw_sp_from_sa_baltic(double sa, double lon, double lat);
		double gsw_sp_from_c(double c, double t, double p);
		double gsw_specvol(double sa, double ct, double p);
		double gsw_specvol_ice(double t, double p);
		double gsw_specvol_t_exact(double sa, double t, double p);
		double gsw_specvol_sso_0(double p);
		double gsw_sr_from_sp(double sp);
		double gsw_sstar_from_sa(double sa, double p, double lon, double lat);
		double gsw_sstar_from_sp(double sp, double p, double lon, double lat);
		double gsw_sum(double *x, int n);
		/** 'T' */
		double gsw_t_deriv_chem_potential_water_t_exact(double sa, double t, double p);
		double gsw_t_freezing(double sa, double p, double saturation_fraction);

		void gsw_t_freezing_first_derivatives(double sa, double p, double saturation_fraction,
														  double *tfreezing_sa, double *tfreezing_p);

		void gsw_t_freezing_first_derivatives_poly(double sa, double p,
				double saturation_fraction,
				double *tfreezing_sa, double *tfreezing_p);

		double gsw_t_freezing_poly(double sa, double p,
											double saturation_fraction);

		double gsw_t_freezing_poly(double sa, double p,
											double saturation_fraction,
											int polynomial);

		double gsw_t_from_ct(double sa, double ct, double p);
		/** 'U' */
		int    gsw_util_indx(double *x, int n, double z);

		void   gsw_util_sort_real(double *rarray, int nx, int *iarray);
		double gsw_util_xinterp1(double *x, double *y, int n, double x0);

		/** 'Z' */
		double gsw_z_from_p(double p, double lat,
								  double geo_strf_dyn_height=0.0, /** default value */
								  double sea_surface_geopotential=0.0); /** default value */

	protected:

		/** method/member function to simplify code */
		int gsw_util_indxPref(int n, double z);

		/** process values for various Teos methods/member functions */
		const int cppINTERP_METHOD_LINEAR = 1;
		const int cppINTERP_METHOD_PCHIP = 2;

		/********************************
		'INLINE' methods/member functions
		*********************************/
		inline double gsw_max(double a, double b)
		{
			if (a > b) return a;
			else return b;
		};
		inline unsigned gsw_max(unsigned a, unsigned b)
		{
			if (a > b) return a;
			else return b;
		};
		inline int gsw_max(int a, int b)
		{
			if (a > b) return a;
			else return b;
		};

		inline double gsw_min(double a, double b)
		{
			if (a < b) return a;
			else return b;
		};

		inline unsigned gsw_min(unsigned a, unsigned b)
		{
			if (a < b) return a;
			else return b;
		};

		inline int gsw_min(int a, int b)
		{
			if (a < b) return a;
			else return b;
		};

	private: /** method/member functions used only by TeosBase */

		/** supplements gsw_alpha clarity */
		double gsw_get_v_ct_part(double xs, double ys, double z);

		/** supplements gsw_ct_from_pt clarity */
		double gsw_get_pot_ent(double x2, double x, double y);

		/** supplements gsw_dynamic_enthalpy clarity */
		double gsw_get_dyn_enth_prt(double xs, double ys, double z);

		/** supplements gsw_entropy_part clarity */
		double gsw_get_g03plusg08(double z, double y, double x2, double x);

		/** supplements gsw_dilution_coefficient_t_exact */
      double gsw_g08_d(double x2, double z, double y, double x);

      /** supplement to gsw_chem_potential_water_t_exact (sa, t, p) */
      double gsw_get_g03plusg08_h(double z, double y, double x2, double x);

		/**
		   these next ten methods/member functions supplement
		   the gsw_gibbs method/member function's clarity
      */
		double gsw_gibbs_g03plusg08_a(double sa, double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_b(double sa, double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_c(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_d(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_e(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_f(double z, double y, double x2, double x);
		double gsw_gibbs_g03plusg08_g(double z, double y, double x2, double x);

		double gsw_gibbs_g08_a(double sa, double z, double y, double x);
		double gsw_gibbs_g08_b(double z, double y, double x);
		double gsw_gibbs_g08_c(double sa, double z, double y, double x);
};

#endif // TEOSBASE_H

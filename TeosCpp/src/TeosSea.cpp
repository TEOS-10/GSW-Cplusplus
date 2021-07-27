#include "TeosSea.h"

/*****************************************
  "TeosSea.cpp" Version 1.05
  by Randall Kent Whited
  rkwhited@gmail.com
  ---------------------------------------
  All copyrights and all license issues
  are the same as for previous versions
******************************************/

/*****************************************
  Most C function calls in the original
  TEOS-10 "gsw_oceanographic_toolbox.c"
  and "gsw_saar.c" source code file that
  have been used are in the base-class
  "TeosBase.cpp" or in the descendant
  classes "TeosIce.cpp" or "TeosSea.cpp".
  ---------------------------------------
  Methods used by both "TeosIce",
  "TeosSea" and "Teos_Matlab" are
  implemented in the "TeosBase" class.
*****************************************/

/** descendant class constructor : base class constructor*/
TeosSea::TeosSea(): TeosBase()
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
TeosSea::~TeosSea()
{
}


/**************************************************************************
==========================================================================
method: gsw_ct_second_derivatives (sa, pt, ct_sa_sa, ct_sa_pt, &
                                                ct_pt_pt)
==========================================================================
   Calculates the following three, second-order derivatives of Conservative
   Temperature
    (1) CT_SA_SA, the second derivative with respect to Absolute Salinity
        at constant potential temperature (with p_ref = 0 dbar),
    (2) CT_SA_pt, the derivative with respect to potential temperature
        (the regular potential temperature which is referenced to 0 dbar)
        and Absolute Salinity, and
    (3) CT_pt_pt, the second derivative with respect to potential
        temperature (the regular potential temperature which is referenced
        to 0 dbar) at constant Absolute Salinity.
   SA    :  Absolute Salinity                                        [ g/kg ]
   pt  =  potential temperature (ITS-90)                          [ deg C ]
          (whose reference pressure is 0 dbar)
   CT_SA_SA  =  The second derivative of Conservative Temperature with
                respect to Absolute Salinity at constant potential
                temperature (the regular potential temperature which
                has reference sea pressure of 0 dbar).
                CT_SA_SA has units of:                     [ K/((g/kg)^2) ]
   CT_SA_pt  =  The derivative of Conservative Temperature with
                respect to potential temperature (the regular one with
                p_ref = 0 dbar) and Absolute Salinity.
                CT_SA_pt has units of:                        [ 1/(g/kg) ]
   CT_pt_pt  =  The second derivative of Conservative Temperature with
                respect to potential temperature (the regular one with
                p_ref = 0 dbar) at constant SA.
                CT_pt_pt has units of:                              [ 1/K ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosSea::gsw_ct_second_derivatives(double sa,double pt,double *ct_sa_sa,double *ct_sa_pt,double *ct_pt_pt)
{
	double  ct_pt_l, ct_pt_u, ct_sa_l, ct_sa_u, pt_l, pt_u, sa_l, sa_u,
			  dsa = 1e-3, dpt = 1e-2;
	if ((ct_sa_sa != NULL))
	{
		if ((sa_l = sa - dsa) < 0.0)
			sa_l = 0.0;
		sa_u = sa + dsa;
		gsw_ct_first_derivatives(sa_l,pt,&ct_sa_l,NULL);
		gsw_ct_first_derivatives(sa_u,pt,&ct_sa_u,NULL);
		*ct_sa_sa = (ct_sa_u - ct_sa_l)/(sa_u - sa_l);
	}
	if ((ct_sa_pt != NULL) || (ct_pt_pt != NULL))
	{
		pt_l = pt - dpt;
		pt_u = pt + dpt;
		if ((ct_sa_pt != NULL) && (ct_pt_pt != NULL))
		{
			gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,&ct_pt_l);
			gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,&ct_pt_u);
			*ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);
			*ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);
		}
		else if ((ct_sa_pt != NULL) && (ct_pt_pt == NULL))
		{
			gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,NULL);
			gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,NULL);
			*ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);
		}
		else if ((ct_sa_pt == NULL) && (ct_pt_pt != NULL))
		{
			gsw_ct_first_derivatives(sa,pt_l,NULL,&ct_pt_l);
			gsw_ct_first_derivatives(sa,pt_u,NULL,&ct_pt_u);
			*ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);
		}
	}
}

/**************************************************************************
==========================================================================
method: gsw_geo_strf_dyn_height (sa, ct, p, p_ref)
==========================================================================
   Calculates dynamic height anomaly as the integral of specific volume
   anomaly from the pressure p of the bottle to the reference pressure
   p_ref.
   Hence, geo_strf_dyn_height is the dynamic height anomaly with respect
   to a given reference pressure.  This is the geostrophic streammethod
   for the difference between the horizontal velocity at the pressure
   concerned, p, and the horizontal velocity at p_ref.  Dynamic height
   anomaly is the geostrophic streammethod in an isobaric surface.  The
   reference values used for the specific volume anomaly are
   SSO = 35.16504 g/kg and CT = 0 deg C.  This method calculates
   specific volume anomaly using the computationally efficient
   expression for specific volume of Roquet et al. (2015).
   This method evaluates the pressure integral of specific volume using
   SA and CT interpolated with respect to pressure using the method of
   Reiniger and Ross (1968).  It uses a weighted mean of (i) values
   obtained from linear interpolation of the two nearest data points, and
   (ii) a linear extrapolation of the pairs of data above and below.  This
   "curve fitting" method resembles the use of cubic splines.
   SA    :  Absolute Salinity                                      [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                     [ deg C ]
   p     :  sea pressure                                           [ dbar ]
            ( i.e. absolute pressure - 10.1325 dbar )
   p_ref :  reference pressure                                     [ dbar ]
            ( i.e. reference absolute pressure - 10.1325 dbar )
   geo_strf_dyn_height  :  dynamic height anomaly               [ m^2/s^2 ]
   Note. If p_ref exceeds the pressure of the deepest bottle on a
   vertical profile, the dynamic height anomaly for each bottle
   on the whole vertical profile is returned as NaN.
---------------------------------------------------------------------------
*************************************************************************/

double *TeosSea::gsw_geo_strf_dyn_height(double *sa,double *ct,double *p,double p_ref,int n_levels,double *dyn_height)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int m_levels = (n_levels <= 0) ? 1 : n_levels,
		 p_cnt, top_pad, i, nz, ibottle, ipref,
		 np_max, np, ibpr=0, *iidata=NULL;
	double dp_min, dp_max, p_min, p_max, max_dp_i,
			 *b=NULL, *b_av=NULL, *dp=NULL, *dp_i=NULL, *sa_i=NULL, *ct_i=NULL,
			  *p_i=NULL, *geo_strf_dyn_height0=NULL;
	/**
	----------------------------------------------------------------------------
	   This max_dp_i is the limit we choose for the evaluation of specific
	   volume in the pressure integration.  That is, the vertical integration
	   of specific volume with respect to pressure is perfomed with the pressure
	   increment being no more than max_dp_i (the default value being 1 dbar).
	----------------------------------------------------------------------------
	*/
	max_dp_i = 1.0;
	if ((nz = m_levels) <= 1)
		return (NULL);
	dp = new double[nz];
	dp_min = 11000.0;
	dp_max = -11000.0;
	for (i=0; i<nz-1; i++)
	{
		if ((dp[i] = p[i+1] - p[i]) < dp_min)
			dp_min = dp[i];
		if (dp[i] > dp_max)
			dp_max = dp[i];
	}
	if (dp_min <= 0.0)
	{
		/** pressure must be monotonic */
		if (dp) delete []dp;
		return (NULL);
	}
	p_min = p[0];
	p_max = p[nz-1];
	if (p_ref > p_max)
	{
		/**the reference pressure p_ref is deeper than all bottles*/
		if (dp) delete []dp;
		return (NULL);
	}
	/** Determine if there is a "bottle" at exactly p_ref */
	ipref = -1;
	for (ibottle = 0; ibottle < nz; ibottle++)
	{
		if (p[ibottle] == p_ref)
		{
			ipref = ibottle;
			break;
		}
	}
	if ((dp_max <= max_dp_i) && (p[0] == 0.0) && (ipref >= 0))
	{
		/**
		  vertical resolution is good (bottle gap is no larger than max_dp_i)
		  & the vertical profile begins at the surface (i.e. at p = 0 dbar)
		  & the profile contains a "bottle" at exactly p_ref.
		*/
		b = new double[3*nz];
		b_av = b+nz;
		geo_strf_dyn_height0 = b_av+nz;
		for (i=0; i<nz; i++)
		{
			b[i] = gsw_specvol_anom_standard(sa[i],ct[i],p[i]);
			if (i > 0)
				b_av[i-1] = 0.5*(b[i] + b[i-1]);
		}
		/**
		  "geo_strf_dyn_height0" is the dynamic height anomaly with respect
		  to p_ref = 0 (the surface).
		*/
		geo_strf_dyn_height0[0] = 0.0;
		for (i=1; i<nz; i++)
			geo_strf_dyn_height0[i] = b_av[i]*dp[i]*gtc.db2pa;
		for (i=1; i<nz; i++) /** cumulative sum */
			geo_strf_dyn_height0[i] = geo_strf_dyn_height0[i-1]
											  - geo_strf_dyn_height0[i];
		for (i=0; i<nz; i++)
			dyn_height[i] = geo_strf_dyn_height0[i]
								 - geo_strf_dyn_height0[ipref];
		if (b) delete []b;
	}
	else
	{
		/**
		  Test if there are vertical gaps between adjacent "bottles" which are
		  greater than max_dp_i, and that there is a "bottle" exactly at the
		  reference pressure.
		*/
		iidata = new int[nz+1];
		if ((dp_max <= max_dp_i) && (ipref >= 0))
		{
			/**
			  Vertical resolution is already good (no larger than max_dp_i), and
			  there is a "bottle" at exactly p_ref.
			*/
			sa_i = new double[2*(nz+1)];
			ct_i = sa_i+nz+1;
			p_i = new double[nz+1];
			if (p_min > 0.0)
			{
				/**
				  resolution is fine and there is a bottle at p_ref, but
				  there is not a bottle at p = 0. So add an extra bottle.
				*/
				for (i=0; i<nz; i++)
				{
					sa_i[i+1]       = sa[i];
					ct_i[i+1]       = ct[i];
					p_i[i+1]        = p[i];
				}
				sa_i[0] = sa[0];
				ct_i[0] = ct[0];
				p_i[0] = 0.0;
				ibpr = ipref+1;
				p_cnt = nz+1;
				for (i=0; i<p_cnt; i++)
					iidata[i] = i;
			}
			else
			{
				/**
				  resolution is fine, there is a bottle at p_ref, and
				  there is a bottle at p = 0
				*/
				memmove(sa_i, sa, nz*sizeof (double));
				memmove(ct_i, ct, nz*sizeof (double));
				memmove(p_i, p, nz*sizeof (double));
				ibpr = ipref;
				for (i=0; i<nz; i++)
					iidata[i] = i;
				p_cnt = nz;
			}
		}
		else
		{
			/**
			  interpolation is needed.
			*/
			np_max = 2*rint(p[nz-1]/max_dp_i+0.5);
			if (p_i) delete []p_i;
			p_i = new double[np_max];
			/** sa_i is allocated below, when its size is known */
			if (p_min > 0.0)
			{
				/**
				  there is not a bottle at p = 0.
				*/
				if (p_ref < p_min)
				{
					/**
					  p_ref is shallower than the minimum bottle pressure.
					*/
					p_i[0] = 0.0;
					gsw_p_sequence(p_i[0],p_ref,max_dp_i, p_i+1,&np);
					ibpr = p_cnt = np;
					p_cnt++;
					gsw_p_sequence(p_ref,p_min,max_dp_i, p_i+p_cnt,&np);
					p_cnt += np;
					top_pad = p_cnt;
				}
				else
				{
					/**
					  p_ref is deeper than the minimum bottle pressure.
					*/
					p_i[0] = 0.0;
					p_i[1] = p_min;
					top_pad = 2;
					p_cnt = 2;
				}
			}
			else
			{
				/**
				  there is a bottle at p = 0.
				*/
				p_i[0] = p_min;
				top_pad = 1;
				p_cnt = 1;
			}
			for (ibottle=0; ibottle < nz-1; ibottle++)
			{
				iidata[ibottle] = p_cnt-1;
				if (p[ibottle] == p_ref) ibpr = p_cnt-1;
				if (p[ibottle] < p_ref && p[ibottle+1] > p_ref)
				{
					/**
					  ... reference pressure is spanned by bottle pairs -
					  need to include p_ref as an interpolated pressure.
					*/
					gsw_p_sequence(p[ibottle],p_ref,max_dp_i, p_i+p_cnt,&np);
					p_cnt += np;
					ibpr = p_cnt-1;
					gsw_p_sequence(p_ref,p[ibottle+1],max_dp_i,p_i+p_cnt,&np);
					p_cnt += np;
				}
				else
				{
					/**
					  ... reference pressure is not spanned by bottle pairs.
					*/
					gsw_p_sequence(p[ibottle],p[ibottle+1],max_dp_i,
										p_i+p_cnt,&np);
					p_cnt += np;
				}
			}
			iidata[nz-1] = p_cnt-1;
			if (p[nz-1] == p_ref) ibpr = p_cnt-1;
			sa_i = new double[2*p_cnt];
			ct_i = sa_i+p_cnt;
			if (top_pad > 1)
			{
				gsw_linear_interp_sa_ct(sa,ct,p,nz,
												p_i,top_pad-1,sa_i,ct_i);
			}
			gsw_rr68_interp_sa_ct(sa,ct,p,nz,p_i+top_pad-1,p_cnt-top_pad+1,
										 sa_i+top_pad-1,ct_i+top_pad-1);
		}
		b = new double[4*p_cnt];
		b_av = b+p_cnt;
		dp_i = b_av+p_cnt;
		geo_strf_dyn_height0 = dp_i+p_cnt;
		for (i=0; i<p_cnt; i++)
		{
			b[i] = gsw_specvol_anom_standard(sa_i[i],ct_i[i],p_i[i]);
			if (i > 0)
			{
				dp_i[i-1] = p_i[i]-p_i[i-1];
				b_av[i-1] = 0.5*(b[i] + b[i-1]);
			}
		}
		/**
		  "geo_strf_dyn_height0" is the dynamic height anomaly with respect
		  to p_ref = 0 (the surface).
		*/
		geo_strf_dyn_height0[0] = 0.0;
		for (i=1; i<p_cnt; i++)
			geo_strf_dyn_height0[i] = b_av[i-1]*dp_i[i-1];
		for (i=1; i<p_cnt; i++) /** cumulative sum */
			geo_strf_dyn_height0[i] = geo_strf_dyn_height0[i-1]
											  - geo_strf_dyn_height0[i];
		for (i=0; i<nz; i++)
			dyn_height[i] = (geo_strf_dyn_height0[iidata[i]]
								  - geo_strf_dyn_height0[ibpr])*gtc.db2pa;
		if (b) delete []b;
		if (iidata) delete []iidata;
		if (sa_i) delete []sa_i;
		if (p_i)	delete []p_i;
	}
	if (dp) delete []dp;
	return (dyn_height);
}

/**************************************************************************
=======
method: gsw_p_sequence(p1, p2, max_dp_i, *pseq, *nps)
=======
*************************************************************************/

void TeosSea::gsw_p_sequence(double p1,double p2,double max_dp_i,double *pseq,int *nps)
{
	double dp, pstep;
	int n, i;
	dp = p2 - p1;
	n = ceil(dp/max_dp_i);
	pstep = dp/n;
	if (nps != NULL) *nps = n;
	/**
	  Generate the sequence ensuring that the value of p2 is exact to
	  avoid round-off issues, ie. dont do "pseq = p1+pstep*(i+1)".
	*/
	for (i=0; i<n; i++)
		pseq[i] = p2-pstep*(n-1-i);
}



/**************************************************************************
==========================================================================
method: gsw_pt_second_derivatives (sa, ct, pt_sa_sa, &
                                                pt_sa_ct, pt_ct_ct)
=========================================================================
   Calculates the following three second-order derivatives of potential
   temperature (the regular potential temperature which has a reference
   sea pressure of 0 dbar),
    (1) pt_SA_SA, the second derivative with respect to Absolute Salinity
        at constant Conservative Temperature,
    (2) pt_SA_CT, the derivative with respect to Conservative Temperature
        and Absolute Salinity, and
    (3) pt_CT_CT, the second derivative with respect to Conservative
        Temperature at constant Absolute Salinity.
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   pt_SA_SA  =  The second derivative of potential temperature (the
                regular potential temperature which has reference sea
                pressure of 0 dbar) with respect to Absolute Salinity
                at constant Conservative Temperature.
                pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
   pt_SA_CT  =  The derivative of potential temperature with respect
                to Absolute Salinity and Conservative Temperature.
                pt_SA_CT has units of:                         [ 1/(g/kg) ]
   pt_CT_CT  =  The second derivative of potential temperature (the
                regular one with p_ref = 0 dbar) with respect to
                Conservative Temperature at constant SA.
                pt_CT_CT has units of:                              [ 1/K ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosSea::gsw_pt_second_derivatives(double sa,double ct,double *pt_sa_sa,double *pt_sa_ct,double *pt_ct_ct)
{
	double ct_l, ct_u, pt_ct_l, pt_ct_u, pt_sa_l, pt_sa_u, sa_l, sa_u,
			 dct = 1e-2, dsa = 1e-3;
	if (pt_sa_sa != NULL)
	{
		if ((sa_l = sa - dsa) < 0.0)
			sa_l = 0.0;
		sa_u = sa + dsa;
		gsw_pt_first_derivatives(sa_l,ct,&pt_sa_l,NULL);
		gsw_pt_first_derivatives(sa_u,ct,&pt_sa_u,NULL);
		*pt_sa_sa = (pt_sa_u - pt_sa_l)/(sa_u - sa_l);
	}
	if (pt_sa_ct != NULL || pt_ct_ct != NULL)
	{
		ct_l = ct - dct;
		ct_u = ct + dct;
		if ((pt_sa_ct != NULL) && (pt_ct_ct != NULL))
		{
			gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,&pt_ct_l);
			gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,&pt_ct_u);
			*pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);
			*pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);
		}
		else if ((pt_sa_ct != NULL) && (pt_ct_ct == NULL))
		{
			gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,NULL);
			gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,NULL);
			*pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);
		}
		else if ((pt_sa_ct == NULL) && (pt_ct_ct != NULL))
		{
			gsw_pt_first_derivatives(sa,ct_l,NULL,&pt_ct_l);
			gsw_pt_first_derivatives(sa,ct_u,NULL,&pt_ct_u);
			*pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);
		}
	}
}

/**************************************************************************
===================================================================================
method:gsw_refine_grid_for_dh(*p,p_ref,nz,dp,*p_i,ni_max,*p_indices,*p_ref_ind_ptr)
===================================================================================
*************************************************************************/

int TeosSea::gsw_refine_grid_for_dh(double *p,double p_ref,int nz,double dp,double *p_i,int ni_max,int *p_indices,int *p_ref_ind_ptr)
{
	int i, iuniform, iorig;
	double p_next;
	/** Dont add a new point if it is within p_tol of an original. */
	double p_tol = 0.001 * dp;
	p_i[0] = p[0];
	p_indices[0] = 0;
	*p_ref_ind_ptr = -1;  /** initialize to a flag value */
	if (p_ref <= p[0] + p_tol)
	{
		*p_ref_ind_ptr = 0;
	}
	for (i=1, iuniform=1, iorig=1; i<ni_max && iorig<nz; i++)
	{
		/** Candidate insertion based on uniform grid: */
		p_next = p[0] + dp * iuniform;
		/** See if we need to insert p_ref: */
		if (*p_ref_ind_ptr == -1 && p_ref <= p_next && p_ref <= p[iorig])
		{
			p_i[i] = p_ref;
			*p_ref_ind_ptr = i;
			if (p_ref == p[iorig])
			{
				p_indices[iorig] = i;
				iorig++;
			}
			if (p_ref > p_next - p_tol)
			{
				iuniform++;
			}
			continue;
		}
		/** We did not insert p_ref, so insert either p_next or p[iorig]. */
		if (p_next < p[iorig] - p_tol)
		{
			p_i[i] = p_next;
			iuniform++;
		}
		else
		{
			p_i[i] = p[iorig];
			p_indices[iorig] = i;
			/** Skip this p_next if it is close to the point we just added. */
			if (p_next < p[iorig] + p_tol)
			{
				iuniform++;
			}
			iorig++;
		}
	}
	if (i == ni_max)
	{
		return (-1);  /** error  */
	}
	return (i);  /** number of elements in p_i */
}


/**************************************************************************
==========================================================================
method: gsw_rho_first_derivatives(sa,ct,p,*drho_dsa,*drho_dct,*drho_dp)
==========================================================================
   Calculates the three (3) partial derivatives of in situ density with
   respect to Absolute Salinity, Conservative Temperature and pressure.
   Note that the pressure derivative is done with respect to pressure in
   Pa, not dbar.  This method uses the computationally-efficient expression
   for specific volume in terms of SA, CT and p (Roquet et al., 2014).
  sa        : Absolute Salinity                               [g/kg]
  ct        : Conservative Temperature                        [deg C]
  p         : sea pressure                                    [dbar]
  drho_dsa  : partial derivatives of density                  [kg^2/(g m^3)]
              with respect to Absolute Salinity
  drho_dct  : partial derivatives of density                  [kg/(K m^3)]
              with respect to Conservative Temperature
  drho_dp   : partial derivatives of density                  [kg/(Pa m^3)]
              with respect to pressure in Pa
*************************************************************************/

void TeosSea::gsw_rho_first_derivatives(double sa,double ct,double p,double *drho_dsa,double *drho_dct,double *drho_dp)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double rho2, v_ct, v_p, v_sa, xs, ys, z, v;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v = gsw_gsvco_v(xs, ys, z);
	rho2 = pow(1.0/v, 2.0);
	if (drho_dsa != NULL)
	{
		v_sa = gsw_gsvco_b(xs, ys, z);
		*drho_dsa = -rho2*0.5*gtc.gsw_sfac*v_sa/xs;
	}
	if (drho_dct != NULL)
	{
		v_ct = gsw_gsvco_a(xs, ys, z);
		*drho_dct = -rho2*0.025*v_ct;
	}
	if (drho_dp != NULL)
	{
		v_p = gsw_gsvco_c(xs, ys, z);
		*drho_dp = 1e-4*gtc.pa2db*-rho2*v_p;
	}
	return;
}

/**************************************************************************
======================================================================================
method: gsw_rho_second_derivatives_wrt_enthalpy(sa,ct,p,*rho_sa_sa,*rho_sa_h,*rho_h_h)
======================================================================================
   Calculates three second-order derivatives of rho with respect to enthalpy.
   Note that this method uses the using the computationally-efficient
   expression for specific volume (Roquet et al., 2014).
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   rho_SA_SA : The second-order derivative of rho with respect to
               Absolute Salinity at constant h & p.     [ J/(kg (g/kg)^2) ]
   rho_SA_h  : The second-order derivative of rho with respect to
               SA and h at constant p.                   [ J/(kg K(g/kg)) ]
   rho_h_h   : The second-order derivative of rho with respect to h at
               constant SA & p[J/(kg K(g/kg))]
---------------------------------------------------------------------------
*************************************************************************/

void TeosSea::gsw_rho_second_derivatives_wrt_enthalpy(double sa,double ct,double p,
		double *rho_sa_sa,double *rho_sa_h,double *rho_h_h)
{
	double  rec_v, rec_v2, rec_v3, v_h, v_h_h, v_sa, v_sa_h, v_sa_sa,
			  *pv_sa=NULL, *pv_h=NULL, *pv_sa_sa=NULL, *pv_sa_h=NULL, *pv_h_h=NULL;

	pv_sa   = ((rho_sa_sa != NULL) || (rho_sa_h != NULL)) ? &v_sa : NULL;
	pv_h    = ((rho_sa_h != NULL) || (rho_h_h != NULL)) ?  &v_h : NULL;

	gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,pv_sa,pv_h);

	pv_sa_sa = ((rho_sa_sa != NULL)) ? &v_sa_sa : NULL;
	pv_sa_h  = ((rho_sa_h != NULL)) ? &v_sa_h : NULL;
	pv_h_h   = ((rho_h_h != NULL)) ? &v_h_h : NULL;

	gsw_specvol_second_derivatives_wrt_enthalpy(sa,ct,p,pv_sa_sa,pv_sa_h,pv_h_h);

	rec_v = 1.0/gsw_specvol(sa,ct,p);
	rec_v2 = rec_v*rec_v;
	rec_v3 = rec_v2*rec_v;

	if (rho_sa_sa != NULL)
		*rho_sa_sa = -v_sa_sa*rec_v2 + 2.0*v_sa*v_sa*rec_v3;

	if (rho_sa_h != NULL)
		*rho_sa_h = -v_sa_h*rec_v2 + 2.0*v_sa*v_h*rec_v3;

	if (rho_h_h != NULL)
		*rho_h_h = -v_h_h*rec_v2 + 2.0*v_h*v_h*rec_v3;

}







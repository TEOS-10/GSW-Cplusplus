/******************************************
  TeosBase.cpp Version 1.05
  Base class for TeosIce, TeosSea, and
  Teos_Matlab classes
  ---------------------------------------
  by Randall Kent Whited
  rkwhited@gmail.com
  ---------------------------------------
  All copyrights and all license issues
  are the same as for previous versions
*******************************************/
#include "TeosBase.h"
#include "TeosCppSupp.h"


TeosBase::TeosBase()
	: /** initialize the complex numbers */
	  t1(3.68017112855051e-2, 5.10878114959572e-2),
	  t2(3.37315741065416e-1, 3.35449415919309e-1),
	  r1(4.47050716285388e1,  6.56876847463481e1),
	  r20(-7.25974574329220e1, -7.81008427112870e1),
	  r21(-5.57107698030123e-5, 4.64578634580806e-5),
	  r22(2.34801409215913e-11, -2.85651142904972e-11)
{
	try
	{
		/** instantiate a SaarDataHandler object */
		pSaarData = new SaarDataHandler();

		/** initialize GSW_SPECVOL_COEFFICIENTS */
		gsvco.init_GSW_SPECVOL_COEFFICIENTS();

		/** initialize GSW_TEOS10_CONSTANTS */
		gtc.init_GSW_TEOS10_CONSTANTS();

		/** initialize GSW_BALTIC_DATA */
		gbalt.init_GSW_BALTIC_DATA();

		/** initialize GSW_SAAR_DATA */
		gsaar.init_GSW_SAAR_DATA();

		/** initialize GSW_SP_COEFFICIENTS */
		gspc.init_GSW_SP_COEFFICIENTS();

		/** initialize GSW_FREEZING_POLY_COEFFICIENTS */
		gfpc.init_GSW_FREEZING_POLY_COEFFICIENTS();

		/** initialize GSW_GIBBS_ICE_COEFFICIENTS */
		gic.init_GSW_GIBBS_ICE_COEFFICIENTS();
	}
	catch(bad_alloc &e)
	{
		throw badBaseAlloc;
	}
}
/***************************************
    This destructor deletes "pSaarData"
    ----------------------------------
    (a pointer to an object of the
    "SaarDataHandler" class)
****************************************/
TeosBase::~TeosBase()
{
	if (pSaarData)
	{
		delete pSaarData;
		pSaarData = NULL;
	}
}


/*************************************************************************
==========================================================================
method: gsw_add_barrier(input_data,lon,lat,long_grid,lat_grid,dlong_grid,
                        dlat_grid,output_data)
==========================================================================
  Adds a barrier through Central America (Panama) and then averages
  over the appropriate side of the barrier
  data_in      :  data                                         [unitless]
  lon          :  Longitudes of data degrees east              [0 ... +360]
  lat          :  Latitudes of data degrees north              [-90 ... +90]
  longs_grid   :  Longitudes of regular grid degrees east      [0 ... +360]
  lats_grid    :  Latitudes of regular grid degrees north      [-90 ... +90]
  dlongs_grid  :  Longitude difference of regular grid degrees [deg longitude]
  dlats_grid   :  Latitude difference of regular grid degrees  [deg latitude]
  output_data  : average of data depending on which side of the
                 Panama canal it is on                         [unitless]
**************************************************************************/

void TeosBase::gsw_add_barrier(double *input_data,double lon,double lat,double long_grid,double lat_grid,double dlong_grid,double dlat_grid,double *output_data)
{
	/** for GSW_SAAR_DATA use gsaar */
	int above_line[4], k, nmean, above_line0, kk;
	double r, lats_line, data_mean;
	k = gsw_util_indx(gsaar.longs_pan,gsaar.npan,lon);
	/**   the lon/lat point */
	r = (lon-gsaar.longs_pan[k])/(gsaar.longs_pan[k+1]-gsaar.longs_pan[k]);
	lats_line	= gsaar.lats_pan[k] + r*(gsaar.lats_pan[k+1]-gsaar.lats_pan[k]);
	above_line0	= (lats_line <= lat);
	k = gsw_util_indx(gsaar.longs_pan,gsaar.npan,long_grid);
	/**the 1 & 4 lon/lat points*/
	r = (long_grid-gsaar.longs_pan[k])/
		 (gsaar.longs_pan[k+1]-gsaar.longs_pan[k]);
	lats_line	= gsaar.lats_pan[k] + r*(gsaar.lats_pan[k+1]-gsaar.lats_pan[k]);
	above_line[0] = (lats_line <= lat_grid);
	above_line[3] = (lats_line <= lat_grid+dlat_grid);
	k = gsw_util_indx(gsaar.longs_pan,6,long_grid+dlong_grid);
	/**the 2 & 3 lon/lat points */
	r = (long_grid+dlong_grid-gsaar.longs_pan[k])/
		 (gsaar.longs_pan[k+1]-gsaar.longs_pan[k]);
	lats_line = gsaar.lats_pan[k] + r*(gsaar.lats_pan[k+1]-gsaar.lats_pan[k]);
	above_line[1] = (lats_line <= lat_grid);
	above_line[2] = (lats_line <= lat_grid+dlat_grid);
	nmean = 0;
	data_mean = 0.0;
	for (kk=0; kk<4; kk++)
	{
		if ((fabs(input_data[kk]) <= 100.0) && above_line0 == above_line[kk])
		{
			nmean	= nmean+1;
			data_mean = data_mean+input_data[kk];
		}
	}
	if (nmean == 0)
		data_mean = 0.0;	/**error return*/
	else
		data_mean = data_mean/nmean;
	for (kk=0; kk<4; kk++)
	{
		if ((fabs(input_data[kk]) >= 1.0e10) || above_line0 != above_line[kk])
		{
			output_data[kk] = data_mean;
		}
		else
			output_data[kk] = input_data[kk];
	}
}

/*************************************************************************
==========================================================================
method: gsw_add_mean(data_in,data_out)
==========================================================================
 Replaces NaNs with non-nan mean of the 4 adjacent neighbours
 data_in   : data set of the 4 adjacent neighbours[double]
 data_out : non-nan mean of the 4 adjacent neighbours     [unitless]
*************************************************************************/

void TeosBase::gsw_add_mean(double *data_in, double *data_out)
{
	int k, nmean = 0;
	double data_mean = 0.0;
	for (k=0; k<4; k++)
	{
		if (fabs(data_in[k]) <= 100.0)
		{
			nmean++;
			data_mean = data_mean+data_in[k];
		}
	}
	if (nmean == 0.0)
		data_mean = 0.0;    /**error return*/
	else
		data_mean = data_mean/nmean;
	for (k=0; k<4; k++)
	{
		if (fabs(data_in[k]) >= 100.0)
		{
			data_out[k]	= data_mean;
		}
		else
			data_out[k]	= data_in[k];
	}
}

/**************************************************************************
==========================================================================
method: gsw_adiabatic_lapse_rate_ice (t, p)
==========================================================================
   Calculates the adiabatic lapse rate of ice.
   t  :  in-situ temperature (ITS-90)                             [ deg C ]
   p  :  sea pressure                                              [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   Note:  The output is in unit of degress Celsius per Pa,
       (or equivilently K/Pa) not in units of K/dbar.
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_adiabatic_lapse_rate_ice(double t,double p)
{
	return (-gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(2,0,t,p));
}

/*************************************************************************
==========================================================================
method: gsw_alpha(sa,ct,p)
==========================================================================
  Calculates the thermal expansion coefficient of seawater with respect to
  Conservative Temperature using the computationally-efficient 48-term
  expression for density in terms of SA, CT and p (IOC et al., 2010)
 sa     : Absolute Salinity                               [g/kg]
 ct     : Conservative Temperature                        [deg C]
 p      : sea pressure                                    [dbar]
 gsw_alpha : thermal expansion coefficient of seawater (48 term equation)
**************************************************************************/

double TeosBase::gsw_alpha(double sa, double ct, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double xs, ys, z, v_ct_part, rval;
	xs	= sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys	= ct*0.025;
	z	= p*1e-4;
	v_ct_part = gsw_get_v_ct_part(xs, ys, z);
	rval = 0.025*v_ct_part/gsw_specvol(sa,ct,p);
	return rval;
}

/**************************************************************************
==========================================================================
method: gsw_beta(sa,ct,p)
==========================================================================
   Calculates the saline (i.e. haline) contraction coefficient of seawater
   at constant Conservative Temperature using the computationally-efficient
   expression for specific volume in terms of SA, CT and p
   (Roquet et al., 2014).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
          ( i.e. absolute pressure - 10.1325 dbar )
  beta   : saline contraction coefficient of seawater      [kg/g]
           at constant Conservative Temperature
*************************************************************************/

double TeosBase::gsw_beta(double sa,double ct,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double xs, ys, z, v_sa_part, rval;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v_sa_part = gsw_gsvco_b(xs, ys, z);
	rval = (-v_sa_part*0.5*gtc.gsw_sfac/(gsw_specvol(sa,ct,p)*xs));
	return rval;
}

/**************************************************************************
================================================================================
method: gsw_chem_potential_water_t_exact (sa, t, p)
================================================================================
   Calculates the chemical potential of water in seawater.
   SA  :  Absolute Salinity                                        [ g/kg ]
   t   :  in-situ temperature (ITS-90)                            [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   chem_potential_water_t_exact  :  chemical potential of water in seawater[ J/g ]
----------------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_chem_potential_water_t_exact(double sa,double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double g03g_plus_g08g, g_sa_part, x, x2, y, z, kg2g = 1e-3, rval;
	x2 = gtc.gsw_sfac*sa;
	x = sqrt(x2);
	y = t*0.025;
	z = p*1e-4;
	g03g_plus_g08g = gsw_get_g03plusg08_h(z, y, x2, x);
	g_sa_part = 8645.36753595126 +
					x*(-7296.43987145382 + x*(8103.20462414788 +
													  y*(2175.341332000392 + y*(-274.2290036817964 +
															  y*(197.4670779425016 + y*(-68.5590309679152 + 9.98788038278032*y))) - 90.6734234051316*z) +
													  x*(-5458.34205214835 - 980.14153344888*y +
															  x*(2247.60742726704 - 340.1237483177863*x + 220.542973797483*y) + 180.142097805543*z) +
													  z*(-219.1676534131548 + (-16.32775915649044 - 120.7020447884644*z)*z)) +
						z*(598.378809221703 + z*(-156.8822727844005 + (204.1334828179377 - 10.23755797323846*z)*z)) +
						y*(-1480.222530425046 + z*(-525.876123559641 + (249.57717834054571 - 88.449193048287*z)*z) +
							y*(-129.1994027934126 + z*(1149.174198007428 +
									z*(-162.5751787551336 + 76.9195462169742*z)) +
								y*(-30.0682112585625 - 1380.9597954037708*z +
									y*(2.626801985426835 + 703.695562834065*z))))) +
					y*(1187.3715515697959);
	rval = (kg2g*(g03g_plus_g08g - 0.5*x2*g_sa_part));
	return rval;
}

/**************************************************************************
==========================================================================
method: gsw_cp_ice (t,p)
==========================================================================
   Calculates the isobaric heat capacity of seawater.
   t   :  in-situ temperature (ITS-90)                            [ deg C ]
   p   :  sea pressure                                             [ dbar ]
           ( i.e. absolute pressure - 10.1325 dbar )
   gsw_cp_ice  :  heat capacity of ice                       [J kg^-1 K^-1]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_cp_ice(double t, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (-(t + gtc.gsw_t0)*gsw_gibbs_ice(2,0,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_ct_first_derivatives_wrt_t_exact (sa, t, p, &
                                       ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
==========================================================================
   Calculates the following three derivatives of Conservative Temperature.
   These derivatives are done with respect to in-situ temperature t (in the
   case of CT_T_wrt_t) or at constant in-situ tempertature (in the cases of
   CT_SA_wrt_t and CT_P_wrt_t).
    (1) CT_SA_wrt_t, the derivative of CT with respect to Absolute Salinity
        at constant t and p, and
    (2) CT_T_wrt_t, derivative of CT with respect to in-situ temperature t
        at constant SA and p.
    (3) CT_P_wrt_t, derivative of CT with respect to pressure P (in Pa) at
        constant SA and t.
   This method uses the full Gibbs function. Note that this method
   avoids the NaN that would exist in CT_SA_wrt_t at SA = 0 if it were
   evaluated in the straightforward way from the derivatives of the Gibbs
   function method.
   SA  :  Absolute Salinity                                        [ g/kg ]
   t   :  in-situ temperature (ITS-90)                            [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar)
   CT_SA_wrt_t  :  The first derivative of Conservative Temperature with
                   respect to Absolute Salinity at constant t and p. [ K/(g/kg)]  i.e. [ K kg/g ]
   CT_T_wrt_t  :  The first derivative of Conservative Temperature with
                  respect to in-situ temperature, t, at constant SA and p. [ unitless ]
   CT_P_wrt_t  :  The first derivative of Conservative Temperature with
                  respect to pressure P (in Pa) at constant SA and t. [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_ct_first_derivatives_wrt_t_exact(double sa,double t,double p,double *ct_sa_wrt_t,double *ct_t_wrt_t,double *ct_p_wrt_t)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double g_sa_mod, g_sa_t_mod, pt0, x, y, y_pt, z;
	pt0 = gsw_pt0_from_t(sa,t,p);
	if (ct_sa_wrt_t != NULL)
	{
		x = sqrt(gtc.gsw_sfac*sa);
		y = 0.025*t;
		y_pt = 0.025*pt0;
		z = gtc.rec_db2pa*p;
		/** the input pressure (p) is sea pressure in units of dbar */
		g_sa_t_mod = 1187.3715515697959 + z*(1458.233059470092 +
														 z*(-687.913805923122 + z*(249.375342232496 +
																 z*(-63.313928772146 + 14.09317606630898*z)))) +
						 x*(-1480.222530425046 + x*(2175.341332000392 +
															 x*(-980.14153344888 + 220.542973797483*x) +
															 y*(-548.4580073635929 + y*(592.4012338275047 +
																	 y*(-274.2361238716608 + 49.9394019139016*y))) -
															 90.6734234051316*z) + z*(-525.876123559641 +
																	 (249.57717834054571 - 88.449193048287*z)*z) +
							 y*(-258.3988055868252 + z*(2298.348396014856 +
									 z*(-325.1503575102672 + 153.8390924339484*z)) +
								 y*(-90.2046337756875 - 4142.8793862113125*z +
									 y*(10.50720794170734 + 2814.78225133626*z)))) +
						 y*(3520.125411988816 + y*(-1351.605895580406 +
															y*(731.4083582010072 + y*(-216.60324087531103 +
																	25.56203650166196*y) + z*(-2381.829935897496 +
																			(597.809129110048 - 291.8983352012704*z)*z)) +
															z*(4165.4688847996085 + z*(-1229.337851789418 +
																	(681.370187043564 - 66.7696405958478*z)*z))) +
							 z*(-3443.057215135908 + z*(1349.638121077468 +
									 z*(-713.258224830552 +
										 (176.8161433232 - 31.68006188846728*z)*z))));
		g_sa_t_mod = 0.5*gtc.gsw_sfac*0.025*g_sa_t_mod;
		g_sa_mod = 8645.36753595126 +
					  x*(-7296.43987145382 + x*(8103.20462414788 +
														 y_pt*(2175.341332000392 + y_pt*(-274.2290036817964 +
																 y_pt*(197.4670779425016 + y_pt*(-68.5590309679152 +
																		 9.98788038278032*y_pt)))) +
														 x*(-5458.34205214835 - 980.14153344888*y_pt +
															 x*(2247.60742726704 - 340.1237483177863*x +
																	 220.542973797483*y_pt))) +
						  y_pt*(-1480.222530425046 +
								  y_pt*(-129.1994027934126 +
										  y_pt*(-30.0682112585625 + y_pt*(2.626801985426835 ))))) +
					  y_pt*(1187.3715515697959 +
							  y_pt*(1760.062705994408 + y_pt*(-450.535298526802 +
									  y_pt*(182.8520895502518 + y_pt*(-43.3206481750622 +
											  4.26033941694366*y_pt)))));
		g_sa_mod = 0.5*gtc.gsw_sfac*g_sa_mod;
		*ct_sa_wrt_t = (g_sa_mod - (gtc.gsw_t0+pt0)*g_sa_t_mod)/gtc.gsw_cp0;
	}
	if (ct_t_wrt_t != NULL)
	{
		*ct_t_wrt_t = -(gtc.gsw_t0+pt0)*gsw_gibbs(0,2,0,sa,t,p)/gtc.gsw_cp0;
	}
	if (ct_p_wrt_t != NULL)
	{
		*ct_p_wrt_t = -(gtc.gsw_t0+pt0)*gsw_gibbs(0,1,1,sa,t,p)/gtc.gsw_cp0;
	}
}

/**************************************************************************
==========================================================================
method: gsw_ct_freezing (sa, p, saturation_fraction)
==========================================================================
   Calculates the Conservative Temperature at which seawater freezes.  The
   Conservative Temperature freezing point is calculated from the exact
   in-situ freezing temperature which is found by a modified Newton-Raphson
   iteration (McDougall and Wotherspoon, 2013) of the equality of the
   chemical potentials of water in seawater and in ice.
   An alternative GSW method, gsw_CT_freezing_poly, it is based on a
   computationally-efficient polynomial, and is accurate to within -5e-4 K
   and 6e-4 K, when compared with this method.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
                         seawater
   CT_freezing : Conservative Temperature at freezing of seawater [ deg C ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_ct_freezing(double sa, double p, double saturation_fraction)
{
	double  t_freezing = gsw_t_freezing(sa,p,saturation_fraction);
	return (gsw_ct_from_t(sa,t_freezing,p));
}

/**************************************************************************
==========================================================================
method: gsw_ct_freezing_first_derivatives (sa, p, &
                          saturation_fraction, ctfreezing_sa, ctfreezing_p)
==========================================================================
   Calculates the first derivatives of the Conservative Temperature at
   which seawater freezes,with respect to Absolute Salinity SA and
   pressure P (in Pa).
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction = the saturation fraction of dissolved air in
                         seawater
   CTfreezing_SA : the derivative of the Conservative Temperature at
                   freezing (ITS-90) with respect to Absolute Salinity at
                   fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
   CTfreezing_P : the derivative of the Conservative Temperature at
                   freezing (ITS-90) with respect to pressure (in Pa) at
                   fixed Absolute Salinity                         [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_ct_freezing_first_derivatives(double sa,double p,double saturation_fraction,double *ctfreezing_sa,double *ctfreezing_p)
{
	double tf_sa, tf_p, ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t, tf;
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	if (ctfreezing_sa != NULL && ctfreezing_p != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa,&tf_p);
		gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
														 &ct_sa_wrt_t,&ct_t_wrt_t,&ct_p_wrt_t);
		*ctfreezing_sa = ct_sa_wrt_t + ct_t_wrt_t*tf_sa;
		*ctfreezing_p  = ct_p_wrt_t  + ct_t_wrt_t*tf_p;
	}
	else if (ctfreezing_sa != NULL && ctfreezing_p == NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa, NULL);
		gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
														 &ct_sa_wrt_t,&ct_t_wrt_t,NULL);
		*ctfreezing_sa = ct_sa_wrt_t + ct_t_wrt_t*tf_sa;
	}
	else if (ctfreezing_sa == NULL && ctfreezing_p != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													NULL, &tf_p);
		gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
														 NULL,&ct_t_wrt_t,&ct_p_wrt_t);
		*ctfreezing_p  = ct_p_wrt_t  + ct_t_wrt_t*tf_p;
	}
}

/**************************************************************************
==========================================================================
method: gsw_ct_freezing_first_derivatives_poly (sa, p, &
                          saturation_fraction, ctfreezing_sa, ctfreezing_p)
==========================================================================
   Calculates the first derivatives of the Conservative Temperature at
   which seawater freezes, with respect to Absolute Salinity SA and
   pressure P (in Pa) of the comptationally efficient polynomial fit of the
   freezing temperature (McDougall et al., 2014).
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
                         seawater
   CTfreezing_SA : the derivative of the Conservative Temperature at
                   freezing (ITS-90) with respect to Absolute Salinity at
                   fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
   CTfreezing_P  : the derivative of the Conservative Temperature at
                   freezing (ITS-90) with respect to pressure (in Pa) at
                   fixed Absolute Salinity                         [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_ct_freezing_first_derivatives_poly(double sa,double p,double saturation_fraction,double *ctfreezing_sa,double *ctfreezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x,
			 d = -gfpc.a - gfpc.a*gfpc.b - 2.4*gfpc.b/gtc.gsw_sso,
			 e = 2.0*gfpc.a*gfpc.b/gtc.gsw_sso;
	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;
	if (ctfreezing_sa != NULL)
	{
		*ctfreezing_sa =
			(gfpc.c1 + x*(1.5*gfpc.c2 + x*(2.0*gfpc.c3 + x*(2.5*gfpc.c4 + x*(3.0*gfpc.c5
													 + 3.5*gfpc.c6*x)))) + p_r*(gfpc.c10 + x*(1.5*gfpc.c11 + x*(2.0*gfpc.c13
															 + x*(2.5*gfpc.c16 + x*(3.0*gfpc.c19 + 3.5*gfpc.c22*x))))
															 + p_r*(gfpc.c12 + x*(1.5*gfpc.c14 + x*(2.0*gfpc.c17 + 2.5*gfpc.c20*x))
																	 + p_r*(gfpc.c15 + x*(1.5*gfpc.c18 + 2.0*gfpc.c21*x)))))*1e-2
			- saturation_fraction*1e-3*(d - sa*e);
	}
	if (ctfreezing_p != NULL)
	{
		*ctfreezing_p =
			(gfpc.c7 + sa_r*(gfpc.c10 + x*(gfpc.c11 + x*(gfpc.c13 + x*(gfpc.c16
													 + x*(gfpc.c19 + gfpc.c22*x)))))
			 + p_r*(2.0*gfpc.c8 + sa_r*(2.0*gfpc.c12 + x*(2.0*gfpc.c14 + x*(2.0*gfpc.c17
												 + 2.0*gfpc.c20*x))) + p_r*(3.0*gfpc.c9 + sa_r*(3.0*gfpc.c15 + x*(3.0*gfpc.c18
														 + 3.0*gfpc.c21*x)))))*1e-8;
	}
}

/**************************************************************************
==========================================================================
method: gsw_ct_freezing_poly (sa, p, saturation_fraction)
==========================================================================
   Calculates the Conservative Temperature at which seawater freezes.
   The error of this fit ranges between -5e-4 K and 6e-4 K when compared
   with the Conservative Temperature calculated from the exact in-situ
   freezing temperature which is found by a Newton-Raphson iteration of the
   equality of the chemical potentials of water in seawater and in ice.
   Note that the Conservative temperature freezing temperature can be found
   by this exact function using the method gsw_CT_freezing.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in seawater
   CT_freezing : Conservative Temperature at freezing of seawater [ deg C ]
                 That is, the freezing temperature expressed in
                 terms of Conservative Temperature (ITS-90).
----------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_ct_freezing_poly(double sa,double p,double saturation_fraction)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x, return_value;
	sa_r = sa*1.0e-2;
	x = sqrt(sa_r);
	p_r = p*1.0e-4;
	return_value = gfpc.c0
						+ sa_r*(gfpc.c1 + x*(gfpc.c2 + x*(gfpc.c3 + x*(gfpc.c4 + x*(gfpc.c5 + gfpc.c6*x)))))
						+ p_r*(gfpc.c7 + p_r*(gfpc.c8 + gfpc.c9*p_r)) + sa_r*p_r*(gfpc.c10 + p_r*(gfpc.c12
								+ p_r*(gfpc.c15 + gfpc.c21*sa_r)) + sa_r*(gfpc.c13 + gfpc.c17*p_r + gfpc.c19*sa_r)
								+ x*(gfpc.c11 + p_r*(gfpc.c14 + gfpc.c18*p_r) + sa_r*(gfpc.c16 + gfpc.c20*p_r
										+ gfpc.c22*sa_r)));
	/** Adjust for the effects of dissolved air */
	return_value = return_value - saturation_fraction*
						(1e-3)*(2.4 - gfpc.a*sa)*(1.0 + gfpc.b*(1.0 - sa/gtc.gsw_sso));
	return (return_value);
}

/**************************************************************************
==========================================================================
method: gsw_ct_from_enthalpy (sa, h, p)
==========================================================================
   Calculates the Conservative Temperature of seawater, given the Absolute
   Salinity, specific enthalpy, h, and pressure p.
   SA  :  Absolute Salinity                                        [ g/kg ]
   h   :  specific enthalpy                                        [ J/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325d0 dbar )
   CT  :  Conservative Temperature ( ITS-90)                      [ deg C ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_ct_from_enthalpy(double sa,double h,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double ct, ct_freezing, ct_mean, ct_old, f, h_freezing,
			 h_ct, h_40, ct_40 = 40.0;
	ct_freezing = gsw_ct_freezing(sa,p,0.0);
	h_freezing = gsw_enthalpy(sa,ct_freezing,p);
	if (h < (h_freezing - gtc.gsw_cp0))
	{
		/**
		  The input, seawater enthalpy h, is less than the enthalpy at the
		  freezing temperature, i.e. the water is frozen.
		*/
		return (cppGSW_INVALID_VALUE);
	}
	h_40 = gsw_enthalpy(sa,ct_40,p);
	if (h > h_40)
	{
		/**
		  The input seawater enthalpy is greater than the enthalpy
		  when CT is 40C
		*/
		return (cppGSW_INVALID_VALUE);
	}
	/** first guess of ct */
	ct = ct_freezing + (ct_40 - ct_freezing)*
		  (h - h_freezing)/(h_40 - h_freezing);
	gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);
	/**
	   ------------------------------------------------------
	   Begin the modified Newton-Raphson iterative procedure
	   ------------------------------------------------------
	*/
	ct_old = ct;
	f = gsw_enthalpy(sa,ct_old,p) - h;
	ct = ct_old - f/h_ct;
	ct_mean = 0.5*(ct + ct_old);
	gsw_enthalpy_first_derivatives(sa,ct_mean,p,NULL,&h_ct);
	ct = ct_old - f/h_ct;
	ct_old = ct;
	f = gsw_enthalpy(sa,ct_old,p) - h;
	ct = ct_old - f/h_ct;
	/**
	   After 1.5d0 iterations of this modified Newton-Raphson iteration,
	   the error in CT is no larger than 4x10^-13 degrees C, which
	   is machine precision for this calculation.
	*/
	return (ct);
}

/**************************************************************************
==========================================================================
method: gsw_ct_from_enthalpy_exact (sa, h, p)
==========================================================================
   Calculates the Conservative Temperature of seawater, given the Absolute
   Salinity, specific enthalpy, h, and pressure p.
   SA  :  Absolute Salinity                                        [ g/kg ]
   h   :  specific enthalpy                                        [ J/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325d0 dbar )
   CT  :  Conservative Temperature ( ITS-90)                      [ deg C ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_ct_from_enthalpy_exact(double sa, double h, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double ct, ct_freezing, ct_mean, ct_old, f, h_freezing,
			 h_ct, h_40, ct_40 = 40.0;
	ct_freezing = gsw_t_freezing(sa,p,0.0);
	h_freezing = gsw_enthalpy_ct_exact(sa,ct_freezing,p);
	if (h < (h_freezing - gtc.gsw_cp0))
	{
		/**
		   The input, seawater enthalpy h, is less than the enthalpy at the
		   freezing temperature, i.e. the water is frozen.
		*/
		return (cppGSW_INVALID_VALUE);
	}
	h_40 = gsw_enthalpy_ct_exact(sa,ct_40,p);
	if (h > h_40)
	{
		/**
		   The input seawater enthalpy is greater than the enthalpy
		   when CT is 40C
		*/
		return (cppGSW_INVALID_VALUE);
	}
	/** First guess of ct */
	ct = ct_freezing + (ct_40 - ct_freezing)*
		  (h - h_freezing)/(h_40 - h_freezing);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,NULL,&h_ct);
	/**
	   ------------------------------------------------------
	   Begin the modified Newton-Raphson iterative procedure
	   ------------------------------------------------------
	*/
	ct_old = ct;
	f = gsw_enthalpy_ct_exact(sa,ct_old,p) - h;
	ct = ct_old - f/h_ct;
	ct_mean = 0.5*(ct + ct_old);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ct_mean,p,NULL,&h_ct);
	ct = ct_old - f/h_ct;
	/**
	   After 1 iteration of this modified Newton-Raphson iteration,
	   the error in CT is no larger than 5x10^-14 degrees C, which
	   is machine precision for this calculation.
	*/
	return (ct);
}

/*************************************************************************
==========================================================================
method: gsw_ct_from_pt(sa,pt)
==========================================================================
 Calculates Conservative Temperature from potential temperature of seawater
 sa      : Absolute Salinity                              [g/kg]
 pt      : potential temperature with                     [deg C]
           reference pressure of 0 dbar
 gsw_ct_from_pt : Conservative Temperature                [deg C]
*************************************************************************/

double TeosBase::gsw_ct_from_pt(double sa, double pt)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double	x2, x, y, pot_enthalpy, rval;
	x2 = gtc.gsw_sfac*sa;
	x = sqrt(x2);
	y = pt*0.025e0;	/*  normalize for F03 and F08 */
	pot_enthalpy = gsw_get_pot_ent(x2, x, y);
	rval = pot_enthalpy/gtc.gsw_cp0;
	return rval;
}

/**************************************************************************
==========================================================================
method: gsw_ct_from_rho (rho, sa, p, ct, ct_multiple)
=========================================================================
   Calculates the Conservative Temperature of a seawater sample, for given
   values of its density, Absolute Salinity and sea pressure (in dbar).
   rho  :  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
    Note. This input has not had 1000 kg/m^3 subtracted from it.
      That is, it is  density , not  density anomaly .
   SA   :  Absolute Salinity                                       [ g/kg ]
   p    :  sea pressure                                            [ dbar ]
           ( i.e. absolute pressure - 10.1325 dbar )
   CT  :  Conservative Temperature  (ITS-90)                      [ deg C ]
   CT_multiple  :  Conservative Temperature  (ITS-90)             [ deg C ]
     Note that at low salinities, in brackish water, there are two possible
       Conservative Temperatures for a single density.  This programme will
       output both valid solutions.  To see this second solution the user
       must call the programme with two outputs (i.e. [CT,CT_multiple]), if
       there is only one possible solution and the programme has been
       called with two outputs the second variable will be set to NaN.
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_ct_from_rho(double rho,double sa,double p,double *ct,double *ct_multiple)
{
	int number_of_iterations;
	double a, alpha_freezing, alpha_mean, b, c, ct_a, ct_b, ct_diff,
			 ct_freezing, ct_max_rho, ct_mean, ct_old,
			 delta_ct, delta_v, factor, factorqa, factorqb,
			 rho_40, rho_extreme, rho_freezing, rho_max, rho_mean,
			 rho_old, sqrt_disc, top, v_ct, v_lab;
	/**
	   alpha_limit is the positive value of the thermal expansion coefficient
	   which is used at the freezing temperature to distinguish between
	   salty and fresh water.
	*/
	double alpha_limit = 1e-5;
	/**
	   rec_half_rho_TT is a constant representing the reciprocal of half the
	   second derivative of density with respect to temperature near the
	   temperature of maximum density.
	*/
	double rec_half_rho_tt = -110.0;
	rho_40 = gsw_rho(sa,40.0,p);
	if (rho < rho_40)
	{
		*ct = cppGSW_INVALID_VALUE;
		if (ct_multiple != NULL) *ct_multiple = *ct;
		return;
	}
	ct_max_rho = gsw_ct_maxdensity(sa,p);
	rho_max = gsw_rho(sa,ct_max_rho,p);
	rho_extreme = rho_max;
	/** Assumes that the seawater is always unsaturated with air */
	ct_freezing = gsw_ct_freezing_poly(sa,p,0.0);
	gsw_rho_alpha_beta(sa,ct_freezing,p,&rho_freezing,&alpha_freezing,NULL);
	/** reset the extreme values */
	if (ct_freezing > ct_max_rho)
	{
		rho_extreme = rho_freezing;
	}
	if (rho > rho_extreme)
	{
		*ct = cppGSW_INVALID_VALUE;
		if (ct_multiple != NULL) *ct_multiple = *ct;
		return;
	}
	if (alpha_freezing > alpha_limit)
	{
		ct_diff = 40.0 - ct_freezing;
		top = rho_40 - rho_freezing + rho_freezing*alpha_freezing*ct_diff;
		a = top/(ct_diff*ct_diff);
		b = -rho_freezing*alpha_freezing;
		c = rho_freezing - rho;
		sqrt_disc = sqrt(b*b - 4*a*c);
		*ct = ct_freezing + 0.5*(-b - sqrt_disc)/a;
	}
	else
	{
		ct_diff = 40.0 - ct_max_rho;
		factor = (rho_max - rho)/(rho_max - rho_40);
		delta_ct = ct_diff*sqrt(factor);
		if (delta_ct > 5.0)
		{
			*ct = ct_max_rho + delta_ct;
		}
		else
		{
			/** Set the initial value of the quadratic solution roots. */
			ct_a = ct_max_rho + sqrt(rec_half_rho_tt*(rho - rho_max));
			for (number_of_iterations = 1; number_of_iterations <= 7;
					number_of_iterations++)
			{
				ct_old = ct_a;
				rho_old = gsw_rho(sa,ct_old,p);
				factorqa = (rho_max - rho)/(rho_max - rho_old);
				ct_a = ct_max_rho + (ct_old - ct_max_rho)*sqrt(factorqa);
			}
			if ((ct_freezing - ct_a) < 0.0)
			{
				*ct = cppGSW_INVALID_VALUE;
				if (ct_multiple != NULL) *ct_multiple = *ct;
				return;
			}
			*ct = ct_a;
			if (ct_multiple == NULL)
			{
				return;
			}
			/** Set the initial value of the quadratic solution roots. */
			ct_b = ct_max_rho - sqrt(rec_half_rho_tt*(rho - rho_max));
			for (number_of_iterations = 1; number_of_iterations <= 7;
					number_of_iterations++)
			{
				ct_old = ct_b;
				rho_old = gsw_rho(sa,ct_old,p);
				factorqb = (rho_max - rho)/(rho_max - rho_old);
				ct_b = ct_max_rho + (ct_old - ct_max_rho)*sqrt(factorqb);
			}
			/**
			   After seven iterations of this quadratic iterative procedure,
			   the error in rho is no larger than 4.6x10^-13 kg/m^3.
			*/
			if ((ct_freezing - ct_b) < 0.0)
			{
				*ct = cppGSW_INVALID_VALUE;
				*ct_multiple = *ct;
				return;
			}
			*ct_multiple = ct_b;
			return;
		} /** if delta ct > 5 else */
	} /** if alpha_freezing else */
	/** Begin the modified Newton-Raphson iterative function */
	v_lab = 1.0/rho;
	gsw_rho_alpha_beta(sa,*ct,p,&rho_mean,&alpha_mean,NULL);
	v_ct = alpha_mean/rho_mean;
	for (number_of_iterations = 1; number_of_iterations <= 3;
			number_of_iterations++)
	{
		ct_old = *ct;
		delta_v = gsw_specvol(sa,ct_old,p) - v_lab;
		*ct = ct_old - delta_v/v_ct;
		ct_mean = 0.5*(*ct + ct_old);
		gsw_rho_alpha_beta(sa,ct_mean,p,&rho_mean,&alpha_mean,NULL);
		v_ct = alpha_mean/rho_mean;
		*ct = ct_old - delta_v/v_ct ;
	}
	/**
	   After three iterations of this modified Newton-Raphson iteration,
	   the error in rho is no larger than 1.6x10^-12 kg/m^3.
	*/
	if (ct_multiple != NULL)
	{
		*ct_multiple = cppGSW_INVALID_VALUE;
	}
	return;
}

/*************************************************************************
==========================================================================
method: gsw_ct_from_t(sa,t,p)
==========================================================================
 Calculates Conservative Temperature from in-situ temperature
 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]
 gsw_ct_from_t : Conservative Temperature                 [deg C]
*************************************************************************/

double TeosBase::gsw_ct_from_t(double sa, double t, double p)
{
	double	pt0;
	pt0	= gsw_pt0_from_t(sa,t,p);
	return (gsw_ct_from_pt(sa,pt0));
}

/**************************************************************************
==========================================================================
method: gsw_ct_maxdensity (sa, p)
==========================================================================
   Calculates the Conservative Temperature of maximum density of seawater.
   This method returns the Conservative temperature at which the density
   of seawater is a maximum, at given Absolute Salinity, SA, and sea
   pressure, p (in dbar).
   SA :  Absolute Salinity                                         [ g/kg ]
   p  :  sea pressure                                              [ dbar ]
         ( i.e. absolute pressure - 10.1325 dbar )
   CT_maxdensity :  Conservative Temperature at which            [ deg C ]
                     the density of seawater is a maximum for
                     given Absolute Salinity and pressure.
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_ct_maxdensity(double sa, double p)
{
	int number_of_iterations;
	double alpha, ct, ct_mean, ct_old, dalpha_dct,
			 dct = 0.001;
	ct = 3.978 - 0.22072*sa;         /**the initial guess of ct.*/
	dalpha_dct = 1.1e-5;             /**the initial guess for dalpha_dct.*/
	for (number_of_iterations = 1; number_of_iterations <= 3;
			number_of_iterations++)
	{
		ct_old = ct;
		alpha = gsw_alpha(sa,ct_old,p);
		ct = ct_old - alpha/dalpha_dct;
		ct_mean = 0.5*(ct + ct_old);
		dalpha_dct = (gsw_alpha(sa,ct_mean+dct,p)
						  - gsw_alpha(sa,ct_mean-dct,p))/(dct + dct);
		ct = ct_old - alpha/dalpha_dct;
	}
	/**
	  After three iterations of this modified Newton-Raphson (McDougall and
	  Wotherspoon, 2012) iteration, the error in CT_maxdensity is typically
	  no larger than 1x10^-15 degress C.
	*/
	return (ct);
}

/**************************************************************************
==========================================================================
method: gsw_deltasa_atlas(p,lon,lat)
==========================================================================
 Calculates the Absolute Salinity Anomaly atlas value, delta_SA_atlas.
 p      : sea pressure                                    [dbar]
 lon    : longiture                                       [deg E]
 lat    : latitude                                        [deg N]
 deltasa_atlas : Absolute Salinity Anomaly atlas value    [g/kg]
*************************************************************************/

double TeosBase::gsw_deltasa_atlas(double p,double lon,double lat)
{
	/** for GSW_SAAR_DATA use gsaar */
	int nx=gsw_nx, ny=gsw_ny, nz=gsw_nz;
	int indx0, indy0, indz0, k, ndepth_index;
	double dsar[4], dsar_old[4], dlong, dlat;
	double return_value, sa_upper, sa_lower;
	double r1, s1, t1, ndepth_max;
	return_value = cppGSW_INVALID_VALUE;
	if (isnan(lat) || isnan(lon) || isnan(p))
	{
		return (return_value);
	}
	if (lat < -86.0  ||  lat  >  90.0)
	{
		return (return_value);
	}
	if (lon < 0.0)
	{
		lon += 360.0;
	}
	dlong	= pSaarData->get_longs_ref(1)-pSaarData->get_longs_ref(0);
	dlat	= pSaarData->get_lats_ref(1)-pSaarData->get_lats_ref(0);
	/** these were set apart to resolve a C++ compiler problem calculating indx0 */
	double longsrefzero = pSaarData->get_longs_ref(0);
	double nxminusone = (double)nx-1.0;
	double longsrefnxminusone = pSaarData->get_longs_ref((unsigned)nxminusone);
	double lonminuslongsrefzero = lon-longsrefzero;
	indx0	= floor(nxminusone*(lonminuslongsrefzero/(longsrefnxminusone-longsrefzero)));
	if (indx0 == nx-1)
	{
		indx0	= nx-2;
	}
	/** these were set apart to resolve a C++ compiler problem calculating indy0 */
	double latsrefzero = pSaarData->get_lats_ref(0);
	double nyminusone = (double)ny-1.0;
	double latsrefnyminusone = pSaarData->get_lats_ref((unsigned)nyminusone);
	double latminuslatsrefzero = lat-latsrefzero;
	indy0	= floor(nyminusone*(latminuslatsrefzero/(latsrefnyminusone-latsrefzero)));
	if (indy0 == ny-1)
	{
		indy0	= ny-2;
	}
	ndepth_max	= -1;
	for (k=0; k<4; k++)
	{
		ndepth_index	= indy0+gsaar.delj[k]+(indx0+gsaar.deli[k])*ny;
		if (pSaarData->get_ndepth_ref(ndepth_index) > 0.0)
		{
			ndepth_max = gsw_max_d(ndepth_max, pSaarData->get_ndepth_ref(ndepth_index));
			if (ndepth_max < pSaarData->get_ndepth_ref(ndepth_index))
			{
				ndepth_max = pSaarData->get_ndepth_ref(ndepth_index);
			}
		}
	} /** for */
	if (ndepth_max == -1.0)
	{
		return (0.0);
	}
	if (p > pSaarData->get_p_ref((int)(ndepth_max)-1))
	{
		p	= pSaarData->get_p_ref((int)(ndepth_max)-1);
	}
	/*************************************************
	   gsw_util_indxPref is relocated to TeosBase.cpp
	   for a more clear continuity with gsw_util_indx
	*************************************************/
	indz0 = gsw_util_indxPref(nz,p);
	r1	= (lon-pSaarData->get_longs_ref(indx0))/
		  (pSaarData->get_longs_ref(indx0+1)-pSaarData->get_longs_ref(indx0));
	s1	= (lat-pSaarData->get_lats_ref(indy0))/
		  (pSaarData->get_lats_ref(indy0+1)-pSaarData->get_lats_ref(indy0));
	t1	= (p-pSaarData->get_p_ref(indz0))/
		  (pSaarData->get_p_ref(indz0+1)-pSaarData->get_p_ref(indz0));
	for (k=0; k < 4; k++)
	{
		dsar[k] = pSaarData->get_delta_sa_ref(indz0+nz*(indy0+gsaar.delj[k]+
														  (indx0+gsaar.deli[k])*ny));
	}
	if (gsaar.longs_pan[0] <= lon && lon <= gsaar.longs_pan[gsaar.npan-1]-0.001 &&
			gsaar.lats_pan[gsaar.npan-1] <= lat && lat <= gsaar.lats_pan[0])
	{
		memmove(dsar_old,dsar,4*sizeof (double));
		gsw_add_barrier(dsar_old,lon,lat,pSaarData->get_longs_ref(indx0),
							 pSaarData->get_lats_ref(indy0),dlong,dlat,dsar);
	}
	else if (fabs(gsw_sum(dsar, sizeof (dsar)/sizeof (double))) >= cppGSW_ERROR_LIMIT)
	{
		memmove(dsar_old,dsar,4*sizeof (double));
		gsw_add_mean(dsar_old,dsar);
	}
	sa_upper	= (1.0-s1)*(dsar[0] + r1*(dsar[1]-dsar[0])) +
				  s1*(dsar[3] + r1*(dsar[2]-dsar[3]));
	for (k=0; k<4; k++)
	{
		dsar[k] = pSaarData->get_delta_sa_ref(indz0+1+nz*(indy0+gsaar.delj[k]+
														  (indx0+gsaar.deli[k])*ny));
	}
	if (gsaar.longs_pan[0] <= lon && lon <= gsaar.longs_pan[gsaar.npan-1] &&
			gsaar.lats_pan[gsaar.npan-1] <= lat && lat <= gsaar.lats_pan[0])
	{
		memmove(dsar_old,dsar,4*sizeof (double));
		gsw_add_barrier(dsar_old,lon,lat,pSaarData->get_longs_ref(indx0),
							 pSaarData->get_lats_ref(indy0),dlong,dlat,dsar);
	}
	else if (fabs(gsw_sum(dsar, sizeof (dsar)/sizeof (double))) >= cppGSW_ERROR_LIMIT)
	{
		memmove(dsar_old,dsar,4*sizeof (double));
		gsw_add_mean(dsar_old,dsar);
	}
	sa_lower	= (1.0-s1)*(dsar[0] + r1*(dsar[1]-dsar[0])) +
				  s1*(dsar[3] + r1*(dsar[2]-dsar[3]));
	if (fabs(sa_lower) >= cppGSW_ERROR_LIMIT)
	{
		sa_lower	= sa_upper;
	}
	return_value = sa_upper + t1*(sa_lower-sa_upper);
	if (fabs(return_value) >= cppGSW_ERROR_LIMIT)
	{
		return (cppGSW_INVALID_VALUE);
	}
	return (return_value);
}

/**************************************************************************
==========================================================================
method: gsw_dilution_coefficient_t_exact (sa, t, p)
==========================================================================
   Calculates the dilution coefficient of seawater.  The dilution
   coefficient of seawater is defined as the Absolute Salinity times the
   second derivative of the Gibbs function with respect to Absolute
   Salinity, that is, SA.*g_SA_SA.
   SA :  Absolute Salinity                                         [ g/kg ]
   t  =  in-situ temperature (ITS-90)                             [ deg C ]
   p  =  sea pressure                                              [ dbar ]
         ( i.e. absolute pressure - 10.1325 dbar )
   dilution_coefficient_t_exact  =  dilution coefficient   [ (J/kg)(kg/g) ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_dilution_coefficient_t_exact(double sa,double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  g08, x, x2, y, z, rval;
	x2 = gtc.gsw_sfac*sa;
	x = sqrt(x2);
	y = t*0.025;
	z = p*1e-4;
	/**note. the input pressure (p) is sea pressure in units of dbar.*/
	g08 = gsw_g08_d(x2, z, y, x);
	rval = 0.25*gtc.gsw_sfac*g08;
	return rval;
	/**
	   Note that this method avoids the singularity that occurs at SA = 0 if
	   the straightforward expression for the dilution coefficient of seawater,
	   SA*g_SA_SA is simply evaluated as SA.*gsw_gibbs(2,0,0,SA,t,p).
	*/
}

/**************************************************************************
==========================================================================
method: gsw_dynamic_enthalpy(sa,ct,p)
==========================================================================
  Calculates dynamic enthalpy of seawater using the computationally-
  efficient expression for specific volume in terms of SA, CT and p
  (Roquet et al., 2014).  Dynamic enthalpy is defined as enthalpy minus
  potential enthalpy (Young, 2010).
 sa     : Absolute Salinity                               [g/kg]
 ct     : Conservative Temperature (ITS-90)               [deg C]
 p      : sea pressure                                    [dbar]
         ( i.e. absolute pressure - 10.1325 dbar )
 dynamic_enthalpy  :  dynamic enthalpy                    [J/kg]
*************************************************************************/

double TeosBase::gsw_dynamic_enthalpy(double sa, double ct, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double dynamic_enthalpy_part, xs, ys, z, rval;
	xs	= sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys	= ct*0.025;
	z	= p*1e-4;
	dynamic_enthalpy_part = gsw_get_dyn_enth_prt(xs, ys, z);
	rval = dynamic_enthalpy_part*gtc.db2pa*1e4;
	return rval;
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy(sa,ct,p)
==========================================================================
   Calculates specific enthalpy of seawater using the computationally-
   efficient expression for specific volume in terms of SA, CT and p
   (Roquet et al., 2014).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
          ( i.e. absolute pressure - 10.1325 dbar )
  enthalpy  :  specific enthalpy of seawater               [J/kg]
*************************************************************************/

double TeosBase::gsw_enthalpy(double sa, double ct, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (gtc.gsw_cp0*ct + gsw_dynamic_enthalpy(sa,ct,p));
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy_ct_exact (sa, ct, p)
==========================================================================
  Calculates specific enthalpy of seawater from Absolute Salinity and
  Conservative Temperature and pressure.
  Note that this method uses the full Gibbs function.
  SA  :  Absolute Salinity                                        [ g/kg ]
  CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
  p   :  sea pressure                                             [ dbar ]
         ( i.e. absolute pressure - 10.1325 dbar )
  enthalpy_CT_exact :  specific enthalpy                         [ J/kg ]
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_enthalpy_ct_exact(double sa, double ct, double p)
{
	double t = gsw_t_from_ct(sa,ct,p);
	return (gsw_enthalpy_t_exact(sa,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy_first_derivatives (sa, ct, p, h_sa, h_ct)
==========================================================================
   Calculates the following two derivatives of specific enthalpy (h) of
   seawater using the computationally-efficient expression for
   specific volume in terms of SA, CT and p (Roquet et al., 2014).
    (1) h_SA, the derivative with respect to Absolute Salinity at
        constant CT and p, and
    (2) h_CT, derivative with respect to CT at constant SA and p.
   Note that h_P is specific volume (1/rho) it can be caclulated by calling
   gsw_specvol(SA,CT,p).
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   h_SA  :  The first derivative of specific enthalpy with respect to
            Absolute Salinity at constant CT and p.[ J/(kg (g/kg))]  i.e. [ J/g ]
   h_CT  :  The first derivative of specific enthalpy with respect to
            CT at constant SA and p.                           [ J/(kg K) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_enthalpy_first_derivatives(double sa,double ct,double p,double *h_sa,double *h_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double dynamic_h_ct_part, dynamic_h_sa_part, xs, ys, z;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*gtc.rec_db2pa;
	if (h_sa != NULL)
	{
		dynamic_h_sa_part = z*(gsvco.h101
									  + xs*(2.0*gsvco.h201 + xs*(3.0*gsvco.h301
											  + xs*(4.0*gsvco.h401 + xs*(5.0*gsvco.h501
													  + 6.0*gsvco.h601*xs)))) + ys*(gsvco.h111
															  + xs*(2.0*gsvco.h211 + xs*(3.0*gsvco.h311
																	  + xs*(4.0*gsvco.h411 + 5.0*gsvco.h511*xs)))
															  + ys*(gsvco.h121 + xs*(2.0*gsvco.h221
																	  + xs*(3.0*gsvco.h321 + 4.0*gsvco.h421*xs))
																	  + ys*(gsvco.h131 + xs*(2.0*gsvco.h231
																			  + 3.0*gsvco.h331*xs) + ys*(gsvco.h141
																					  + 2.0*gsvco.h241*xs + gsvco.h151*ys))))
									  + z*(gsvco.h102 + xs*(2.0*gsvco.h202
											  + xs*(3.0*gsvco.h302 + xs*(4.0*gsvco.h402
													  + 5.0*gsvco.h502*xs))) + ys*(gsvco.h112
															  + xs*(2.0*gsvco.h212 + xs*(3.0*gsvco.h312
																	  + 4.0*gsvco.h412*xs)) + ys*(gsvco.h122
																			  + xs*(2.0*gsvco.h222 + 3.0*gsvco.h322*xs)
																			  + ys*(gsvco.h132 + 2.0*gsvco.h232*xs
																					  + gsvco.h142*ys ))) + z*(gsvco.h103
																							  + xs*(2.0*gsvco.h203 + xs*(3.0*gsvco.h303
																									  + 4.0*gsvco.h403*xs)) + ys*(gsvco.h113
																											  + xs*(2.0*gsvco.h213 + 3.0*gsvco.h313*xs)
																											  + ys*(gsvco.h123 + 2.0*gsvco.h223*xs
																													  + gsvco.h133*ys)) + z*(gsvco.h104
																															  + 2.0*gsvco.h204*xs + gsvco.h114*ys
																															  + gsvco.h105*z))));
		*h_sa = 1e8*0.5*gtc.gsw_sfac*dynamic_h_sa_part/xs;
	}
	if (h_ct != NULL)
	{
		dynamic_h_ct_part =
			z*(gsvco.h011 + xs*(gsvco.h111
									  + xs*(gsvco.h211 + xs*(gsvco.h311
											  + xs*(gsvco.h411 + gsvco.h511*xs))))
				+ ys*(2.0*(gsvco.h021 + xs*(gsvco.h121
													 + xs*(gsvco.h221 + xs*(gsvco.h321
															 + gsvco.h421*xs)))) + ys*(3.0*(gsvco.h031
																	 + xs*(gsvco.h131 + xs*(gsvco.h231
																			 + gsvco.h331*xs))) + ys*(4.0*(gsvco.h041
																					 + xs*(gsvco.h141 + gsvco.h241*xs))
																					 + ys*(5.0*(gsvco.h051 + gsvco.h151*xs)
																							 + 6.0*gsvco.h061*ys)))) + z*(gsvco.h012
																									 + xs*(gsvco.h112 + xs*(gsvco.h212
																											 + xs*(gsvco.h312 + gsvco.h412*xs)))
																									 + ys*(2.0*(gsvco.h022 + xs*(gsvco.h122
																											 + xs*(gsvco.h222 + gsvco.h322*xs)))
																											 + ys*(3.0*(gsvco.h032 + xs*(gsvco.h132
																													 + gsvco.h232*xs)) + ys*(4.0*(gsvco.h042
																															 + gsvco.h142*xs) + 5.0*gsvco.h052*ys)))
																									 + z*(gsvco.h013 + xs*(gsvco.h113
																											 + xs*(gsvco.h213 + gsvco.h313*xs))
																											 + ys*(2.0*(gsvco.h023 + xs*(gsvco.h123
																													 + gsvco.h223*xs)) + ys*(3.0*(gsvco.h033
																															 + gsvco.h133*xs) + 4.0*gsvco.h043*ys))
																											 + z*(gsvco.h014 + gsvco.h114*xs
																													 + 2.0*gsvco.h024*ys + gsvco.h015*z))));
		*h_ct = gtc.gsw_cp0 + 1e8*0.025*dynamic_h_ct_part;
	}
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy_first_derivatives_ct_exact (sa,ct,p,*h_sa,*h_ct)
==========================================================================
   Calculates the following two derivatives of specific enthalpy (h)
    (1) h_SA, the derivative with respect to Absolute Salinity at
        constant CT and p, and
    (2) h_CT, derivative with respect to CT at constant SA and p.
   Note that h_P is specific volume (1/rho) it can be calulated by calling
   gsw_specvol_CT_exact(SA,CT,p). This method uses the full Gibbs method.
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   h_SA  :  The first derivative of specific enthalpy with respect to
            Absolute Salinity at constant CT and p.[ J/(kg (g/kg))]  i.e. [ J/g ]
   h_CT  :  The first derivative of specific enthalpy with respect to
            CT at constant SA and p.                           [ J/(kg K) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_enthalpy_first_derivatives_ct_exact(double sa,double ct,double p,double *h_sa,double *h_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double g_sa_mod_pt, g_sa_mod_t, pt0, t, temp_ratio, x, y, y_pt, z;
	t = gsw_t_from_ct(sa,ct,p);
	pt0 = gsw_pt_from_ct(sa,ct);
	temp_ratio = (gtc.gsw_t0 + t)/(gtc.gsw_t0 + pt0);
	if (h_ct != NULL)
	{
		*h_ct = gtc.gsw_cp0*temp_ratio;
	}
	if (h_sa == NULL)
	{
		return;
	}
	x = sqrt(gtc.gsw_sfac*sa);
	y = 0.025*t;
	z = gtc.rec_db2pa*p;
	/**note.the input pressure (p) is sea pressure in units of dbar.*/
	g_sa_mod_t = 8645.36753595126 + z*(-6620.98308089678 +
												  z*(769.588305957198 + z*(-193.0648640214916 +
														  (31.6816345533648 - 5.24960313181984*z)*z))) +
					 x*(-7296.43987145382 + x*(8103.20462414788 +
														y*(2175.341332000392 + y*(-274.2290036817964 +
																y*(197.4670779425016 + y*(-68.5590309679152 + 9.98788038278032*y)))
																- 90.6734234051316*z) +
														x*(-5458.34205214835 - 980.14153344888*y +
																x*(2247.60742726704 - 340.1237483177863*x + 220.542973797483*y)
																+ 180.142097805543*z) +
														z*(-219.1676534131548 + (-16.32775915649044
																- 120.7020447884644*z)*z)) +
						 z*(598.378809221703 + z*(-156.8822727844005 + (204.1334828179377
								 - 10.23755797323846*z)*z)) +
						 y*(-1480.222530425046 + z*(-525.876123559641 + (249.57717834054571
								 - 88.449193048287*z)*z) +
							 y*(-129.1994027934126 + z*(1149.174198007428 +
									 z*(-162.5751787551336 + 76.9195462169742*z)) +
								 y*(-30.0682112585625 - 1380.9597954037708*z + y*(2.626801985426835
										 + 703.695562834065*z))))) +
					 y*(1187.3715515697959 + z*(1458.233059470092 +
														 z*(-687.913805923122 + z*(249.375342232496 + z*(-63.313928772146
																 + 14.09317606630898*z)))) +
						 y*(1760.062705994408 + y*(-450.535298526802 +
								 y*(182.8520895502518 + y*(-43.3206481750622 + 4.26033941694366*y) +
									 z*(-595.457483974374 + (149.452282277512 - 72.9745838003176*z)*z)) +
								 z*(1388.489628266536 + z*(-409.779283929806 + (227.123395681188
										 - 22.2565468652826*z)*z))) +
							 z*(-1721.528607567954 + z*(674.819060538734 +
									 z*(-356.629112415276 + (88.4080716616 - 15.84003094423364*z)*z)))));
	g_sa_mod_t = 0.5*gtc.gsw_sfac*g_sa_mod_t;
	y_pt = 0.025*pt0;
	g_sa_mod_pt = 8645.36753595126 +
					  x*(-7296.43987145382 + x*(8103.20462414788 +
														 y_pt*(2175.341332000392 + y_pt*(-274.2290036817964 +
																 y_pt*(197.4670779425016 + y_pt*(-68.5590309679152
																		 + 9.98788038278032*y_pt)))) +
														 x*(-5458.34205214835 - 980.14153344888*y_pt +
																 x*(2247.60742726704 - 340.1237483177863*x
																		 + 220.542973797483*y_pt))) +
						  y_pt*(-1480.222530425046 + y_pt*(-129.1994027934126 +
								  y_pt*(-30.0682112585625 + y_pt*2.626801985426835)))) +
					  y_pt*(1187.3715515697959 + y_pt*(1760.062705994408
							  + y_pt*(-450.535298526802 +
										 y_pt*(182.8520895502518 + y_pt*(-43.3206481750622
												 + 4.26033941694366*y_pt)))));
	g_sa_mod_pt = 0.5*gtc.gsw_sfac*g_sa_mod_pt;
	*h_sa = g_sa_mod_t - temp_ratio*g_sa_mod_pt;
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy_ice (t, p)
==========================================================================
  Calculates the specific enthalpy of ice (h_Ih).
   t  =  in-situ temperature (ITS-90)                             [ deg C ]
   p  =  sea pressure                                              [ dbar ]
         ( i.e. absolute pressure - 10.1325 dbar )
   gsw_enthalpy_ice  :  specific enthalpy of ice                   [ J/kg ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_enthalpy_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double tau, dzi, g0;
	complex<double> r2, sqtau_t1, sqtau_t2, g;
	tau = (t + gtc.gsw_t0)*gic.rec_tt;
	dzi = gtc.db2pa*p*gic.rec_pt;
	g0 = gic.g00 + dzi*(gic.g01 + dzi*(gic.g02 + dzi*(gic.g03 + gic.g04*dzi)));
	r2 = r20 + dzi*(r21 + r22*dzi);
	sqtau_t1 = (tau*tau)/(t1*t1);
	sqtau_t2 = (tau*tau)/(t2*t2);
	g = r1*t1*(log(1.0 - sqtau_t1) + sqtau_t1)
		 + r2*t2*(log(1.0 - sqtau_t2) + sqtau_t2);
	return (g0 + gic.tt*real(g));
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy_second_derivatives_ct_exact (sa, ct, p, &
                                                 h_sa_sa, h_sa_ct, h_ct_ct)
==========================================================================
   Calculates three second-order derivatives of specific enthalpy (h).
   Note that this method uses the full Gibbs function.
   sa  :  Absolute Salinity                                        [ g/kg ]
   ct  :  Conservative Temperature (ITS-90)                       [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   h_sa_sa  :  The second derivative of specific enthalpy with respect to
               Absolute Salinity at constant ct & p.    [ J/(kg (g/kg)^2) ]
   h_sa_ct  :  The second derivative of specific enthalpy with respect to
               sa and ct at constant p.                  [ J/(kg K(g/kg)) ]
   h_ct_ct  :  The second derivative of specific enthalpy with respect to
               ct at constant sa and p.                      [ J/(kg K^2) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_enthalpy_second_derivatives_ct_exact(double sa,double ct,double p,
		double &h_sa_sa,double &h_sa_ct,double &h_ct_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double factor, gsa_pt0, gsat_pt0, gsat, part_b, pt0, h_ct_ct_val,
			 rec_abs_pt0, rec_gtt_pt0, rec_gtt, t, temp_ratio,
			 gsasa, gsasa_pt0, pr0 = 0.0, sa_small = 1e-100;
	int n0=0, n1=1, n2=2;
	pt0 = gsw_pt_from_ct(sa,ct);
	rec_abs_pt0 = 1.0/(gtc.gsw_t0 + pt0);
	t = gsw_pt_from_t(sa,pt0,pr0,p);
	temp_ratio = (gtc.gsw_t0 + t)*rec_abs_pt0;
	rec_gtt_pt0 = 1.0/gsw_gibbs(n0,n2,n0,sa,pt0,pr0);
	rec_gtt = 1.0/gsw_gibbs(n0,n2,n0,sa,t,p);
	/**h_ct_ct is naturally well-behaved as sa approaches zero.*/
	h_ct_ct_val = gtc.gsw_cp0*gtc.gsw_cp0*
					  (temp_ratio*rec_gtt_pt0 - rec_gtt)*(rec_abs_pt0*rec_abs_pt0);
	//if (h_ct_ct != NULL)
	//{
	h_ct_ct = h_ct_ct_val;
	//}
	//if ((h_sa_sa == NULL) && (h_sa_ct == NULL))
	//{
	//return;
	//}
	gsat_pt0 = gsw_gibbs(n1,n1,n0,sa,pt0,pr0);
	gsat = gsw_gibbs(n1,n1,n0,sa,t,p);
	gsa_pt0 = gsw_gibbs(n1,n0,n0,sa,pt0,pr0);
	part_b = (temp_ratio*gsat_pt0*rec_gtt_pt0 - gsat*rec_gtt)*rec_abs_pt0;
	factor = gsa_pt0/gtc.gsw_cp0;
	//if (h_sa_sa != NULL)
	//{
	gsasa = gsw_gibbs(n2,n0,n0,sa,t,p);
	gsasa_pt0 = gsw_gibbs(n2,n0,n0,sa,pt0,pr0);
	/**
	    h_sa_sa has a singularity at sa = 0, and blows up as sa
	    approaches zero.
	*/
	h_sa_sa = gsasa - temp_ratio*gsasa_pt0
				 + temp_ratio*gsat_pt0*gsat_pt0*rec_gtt_pt0
				 - gsat*gsat*rec_gtt
				 - 2.0*gsa_pt0*part_b + (factor*factor)*h_ct_ct_val;
	//}
	//if (h_sa_ct == NULL)
	//{
	//return;
	//}
	/**
	   h_sa_ct should not blow up as sa approaches zero.
	   The following lines of code ensure that the h_sa_ct output
	   of this method does not blow up in this limit.
	   That is, when sa < 1e-100 g/kg, we force the h_sa_ct
	   output to be the same as if sa = 1e-100 g/kg.
	*/
	if (sa < sa_small)
	{
		rec_gtt_pt0 = 1.0/gsw_gibbs(n0,n2,n0,sa_small,pt0,pr0);
		rec_gtt = 1.0/gsw_gibbs(n0,n2,n0,sa_small,t,p);
		gsat_pt0 = gsw_gibbs(n1,n1,n0,sa_small,pt0,pr0);
		gsat = gsw_gibbs(n1,n1,n0,sa_small,t,p);
		gsa_pt0 = gsw_gibbs(n1,n0,n0,sa_small,pt0,pr0);
		part_b = (temp_ratio*gsat_pt0*rec_gtt_pt0
					 - gsat*rec_gtt)*rec_abs_pt0;
		factor = gsa_pt0/gtc.gsw_cp0;
	}
	h_sa_ct  = gtc.gsw_cp0 * part_b - factor * h_ct_ct_val;
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy_sso_0(p)
==========================================================================
   This method calculates enthalpy at the Standard Ocean Salinity, SSO,
   and at a Conservative Temperature of zero degrees C, as a function of
   pressure, p, in dbar, using a streamlined version of the
   computationally-efficient expression for specific volume, that is, a
   streamlined version of the code "gsw_enthalpy(SA,CT,p)".
  p      : sea pressure                                    [dbar]
  enthalpy_sso_0 : enthalpy(sso,0,p)
*************************************************************************/

double TeosBase::gsw_enthalpy_sso_0(double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double dynamic_enthalpy_sso_0_p, z;
	z = p*1.0e-4;
	dynamic_enthalpy_sso_0_p =
		z*( 9.726613854843870e-4 + z*(-2.252956605630465e-5
												+ z*( 2.376909655387404e-6 + z*(-1.664294869986011e-7
														+ z*(-5.988108894465758e-9 + z*(gsvco.h006 + gsvco.h007*z))))));
	return (dynamic_enthalpy_sso_0_p*gtc.db2pa*1.0e4);
}

/**************************************************************************
==========================================================================
method: gsw_enthalpy_t_exact(sa,t,p)
==========================================================================
 Calculates the specific enthalpy of seawater
 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]
 gsw_enthalpy_t_exact : specific enthalpy                 [J/kg]
*************************************************************************/

double TeosBase::gsw_enthalpy_t_exact(double sa, double t, double p)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	int n0=0, n1=1;
	return (gsw_gibbs(n0,n0,n0,sa,t,p) -
			  (t+gtc.gsw_t0)*gsw_gibbs(n0,n1,n0,sa,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_entropy_first_derivatives (sa, ct, eta_sa, eta_ct)
=========================================================================
   Calculates the following two partial derivatives of specific entropy
   (eta)
    (1) eta_SA, the derivative with respect to Absolute Salinity at
        constant Conservative Temperature, and
    (2) eta_CT, the derivative with respect to Conservative Temperature at
        constant Absolute Salinity.
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   eta_SA :  The derivative of specific entropy with respect to
             Absolute Salinity (in units of g kg^-1) at constant
             Conservative Temperature.
             eta_SA has units of:         [ J/(kg K(g/kg))]  or [ J/(g K) ]
   eta_CT :  The derivative of specific entropy with respect to
             Conservative Temperature at constant Absolute Salinity.
             eta_CT has units of:                            [ J/(kg K^2) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_entropy_first_derivatives(double sa,double ct,double *eta_sa,double *eta_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double pt, pr0 = 0.0;
	int n0=0, n1=1;
	pt = gsw_pt_from_ct(sa,ct);
	if (eta_sa != NULL)
	{
		*eta_sa = -(gsw_gibbs(n1,n0,n0,sa,pt,pr0))/(gtc.gsw_t0 + pt);
	}
	if (eta_ct != NULL)
	{
		*eta_ct = gtc.gsw_cp0/(gtc.gsw_t0 + pt);
	}
}

/**************************************************************************
=========================================================================
method: gsw_entropy_from_ct (sa, ct)
=========================================================================
   Calculates specific entropy of seawater from Conservative Temperature.
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   entropy  =  specific entropy                                   [ deg C ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_entropy_from_ct(double sa,double ct)
{
	double pt0 = gsw_pt_from_ct(sa, ct);
	return (-gsw_gibbs(0,1,0,sa,pt0,0));
}

/**************************************************************************
==========================================================================
method: gsw_entropy_ice (t, p)
==========================================================================
   Calculates specific entropy of ice.
   t  :  in-situ temperature (ITS-90)                             [ deg C ]
   p  :  sea pressure                                              [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   ice_entropy  :  specific entropy of ice                 [ J kg^-1 K^-1 ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_entropy_ice(double t,double p)
{
	return (-gsw_gibbs_ice(1,0,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_entropy_part(sa,t,p)
==========================================================================
 entropy minus the terms that are a function of only SA
 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]
 entropy_part : entropy part
*************************************************************************/

double TeosBase::gsw_entropy_part(double sa, double t, double p)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double x2, x, y, z, rval;
	x2	= gtc.gsw_sfac*sa;
	x	= sqrt(x2);
	y	= t*0.025;
	z	= p*1e-4;
	rval = gsw_get_g03plusg08(z, y, x2, x);
	return rval;
}

/**************************************************************************
==========================================================================
method: gsw_entropy_part_zerop(sa,pt0)
==========================================================================
  entropy part evaluated at the sea surface
  sa     : Absolute Salinity                               [g/kg]
  pt0    : insitu temperature                              [deg C]
  entropy_part_zerop : entropy part at the sea surface
*************************************************************************/

double TeosBase::gsw_entropy_part_zerop(double sa, double pt0)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double x2, x, y, g03, g08;
	x2 = gtc.gsw_sfac*sa;
	x = sqrt(x2);
	y = pt0*0.025;
	g03 = y*(-24715.571866078 + y*(2210.2236124548363 +
											 y*(-592.743745734632 + y*(290.12956292128547 +
													 y*(-113.90630790850321 + y*21.35571525415769)))));
	g08 = x2*(x*(x*(y*(-137.1145018408982 + y*(148.10030845687618 +
							 y*(-68.5590309679152 + 12.4848504784754*y)))) +
					 y*(-86.1329351956084 + y*(-30.0682112585625 +
														y*3.50240264723578))) +
				 y*(1760.062705994408 + y*(-675.802947790203 +
													y*(365.7041791005036 + y*(-108.30162043765552 +
															12.78101825083098*y)))));
	return (-(g03 + g08)*0.025);
}

/**************************************************************************
==========================================================================
method: gsw_entropy_second_derivatives (sa, ct, eta_sa_sa, &
                                                     eta_sa_ct, eta_ct_ct)
=========================================================================
   Calculates the following three second-order partial derivatives of
   specific entropy (eta)
    (1) eta_SA_SA, the second derivative with respect to Absolute
        Salinity at constant Conservative Temperature, and
    (2) eta_SA_CT, the derivative with respect to Absolute Salinity and
        Conservative Temperature.
    (3) eta_CT_CT, the second derivative with respect to Conservative
        Temperature at constant Absolute Salinity.
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   eta_SA_SA :  The second derivative of specific entropy with respect
                to Absolute Salinity (in units of g kg^-1) at constant
                Conservative Temperature.
                eta_SA_SA has units of:                 [ J/(kg K(g/kg)^2)]
   eta_SA_CT :  The second derivative of specific entropy with respect
                to Conservative Temperature at constant Absolute
                Salinity. eta_SA_CT has units of:     [ J/(kg (g/kg) K^2) ]
   eta_CT_CT :  The second derivative of specific entropy with respect
                to Conservative Temperature at constant Absolute
                Salinity.  eta_CT_CT has units of:           [ J/(kg K^3) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_entropy_second_derivatives(double sa,double ct,double *eta_sa_sa,double *eta_sa_ct,double *eta_ct_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double abs_pt, ct_pt, ct_sa, pt, ct_ct, pr0 = 0.0;
	int n0=0, n1=1, n2=2;
	pt = gsw_pt_from_ct(sa,ct);
	abs_pt = gtc.gsw_t0 + pt;
	ct_pt = -(abs_pt*gsw_gibbs(n0,n2,n0,sa,pt,pr0))/gtc.gsw_cp0;
	ct_ct = -gtc.gsw_cp0/(ct_pt*abs_pt*abs_pt);
	if ((eta_sa_ct != NULL) || (eta_sa_sa != NULL))
	{
		ct_sa = (gsw_gibbs(n1,n0,n0,sa,pt,pr0) -
					(abs_pt*gsw_gibbs(n1,n1,n0,sa,pt,pr0)))/gtc.gsw_cp0;
		if (eta_sa_ct != NULL)
		{
			*eta_sa_ct = -ct_sa*ct_ct;
		}
		if (eta_sa_sa != NULL)
		{
			*eta_sa_sa = -gsw_gibbs(n2,n0,n0,sa,pt,pr0)/abs_pt +
							 ct_sa*ct_sa*ct_ct;
		}
	}
	if (eta_ct_ct != NULL)
	{
		*eta_ct_ct = ct_ct;
	}
}

/**************************************************************************
==========================================================================
method: gsw_fdelta(p,lon,lat)
==========================================================================
  Calculates fdelta.
  p      : sea pressure                                    [dbar]
  lon    : longitude                                       [deg E]
  lat    : latitude                                        [deg N]
  gsw_fdelta : Absolute Salinty Anomaly                    [unitless]
*************************************************************************/

double TeosBase::gsw_fdelta(double p, double lon, double lat)
{
	double sa, saar;
	saar = gsw_saar(p,lon,lat);
	if (saar >= cppGSW_ERROR_LIMIT)
	{
		sa  = cppGSW_INVALID_VALUE;
	}
	else
	{
		sa  = ((1.0 + 0.35)*saar)/(1.0 - 0.35*saar);
	}
	return (sa);
}

/**************************************************************************
==========================================================================
method: gsw_frazil_properties (sa_bulk, h_bulk, p, &
                                            sa_final, ct_final, w_ih_final)
==========================================================================
   Calculates the mass fraction of ice (mass of ice divided by mass of ice
   plus seawater), w_Ih_final, which results from given values of the bulk
   Absolute Salinity, SA_bulk, bulk enthalpy, h_bulk, occuring at pressure
   p.  The final values of Absolute Salinity, SA_final, and Conservative
   Temperature, CT_final, of the interstitial seawater phase are also
   returned.  This code assumes that there is no dissolved air in the
   seawater (that is, saturation_fraction is assumed to be zero
   throughout the code).
   When the mass fraction w_Ih_final is calculated as being a positive
   value, the seawater-ice mixture is at thermodynamic equlibrium.
   This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
   is sufficiently large (i.e. sufficiently "warm") so that there is no ice
   present in the final state.  In this case the final state consists of
   only seawater rather than being an equlibrium mixture of seawater and
   ice which occurs when w_Ih_final is positive.  Note that when
   w_Ih_final = 0, the final seawater is not at the freezing temperature.
   SA_bulk :  bulk Absolute Salinity of the seawater and ice mixture[ g/kg ]
   h_bulk  :  bulk enthalpy of the seawater and ice mixture        [ J/kg ]
   p       :  sea pressure                                         [ dbar ]
              ( i.e. absolute pressure - 10.1325 dbar )
   SA_final    :  Absolute Salinity of the seawater in the final state,
                  whether or not any ice is present.               [ g/kg ]
   CT_final    :  Conservative Temperature of the seawater in the the final
                  state, whether or not any ice is present.       [ deg C ]
   w_Ih_final  :  mass fraction of ice in the final seawater-ice mixture.
                  If this ice mass fraction is positive, the system is at
                  thermodynamic equilibrium.  If this ice mass fraction is
                  zero there is no ice in the final state which consists
                  only of seawater which is warmer than the freezing
                  temperature.                                   [unitless]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_frazil_properties(double sa_bulk,double h_bulk,double p,double *sa_final,double *ct_final,double *w_ih_final)
{
	int number_of_iterations;
	double cp_ih, ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly,
			 func, func0, hf, h_hat_ct, h_hat_sa,
			 h_ihf, sa, tf_sa, tf, w_ih_mean, w_ih_old, w_ih,
			 /**
			    Throughout this code seawater is taken to contain
			    no dissolved air.
			 */
			 saturation_fraction = 0.0,
			 num_f = 5.0e-2, num_f2 = 6.9e-7, num_p = 2.21;
	/**
	    ---------------
	      Finding func0
	    --------------
	*/
	ctf = gsw_ct_freezing(sa_bulk,p,saturation_fraction);
	func0 = h_bulk - gsw_enthalpy_ct_exact(sa_bulk,ctf,p);
	/**
	    -----------------------------------------------------------------------
	      When func0 is zero or positive we can immediately calculate the three
	      outputs, as the bulk enthalpy, h_bulk, is too large to allow any ice
	      at thermodynamic equilibrium. The result will be (warm) seawater with
	      no frazil ice being present. The three outputs can be set and the rest
	      of this code does not need to be performed.
	    -----------------------------------------------------------------------
	*/
	if (func0 >= 0.0)
	{
		*sa_final = sa_bulk;
		*ct_final = gsw_ct_from_enthalpy_exact(sa_bulk,h_bulk,p);
		*w_ih_final = 0.0;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	      Begin to find the solution for those data points that have func0 < 0,
	      implying that the output will be a positive ice mass fraction
	      w_Ih_final.
	      Do a quasi-Newton step with a separate polynomial estimate of the
	      derivative of func with respect to the ice mass fraction.  This
	      section of the code delivers initial values of both w_Ih and SA to
	      the rest of the more formal modified Newtons Method approach of
	      McDougall and Wotherspoon (2014).
	    -----------------------------------------------------------------------
	*/
	dfunc_dw_ih_mean_poly = 3.347814e+05
									- num_f*func0*(1.0 + num_f2*func0) - num_p*p;
	w_ih = gsw_min_d(-func0/dfunc_dw_ih_mean_poly, 0.95);
	sa = sa_bulk/(1.0 - w_ih);
	if (sa < 0.0 || sa > 120.0)
	{
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		*w_ih_final = *sa_final;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	      Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
	      fed into the iterative Newton s Method.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	hf = gsw_enthalpy_ct_exact(sa,ctf,p);
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	h_ihf = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p, &h_hat_sa,&h_hat_ct);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa, NULL);
	gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa,NULL);
	dfunc_dw_ih = hf - h_ihf - sa*(h_hat_sa + h_hat_ct*ctf_sa +
											 w_ih*cp_ih*tf_sa/(1.0 - w_ih));
	/**
	    -----------------------------------------------------------------------
	     Enter the main McDougall-Wotherspoon (2014) modified Newton-Raphson
	     loop
	    -----------------------------------------------------------------------
	*/
	for (number_of_iterations = 1; number_of_iterations <= 3;
			number_of_iterations++)
	{
		if (number_of_iterations > 1)
		{
			/** on the first iteration these values are already known */
			ctf = gsw_ct_freezing(sa,p,saturation_fraction);
			hf = gsw_enthalpy_ct_exact(sa,ctf,p);
			tf = gsw_t_freezing(sa,p,saturation_fraction);
			h_ihf = gsw_enthalpy_ice(tf,p);
		}
		func = h_bulk - (1.0 - w_ih)*hf - w_ih*h_ihf;
		w_ih_old = w_ih;
		w_ih = w_ih_old - func/dfunc_dw_ih;
		w_ih_mean = 0.5*(w_ih + w_ih_old);
		if (w_ih_mean > 0.9)
		{
			/**This ensures that the mass fraction of ice never exceeds 0.9*/
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih_mean);
		ctf = gsw_ct_freezing(sa,p,saturation_fraction);
		hf = gsw_enthalpy_ct_exact(sa,ctf,p);
		tf = gsw_t_freezing(sa,p,saturation_fraction);
		h_ihf = gsw_enthalpy_ice(tf,p);
		cp_ih = gsw_cp_ice(tf,p);
		gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p, &h_hat_sa, &h_hat_ct);
		gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa, NULL);
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa, NULL);
		dfunc_dw_ih = hf - h_ihf - sa*(h_hat_sa + h_hat_ct*ctf_sa
												 + w_ih_mean*cp_ih*tf_sa/
												 (1.0 - w_ih_mean));
		w_ih = w_ih_old - func/dfunc_dw_ih;
		if (w_ih > 0.9)
		{
			/**This ensures that the mass fraction of ice never exceeds 0.9*/
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih);
	} /** for */
	*sa_final = sa;
	*ct_final = gsw_ct_freezing(sa,p,saturation_fraction);
	*w_ih_final = w_ih;
	if (*w_ih_final < 0.0)
	{
		/**
		      This will only trap cases that are smaller than zero by just
		      machine precision
		*/
		*sa_final = sa_bulk;
		*ct_final = gsw_ct_from_enthalpy_exact(*sa_final,h_bulk,p);
		*w_ih_final = 0.0;
	}
}

/**************************************************************************
==========================================================================
method: gsw_frazil_properties_potential (sa_bulk, h_pot_bulk,&
                                         p, sa_final, ct_final, w_ih_final)
==========================================================================
   Calculates the mass fraction of ice (mass of ice divided by mass of ice
   plus seawater), w_Ih_final, which results from given values of the bulk
   Absolute Salinity, SA_bulk, bulk potential enthalpy, h_pot_bulk,
   occuring at pressure p.  The final equilibrium values of Absolute
   Salinity, SA_final, and Conservative Temperature, CT_final, of the
   interstitial seawater phase are also returned.  This code assumes that
   there is no dissolved air in the seawater (that is, saturation_fraction
   is assumed to be zero thoughout the code).
   When the mass fraction w_Ih_final is calculated as being a positive
   value, the seawater-ice mixture is at thermodynamic equlibrium.
   This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
   is sufficiently large (i.e. sufficiently "warm") so that there is no ice
   present in the final state.  In this case the final state consists of
   only seawater rather than being an equlibrium mixture of seawater and
   ice which occurs when w_Ih_final is positive.  Note that when
   w_Ih_final = 0, the final seawater is not at the freezing temperature.
   Note that this code uses the exact forms of CT_freezing and
   pot_enthalpy_ice_freezing.
   SA_bulk     :  bulk Absolute Salinity of the seawater and ice mixture[ g/kg ]
   h_pot_bulk  :  bulk potential enthalpy of the seawater and ice mixture[ J/kg ]
   p           :  sea pressure                                  [ dbar ]
                   ( i.e. absolute pressure - 10.1325 dbar )
   SA_final    :  Absolute Salinity of the seawater in the final state,
                  whether or not any ice is present.               [ g/kg ]
   CT_final    :  Conservative Temperature of the seawater in the the final
                  state, whether or not any ice is present.       [ deg C ]
   w_Ih_final  :  mass fraction of ice in the final seawater-ice mixture.
                  If this ice mass fraction is positive, the system is at
                  thermodynamic equilibrium.  If this ice mass fraction is
                  zero there is no ice in the final state which consists
                  only of seawater which is warmer than the freezing
                  temperature.                                   [unitless]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_frazil_properties_potential(double sa_bulk,double h_pot_bulk,double p,double *sa_final,double *ct_final,double *w_ih_final)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int iterations, max_iterations;
	double ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly,
			 dpot_h_ihf_dsa, func, func0, h_pot_ihf, sa, w_ih_old, w_ih,
			 x, xa, y, z;
	double f01 = -9.041191886754806e-1,
			 f02 =  4.169608567309818e-2,
			 f03 = -9.325971761333677e-3,
			 f04 =  4.699055851002199e-2,
			 f05 = -3.086923404061666e-2,
			 f06 =  1.057761186019000e-2,
			 f07 = -7.349302346007727e-2,
			 f08 =  1.444842576424337e-1,
			 f09 = -1.408425967872030e-1,
			 f10 =  1.070981398567760e-1,
			 f11 = -1.768451760854797e-2,
			 f12 = -4.013688314067293e-1,
			 f13 =  7.209753205388577e-1,
			 f14 = -1.807444462285120e-1,
			 f15 =  1.362305015808993e-1,
			 f16 = -9.500974920072897e-1,
			 f17 =  1.192134856624248,
			 f18 = -9.191161283559850e-2,
			 f19 = -1.008594411490973,
			 f20 =  8.020279271484482e-1,
			 f21 = -3.930534388853466e-1,
			 f22 = -2.026853316399942e-2,
			 f23 = -2.722731069001690e-2,
			 f24 =  5.032098120548072e-2,
			 f25 = -2.354888890484222e-2,
			 f26 = -2.454090179215001e-2,
			 f27 =  4.125987229048937e-2,
			 f28 = -3.533404753585094e-2,
			 f29 =  3.766063025852511e-2,
			 f30 = -3.358409746243470e-2,
			 f31 = -2.242158862056258e-2,
			 f32 =  2.102254738058931e-2,
			 f33 = -3.048635435546108e-2,
			 f34 = -1.996293091714222e-2,
			 f35 =  2.577703068234217e-2,
			 f36 = -1.292053030649309e-2,
			 g01 =  3.332286683867741e5,
			 g02 =  1.416532517833479e4,
			 g03 = -1.021129089258645e4,
			 g04 =  2.356370992641009e4,
			 g05 = -8.483432350173174e3,
			 g06 =  2.279927781684362e4,
			 g07 =  1.506238790315354e4,
			 g08 =  4.194030718568807e3,
			 g09 = -3.146939594885272e5,
			 g10 = -7.549939721380912e4,
			 g11 =  2.790535212869292e6,
			 g12 =  1.078851928118102e5,
			 g13 = -1.062493860205067e7,
			 g14 =  2.082909703458225e7,
			 g15 = -2.046810820868635e7,
			 g16 =  8.039606992745191e6,
			 g17 = -2.023984705844567e4,
			 g18 =  2.871769638352535e4,
			 g19 = -1.444841553038544e4,
			 g20 =  2.261532522236573e4,
			 g21 = -2.090579366221046e4,
			 g22 = -1.128417003723530e4,
			 g23 =  3.222965226084112e3,
			 g24 = -1.226388046175992e4,
			 g25 =  1.506847628109789e4,
			 g26 = -4.584670946447444e4,
			 g27 =  1.596119496322347e4,
			 g28 = -6.338852410446789e4,
			 g29 =  8.951570926106525e4,
			 saturation_fraction = 0.0;
	/**
	-----------------------------------------------------------------------
	   Finding func0.  This is the value of the function, func, that would
	   result in the output w_Ih_final being exactly zero.
	-----------------------------------------------------------------------
	*/
	func0 = h_pot_bulk - gtc.gsw_cp0
			  *gsw_ct_freezing(sa_bulk,p,saturation_fraction);
	/**
	-------------------------------------------------------------------------
	   Setting the three outputs for data points that have func0 non-negative
	-------------------------------------------------------------------------
	*/
	if (func0 >= 0.0)
	{
		/**
		    When func0 is zero or positive then the final answer will contain
		    no frazil ice; that is, it will be pure seawater that is warmer
		    than the freezing temperature. If func0 >= 0 we do not need to go
		    through the modified Newton-Raphson procedure and we can simply
		    write down the answer, as in the following 4 lines of code.
		*/
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	    Begin finding the solution for data points that have func0 < 0, so that
	    the output will have a positive ice mass fraction w_Ih_final.
	    -----------------------------------------------------------------------
	*/
	/**Evalaute a polynomial for w_Ih in terms of SA_bulk, func0 and p*/
	x = sa_bulk*1e-2;
	y = func0/3e5;
	z = p*1e-4;
	w_ih = y*(f01 + x*(f02 + x*(f03 + x*(f04 + x*(f05 + f06*x))))
				 + y*(f07 + x*(f08 + x*(f09 + x*(f10 + f11*x))) + y*(f12 + x*(f13
						 + x*(f14 + f15*x)) + y*(f16 + x*(f17 + f18*x) + y*(f19 + f20*x
														 + f21*y)))) + z*(f22 + x*(f23 + x*(f24 + f25*x)) + y*(x*(f26
																 + f27*x)
																 + y*(f28 + f29*x + f30*y)) + z*(f31 + x*(f32 + f33*x) + y*(f34
																		 + f35*x + f36*y))));
	if (w_ih > 0.9)
	{
		/**
		   The ice mass fraction out of this code is restricted to be
		   less than 0.9.
		*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		*w_ih_final = *sa_final;
		return;
	}
	/**
	    The initial guess at the absolute salinity of the interstitial seawater
	*/
	sa = sa_bulk/(1.0 - w_ih);
	/**
	    -----------------------------------------------------------------------
	      Doing a Newton step with a separate polynomial estimate of the mean
	      derivative dfunc_dw_Ih_mean_poly.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing(sa,p);
	func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
	xa = sa*1e-2;
	dfunc_dw_ih_mean_poly = g01 + xa*(g02 + xa*(g03 + xa*(g04 + g05*xa)))
									+ w_ih*(xa*(g06 + xa*(g07 + g08*xa)) + w_ih*(xa*(g09 + g10*xa)
											  + w_ih*xa*(g11 + g12*xa + w_ih*(g13 + w_ih*(g14 + w_ih*(g15
													  + g16*w_ih)))))) + z*(g17 + xa*(g18 + g19*xa) + w_ih*(g20
															  + w_ih*(g21 + g22*w_ih) + xa*(g23 + g24*xa*w_ih))
															  + z*(g25 + xa*(g26 + g27*xa) + w_ih*(g28 + g29*w_ih)));
	w_ih_old = w_ih;
	w_ih = w_ih_old - func/dfunc_dw_ih_mean_poly;
	sa = sa_bulk/(1.0 - w_ih);
	/**
	  -----------------------------------------------------------------------
	     Calculating the estimate of the derivative of func, dfunc_dw_Ih, to
	     be fed into Newton s Method.
	  -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing(sa,p);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa, NULL);
	gsw_pot_enthalpy_ice_freezing_first_derivatives(sa,p,&dpot_h_ihf_dsa, NULL);
	dfunc_dw_ih = gtc.gsw_cp0*ctf - h_pot_ihf - sa*
					  (gtc.gsw_cp0*ctf_sa + w_ih*dpot_h_ihf_dsa/(1.0 - w_ih));
	if (w_ih >= 0.0 && w_ih <= 0.20 && sa > 15.0
			&& sa < 60.0 && p <= 3000.0)
	{
		max_iterations = 1;
	}
	else if (w_ih >= 0.0 && w_ih <= 0.85 && sa > 0.0
				&& sa < 120.0 && p <= 3500.0)
	{
		max_iterations = 2;
	}
	else
	{
		max_iterations = 3;
	}
	for (iterations = 1; iterations <= max_iterations; iterations++)
	{
		if (iterations > 1)
		{
			/**On the first iteration ctf and h_pot_ihf are both known*/
			ctf = gsw_ct_freezing(sa,p,saturation_fraction);
			h_pot_ihf = gsw_pot_enthalpy_ice_freezing(sa,p);
		}
		/**This is the method, func, whose zero we seek ...*/
		func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
		w_ih_old = w_ih;
		w_ih = w_ih_old - func/dfunc_dw_ih;
		if (w_ih > 0.9)
		{
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih);
	} /** for */
	if (w_ih < 0.0)
	{
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
	}
	else
	{
		*sa_final = sa;
		*ct_final = gsw_ct_freezing(sa,p,saturation_fraction);
		*w_ih_final = w_ih;
	}
}

/**************************************************************************
==========================================================================
method: gsw_frazil_properties_potential_poly (sa_bulk, &
                             h_pot_bulk, p, sa_final, ct_final, w_ih_final)
==========================================================================
   Calculates the mass fraction of ice (mass of ice divided by mass of ice
   plus seawater), w_Ih_final, which results from given values of the bulk
   Absolute Salinity, SA_bulk, bulk potential enthalpy, h_pot_bulk,
   occuring at pressure p.  The final equilibrium values of Absolute
   Salinity, SA_final, and Conservative Temperature, CT_final, of the
   interstitial seawater phase are also returned.  This code assumes that
   there is no dissolved air in the seawater (that is, saturation_fraction
   is assumed to be zero thoughout the code).
   When the mass fraction w_Ih_final is calculated as being a positive
   value, the seawater-ice mixture is at thermodynamic equlibrium.
   This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
   is sufficiently large (i.e. sufficiently "warm") so that there is no ice
   present in the final state.  In this case the final state consists of
   only seawater rather than being an equlibrium mixture of seawater and
   ice which occurs when w_Ih_final is positive.  Note that when
   w_Ih_final = 0, the final seawater is not at the freezing temperature.
   Note that this code uses the polynomial forms of CT_freezing and
   pot_enthalpy_ice_freezing. This code is intended to be used in ocean
   models where the model prognostic variables are SA_bulk and h_pot_bulk.
   SA_bulk     :  bulk Absolute Salinity of the seawater and ice mixture[ g/kg ]
   h_pot_bulk  :  bulk potential enthalpy of the seawater and ice mixture[ J/kg ]
   p           :  sea pressure                                  [ dbar ]
                   ( i.e. absolute pressure - 10.1325 dbar )
   SA_final    :  Absolute Salinity of the seawater in the final state,
                  whether or not any ice is present.               [ g/kg ]
   CT_final    :  Conservative Temperature of the seawater in the the final
                  state, whether or not any ice is present.       [ deg C ]
   w_Ih_final  :  mass fraction of ice in the final seawater-ice mixture.
                  If this ice mass fraction is positive, the system is at
                  thermodynamic equilibrium.  If this ice mass fraction is
                  zero there is no ice in the final state which consists
                  only of seawater which is warmer than the freezing
                  temperature.                                   [unitless]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_frazil_properties_potential_poly(double sa_bulk,double h_pot_bulk,double p,double *sa_final,double *ct_final,double *w_ih_final)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int     iterations, max_iterations;
	double  ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly, dpot_h_ihf_dsa,
			  func, func0, h_pot_ihf, sa, w_ih_old, w_ih, x, xa, y, z;
	double  f01 = -9.041191886754806e-1,
			  f02 =  4.169608567309818e-2,
			  f03 = -9.325971761333677e-3,
			  f04 =  4.699055851002199e-2,
			  f05 = -3.086923404061666e-2,
			  f06 =  1.057761186019000e-2,
			  f07 = -7.349302346007727e-2,
			  f08 =  1.444842576424337e-1,
			  f09 = -1.408425967872030e-1,
			  f10 =  1.070981398567760e-1,
			  f11 = -1.768451760854797e-2,
			  f12 = -4.013688314067293e-1,
			  f13 =  7.209753205388577e-1,
			  f14 = -1.807444462285120e-1,
			  f15 =  1.362305015808993e-1,
			  f16 = -9.500974920072897e-1,
			  f17 =  1.192134856624248,
			  f18 = -9.191161283559850e-2,
			  f19 = -1.008594411490973,
			  f20 =  8.020279271484482e-1,
			  f21 = -3.930534388853466e-1,
			  f22 = -2.026853316399942e-2,
			  f23 = -2.722731069001690e-2,
			  f24 =  5.032098120548072e-2,
			  f25 = -2.354888890484222e-2,
			  f26 = -2.454090179215001e-2,
			  f27 =  4.125987229048937e-2,
			  f28 = -3.533404753585094e-2,
			  f29 =  3.766063025852511e-2,
			  f30 = -3.358409746243470e-2,
			  f31 = -2.242158862056258e-2,
			  f32 =  2.102254738058931e-2,
			  f33 = -3.048635435546108e-2,
			  f34 = -1.996293091714222e-2,
			  f35 =  2.577703068234217e-2,
			  f36 = -1.292053030649309e-2,
			  g01 =  3.332286683867741e5,
			  g02 =  1.416532517833479e4,
			  g03 = -1.021129089258645e4,
			  g04 =  2.356370992641009e4,
			  g05 = -8.483432350173174e3,
			  g06 =  2.279927781684362e4,
			  g07 =  1.506238790315354e4,
			  g08 =  4.194030718568807e3,
			  g09 = -3.146939594885272e5,
			  g10 = -7.549939721380912e4,
			  g11 =  2.790535212869292e6,
			  g12 =  1.078851928118102e5,
			  g13 = -1.062493860205067e7,
			  g14 =  2.082909703458225e7,
			  g15 = -2.046810820868635e7,
			  g16 =  8.039606992745191e6,
			  g17 = -2.023984705844567e4,
			  g18 =  2.871769638352535e4,
			  g19 = -1.444841553038544e4,
			  g20 =  2.261532522236573e4,
			  g21 = -2.090579366221046e4,
			  g22 = -1.128417003723530e4,
			  g23 =  3.222965226084112e3,
			  g24 = -1.226388046175992e4,
			  g25 =  1.506847628109789e4,
			  g26 = -4.584670946447444e4,
			  g27 =  1.596119496322347e4,
			  g28 = -6.338852410446789e4,
			  g29 =  8.951570926106525e4,
			  saturation_fraction = 0.0;
	/**
	    -----------------------------------------------------------------------
	      Finding func0.  This is the value of the function, func, that would
	      result in the output w_Ih_final being exactly zero.
	    -----------------------------------------------------------------------
	*/
	func0 = h_pot_bulk - gtc.gsw_cp0
			  *gsw_ct_freezing_poly(sa_bulk,p,saturation_fraction);
	/**
	    -----------------------------------------------------------------------
	      Setting the three outputs for data points that have func0 non-negative
	    -----------------------------------------------------------------------
	*/
	if (func0 >= 0.0)
	{
		/**
		      When func0 is zero or positive then the final answer will contain
		      no frazil ice; that is, it will be pure seawater that is warmer
		      han the freezing temperature.  If func0 >= 0 we do not need to go
		      through the modified Newton-Raphson procedure and we can simply
		      write down the answer, as in the following 4 lines of code.
		*/
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
		return;
	}
	/**
	    -----------------------------------------------------------------------
	    Begin finding the solution for data points that have func0 < 0, so that
	    the output will have a positive ice mass fraction w_Ih_final.
	    -----------------------------------------------------------------------
	*/
	/**Evalaute a polynomial for w_Ih in terms of SA_bulk, func0 and p*/
	x = sa_bulk*1e-2;
	y = func0/3e5;
	z = p*1e-4;
	w_ih = y*(f01 + x*(f02 + x*(f03 + x*(f04 + x*(f05 + f06*x))))
				 + y*(f07 + x*(f08 + x*(f09 + x*(f10 + f11*x))) + y*(f12 + x*(f13
						 + x*(f14 + f15*x)) + y*(f16 + x*(f17 + f18*x) + y*(f19 + f20*x
														 + f21*y)))) + z*(f22 + x*(f23 + x*(f24 + f25*x)) + y*(x*(f26
																 + f27*x) + y*(f28 + f29*x + f30*y)) + z*(f31 + x*(f32 + f33*x)
																		 + y*(f34 + f35*x + f36*y))));
	if (w_ih > 0.9)
	{
		/**
		    The ice mass fraction out of this code is restricted to be
		    less than 0.9.
		*/
		*sa_final = cppGSW_INVALID_VALUE;
		*ct_final = *sa_final;
		*w_ih_final = *sa_final;
		return;
	}
	/**
	     The initial guess at the absolute salinity of the interstitial
	     seawater
	*/
	sa = sa_bulk/(1.0 - w_ih);
	/**
	    -----------------------------------------------------------------------
	      Doing a Newton step with a separate polynomial estimate of the mean
	      derivative dfunc_dw_Ih_mean_poly.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing_poly(sa,p);
	func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
	xa = sa*1e-2;
	dfunc_dw_ih_mean_poly = g01 + xa*(g02 + xa*(g03 + xa*(g04 + g05*xa)))
									+ w_ih*(xa*(g06 + xa*(g07 + g08*xa)) + w_ih*(xa*(g09 + g10*xa)
											  + w_ih*xa*(g11 + g12*xa + w_ih*(g13 + w_ih*(g14 + w_ih*(g15
													  + g16*w_ih)))))) + z*(g17 + xa*(g18 + g19*xa) + w_ih*(g20
															  + w_ih*(g21 + g22*w_ih) + xa*(g23 + g24*xa*w_ih))
															  + z*(g25 + xa*(g26 + g27*xa) + w_ih*(g28 + g29*w_ih)));
	w_ih_old = w_ih;
	w_ih = w_ih_old - func/dfunc_dw_ih_mean_poly;
	sa = sa_bulk/(1.0 - w_ih);
	/**
	    -----------------------------------------------------------------------
	      Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
	      fed into Newton s Method.
	    -----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	h_pot_ihf = gsw_pot_enthalpy_ice_freezing_poly(sa,p);
	gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,&ctf_sa,
														NULL);
	gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(sa,p,
			&dpot_h_ihf_dsa, NULL);
	dfunc_dw_ih = gtc.gsw_cp0*ctf - h_pot_ihf -
					  sa*(gtc.gsw_cp0*ctf_sa + w_ih*dpot_h_ihf_dsa/(1.0 - w_ih));
	if (w_ih >= 0.0 && w_ih <= 0.20 && sa > 15.0
			&& sa < 60.0 && p <= 3000.0)
	{
		max_iterations = 1;
	}
	else if (w_ih >= 0.0 && w_ih <= 0.85 && sa > 0.0
				&& sa < 120.0 && p <= 3500.0)
	{
		max_iterations = 2;
	}
	else
	{
		max_iterations = 3;
	}
	for (iterations = 1; iterations <= max_iterations; iterations++)
	{
		if (iterations > 1)
		{
			/**On the first iteration ctf and h_pot_ihf are both known*/
			ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
			h_pot_ihf = gsw_pot_enthalpy_ice_freezing_poly(sa,p);
		}
		/**This is the function, func, whose zero we seek ...*/
		func = h_pot_bulk - (1.0 - w_ih)*gtc.gsw_cp0*ctf - w_ih*h_pot_ihf;
		w_ih_old = w_ih;
		w_ih = w_ih_old - func/dfunc_dw_ih;
		if (w_ih > 0.9)
		{
			*sa_final = cppGSW_INVALID_VALUE;
			*ct_final = *sa_final;
			*w_ih_final = *sa_final;
			return;
		}
		sa = sa_bulk/(1.0 - w_ih);
	} /** for */
	if (w_ih < 0.0)
	{
		*sa_final = sa_bulk;
		*ct_final = h_pot_bulk/gtc.gsw_cp0;
		*w_ih_final = 0.0;
	}
	else
	{
		*sa_final = sa;
		*ct_final = gsw_ct_freezing_poly(sa,p,saturation_fraction);
		*w_ih_final = w_ih;
	}
}

/**************************************************************************
==========================================================================
method: gsw_frazil_ratios_adiabatic (sa, p, w_ih, &
                              dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
==========================================================================
   Calculates the ratios of SA, CT and P changes when frazil ice forms or
   melts in response to an adiabatic change in pressure of a mixture of
   seawater and frazil ice crystals.
   Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than
   dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero
   whereas dCT/dSA would then be infinite.
   Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
   derivatives with the pressure measured in Pa not dbar.
   SA  :  Absolute Salinity of seawater                            [ g/kg ]
   p   :  sea pressure of seawater at which melting occurs         [ dbar ]
          ( i.e. absolute pressure - 10.1325d0 dbar )
   w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
            sum of the masses of ice and seawater.  That is, the mass of
            ice divided by the mass of the final mixed fluid.
            w_Ih must be between 0 and 1.                      [ unitless ]
   dSA_dCT_frazil :  the ratio of the changes in Absolute Salinity
                     to that of Conservative Temperature       [ g/(kg K) ]
   dSA_dP_frazil  :  the ratio of the changes in Absolute Salinity
                     to that of pressure (in Pa)              [ g/(kg Pa) ]
   dCT_dP_frazil  :  the ratio of the changes in Conservative Temperature
                     to that of pressure (in Pa)                   [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_frazil_ratios_adiabatic(double sa,double p,double w_ih,double *dsa_dct_frazil,double *dsa_dp_frazil,double *dct_dp_frazil)
{
	double bracket1, bracket2, cp_ih, gamma_ih, h, h_ih, part,
			 rec_bracket3, tf, wcp, h_hat_sa, h_hat_ct, tf_sa, tf_p,
			 ctf, ctf_sa, ctf_p, saturation_fraction = 0.0;
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	h = gsw_enthalpy_ct_exact(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gamma_ih = gsw_adiabatic_lapse_rate_ice(tf,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa,&tf_p);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction, &ctf_sa,&ctf_p);
	wcp = cp_ih*w_ih/(1.0 - w_ih);
	part = (tf_p - gamma_ih)/ctf_p;
	bracket1 = h_hat_ct + wcp*part;
	bracket2 = h - h_ih - sa*(h_hat_sa + wcp*(tf_sa - part*ctf_sa));
	rec_bracket3 = 1.0/(h - h_ih - sa*(h_hat_sa + h_hat_ct*ctf_sa + wcp*tf_sa));
	*dsa_dct_frazil = sa*(bracket1/bracket2);
	*dsa_dp_frazil = sa*ctf_p*bracket1*rec_bracket3;
	*dct_dp_frazil = ctf_p*bracket2*rec_bracket3;
}

/**************************************************************************
==========================================================================
method: gsw_frazil_ratios_adiabatic_poly (sa, p, w_ih, &
                              dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
==========================================================================
   Calculates the ratios of SA, CT and P changes when frazil ice forms or
   melts in response to an adiabatic change in pressure of a mixture of
   seawater and frazil ice crystals.
   Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than
   dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero
   whereas dCT/dSA would then be infinite.
   Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
   derivatives with the pressure measured in Pa not dbar.
   SA  :  Absolute Salinity of seawater                            [ g/kg ]
   p   :  sea pressure of seawater at which melting occurs         [ dbar ]
          ( i.e. absolute pressure - 10.1325d0 dbar )
   w_Ih  :  mass fraction of ice, that is the mass of ice divided by the
            sum of the masses of ice and seawater.  That is, the mass of
            ice divided by the mass of the final mixed fluid.
            w_Ih must be between 0 and 1.                      [ unitless ]
   dSA_dCT_frazil :  the ratio of the changes in Absolute Salinity
                     to that of Conservative Temperature       [ g/(kg K) ]
   dSA_dP_frazil  :  the ratio of the changes in Absolute Salinity
                     to that of pressure (in Pa)              [ g/(kg Pa) ]
   dCT_dP_frazil  :  the ratio of the changes in Conservative Temperature
                     to that of pressure (in Pa)                   [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_frazil_ratios_adiabatic_poly(double sa,double p,double w_ih,double *dsa_dct_frazil,double *dsa_dp_frazil,double *dct_dp_frazil)
{
	double bracket1, bracket2, cp_ih, gamma_ih, h, h_ih, part,
			 rec_bracket3, tf, wcp, h_hat_sa, h_hat_ct, tf_sa, tf_p,
			 ctf, ctf_sa, ctf_p, saturation_fraction = 0.0;
	tf = gsw_t_freezing_poly(sa,p,saturation_fraction);
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	h = gsw_enthalpy(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gamma_ih = gsw_adiabatic_lapse_rate_ice(tf,p);
	gsw_enthalpy_first_derivatives(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction, &tf_sa,&tf_p);
	gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction, &ctf_sa,&ctf_p);
	wcp = cp_ih*w_ih/(1.0 - w_ih);
	part = (tf_p - gamma_ih)/ctf_p;
	bracket1 = h_hat_ct + wcp*part;
	bracket2 = h - h_ih - sa*(h_hat_sa + wcp*(tf_sa - part*ctf_sa));
	rec_bracket3 = 1.0/(h - h_ih - sa*(h_hat_sa + h_hat_ct*ctf_sa + wcp*tf_sa));
	*dsa_dct_frazil = sa*(bracket1/bracket2);
	*dsa_dp_frazil = sa*ctf_p*bracket1*rec_bracket3;
	*dct_dp_frazil = ctf_p*bracket2*rec_bracket3;
}

/**************************************************************************
==========================================================================
method: gsw_gibbs(ns,nt,np,sa,t,p)
==========================================================================
 seawater specific Gibbs free energy and derivatives up to order 2
 ns     : order of s derivative
 nt     : order of t derivative
 np     : order of p derivative
 sa     : Absolute Salinity                               [g/kg]
 t      : temperature                                     [deg C]
 p      : sea pressure                                    [dbar]-1
 gsw_gibbs  : specific Gibbs energy or its derivative	   [J kg  ]
*************************************************************************/

double TeosBase::gsw_gibbs(int ns, int nt, int np, double sa, double t, double p)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double x2, x, y, z, return_value = cppGSW_INVALID_VALUE;
	x2	= gtc.gsw_sfac * sa;
	x	= sqrt(x2);
	y	= t * 0.025;
	z	= p * gtc.rec_db2pa;

   if (ns == 0  && nt == 0  && np == 0)
	{
		return_value = gsw_gibbs_g03plusg08_a(sa, z, y, x2, x);
	}
	else if (ns == 1  && nt == 0  && np == 0)
	{
		return_value = gsw_gibbs_g08_a(sa, z, y, x);
	}
	else if (ns == 0  && nt == 1  && np == 0)
	{
		return_value = gsw_gibbs_g03plusg08_b(sa, z, y, x2, x);
	}
	else if (ns == 0  && nt == 0  && np == 1)
	{
		return_value = gsw_gibbs_g03plusg08_c(z, y, x2, x);
	}
	else if (ns == 0  && nt == 2  && np == 0)
	{
		return_value	= gsw_gibbs_g03plusg08_d(z, y, x2, x);
	}
	else if (ns == 1  && nt == 0  && np == 1)
	{
		return_value = gsw_gibbs_g08_b(z, y, x);
	}
	else if (ns == 0  && nt == 1  && np == 1)
	{
		return_value = gsw_gibbs_g03plusg08_e(z, y, x2, x);
	}
	else if (ns == 1  && nt == 1  && np == 0)
	{
		return_value = gsw_gibbs_g08_c(sa, z, y, x);
	}
	else if (ns == 2  && nt == 0  && np == 0)
	{
		return_value = gsw_gibbs_g03plusg08_f(z, y, x2, x);
	}
	else if (ns == 0  && nt == 0  && np == 2)
	{
		return_value = gsw_gibbs_g03plusg08_g(z, y, x2, x);
	}


	return return_value;
}

/**************************************************************************
=========================================================================
method: gsw_gibbs_ice (nt, np, t, p)
=========================================================================
   Ice specific Gibbs energy and derivatives up to order 2.
   nt  :  order of t derivative                      [ integers 0, 1 or 2 ]
   np  :  order of p derivative                      [ integers 0, 1 or 2 ]
   t   :  in-situ temperature (ITS-90)                            [ deg C ]
   p   :  sea pressure                                             [ dbar ]
   gibbs_ice : Specific Gibbs energy of ice or its derivatives.
   The Gibbs energy (when nt = np = 0) has units of:     [ J/kg ]
   The temperature derivatives are output in units of:[ (J/kg) (K)^(-nt) ]
   The pressure derivatives are output in units of:[ (J/kg) (Pa)^(-np) ]
   The mixed derivatives are output in units of:[ (J/kg) (K)^(-nt) (Pa)^(-np) ]
   Note: The derivatives are taken with respect to pressure in Pa, not
     withstanding that the pressure input into this routine is in dbar.
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_gibbs_ice(int nt,int np,double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  dzi, g0, g0p, g0pp, sqrec_pt;
	complex<double> r2, r2p, r2pp, g, sqtau_t1, sqtau_t2, tau, tau_t1, tau_t2;
	double  s0 = -3.32733756492168e3;
	tau = (t + gtc.gsw_t0)*gic.rec_tt;
	dzi = gtc.db2pa*p*gic.rec_pt;
	if (nt == 0 && np == 0)
	{
		tau_t1 = tau/t1;
		sqtau_t1 = tau_t1*tau_t1;
		tau_t2 = tau/t2;
		sqtau_t2 = tau_t2*tau_t2;
		g0 = gic.g00 + dzi*(gic.g01 + dzi*(gic.g02 + dzi*(gic.g03 + gic.g04*dzi)));
		r2 = r20 + dzi*(r21 + r22*dzi);
		g = r1*(tau*log((1.0 + tau_t1)/(1.0 - tau_t1))
				  + t1*(log(1.0 - sqtau_t1) - sqtau_t1))
			 + r2*(tau*log((1.0 + tau_t2)/(1.0 - tau_t2))
					 + t2*(log(1.0 - sqtau_t2) - sqtau_t2));
		return real(g0 - gic.tt*(s0*tau - real(g)));  // EF: bug in original.
	}
	else if (nt == 1 && np == 0)
	{
		tau_t1 = tau/t1;
		tau_t2 = tau/t2;
		r2 = r20 + dzi*(r21 + r22*dzi);
		g = r1*(log((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
			 + r2*(log((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);
		return (-s0 + real(g));
	}
	else if (nt == 0 && np == 1)
	{
		tau_t2 = tau/t2;
		sqtau_t2 = tau_t2*tau_t2;
		g0p = gic.rec_pt*(gic.g01 + dzi*(2.0*gic.g02 + dzi*(3.0*gic.g03 + 4.0*gic.g04*dzi)));
		r2p = gic.rec_pt*(r21 + 2.0*r22*dzi);
		g = r2p*(tau*log((1.0 + tau_t2)/(1.0 - tau_t2))
					+ t2*(log(1.0 - sqtau_t2) - sqtau_t2));
		return (g0p + gic.tt*real(g));
	}
	else if (nt == 1 && np == 1)
	{
		tau_t2 = tau/t2;
		r2p = gic.rec_pt*(r21 + 2.0*r22*dzi) ;
		g = r2p*(log((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);
		return (real(g));
	}
	else if (nt == 2 && np == 0)
	{
		r2 = r20 + dzi*(r21 + r22*dzi);
		g = r1*(1.0/(t1 - tau) + 1.0/(t1 + tau) - 2.0/t1)
			 + r2*(1.0/(t2 - tau) + 1.0/(t2 + tau) - 2.0/t2);
		return (gic.rec_tt*real(g));
	}
	else if (nt == 0 && np == 2)
	{
		sqrec_pt = gic.rec_pt*gic.rec_pt;
		tau_t2 = tau/t2;
		sqtau_t2 = tau_t2*tau_t2;
		g0pp = sqrec_pt*(2.0*gic.g02 + dzi*(6.0*gic.g03 + 12.0*gic.g04*dzi));
		r2pp = 2.0*r22*sqrec_pt;
		g = r2pp*(tau*log((1.0 + tau_t2)/(1.0 - tau_t2))
					 + t2*(log(1.0 - sqtau_t2) - sqtau_t2));
		return (g0pp + gic.tt*real(g));
	}
	else
		return (cppGSW_INVALID_VALUE);
}

/**************************************************************************
=========================================================================
method: gsw_gibbs_ice_part_t (t, p)
=========================================================================
   part of the the first temperature derivative of Gibbs energy of ice
   that is the outout is gibbs_ice(1,0,t,p) + S0
   t   :  in-situ temperature (ITS-90)                            [ deg C ]
   p   :  sea pressure                                             [ dbar ]
   gibbs_ice_part_t = part of temperature derivative       [ J kg^-1 K^-1 ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_gibbs_ice_part_t(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  dzi, tau;
	complex<double> g, tau_t1, tau_t2, r2;
	tau = (t + gtc.gsw_t0)*gic.rec_tt;
	dzi = gtc.db2pa*p*gic.rec_pt;
	tau_t1 = tau/t1;
	tau_t2 = tau/t2;
	r2 = r20 + dzi*(r21 + r22*dzi);
	g = r1*(log((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
		 + r2*(log((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);
	return (real(g));
}

/**************************************************************************
==========================================================================
method: gsw_gibbs_ice_pt0 (pt0)
==========================================================================
   Part of the the first temperature derivative of Gibbs energy of ice
   that is the outout is "gibbs_ice(1,0,pt0,0) + s0"
   pt0  :  potential temperature with reference sea pressure of zero dbar
                                                                  [ deg C ]
   gsw_gibbs_ice_pt0 : part of temperature derivative     [ J kg^-1 K^-1 ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_gibbs_ice_pt0(double pt0)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  tau;
	complex<double> g, tau_t1, tau_t2;
	tau = (pt0 + gtc.gsw_t0)*gic.rec_tt;
	tau_t1 = tau/t1;
	tau_t2 = tau/t2;
	g = r1*(log((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
		 + r20*(log((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);
	return (real(g));
}

/**************************************************************************
==========================================================================
method: gsw_gibbs_ice_pt0_pt0 (pt0)
==========================================================================
   The second temperature derivative of Gibbs energy of ice at the
   potential temperature with reference sea pressure of zero dbar.  That is
   the output is gibbs_ice(2,0,pt0,0).
   pt0  :  potential temperature with reference sea pressure of zero dbar
                                                                  [ deg C ]
   gsw_gibbs_ice_pt0_pt0 : temperature second derivative at pt0
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_gibbs_ice_pt0_pt0(double pt0)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  tau;
	complex<double> g;
	tau = (pt0 + gtc.gsw_t0)*gic.rec_tt;
	g = r1*(1.0/(t1 - tau) + 1.0/(t1 + tau) - 2.0/t1)
		 + r20*(1.0/(t2 - tau) + 1.0/(t2 + tau) - 2.0/t2);
	return (gic.rec_tt*real(g));
}

/**************************************************************************
==========================================================================
method: gsw_gibbs_pt0_pt0(sa,pt0)
==========================================================================
 gibbs_tt at (sa,pt,0)
 sa     : Absolute Salinity                            [g/kg]
 pt0    : potential temperature                        [deg C]
 gibbs_pt0_pt0 : gibbs_tt at (sa,pt,0)
*************************************************************************/

double TeosBase::gsw_gibbs_pt0_pt0(double sa, double pt0)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double x2, x, y, g03, g08;
	x2	= gtc.gsw_sfac*sa;
	x = sqrt(x2);
	y = pt0*0.025;
	g03 = -24715.571866078 +
			y*(4420.4472249096725 +
				y*(-1778.231237203896 +
					y*(1160.5182516851419 +
						y*(-569.531539542516 + y*128.13429152494615))));
	g08 = x2*(1760.062705994408 + x*(-86.1329351956084 +
												x*(-137.1145018408982 + y*(296.20061691375236 +
														y*(-205.67709290374563 + 49.9394019139016*y))) +
												y*(-60.136422517125 + y*10.50720794170734)) +
				 y*(-1351.605895580406 + y*(1097.1125373015109 +
													 y*(-433.20648175062206 + 63.905091254154904*y))));
	return ((g03 + g08)*0.000625);
}

/**************************************************************************
==========================================================================
method: gsw_grav(lat,p)
==========================================================================
  Calculates acceleration due to gravity as a function of latitude and as
   a function of pressure in the ocean.
  lat  :  latitude in decimal degress north                [ -90 ... +90 ]
  p    :  sea pressure                                     [ dbar ]
  grav :  gravitational acceleration               [ m s^-2 ]
*************************************************************************/

double TeosBase::gsw_grav(double lat, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double x, sin2, gs, z;
	x = sin(lat*gtc.deg2rad);  /** convert to radians */
	sin2 = x*x;
	gs = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);
	z = gsw_z_from_p(p,lat,0.0,0.0);
	return (gs*(1.0 - gtc.gamma*z));    /** z is the height corresponding to p.
                                           Note. In the ocean z is negative. */
}

/**************************************************************************
==========================================================================
method:  gsw_hill_ratio_at_sp2(t)
==========================================================================
   Calculates the Hill ratio, which is the adjustment needed to apply for
   Practical Salinities smaller than 2.  This ratio is defined at a
   Practical Salinity = 2 and in-situ temperature, t using PSS-78. The Hill
   ratio is the ratio of 2 to the output of the Hill et al. (1986) formula
   for Practical Salinity at the conductivity ratio, Rt, at which Practical
   Salinity on the PSS-78 scale is exactly 2.
   t : in-situ temperature (ITS-90)              [deg C]
   hill_ratio_at_sp2 : Hill ratio                [dimensionless]
*************************************************************************/

double TeosBase::gsw_hill_ratio_at_sp2(double t)
{
	/** for GSW_SP_COEFFICIENTS use gspc */
	double  g0 = 2.641463563366498e-1, g1 = 2.007883247811176e-4,
			  g2 = -4.107694432853053e-6, g3 = 8.401670882091225e-8,
			  g4 = -1.711392021989210e-9, g5 = 3.374193893377380e-11,
			  g6 = -5.923731174730784e-13, g7 = 8.057771569962299e-15,
			  g8 = -7.054313817447962e-17, g9 = 2.859992717347235e-19,
			  sp2 = 2.0;
	double t68, ft68, rtx0, dsp_drtx, sp_est, rtx, rtxm, x, part1, part2;
	double sqrty, sp_hill_raw_at_sp2;
	t68 = t*1.00024;
	ft68 = (t68 - 15.0)/(1.0 + gspc.k*(t68 - 15.0));
	/**
	------------------------------------------------------------------------
	Find the initial estimates of Rtx (Rtx0) and of the derivative dSP_dRtx
	at SP = 2.
	------------------------------------------------------------------------
	*/
	rtx0 = g0 + t68*(g1 + t68*(g2 + t68*(g3 + t68*(g4 + t68*(g5
													 + t68*(g6 + t68*(g7 + t68*(g8 + t68*g9))))))));
	dsp_drtx = gspc.a1
				  + (2*gspc.a2 + (3*gspc.a3 + (4*gspc.a4
														 + 5*gspc.a5*rtx0)*rtx0)*rtx0)*rtx0
				  + ft68*(gspc.b1 + (2*gspc.b2 + (3*gspc.b3
											+ (4*gspc.b4 + 5*gspc.b5*rtx0)*rtx0)*rtx0)*rtx0);
	/**
	-------------------------------------------------------------------------
	  Begin a single modified Newton-Raphson iteration to find Rt at SP = 2.
	-------------------------------------------------------------------------
	*/
	sp_est = gspc.a0
				+ (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4
												 + gspc.a5*rtx0)*rtx0)*rtx0)*rtx0)*rtx0
				+ ft68*(gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3
														+ (gspc.b4 + gspc.b5*rtx0)*rtx0)*rtx0)* rtx0)*rtx0);
	rtx = rtx0 - (sp_est - sp2)/dsp_drtx;
	rtxm = 0.5*(rtx + rtx0);
	dsp_drtx = gspc.a1
				  + (2*gspc.a2 + (3*gspc.a3 + (4*gspc.a4
														 + 5*gspc.a5*rtxm)*rtxm)*rtxm)*rtxm
				  + ft68*(gspc.b1 + (2*gspc.b2 + (3*gspc.b3
											+ (4*gspc.b4 + 5*gspc.b5*rtxm)* rtxm)*rtxm)*rtxm);
	rtx = rtx0 - (sp_est - sp2)/dsp_drtx;
	/**
	  This is the end of one full iteration of the modified Newton-Raphson
	  iterative equation solver. The error in Rtx at this point is equivalent
	  to an error in SP of 9e-16 psu.
	*/
	x = 400.0*rtx*rtx;
	sqrty = 10.0*rtx;
	part1 = 1.0 + x*(1.5 + x);
	part2 = 1.0 + sqrty*(1.0 + sqrty*(1.0 + sqrty));
	sp_hill_raw_at_sp2 = sp2 - gspc.a0/part1 - gspc.b0*ft68/part2;
	return (2.0/sp_hill_raw_at_sp2);
}

/**************************************************************************
==========================================================================
method: gsw_internal_energy(sa,ct,p)
==========================================================================
   Calculates internal energy of seawater.
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
  internal_energy  :  internal_energy of seawater          [J/kg]
*************************************************************************/

double TeosBase::gsw_internal_energy(double sa,double ct,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (gsw_enthalpy(sa,ct,p) - (gtc.gsw_p0 + gtc.db2pa*p)
			  *gsw_specvol(sa,ct,p));
}

/**************************************************************************
==========================================================================
method: gsw_kappa(sa,ct,p)
==========================================================================
   Calculates isentropic compressibility of seawater.  This method
   has inputs of Absolute Salinity and Conservative Temperature.  This
   method uses the computationally-efficient expression for
   specific volume in terms of SA, CT and p (Roquet et al., 2014).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
  kappa  :  isentropic compressibility                     [1.0/Pa]
*************************************************************************/

double TeosBase::gsw_kappa(double sa, double ct, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double v, v_p, xs, ys, z;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v = gsw_gsvco_v(xs, ys, z);
	v_p = gsw_gsvco_c(xs, ys, z);
	return (-1e-8*v_p/v);
}

/**************************************************************************
==========================================================================
method: gsw_kappa_t_exact(sa,t,p)
==========================================================================
  isentropic compressibility of seawater
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  gsw_kappa_t_exact : isentropic compressibility           [1/Pa]
*************************************************************************/

double TeosBase::gsw_kappa_t_exact(double sa,double t,double p)
{
	int n0=0, n1=1, n2=2;
	double g_tt, g_tp;
	g_tt = gsw_gibbs(n0,n2,n0,sa,t,p);
	g_tp = gsw_gibbs(n0,n1,n1,sa,t,p);
	return ((g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p)) /
			  (gsw_gibbs(n0,n0,n1,sa,t,p)*g_tt));
}

/**************************************************************************
--------------------------------------------------------------------------
  isobaric melting enthalpy and isobaric evaporation enthalpy
--------------------------------------------------------------------------
==========================================================================
method: gsw_latentheat_melting(sa,p)
==========================================================================
  Calculates latent heat, or enthalpy, of melting.
  sa     : Absolute Salinity                               [g/kg]
  p      : sea pressure                                    [dbar]
  latentheat_melting : latent heat of melting              [kg/m^3]
*************************************************************************/

double TeosBase::gsw_latentheat_melting(double sa,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double tf = gsw_t_freezing(sa,p,0.0);
	return (1000.0*(gsw_chem_potential_water_t_exact(sa,tf,p)
						 - (gtc.gsw_t0 + tf)*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p))
			  - gsw_enthalpy_ice(tf,p));
}

/**************************************************************************
==========================================================================
method: gsw_linear_interp_sa_ct (sa, ct, p, p_i, sa_i, ct_i)
==========================================================================
  This method interpolates the cast with respect to the interpolating
  variable p. This method finds the values of SA, CT at p_i on this cast.
  VERSION NUMBER: 3.05 (27th January 2015)
  This method was adapted from Matlab s interp1q.
==========================================================================
*************************************************************************/

void TeosBase::gsw_linear_interp_sa_ct(double *sa,double *ct,double *p,int np,double *p_i,int npi,double *sa_i,double *ct_i)
{
	char *in_rng=NULL;
	int *j=NULL, *k=NULL, *r=NULL, *jrev=NULL, *ki=NULL;
	int imax_p, imin_p, i, n, m, ii;
	double *xi=NULL, *xxi=NULL, u, max_p, min_p;
	min_p = max_p = p[0];
	imin_p = imax_p = 0;
	for (i=1; i<np; i++)
	{
		if (p[i] < min_p)
		{
			min_p = p[i];
			imin_p = i;
		}
		else if (p[i] > max_p)
		{
			max_p = p[i];
			imax_p = i;
		}
	}
	in_rng = new char[npi];
	memset(in_rng, 0, npi*sizeof (char));
	for (i=n=0; i<npi; i++)
	{
		if (p_i[i] <= min_p)
		{
			sa_i[i] = sa[imin_p];
			ct_i[i] = ct[imin_p];
		}
		else if (p_i[i] >= max_p)
		{
			sa_i[i] = sa[imax_p];
			ct_i[i] = ct[imax_p];
		}
		else
		{
			in_rng[i] = 1;
			n++;
		}
	}
	if (n==0)
	{
		if (in_rng) delete []in_rng;
		return;
	}
	xi = new double[n];
	k = new int[3*n];
	ki = k+n;
	r = ki+n;
	m = np + n;
	xxi = new double[m];
	j = new int[2*m];
	jrev = j+m;
	ii = 0;
	for (i = 0; i<npi; i++)
	{
		if (in_rng[i])
		{
			xi[ii] = p_i[i];
			ki[ii] = i;
			ii++;
		}
	}
	if (in_rng)
	{
		delete []in_rng;
		in_rng = NULL;
	}
	/**
	  Note that the following operations on the index
	  vectors jrev and r depend on the sort utility
	  gsw_util_sort_dbl() consistently ordering the
	  sorting indexes either in ascending or descending
	  sequence for replicate values in the real  vector .
	*/
	gsw_util_sort_dbl(xi, n, k);
	for (i = 0; i<np; i++)
		xxi[i] = p[i];
	for (i = 0; i<n; i++)
		xxi[np+i] = xi[k[i]];
	gsw_util_sort_dbl(xxi, m, j);
	for (i = 0; i<m; i++)
		jrev[j[i]] = i;
	for (i = 0; i<n; i++)
		r[k[i]] = jrev[np+i]-i-1;
	for (i = 0; i<n; i++)
	{
		u = (xi[i]-p[r[i]])/(p[r[i]+1]-p[r[i]]);
		sa_i[ki[i]] = sa[r[i]] + (sa[r[i]+1]-sa[r[i]])*u;
		ct_i[ki[i]] = ct[r[i]] + (ct[r[i]+1]-ct[r[i]])*u;
	}
	if (xi) delete []xi;
	if (k) delete []k;
	if (xxi) delete []xxi;
	if (j) delete []j;
}

/**************************************************************************
==========================================================================
method: gsw_o2sol(sa, ct, p, lon, lat)
==========================================================================
  Calculates the oxygen concentration expected at equilibrium with air at
  an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
  saturated water vapor.  This method uses the solubility coefficients
  derived from the data of Benson and Krause (1984), as fitted by Garcia
  and Gordon (1992, 1993).
  Note that this algorithm has not been approved by IOC and is not work
  from SCOR/IAPSO Working Group 127. It is included in the GSW
  Oceanographic Toolbox as it seems to be oceanographic best practice.
  SA  :  Absolute Salinity of seawater                           [ g/kg ]
  CT  :  Conservative Temperature of seawater (ITS-90)           [ deg C ]
  p   :  sea pressure at which the melting occurs                [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
  lat : latitude                                                 [deg]
  lon : longitude                                                [deg]
  gsw_o2sol : solubility of oxygen in micro-moles per kg          [umol/kg]
*************************************************************************/

double TeosBase::gsw_o2sol(double sa, double ct, double p, double lon, double lat)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double sp, pt, pt68, x, y, o2sol,
			 a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, c0;
	sp = gsw_sp_from_sa(sa, p, lon, lat);
	x = sp;
	pt = gsw_pt_from_ct(sa, ct);
	pt68 = pt*1.00024;
	y = log((298.15 - pt68)/(gtc.gsw_t0 + pt68));
	a0 =  5.80871;
	a1 =  3.20291;
	a2 =  4.17887;
	a3 =  5.10006;
	a4 = -9.86643e-2;
	a5 =  3.80369;
	b0 = -7.01577e-3;
	b1 = -7.70028e-3;
	b2 = -1.13864e-2;
	b3 = -9.51519e-3;
	c0 = -2.75915e-7;
	o2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y))))
					+ x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x));
	return o2sol;
}

/**************************************************************************
==========================================================================
method: gsw_o2sol_sp_pt(sp, pt)
==========================================================================
  Calculates the oxygen concentration expected at equilibrium with air at
  an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
  saturated water vapor.  This method uses the solubility coefficients
  derived from the data of Benson and Krause (1984), as fitted by Garcia
  and Gordon (1992, 1993).
  Note that this algorithm has not been approved by IOC and is not work
  from SCOR/IAPSO Working Group 127. It is included in the GSW
  Oceanographic Toolbox as it seems to be oceanographic best practice.
  SP  :  Practical Salinity  (PSS-78)                         [ unitless ]
  pt  :  potential temperature (ITS-90) referenced               [ dbar ]
          to one standard atmosphere (0 dbar).
  gsw_o2sol_sp_pt : solubility of oxygen in micro-moles per kg     [umol/kg]
*************************************************************************/

double TeosBase::gsw_o2sol_sp_pt(double sp, double pt)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double pt68, x, y, o2sol,
			 a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, c0;
	x = sp;
	pt68 = pt*1.00024;
	y = log((298.15 - pt68)/(gtc.gsw_t0 + pt68));
	a0 =  5.80871;
	a1 =  3.20291;
	a2 =  4.17887;
	a3 =  5.10006;
	a4 = -9.86643e-2;
	a5 =  3.80369;
	b0 = -7.01577e-3;
	b1 = -7.70028e-3;
	b2 = -1.13864e-2;
	b3 = -9.51519e-3;
	c0 = -2.75915e-7;
	o2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y))))
					+ x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x));
	return o2sol;
}

/**************************************************************************
==========================================================================
method: gsw_p_from_z(z,lat, double geo_strf_dyn_height,
		               double sea_surface_geopotential)
==========================================================================
 Calculates the pressure p from height z
 z      : height                                          [m]
 lat    : latitude                                        [deg]
 gsw_p_from_z : pressure                                  [dbar]
*************************************************************************/

double TeosBase::gsw_p_from_z(double z, double lat,double geo_strf_dyn_height,double sea_surface_geopotential)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double sinlat, sin2, gs, c1, p, df_dp, f, p_old, p_mid;
	if (z > 5) return cppGSW_INVALID_VALUE;
	sinlat = sin(lat*gtc.deg2rad);
	sin2 = sinlat*sinlat;
	gs = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);
	/** get the first estimate of p from Saunders (1981) */
	c1 =  5.25e-3*sin2 + 5.92e-3;
	p  = -2.0*z/((1-c1) + sqrt((1-c1)*(1-c1) + 8.84e-6*z)) ;
	/** end of the first estimate from Saunders (1981) */
	df_dp = gtc.db2pa*gsw_specvol_sso_0(p); /** initial value of the derivative of f */
	f = gsw_enthalpy_sso_0(p) + gs*(z - 0.5*gtc.gamma*(z*z))
		 - (geo_strf_dyn_height + sea_surface_geopotential);
	p_old = p;
	p = p_old - f/df_dp;
	p_mid = 0.5*(p + p_old);
	df_dp = gtc.db2pa*gsw_specvol_sso_0(p_mid);
	p = p_old - f/df_dp;
	return p;
}

/**************************************************************************
=======
method: gsw_pchip_derivs(*x,*y,n,*d)
=======
    Calculation of the derivatives is the key to the shape-preservation.
    There are other algorithms that could be used, but this appears to be
    the simplest, and adequate for our purposes.
    At minimal computational cost, we include here a check for increasing x.
    Returns 0 on success, 1 if x is not strictly increasing.
****************************************************************************/

int TeosBase::gsw_pchip_derivs(double *x,double *y,int n,double *d)
{
	double mm, mp;      /** slopes bracketing x */
	double hm, hp;      /** bracketing delta-x values */
	int smm, smp;       /** slope signs */
	double w1, w2;
	int i;
	if (n == 2)
	{
		d[0] = d[1] = (y[1] - y[0]) / (x[1] - x[0]);
		return 0;
	}
	hm = x[1] - x[0];
	hp = x[2] - x[1];
	mm = (y[1] - y[0]) / hm;
	mp = (y[2] - y[1]) / hp;
	d[0] = gsw_pchip_edge_case(hm, hp, mm, mp);
	smm = gsw_sgn(mm);
	smp = gsw_sgn(mp);
	for (i=1; i<n-1; i++)
	{
		if (hm <= 0)
		{
			return 1;
		}
		/** change of sign, or either slope is zero */
		if ((smm != smp) || mp == 0 || mm == 0)
		{
			d[i] = 0.0;
		}
		else
		{
			w1 = 2*hp + hm;
			w2 = hp + 2*hm;
			d[i] = (w1 + w2) / (w1/mm + w2/mp);
		}
		if (i < n-2)
		{
			hm = hp;
			hp = x[i+2] - x[i+1];
			mm = mp;
			mp = (y[i+2] - y[i+1]) / hp;
			smm = smp;
			smp = gsw_sgn(mp);
		}
	}
	if (hp <= 0)
	{
		return 1;
	}
	d[n-1] = gsw_pchip_edge_case(hp, hm, mp, mm);
	return 0;
}

/**************************************************************************
=======
method: gsw_pchip_edge_case(h0,h1,m0,m1)
=======
****************************************************************************/

double TeosBase::gsw_pchip_edge_case(double h0,double h1,double m0,double m1)
{
	double d;
	int mask, mask2;
	d = ((2*h0 + h1)*m0 - h0*m1) / (h0 + h1);
	mask = gsw_sgn(d) != gsw_sgn(m0);
	mask2 = (gsw_sgn(m0) != gsw_sgn(m1)) && (fabs(d) > 3.0*fabs(m0));
	if (mask)
	{
		return 0.0;
	}
	if (!mask && mask2)
	{
		return 3.0*m0;
	}
	return d;
}

/**************************************************************************
==========================================================================
method: gsw_pot_enthalpy_from_pt_ice (pt0_ice)
==========================================================================
   Calculates the potential enthalpy of ice from potential temperature of
   ice (whose reference sea pressure is zero dbar).
   pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
   gsw_pot_enthalpy_ice  =  potential enthalpy of ice              [ J/kg ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_pot_enthalpy_from_pt_ice(double pt0_ice)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_GIBBS_ICE_COEFFICIENTS use gic */
	double  tau;
	complex<double> h0_part, sqtau_t1, sqtau_t2;
	tau = (pt0_ice + gtc.gsw_t0)*gic.rec_tt;
	sqtau_t1 = (tau/t1)*(tau/t1);
	sqtau_t2 = (tau/t2)*(tau/t2);
	h0_part = r1*t1*(log(1.0 - sqtau_t1) + sqtau_t1)
				 + r20*t2*(log(1.0 - sqtau_t2) + sqtau_t2);
	return (gic.g00 + gic.tt*real(h0_part));
}

/**************************************************************************
==========================================================================
method: gsw_pot_enthalpy_ice_freezing (sa, p)
==========================================================================
   Calculates the potential enthalpy of ice at which seawater freezes.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   pot_enthalpy_ice_freezing : potential enthalpy of ice at freezing
                               of seawater                        [ deg C ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_pot_enthalpy_ice_freezing(double sa,double p)
{
	double pt0_ice, t_freezing;
	t_freezing = gsw_t_freezing(sa,p,0.0);
	pt0_ice = gsw_pt0_from_t_ice(t_freezing,p);
	return (gsw_pot_enthalpy_from_pt_ice(pt0_ice));
}

/**************************************************************************
==========================================================================
method: gsw_pot_enthalpy_ice_freezing_first_derivatives (sa, &
              p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
==========================================================================
   Calculates the first derivatives of the potential enthalpy of ice at
   which seawater freezes, with respect to Absolute Salinity SA and
   pressure P (in Pa).
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   pot_enthalpy_ice_freezing_SA : the derivative of the potential enthalpy
                   of ice at freezing (ITS-90) with respect to Absolute
                   salinity at fixed pressure  [ K/(g/kg) ] i.e. [ K kg/g ]
   pot_enthalpy_ice_freezing_P  : the derivative of the potential enthalpy
                   of ice at freezing (ITS-90) with respect to pressure
                   (in Pa) at fixed Absolute Salinity              [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa,double p,double *pot_enthalpy_ice_freezing_sa,double *pot_enthalpy_ice_freezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  cp_ihf, pt_icef, ratio_temp, tf, tf_p, tf_sa;
	double  saturation_fraction = 0.0;
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	pt_icef = gsw_pt0_from_t_ice(tf,p);
	ratio_temp = (gtc.gsw_t0 + pt_icef)/(gtc.gsw_t0 + tf);
	cp_ihf = gsw_cp_ice(tf,p);
	if ((pot_enthalpy_ice_freezing_sa != NULL) &&
			(pot_enthalpy_ice_freezing_p != NULL))
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa,&tf_p);
	}
	else if (pot_enthalpy_ice_freezing_sa != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													&tf_sa, NULL);
	}
	else if (pot_enthalpy_ice_freezing_p != NULL)
	{
		gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
													NULL,&tf_p);
	}
	if (pot_enthalpy_ice_freezing_sa != NULL)
		*pot_enthalpy_ice_freezing_sa = ratio_temp*cp_ihf*tf_sa;
	if (pot_enthalpy_ice_freezing_p != NULL)
		*pot_enthalpy_ice_freezing_p = ratio_temp*cp_ihf*tf_p
												 - (gtc.gsw_t0 + pt_icef)*gsw_gibbs_ice(1,1,tf,p);
}

/**************************************************************************
==========================================================================
method: gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(&
         sa, p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
==========================================================================
   Calculates the first derivatives of the potential enthalpy of ice Ih at
   which ice melts into seawater with Absolute Salinity SA and at pressure
   p.  This code uses the comptationally efficient polynomial fit of the
   freezing potential enthalpy of ice Ih (McDougall et al., 2015).
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   pot_enthalpy_ice_freezing_SA : the derivative of the potential enthalpy
                 of ice at freezing (ITS-90) with respect to Absolute
                 salinity at fixed pressure  [ (J/kg)/(g/kg) ] i.e. [ J/g ]
   pot_enthalpy_ice_freezing_P  : the derivative of the potential enthalpy
                 of ice at freezing (ITS-90) with respect to pressure
                 (in Pa) at fixed Absolute Salinity           [ (J/kg)/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(double sa,double p,double *pot_enthalpy_ice_freezing_sa,double *pot_enthalpy_ice_freezing_p)
{
	double  p_r, sa_r, x,
			  d1 =  -1.249490228128056e4,
			  d2 =   1.336783910789822e4,
			  d3 =  -4.811989517774642e4,
			  d4 =   8.044864276240987e4,
			  d5 =  -7.124452125071862e4,
			  d6 =   2.280706828014839e4,
			  d7 =   0.315423710959628e3,
			  d8 =  -3.592775732074710e2,
			  d9 =   1.644828513129230e3,
			  d10 = -4.809640968940840e3,
			  d11 =  2.901071777977272e3,
			  d12 = -9.218459682855746e2,
			  d13 =  0.379377450285737e3,
			  d14 = -2.672164989849465e3,
			  d15 =  5.044317489422632e3,
			  d16 = -2.631711865886377e3,
			  d17 = -0.160245473297112e3,
			  d18 =  4.029061696035465e2,
			  d19 = -3.682950019675760e2,
			  f1 =  -2.034535061416256e4,
			  f2 =   0.315423710959628e3,
			  f3 =  -0.239518382138314e3,
			  f4 =   0.822414256564615e3,
			  f5 =  -1.923856387576336e3,
			  f6 =   0.967023925992424e3,
			  f7 =  -0.263384562367307e3,
			  f8 =  -5.051613740291480e3,
			  f9 =   7.587549005714740e2,
			  f10 = -3.562886653132620e3,
			  f11 =  5.044317489422632e3,
			  f12 = -2.105369492709102e3,
			  f13 =  6.387082316647800e2,
			  f14 = -4.807364198913360e2,
			  f15 =  8.058123392070929e2,
			  f16 = -5.524425029513641e2;
	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;
	if (pot_enthalpy_ice_freezing_sa != NULL)
		*pot_enthalpy_ice_freezing_sa =
			(d1 + x*(d2  + x*(d3  + x*(d4  + x*(d5  + d6*x))))
			 + p_r*(d7 + x*(d8 + x*(d9 + x*(d10 + x*(d11 + d12*x))))
					  + p_r*(d13 + x*(d14 + x*(d15 + d16*x))
								+ p_r*(d17 + x*(d18 + d19*x)))))*1e-2;
	if (pot_enthalpy_ice_freezing_p != NULL)
		*pot_enthalpy_ice_freezing_p =
			(f1 + sa_r*(f2 + x*(f3 + x*(f4 + x*(f5 + x*(f6 + f7*x)))))
			 + p_r*(f8 + sa_r*(f9 + x*(f10 + x*(f11 + f12*x)))
					  + p_r*(f13 + sa_r*(f14 + x*(f15 + f16*x)))))*1e-8;
}

/**************************************************************************
==========================================================================
method: gsw_pot_enthalpy_ice_freezing_poly (sa, p)
==========================================================================
   Calculates the potential enthalpy of ice at which seawater freezes.
   The error of this fit ranges between -2.5 and 1 J/kg with an rms of
   1.07,between SA of 0 and 120 g/kg and p between 0 and 10,000 dbar (the
   error in the fit is between -0.7 and 0.7 with an rms of
   0.3, between SA of 0 and 120 g/kg and p between 0 and 5,000 dbar) when
   compared with the potential enthalpy calculated from the exact in-situ
   freezing temperature which is found by a Newton-Raphson iteration of the
   equality of the chemical potentials of water in seawater and in ice.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   pot_enthalpy_ice_freezing : potential enthalpy of ice at freezing
                               of seawater                         [ J/kg ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_pot_enthalpy_ice_freezing_poly(double sa,double p)
{
	double p_r, sa_r, x,
			 c0  = -3.333548730778702e5,
			 c1  = -1.249490228128056e4,
			 c2  =  0.891189273859881e4,
			 c3  = -2.405994758887321e4,
			 c4  =  3.217945710496395e4,
			 c5  = -2.374817375023954e4,
			 c6  =  0.651630522289954e4,
			 c7  = -2.034535061416256e4,
			 c8  = -0.252580687014574e4,
			 c9  =  0.021290274388826e4,
			 c10 =  0.315423710959628e3,
			 c11 = -0.239518382138314e3,
			 c12 =  0.379377450285737e3,
			 c13 =  0.822414256564615e3,
			 c14 = -1.781443326566310e3,
			 c15 = -0.160245473297112e3,
			 c16 = -1.923856387576336e3,
			 c17 =  2.522158744711316e3,
			 c18 =  0.268604113069031e3,
			 c19 =  0.967023925992424e3,
			 c20 = -1.052684746354551e3,
			 c21 = -0.184147500983788e3,
			 c22 = -0.263384562367307e3;
	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;
	return (c0 + sa_r*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + c6*x)))))
			  + p_r*(c7 + p_r*(c8 + c9*p_r)) + sa_r*p_r*(c10 + p_r*(c12
					  + p_r*(c15 + c21*sa_r)) + sa_r*(c13 + c17*p_r + c19*sa_r)
					  + x*(c11 + p_r*(c14 + c18*p_r) + sa_r*(c16 + c20*p_r + c22*sa_r))));
}

/**************************************************************************
==========================================================================
method: gsw_pt_from_ct(sa,ct)
==========================================================================
 potential temperature of seawater from conservative temperature
 sa     : Absolute Salinity                               [g/kg]
 ct     : Conservative Temperature                        [deg C]
 p      : sea pressure                                    [dbar]
 gsw_pt_from_ct : potential temperature with              [deg C]
                  reference pressure of  0 dbar
*************************************************************************/

double TeosBase::gsw_pt_from_ct(double sa, double ct)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double a5ct, b3ct, ct_factor, pt_num, pt_recden, ct_diff;
	double pt, pt_old, ptm, dpt_dct, s1;
	double a0 = -1.446013646344788e-2,
			 a1	= -3.305308995852924e-3,
			 a2	=  1.062415929128982e-4,
			 a3	=  9.477566673794488e-1,
			 a4	=  2.166591947736613e-3,
			 a5	=  3.828842955039902e-3,
			 b0	=  1.000000000000000e0,
			 b1	=  6.506097115635800e-4,
			 b2	=  3.830289486850898e-3,
			 b3	=  1.247811760368034e-6;
	s1	= sa/gtc.gsw_ups;
	a5ct	= a5*ct;
	b3ct	= b3*ct;
	ct_factor	= (a3 + a4*s1 + a5ct);
	pt_num = a0 + s1*(a1 + a2*s1) + ct*ct_factor;
	pt_recden	= 1.0/(b0 + b1*s1 + ct*(b2 + b3ct));
	pt = pt_num*pt_recden;
	dpt_dct	= pt_recden*(ct_factor + a5ct - (b2 + b3ct + b3ct)*pt);
	/**
	 Start the 1.5 iterations through the modified Newton-Raphson
	 iterative method.
	*/
	ct_diff	= gsw_ct_from_pt(sa,pt) - ct;
	pt_old	= pt;
	pt	= pt_old - ct_diff*dpt_dct;
	ptm	= 0.5*(pt + pt_old);
	dpt_dct	= -gtc.gsw_cp0/((ptm + gtc.gsw_t0)*gsw_gibbs_pt0_pt0(sa,ptm));
	pt	= pt_old - ct_diff*dpt_dct;
	ct_diff	= gsw_ct_from_pt(sa,pt) - ct;
	pt_old	= pt;
	return (pt_old - ct_diff*dpt_dct);
}

/**************************************************************************
==========================================================================
method: gsw_pt_from_t(sa,t,p,p_ref)
==========================================================================
 Calculates potential temperature of seawater from in-situ temperature
 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]
 p_ref  : reference sea pressure                          [dbar]
 gsw_pt_from_t : potential temperature                    [deg C]
*************************************************************************/

double TeosBase::gsw_pt_from_t(double sa, double t, double p, double p_ref)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	int n0=0, n2=2, no_iter;
	double s1, pt, ptm, pt_old, dentropy, dentropy_dt,
			 true_entropy_part;
	s1	= sa/gtc.gsw_ups;
	pt	= t+(p-p_ref)*( 8.65483913395442e-6  -
							 s1 *  1.41636299744881e-6  -
							 (p+p_ref) *  7.38286467135737e-9  +
							 t  *(-8.38241357039698e-6  +
									s1 *  2.83933368585534e-8  +
									t  *  1.77803965218656e-8  +
									(p+p_ref) *  1.71155619208233e-10));
	dentropy_dt	= gtc.gsw_cp0/((gtc.gsw_t0 + pt)*(1.0-0.05*(1.0 - sa/gtc.gsw_sso)));
	true_entropy_part	= gsw_entropy_part(sa,t,p);
	for (no_iter=1; no_iter <= 2; no_iter++)
	{
		pt_old	= pt;
		dentropy	= gsw_entropy_part(sa,pt_old,p_ref) - true_entropy_part;
		pt = pt_old - dentropy/dentropy_dt;
		ptm = 0.5*(pt + pt_old);
		dentropy_dt	= -gsw_gibbs(n0,n2,n0,sa,ptm,p_ref);
		pt = pt_old - dentropy/dentropy_dt;
	}
	return (pt);
}

/**************************************************************************
==========================================================================
method: gsw_pt0_from_t(sa,t,p)
==========================================================================
 Calculates potential temperature with reference pressure, p_ref = 0 dbar.
 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]
 gsw_pt0_from_t : potential temperature, p_ref = 0        [deg C]
*************************************************************************/

double TeosBase::gsw_pt0_from_t(double sa, double t, double p)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	int no_iter;
	double pt0, pt0_old, dentropy, dentropy_dt;
	double s1, true_entropy_part, pt0m;
	s1	= sa/gtc.gsw_ups;
	pt0 = t+p*( 8.65483913395442e-6  -
					s1 *  1.41636299744881e-6  -
					p *  7.38286467135737e-9  +
					t *(-8.38241357039698e-6  +
						 s1 *  2.83933368585534e-8  +
						 t *  1.77803965218656e-8  +
						 p *  1.71155619208233e-10));
	dentropy_dt	= gtc.gsw_cp0/((gtc.gsw_t0+pt0)*(1.0-0.05*(1.0-sa/gtc.gsw_sso)));
	true_entropy_part = gsw_entropy_part(sa,t,p);
	for (no_iter=1; no_iter <= 2; no_iter++)
	{
		pt0_old	= pt0;
		dentropy	= gsw_entropy_part_zerop(sa,pt0_old) -
					  true_entropy_part;
		pt0 = pt0_old - dentropy/dentropy_dt;
		pt0m	= 0.5*(pt0 + pt0_old);
		dentropy_dt	= -gsw_gibbs_pt0_pt0(sa,pt0m);
		pt0 = pt0_old - dentropy/dentropy_dt;
	}
	return (pt0);
}

/**************************************************************************
=========================================================================
method: gsw_pt0_from_t_ice (t, p)
=========================================================================
   Calculates potential temperature of ice Ih with a reference pressure of
   0 dbar, from in-situ temperature, t.
   t   :  in-situ temperature  (ITS-90)                           [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   pt0_ice  :  potential temperature of ice Ih with reference pressure of
               zero dbar (ITS-90)                                 [ deg C ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_pt0_from_t_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int number_of_iterations;
	double dentropy, dentropy_dt, pt0_ice,
			 pt0_ice_old, ptm_ice, true_entropy,
			 /**This is the starting polynomial for pt0 of ice Ih.*/
			 s1 = -2.256611570832386e-4,
			 s2 = -6.045305921314694e-7,
			 s3 =  5.546699019612661e-9,
			 s4 =  1.795030639186685e-11,
			 s5 =  1.292346094030742e-9,
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
	true_entropy = -gsw_gibbs_ice_part_t(t,p);
	if (t < -45.0 && t > -273.0)
	{
		pt0_ice = t + p*(p1 + p*(p2 + p3*t) + t*(p4 + t*(p5 + p6*t)));
		if (pt0_ice < -gtc.gsw_t0) pt0_ice = -gtc.gsw_t0;
		/**
		   we add 0.05d0 to the initial estimate of pt0_ice at
		   temps less than -273 to ensure that it is never less than -273.15.
		*/
		if (pt0_ice < -273.0) pt0_ice = pt0_ice + 0.05;
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);
		for (number_of_iterations = 1; number_of_iterations <= 3;
				number_of_iterations++)
		{
			pt0_ice_old = pt0_ice;
			dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
			ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
			dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
		}
	}
	else
	{
		pt0_ice = t + p*(s1 + t*(s2 + t*(s3 + t*s4)) + s5*p);
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);
		pt0_ice_old = pt0_ice;
		dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
		pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
		ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
		pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
	}
	if (pt0_ice < -273.0)
	{
		pt0_ice = t + p*(q1 + p*(q2 + q3*t) + t*(q4 + t*(q5 + q6*t)));
		/**
		   add 0.01d0 to the initial estimate of pt_ice used in the
		   derivative to ensure that it is never less than -273.15d0
		   because the derivative approaches zero at absolute zero.
		*/
		dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice+0.01);
		for (number_of_iterations = 1; number_of_iterations <= 3;
				number_of_iterations++)
		{
			pt0_ice_old = pt0_ice;
			dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
			ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
			/**
			   add 0.01d0 to the estimate of ptm_ice for temperatures less
			   than -273 to ensure that they are never less than -273.15d0
			   because the derivative approaches zero at absolute zero and
			   the addition of 0.01d0 degrees c ensures that when we divide
			   by the derivatve in the modified newton routine the method
			   does not blow up.
			*/
			ptm_ice = ptm_ice + 0.01;
			dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
			pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
		}
	}
	/**
	   For temperatures less than -273.1 degsC the maximum error is less
	   than 2x10^-7 degsC. For temperatures between -273.1 and 273 the
	   maximum error is less than 8x10^-8 degsC, and for temperatures
	   greater than -273 degsC the   maximum error is 1.5x10^-12 degsC.
	   These errors are over the whole ocean depths with p varying between
	   0 and 10,000 dbar, while the in-situ temperature varied independently
	   between -273.15 and +2 degsC.
	*/
	return (pt0_ice);
}

/**************************************************************************
==========================================================================
method: gsw_rho(sa,ct,p)
==========================================================================
--------------------------------------------------------------------------
(Density and enthalpy, based on the 48-term expression for density)
--------------------------------------------------------------------------
  Calculates in-situ density from Absolute Salinity and Conservative
  Temperature, using the computationally-efficient expression for
  specific volume in terms of SA, CT and p (Roquet et al., 2014).
 sa     : Absolute Salinity                               [g/kg]
 ct     : Conservative Temperature (ITS-90)               [deg C]
 p      : sea pressure                                    [dbar]
	   ( i.e. absolute pressure - 10.1325 dbar )
 rho    : in-situ density				   [kg/m]
*************************************************************************/

double TeosBase::gsw_rho(double sa, double ct, double p)
{
	return (1.0/gsw_specvol(sa,ct,p));
}

/**************************************************************************
==========================================================================
method: gsw_rho_alpha_beta (sa, ct, p, rho, alpha, beta)
==========================================================================
   Calculates in-situ density, the appropiate thermal expansion coefficient
   and the appropriate saline contraction coefficient of seawater from
   Absolute Salinity and Conservative Temperature.  This method uses the
   computationally-efficient expression for specific volume in terms of
   SA, CT and p (Roquet et al., 2014).
   Note that potential density (pot_rho) with respect to reference pressure
   p_ref is obtained by calling this method with the pressure argument
   being p_ref as in [pot_rho, ~, ~] = gsw_rho_alpha_beta(SA,CT,p_ref).
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   rho    :  in-situ density                                       [ kg/m ]
   alpha  :  thermal expansion coefficient                          [ 1/K ]
             with respect to Conservative Temperature
   beta   :  saline (i.e. haline) contraction                      [ kg/g ]
             coefficient at constant Conservative Temperature
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_rho_alpha_beta(double sa,double ct,double p,double *rho,double *alpha,double *beta)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  v, v_ct_part, v_sa_part, xs, ys, z;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v = gsw_gsvco_v(xs, ys, z);
	if (rho != NULL)
		*rho = 1.0/v;
	if (alpha != NULL)
	{
		v_ct_part = gsw_gsvco_a(xs, ys, z);
		*alpha = 0.025*v_ct_part/v;
	}
	if (beta != NULL)
	{
		v_sa_part = gsw_gsvco_b(xs, ys, z);
		*beta = -v_sa_part*0.5*gtc.gsw_sfac/(v*xs);
	}
}

/***************************************************************************
   method: gsw_rr68_interp_section (sectnum, ip_sect, ip_isect, sa_i, ct_i)
****************************************************************************/

void TeosBase::gsw_rr68_interp_section(int sectnum,double *sa,double *ct,double *p,int mp,int nsect,double *ip_sect,int *ip_isect,double *p_i,double *sa_i,double *ct_i)
{
	int     i, *ip_1=NULL, *ip_2=NULL, *ip_3=NULL, *ip_4=NULL;
	double  m, *ct_12=NULL, *ct_13=NULL, *ct_23=NULL, *ct_34=NULL, ctp1,
					ctp2, *ct_ref=NULL, ctref_denom,
							 ct_ref_minus_ctp1, ct_ref_minus_ctp2,
							 ctref_num,
							 gamma1_23, gamma1_24, gamma2_31,
							 gamma2_34, gamma2_41, gamma3_12,
							 gamma3_42, gamma4_12, gamma4_23,
							 *sa_12=NULL, *sa_13=NULL, *sa_23=NULL, *sa_34=NULL, sap1,
							  sap2, *sa_ref=NULL, saref_denom,
										sa_ref_minus_sap1, sa_ref_minus_sap2,
										saref_num, *p_ii=NULL;
	ip_1 = new int[4*nsect];
	ip_2 = ip_1+nsect;
	ip_3 = ip_2+nsect;
	ip_4 = ip_3+nsect;
	ct_12 = new double[12*nsect];
	sa_12 = ct_12 +  1*nsect;
	sa_13 = ct_12 +  2*nsect;
	sa_23 = ct_12 +  3*nsect;
	sa_34 = ct_12 +  4*nsect;
	sa_ref = ct_12 +  5*nsect;
	ct_13 = ct_12 +  6*nsect;
	ct_23 = ct_12 +  7*nsect;
	ct_34 = ct_12 +  8*nsect;
	ct_ref = ct_12 +  9*nsect;
	p_ii = ct_12 + 10*nsect;
	if (sectnum < 0) /** shallow */
	{
		for (i=0; i<nsect; i++)
		{
			ip_1[i] = floor(ip_sect[i]);
			ip_2[i] = ceil(ip_sect[i]);
			if (ip_1[i] == ip_2[i])
			{
				ip_2[i] = ip_1[i] + 1;
			}
			ip_3[i] = ip_2[i] + 1;
			ip_4[i] = ip_3[i] + 1;
		}
	}
	else if (sectnum == 0) /** central */
	{
		for (i=0; i<nsect; i++)
		{
			ip_2[i] = floor(ip_sect[i]);
			ip_3[i] = ceil(ip_sect[i]);
			if (ip_2[i] == ip_3[i])
			{
				ip_2[i] = ip_3[i] - 1;
			}
			ip_1[i] = ip_2[i] - 1;
			if (ip_1[i] < 0)
			{
				ip_1[i] = 0;
				ip_2[i] = 1;
				ip_3[i] = 2;
			}
			ip_4[i] = ip_3[i] + 1;
		}
	}
	else if (sectnum > 0) /** deep */
	{
		for (i=0; i<nsect; i++)
		{
			ip_1[i] = ceil(ip_sect[i]);
			ip_2[i] = floor(ip_sect[i]);
			if (ip_1[i] == ip_2[i])
			{
				ip_2[i] = ip_1[i] - 1;
			}
			ip_3[i] = ip_2[i] - 1;
			ip_4[i] = ip_3[i] - 1;
		}
	}
	for (i=0; i<nsect; i++)
	{
		p_ii[i] = p_i[ip_isect[i]];
	}
	/** eqn (3d) */
	for (i=0; i<nsect; i++)
	{
		sa_34[i] = sa[ip_3[i]] + ((sa[ip_4[i]] - sa[ip_3[i]])
										  *(p_ii[i] - p[ip_3[i]])/ (p[ip_4[i]] - p[ip_3[i]]));
		ct_34[i] = ct[ip_3[i]] + ((ct[ip_4[i]] - ct[ip_3[i]])
										  *(p_ii[i] - p[ip_3[i]])/ (p[ip_4[i]] - p[ip_3[i]]));
	}
	/**
	  Construct the Reiniger & Ross reference curve equation.
	  m = the power variable
	*/
	m = 1.7;
	if (sectnum == 0)
	{
		gsw_linear_interp_sa_ct(sa,ct,p,mp,p_ii,nsect,sa_23,ct_23);
		/** eqn (3a) */
		for (i=0; i<nsect; i++)
		{
			sa_12[i] = sa[ip_1[i]] + ((sa[ip_2[i]] - sa[ip_1[i]])
											  *(p_ii[i] - p[ip_1[i]])/ (p[ip_2[i]] - p[ip_1[i]]));
			ct_12[i] = ct[ip_1[i]] + ((ct[ip_2[i]] - ct[ip_1[i]])
											  *(p_ii[i] - p[ip_1[i]])/ (p[ip_2[i]] - p[ip_1[i]]));
			saref_num = (pow(fabs(sa_23[i]-sa_34[i]),m))*sa_12[i]
							+ (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i];
			ctref_num = (pow(fabs(ct_23[i]-ct_34[i]),m))*ct_12[i]
							+ (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i];
			saref_denom = pow(fabs(sa_23[i]-sa_34[i]),m)
							  + pow(fabs(sa_12[i]-sa_23[i]),m);
			ctref_denom = pow(fabs(ct_23[i]-ct_34[i]),m)
							  + pow(fabs(ct_12[i]-ct_23[i]),m);
			if (saref_denom == 0.0)
			{
				sa_23[i] = sa_23[i] + 1.0e-6;
				saref_num = (pow(fabs(sa_23[i]-sa_34[i]),m))*sa_12[i]
								+ (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i];
				saref_denom = pow(fabs(sa_23[i]-sa_34[i]),m)
								  + pow(fabs(sa_12[i]-sa_23[i]),m);
			}
			if (ctref_denom == 0.0)
			{
				ct_23[i] = ct_23[i] + 1.0e-6;
				ctref_num = (pow(fabs(ct_23[i]-ct_34[i]),m))*ct_12[i]
								+ (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i];
				ctref_denom = pow(fabs(ct_23[i]-ct_34[i]),m)
								  + pow(fabs(ct_12[i]-ct_23[i]),m);
			}
			sa_ref[i] = 0.5*(sa_23[i] + (saref_num/saref_denom));
			ct_ref[i] = 0.5*(ct_23[i] + (ctref_num/ctref_denom));
		} /** for*/
	}
	else
	{
		gsw_linear_interp_sa_ct(sa,ct,p,mp,p_ii,nsect,sa_12,ct_12);
		for (i=0; i<nsect; i++)
		{
			sa_13[i] = sa[ip_1[i]] + ((sa[ip_3[i]] - sa[ip_1[i]])
											  *(p_ii[i] - p[ip_1[i]])/ (p[ip_3[i]] - p[ip_1[i]]));
			ct_13[i] = ct[ip_1[i]] + ((ct[ip_3[i]] - ct[ip_1[i]])
											  *(p_ii[i] - p[ip_1[i]])/ (p[ip_3[i]] - p[ip_1[i]]));
			sa_23[i] = sa[ip_2[i]] + ((sa[ip_3[i]] - sa[ip_2[i]])
											  *(p_ii[i] - p[ip_2[i]])/ (p[ip_3[i]] - p[ip_2[i]]));
			ct_23[i] = ct[ip_2[i]] + ((ct[ip_3[i]] - ct[ip_2[i]])
											  *(p_ii[i] - p[ip_2[i]])/ (p[ip_3[i]] - p[ip_2[i]]));
			/** eqn (3a ) */
			saref_num = (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i]
							+ (pow(fabs(sa_12[i]-sa_13[i]),m))*sa_23[i];
			ctref_num = (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i]
							+ (pow(fabs(ct_12[i]-ct_13[i]),m))*ct_23[i];
			saref_denom = pow(fabs(sa_12[i]-sa_23[i]),m)
							  + pow(fabs(sa_12[i]-sa_13[i]),m);
			ctref_denom = pow(fabs(ct_12[i]-ct_23[i]),m)
							  + pow(fabs(ct_12[i]-ct_13[i]),m);
			if (saref_denom == 0.0)
			{
				sa_23[i] = sa_23[i] + 1.0e-6;
				saref_num = (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i]
								+ (pow(fabs(sa_12[i]-sa_13[i]),m))*sa_23[i];
				saref_denom = pow(fabs(sa_12[i]-sa_23[i]),m)
								  + pow(fabs(sa_12[i]-sa_13[i]),m);
			}
			if (ctref_denom == 0.0)
			{
				ct_23[i] = ct_23[i] + 1.0e-6;
				ctref_num = (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i]
								+ (pow(fabs(ct_12[i]-ct_13[i]),m))*ct_23[i];
				ctref_denom = pow(fabs(ct_12[i]-ct_23[i]),m)
								  + pow(fabs(ct_12[i]-ct_13[i]),m);
			}
			sa_ref[i] = 0.5*(sa_12[i] + (saref_num/saref_denom));
			ct_ref[i] = 0.5*(ct_12[i] + (ctref_num/ctref_denom));
		} /** for */
	}
	for (i=0; i<nsect; i++)
	{
		/** eqn (3c) */
		gamma1_23 = ((p_ii[i] - p[ip_2[i]])*(p_ii[i] - p[ip_3[i]]))/
						((p[ip_1[i]] - p[ip_2[i]])*(p[ip_1[i]] - p[ip_3[i]]));
		gamma2_31 = ((p_ii[i] - p[ip_3[i]])*(p_ii[i] - p[ip_1[i]]))/
						((p[ip_2[i]] - p[ip_3[i]])*(p[ip_2[i]] - p[ip_1[i]]));
		gamma3_12 = ((p_ii[i] - p[ip_1[i]])*(p_ii[i] - p[ip_2[i]]))/
						((p[ip_3[i]] - p[ip_1[i]])*(p[ip_3[i]] - p[ip_2[i]]));
		if (sectnum == 0)
		{
			gamma2_34 = ((p_ii[i] - p[ip_3[i]])*(p_ii[i] - p[ip_4[i]]))/
							((p[ip_2[i]] - p[ip_3[i]])*(p[ip_2[i]] - p[ip_4[i]]));
			gamma3_42 = ((p_ii[i] - p[ip_4[i]])*(p_ii[i] - p[ip_2[i]]))/
							((p[ip_3[i]] - p[ip_4[i]])*(p[ip_3[i]] - p[ip_2[i]]));
			gamma4_23 = ((p_ii[i] - p[ip_2[i]])*(p_ii[i] - p[ip_3[i]]))/
							((p[ip_4[i]] - p[ip_2[i]])*(p[ip_4[i]] - p[ip_3[i]]));
		}
		else
		{
			gamma1_24 = ((p_ii[i] - p[ip_2[i]])*(p_ii[i] - p[ip_4[i]]))/
							((p[ip_1[i]] - p[ip_2[i]])*(p[ip_1[i]] - p[ip_4[i]]));
			gamma2_41 = ((p_ii[i] - p[ip_4[i]])*(p_ii[i] - p[ip_1[i]]))/
							((p[ip_2[i]] - p[ip_4[i]])*(p[ip_2[i]] - p[ip_1[i]]));
			gamma4_12 = ((p_ii[i] - p[ip_1[i]])*(p_ii[i] - p[ip_2[i]]))/
							((p[ip_4[i]] - p[ip_1[i]])*(p[ip_4[i]] - p[ip_2[i]]));
		}
		/** eqn (3b/3b ) */
		sap1 = gamma1_23*sa[ip_1[i]] + gamma2_31*sa[ip_2[i]]
				 + gamma3_12*sa[ip_3[i]];
		ctp1 = gamma1_23*ct[ip_1[i]] + gamma2_31*ct[ip_2[i]]
				 + gamma3_12*ct[ip_3[i]];
		if (sectnum == 0)
		{
			sap2 = gamma2_34*sa[ip_2[i]] + gamma3_42*sa[ip_3[i]]
					 + gamma4_23*sa[ip_4[i]];
			ctp2 = gamma2_34*ct[ip_2[i]] + gamma3_42*ct[ip_3[i]]
					 + gamma4_23*ct[ip_4[i]];
		}
		else
		{
			sap2 = gamma1_24*sa[ip_1[i]] + gamma2_41*sa[ip_2[i]]
					 + gamma4_12*sa[ip_4[i]];
			ctp2 = gamma1_24*ct[ip_1[i]] + gamma2_41*ct[ip_2[i]]
					 + gamma4_12*ct[ip_4[i]];
		}
		/** eqn (3) */
		sa_ref_minus_sap1 = fabs(sa_ref[i] - sap1);
		sa_ref_minus_sap2 = fabs(sa_ref[i] - sap2);
		if (sa_ref_minus_sap1 == 0.0 && sa_ref_minus_sap2 == 0.0)
		{
			sa_ref[i] = sa_ref[i] + 1.0e-6;
			sa_ref_minus_sap1 = fabs(sa_ref[i] - sap1);
			sa_ref_minus_sap2 = fabs(sa_ref[i] - sap2);
		}
		ct_ref_minus_ctp1 = fabs(ct_ref[i] - ctp1);
		ct_ref_minus_ctp2 = fabs(ct_ref[i] - ctp2);
		if (ct_ref_minus_ctp1 == 0.0 && ct_ref_minus_ctp2 == 0.0)
		{
			ct_ref[i] = ct_ref[i] + 1.0e-6;
			ct_ref_minus_ctp1 = fabs(ct_ref[i] - ctp1);
			ct_ref_minus_ctp2 = fabs(ct_ref[i] - ctp2);
		}
		sa_i[ip_isect[i]] = (sa_ref_minus_sap1*sap2+sa_ref_minus_sap2*sap1)/
								  (sa_ref_minus_sap1 + sa_ref_minus_sap2);
		ct_i[ip_isect[i]] = (ct_ref_minus_ctp1*ctp2+ct_ref_minus_ctp2*ctp1)/
								  (ct_ref_minus_ctp1 + ct_ref_minus_ctp2);
	}
	if (ct_12) delete []ct_12;
	if (ip_1) delete []ip_1;
}

/**************************************************************************
==========================================================================
method: gsw_sa_from_sp(sp,p,lon,lat)
==========================================================================
 Calculates Absolute Salinity, SA, from Practical Salinity, SP
 sp     : Practical Salinity                              [unitless]
 p      : sea pressure                                    [dbar]
 lon	 : longitude                                       [DEG E]
 lat    : latitude                                        [DEG N]
 gsw_sa_from_sp   : Absolute Salinity                     [g/kg]
*************************************************************************/

double TeosBase::gsw_sa_from_sp(double sp, double p, double lon, double lat)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double saar, gsw_sa_baltic;
	gsw_sa_baltic = gsw_sa_from_sp_baltic(sp,lon,lat);
	if (gsw_sa_baltic < cppGSW_ERROR_LIMIT)
		return (gsw_sa_baltic);
	saar = gsw_saar(p,lon,lat);
	if (saar == cppGSW_INVALID_VALUE)
		return (saar);
	return (gtc.gsw_ups*sp*(1.e0 + saar));
}

/**************************************************************************
==========================================================================
method: gsw_sa_from_sp_baltic(sp,lon,lat)
==========================================================================
 For the Baltic Sea, calculates Absolute Salinity with a value
 computed analytically from Practical Salinity
 sp     : Practical Salinity                              [unitless]
 lon    : longitude                                       [deg E]
 lat    : latitude                                        [deg N]
 sa_from_sp_baltic : Absolute Salinity                    [g/kg]
*************************************************************************/

double TeosBase::gsw_sa_from_sp_baltic(double sp, double lon, double lat)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	/**  for GSW_BALTIC_DATA use gbalt */
	double xx_left, xx_right, return_value;
	if (gbalt.xb_left[1] < lon  && lon < gbalt.xb_right[0]
			&& gbalt.yb_left[0] < lat  && lat < gbalt.yb_left[2])
	{
		xx_left	= gsw_util_xinterp1(gbalt.yb_left, gbalt.xb_left, 3, lat);
		xx_right = gsw_util_xinterp1(gbalt.yb_right, gbalt.xb_right, 2, lat);
		if (xx_left <= lon  && lon <= xx_right)
			return_value = ((gtc.gsw_sso - 0.087)/35.0)*sp + 0.087;
		else
			return_value = cppGSW_INVALID_VALUE;
	}
	else
		return_value = cppGSW_INVALID_VALUE;
	return (return_value);
}

/**************************************************************************
==========================================================================
method: gsw_sa_from_sstar(sstar,p,lon,lat)
==========================================================================
  Calculates Absolute Salinity, SA, from Preformed Salinity, Sstar.
  Sstar  : Preformed Salinity                              [g/kg]
  p      : sea pressure                                    [dbar]
  lon   : longitude                                       [deg E]
  lat    : latitude                                        [deg N]
  gsw_sa_from_sstar   : Absolute Salinity                  [g/kg]
*************************************************************************/

double TeosBase::gsw_sa_from_sstar(double sstar,double p,double lon,double lat)
{
	double  saar = gsw_saar(p,lon,lat);
	if (saar == cppGSW_INVALID_VALUE)
		return (saar);
	/**
	  In the Baltic Sea, Sstar = SA, and note that gsw_saar returns zero
	  for SAAR in the Baltic.
	*/
	return (sstar*(1e0 + saar)/(1e0 - 0.35e0*saar));
}

/**************************************************************************
==========================================================================
method: gsw_saar(p,lon,lat)
==========================================================================
 Calculates the Absolute Salinity Anomaly Ratio, SAAR.
 p      : sea pressure                                    [dbar]
 lon    : longitude                                       [deg E]
 lat    : latitude                                        [deg N]
 gsw_saar : Absolute Salinity Anomaly Ratio               [unitless]
*************************************************************************/

double TeosBase::gsw_saar(double p, double lon, double lat)
{
	/** for GSW_SAAR_DATA use gsaar */
	int nx=gsw_nx, ny=gsw_ny, nz=gsw_nz;
	int indx0, indy0, indz0, k, ndepth_index;
	double saar[4], saar_old[4];
	double sa_upper, sa_lower, dlong, dlat;
	double r1, s1, t1, ndepth_max, return_value = cppGSW_INVALID_VALUE;

	if (isnan(lat) || isnan(lon) || isnan(p))
	{
		return (return_value);
	}

	/** further down than -86 is not in the ocean (-86 to -90 is land) */
	if (lat < -86.0 || lat > 90.0)
	{
		return (return_value);
	}

	if (lon < 0.0)
	{
		lon += 360.0;
	}

	dlong = pSaarData->get_longs_ref(1)-pSaarData->get_longs_ref(0);
	dlat = pSaarData->get_lats_ref(1)-pSaarData->get_lats_ref(0);
	indx0 = floor(0 + (nx-1)*(lon-pSaarData->get_longs_ref(0))/
					  (pSaarData->get_longs_ref(nx-1)-pSaarData->get_longs_ref(0)));

	if (indx0 == nx-1)
	{
		indx0 = nx-2;
	}

	indy0 = floor(0 + (ny-1)*(lat-pSaarData->get_lats_ref(0)) /
					  (pSaarData->get_lats_ref(ny-1)-pSaarData->get_lats_ref(0)));

	if(indy0 == ny-1)
	{
		indy0 = ny-2;
	}

	/**
	   Look for the maximum valid "ndepth_ref" value around our point.
	   Note: invalid "ndepth_ref" values are NaNs (a hangover from the codes
	   Matlab origins), but we have replaced the NaNs with a value of "9e90",
	   hence we need an additional upper-limit check in the code below so they
	   will not be recognised as valid values.
	*/
	ndepth_max = -1.0;
	for (k=0; k < 4; k++)
	{
		ndepth_index = indy0+gsaar.delj[k]+(indx0+gsaar.deli[k])*ny;
		if (pSaarData->get_ndepth_ref(ndepth_index) > 0.0 &&
				pSaarData->get_ndepth_ref(ndepth_index) < 1e90)
		{
			ndepth_max = gsw_max_d(ndepth_max, pSaarData->get_ndepth_ref(ndepth_index));
		}
	}
	/**
	   If we are a long way from the ocean then there will be no valid "ndepth_ref"
	   values near the point (ie. surrounded by NaNs) - so just return SAAR = 0.0
	*/
	if (ndepth_max == -1.0)
	{
		return (0.0);
	}
	if (p > pSaarData->get_p_ref((int)(ndepth_max)-1))
	{
		p = pSaarData->get_p_ref((int)(ndepth_max)-1);
	}
	indz0 = gsw_util_indxPref(nz,p);
	r1 = (lon-pSaarData->get_longs_ref(indx0))/
		  (pSaarData->get_longs_ref(indx0+1)-
			pSaarData->get_longs_ref(indx0));
	s1 = (lat-pSaarData->get_lats_ref(indy0))/
		  (pSaarData->get_lats_ref(indy0+1)-
			pSaarData->get_lats_ref(indy0));
	t1 = (p-pSaarData->get_p_ref(indz0))/
		  (pSaarData->get_p_ref(indz0+1)-
			pSaarData->get_p_ref(indz0));
	for (k=0; k<4; k++)
	{
		unsigned kpos = indz0+nz*(indy0+gsaar.delj[k]+(indx0+gsaar.deli[k])*ny);
		saar[k] = pSaarData->get_saar_ref(kpos);
	}
	if (gsaar.longs_pan[0] <= lon && lon <= gsaar.longs_pan[gsaar.npan-1]-0.001 &&
			gsaar.lats_pan[gsaar.npan-1] <= lat && lat <= gsaar.lats_pan[0])
	{
		memmove(saar_old,saar,4*sizeof (double));
		gsw_add_barrier(saar_old,lon,lat,pSaarData->get_longs_ref(indx0),
							 pSaarData->get_lats_ref(indy0),dlong,dlat,saar);
	}
	else if (fabs(gsw_sum(saar, sizeof (saar)/sizeof (double))) >=  cppGSW_ERROR_LIMIT)
	{
		memmove(saar_old,saar,4*sizeof (double));
		gsw_add_mean(saar_old,saar);
	}
	sa_upper = (1.0-s1)*(saar[0] + r1*(saar[1]-saar[0])) +
				  s1*(saar[3] + r1*(saar[2]-saar[3]));
	for (k=0; k<4; k++)
	{
		unsigned kpos = indz0+1+nz*(indy0+gsaar.delj[k]+(indx0+gsaar.deli[k])*ny);
		saar[k] = pSaarData->get_saar_ref(kpos);
	}
	if (gsaar.longs_pan[0] <= lon && lon <= gsaar.longs_pan[gsaar.npan-1]-0.001 &&
			gsaar.lats_pan[gsaar.npan-1] <= lat && lat <= gsaar.lats_pan[0])
	{
		memmove(saar_old,saar,4*sizeof (double));
		gsw_add_barrier(saar_old,lon,lat,pSaarData->get_longs_ref(indx0),
							 pSaarData->get_lats_ref(indy0),dlong,dlat,saar);
	}
	else if (fabs(gsw_sum(saar, sizeof (saar)/sizeof (double)))  >=  cppGSW_ERROR_LIMIT)
	{
		memmove(saar_old,saar,4*sizeof (double));
		gsw_add_mean(saar_old,saar);
	}
	sa_lower = (1.0-s1)*(saar[0] + r1*(saar[1]-saar[0])) +
				  s1*(saar[3] + r1*(saar[2]-saar[3]));
	if (fabs(sa_lower) >= cppGSW_ERROR_LIMIT)
	{
		sa_lower = sa_upper;
	}
	return_value = sa_upper + t1*(sa_lower-sa_upper);
	if (fabs(return_value) >= cppGSW_ERROR_LIMIT)
	{
		return_value = cppGSW_INVALID_VALUE;
	}
	return (return_value);
}

/*************************************************************************************************
    Functions for pchip interpolation
    (Piecewise Cubic Hermite Interpolating Polynomial)
    based on
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html#
    scipy.interpolate.PchipInterpolator
    See references therein.
    This is a shape-preserving algorithm, in the sense that it does not
    overshoot the original points; extrema of the interpolated curve match
    the extrema of the original points.
***************************************************************************************************/

int TeosBase::gsw_sgn(double x)
{
	return (((x) > 0.0) ? 1 : (((x) < 0) ? -1 : 0));
}

/**************************************************************************
==========================================================================
method: gsw_sp_from_sa(sa,p,lon,lat)
==========================================================================
  Calculates Practical salinity, sp, from Absolute salinity, sa
  sa     : Absolute Salinity                               [g/kg]
  p      : sea pressure                                    [dbar]
  lon    : longitude                                       [DEG E]
  lat    : latitude                                        [DEG N]
  gsw_sp_from_sa      : Practical Salinity                 [unitless]
*************************************************************************/

double TeosBase::gsw_sp_from_sa(double sa, double p, double lon, double lat)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double saar, gsw_sp_baltic;
	gsw_sp_baltic = gsw_sp_from_sa_baltic(sa,lon,lat);
	if (gsw_sp_baltic < cppGSW_ERROR_LIMIT)
		return (gsw_sp_baltic);
	saar = gsw_saar(p,lon,lat);
	if (saar == cppGSW_INVALID_VALUE)
		return (saar);
	return ((sa/gtc.gsw_ups)/(1e0 + saar));
}

/**************************************************************************
==========================================================================
method: gsw_sp_from_sa_baltic(sa,lon,lat)
==========================================================================
  For the Baltic Sea, calculates Practical Salinity with a value
  computed analytically from Absolute Salinity
  sa     : Absolute Salinity                               [g/kg]
  lon    : longitude                                       [deg E]
  lat    : latitude                                        [deg N]
  gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]
*************************************************************************/

double TeosBase::gsw_sp_from_sa_baltic(double sa, double lon, double lat)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_BALTIC_DATA use gbalt */
	double xx_left, xx_right, return_value=NAN;
	if (gbalt.xb_left[1] < lon  && lon < gbalt.xb_right[0] &&
			gbalt.yb_left[0] < lat  && lat < gbalt.yb_left[2])
	{
		xx_left = gsw_util_xinterp1(gbalt.yb_left, gbalt.xb_left, 3, lat);
		xx_right = gsw_util_xinterp1(gbalt.yb_right, gbalt.xb_right, 2, lat);
		if (xx_left <= lon  && lon <= xx_right)
			return_value = (35.0/(gtc.gsw_sso - 0.087))*(sa - 0.087);
		else
			return_value = cppGSW_INVALID_VALUE;
	}
	else
		return_value = cppGSW_INVALID_VALUE;
	return (return_value);
}

/**************************************************************************
==========================================================================
method: gsw_sp_from_c(c,t,p)
==========================================================================
--------------------------------------------------------------------------
	Practical Salinity (SP), PSS-78
--------------------------------------------------------------------------
   Calculates Practical Salinity, SP, from conductivity, C, primarily using
   the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical
   Salinity is only valid in the range 2 < SP < 42.  If the PSS-78
   algorithm produces a Practical Salinity that is less than 2 then the
   Practical Salinity is recalculated with a modified form of the Hill et
   al. (1986) formula.  The modification of the Hill et al. (1986)
   expression is to ensure that it is exactly consistent with PSS-78
   at SP = 2.  Note that the input values of conductivity need to be in
   units of mS/cm (not S/m).
  c      : conductivity                                     [ mS/cm ]
  t      : in-situ temperature [ITS-90]                     [deg C]
  p      : sea pressure                                     [dbar]
  sp     : Practical Salinity                               [unitless]
*************************************************************************/

double TeosBase::gsw_sp_from_c(double c, double t, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SP_COEFFICIENTS use gspc */
	double sp, t68, ft68, r, rt_lc, rp, rt, rtx,
			 hill_ratio, x, sqrty, part1, part2,
			 sp_hill_raw;
	t68 = t*1.00024e0;
	ft68 = (t68 - 15e0)/(1e0 + gspc.k*(t68 - 15e0));
	/**
	   The dimensionless conductivity ratio, R, is the conductivity input, C,
	   divided by the present estimate of C(SP=35, t_68=15, p=0) which is
	   42.9140 mS/cm (=4.29140 S/m), (Culkin and Smith, 1980).
	*/
	r = c/gtc.gsw_c3515;        /** 0.023302418791070513 = 1./42.9140 */
	/**rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.*/
	rt_lc = gspc.c0 + (gspc.c1 + (gspc.c2 + (gspc.c3 + gspc.c4*t68)*t68)*t68)*t68;
	rp = 1e0 + (p*(gspc.e1 + gspc.e2*p + gspc.e3*p*p))/(1e0 + gspc.d1*t68 + gspc.d2*t68*t68 +
			(gspc.d3 + gspc.d4*t68)*r);
	rt = r/(rp*rt_lc);
	if (rt < 0.0)
	{
		return (cppGSW_INVALID_VALUE);
	}
	rtx = sqrt(rt);
	sp = gspc.a0
		  + (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4
											+ gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx
		  + ft68*(gspc.b0 + (gspc.b1 + (gspc.b2 + (gspc.b3
												  + (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);
	/**
	   The following section of the code is designed for SP < 2 based on the
	   Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
	   exactly equal to the PSS-78 algorithm at SP = 2.
	*/
	if (sp < 2)
	{
		hill_ratio  = gsw_hill_ratio_at_sp2(t);
		x = 400e0*rt;
		sqrty = 10e0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		sp_hill_raw = sp - gspc.a0/part1 - gspc.b0*ft68/part2;
		sp = hill_ratio*sp_hill_raw;
	}
	/** This line ensures that SP is non-negative. */
	if (sp < 0.0)
	{
		sp  = cppGSW_INVALID_VALUE;
	}
	return (sp);
}

/**************************************************************************
==========================================================================
method: gsw_specvol(sa,ct,p)
==========================================================================
  Calculates specific volume from Absolute Salinity, Conservative
  Temperature and pressure, using the computationally-efficient
  polynomial expression for specific volume (Roquet et al., 2014).
 sa     : Absolute Salinity                               [g/kg]
 ct     : Conservative Temperature (ITS-90)               [deg C]
 p      : sea pressure                                    [dbar]
          ( i.e. absolute pressure - 10.1325 dbar )
 specvol : specific volume                                 [m^3/kg]
*************************************************************************/

double TeosBase::gsw_specvol(double sa, double ct, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double xs, ys, z, value;
	xs	= sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys	= ct*0.025;
	z	= p*1e-4;
	value = gsw_gsvco_v(xs, ys, z);
	return (value);
}

/**************************************************************************
==========================================================================
method: gsw_specvol_ice (t, p)
==========================================================================
   Calculates the specific volume of ice.
   t  :  in-situ temperature (ITS-90)                             [ deg C ]
   p  :  sea pressure                                              [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   specvol_ice  :  specific volume                               [ m^3/kg ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_specvol_ice(double t,double p)
{
	return (gsw_gibbs_ice(0,1,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_specvol_t_exact(sa,t,p)
==========================================================================
  Calculates the specific volume of seawater
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  specvol_t_exact : specific volume                        [kg/m^3]
*************************************************************************/

double TeosBase::gsw_specvol_t_exact(double sa, double t, double p)
{
	return gsw_gibbs(0,0,1,sa,t,p);
}

/**************************************************************************
==========================================================================
method: gsw_specvol_sso_0(p)
==========================================================================
  This method calculates specific volume at the Standard Ocean Salinty,
  SSO, and at a Conservative Temperature of zero degrees C, as a function
  of pressure, p, in dbar, using a streamlined version of the CT version
  of specific volume, that is, a streamlined version of the code
  "gsw_specvol(SA,CT,p)".
 p      : sea pressure                                    [dbar]
 specvol_sso_0 : specvol(sso,0,p)                         [m  kg  ]
*************************************************************************/

double TeosBase::gsw_specvol_sso_0(double p)
{
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double z, return_value;
	z = p*1.0e-4;
	return_value = 9.726613854843870e-04 + z*(-4.505913211160929e-05
						+ z*(7.130728965927127e-06 + z*(-6.657179479768312e-07
								+ z*(-2.994054447232880e-08 + z*(gsvco.v005 + gsvco.v006*z)))));
	return (return_value);
}



/**************************************************************************
==========================================================================
method: gsw_sstar_from_sa(sa,p,lon,lat)
==========================================================================
  Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA.
  sa     : Absolute Salinity                               [g/kg]
  p      : sea pressure                                    [dbar]
  lon   : longitude                                       [deg E]
  lat    : latitude                                        [deg N]
  gsw_sstar_from_sa : Preformed Salinity                   [g/kg]
*************************************************************************/

double TeosBase::gsw_sstar_from_sa(double sa,double p,double lon,double lat)
{
	double  saar = gsw_saar(p,lon,lat);
	/**
	   In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
	   for saar in the Baltic.
	*/
	if (saar == cppGSW_INVALID_VALUE)
	{
		return (saar);
	}
	return (sa*(1e0 - 0.35e0*saar)/(1e0 + saar));
}

/**************************************************************************
==========================================================================
method: gsw_sstar_from_sp(sp,p,lon,lat)
==========================================================================
  Calculates Preformed Salinity, Sstar, from Practical Salinity, SP.
  sp     : Practical Salinity                              [unitless]
  p      : sea pressure                                    [dbar]
  lon    : longitude                                       [deg E]
  lat    : latitude                                        [deg N]
  gsw_sstar_from_sp  : Preformed Salinity                  [g/kg]
*************************************************************************/

double TeosBase::gsw_sstar_from_sp(double sp,double p,double lon,double lat)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  saar, sstar_baltic = gsw_sa_from_sp_baltic(sp,lon,lat);
	/**
	  In the Baltic Sea, Sstar = SA.
	*/
	if (sstar_baltic < cppGSW_ERROR_LIMIT)
	{
		return (sstar_baltic);
	}
	saar = gsw_saar(p,lon,lat);
	if (saar == cppGSW_INVALID_VALUE)
	{
		return (saar);
	}
	return (gtc.gsw_ups*sp*(1 - 0.35e0*saar));
}

/**************************************************************************
==========================================================================
method: gsw_sum(*x,n)
==========================================================================
*************************************************************************/

double TeosBase::gsw_sum(double *x, int n)
{
	int	i;
	double	val;
	for (val=0.0, i=0; i<n; val += x[i], i++);
	return (val);
}

/**************************************************************************
==========================================================================
method: gsw_t_deriv_chem_potential_water_t_exact (sa, t, p)
==========================================================================
   Calculates the temperature derivative of the chemical potential of water
   in seawater so that it is valid at exactly SA = 0.
   SA  :  Absolute Salinity                                        [ g/kg ]
   t   :  in-situ temperature (ITS-90)                            [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   chem_potential_water_dt :  temperature derivative of the chemical
                            potential of water in seawater  [ J g^-1 K^-1 ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_t_deriv_chem_potential_water_t_exact(double sa, double t, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double g03_t, g08_sa_t, x, x2, y, z, g08_t, kg2g = 1e-3;
	/**
	   Note. The kg2g, a factor of 1e-3, is needed to convert the output of this
	   function into units of J/g. See section (2.9) of the TEOS-10 Manual.
	*/
	x2 = gtc.gsw_sfac*sa;
	x = sqrt(x2);
	y = t*0.025;
	z = p*gtc.rec_db2pa;
	/** the input pressure (p) is sea pressure in units of dbar. */
	g03_t = 5.90578347909402
			  + z*(-270.983805184062 + z*(776.153611613101 + z*(-196.51255088122
													+ (28.9796526294175 - 2.13290083518327*z)*z))) + y*(-24715.571866078
															+ z*(2910.0729080936 + z*(-1513.116771538718 + z*(546.959324647056
																	+ z*(-111.1208127634436 + 8.68841343834394*z)))) + y*(2210.2236124548363
																			+ z*(-2017.52334943521 + z*(1498.081172457456 + z*(-718.6359919632359
																					+ (146.4037555781616 - 4.9892131862671505*z)*z))) + y*(-592.743745734632
																							+ z*(1591.873781627888 + z*(-1207.261522487504 + (608.785486935364
																									- 105.4993508931208*z)*z)) + y*(290.12956292128547 + z*(-973.091553087975
																											+ z*(602.603274510125 + z*(-276.361526170076 + 32.40953340386105*z)))
																											+ y*(-113.90630790850321 + y*(21.35571525415769 - 67.41756835751434*z)
																													+ z*(381.06836198507096 + z*(-133.7383902842754 + 49.023632509086724*z)))))));
	g08_t = x2*(168.072408311545 + x*(-493.407510141682 + x*(543.835333000098
												 + x*(-196.028306689776 + 36.7571622995805*x) + y*(-137.1145018408982
														 + y*(148.10030845687618 + y*(-68.5590309679152 + 12.4848504784754*y)))
												 - 22.6683558512829*z) + z*(-175.292041186547
														 + (83.1923927801819 - 29.483064349429*z)*z) + y*(-86.1329351956084
																 + z*(766.116132004952 + z*(-108.3834525034224 + 51.2796974779828*z))
																 + y*(-30.0682112585625 - 1380.9597954037708*z + y*(3.50240264723578
																		 + 938.26075044542*z)))));
	g08_sa_t = 1187.3715515697959
				  + x*(-1480.222530425046 + x*(2175.341332000392 + x*(-980.14153344888
														 + 220.542973797483*x) + y*(-548.4580073635929 + y*(592.4012338275047
																 + y*(-274.2361238716608 + 49.9394019139016*y))) - 90.6734234051316*z)
						 + z*(-525.876123559641 + (249.57717834054571 - 88.449193048287*z)*z)
						 + y*(-258.3988055868252 + z*(2298.348396014856 + z*(-325.1503575102672
								 + 153.8390924339484*z)) + y*(-90.2046337756875 - 4142.8793862113125*z
										 + y*(10.50720794170734 + 2814.78225133626*z))));
	return (kg2g*((g03_t + g08_t)*0.025 - 0.5*gtc.gsw_sfac*0.025*sa*g08_sa_t));
}

/**************************************************************************
==========================================================================
method: gsw_t_freezing(sa,p,saturation_fraction)
==========================================================================
  Calculates the in-situ temperature at which seawater freezes
  sa     : Absolute Salinity                                 [g/kg]
  p      : sea pressure                                      [dbar]
          ( i.e. absolute pressure - 10.1325 dbar )
  saturation_fraction : the saturation fraction of dissolved air
                        in seawater
  t_freezing : in-situ temperature at which seawater freezes.[deg C]
*************************************************************************/

double TeosBase::gsw_t_freezing(double sa, double p, double saturation_fraction)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double sa_r, x, p_r;
	double  df_dt, tf, tfm, tf_old, f, return_value;
	/** The initial value of t_freezing_exact (for air-free seawater) */
	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;
	tf = gfpc.t0
		  + sa_r*(gfpc.t1 + x*(gfpc.t2 + x*(gfpc.t3 + x*(gfpc.t4 + x*(gfpc.t5 + gfpc.t6*x)))))
		  + p_r*(gfpc.t7 + p_r*(gfpc.t8 + gfpc.t9*p_r))
		  + sa_r*p_r*(gfpc.t10 + p_r*(gfpc.t12 + p_r*(gfpc.t15 + gfpc.t21*sa_r))
						  + sa_r*(gfpc.t13 + gfpc.t17*p_r + gfpc.t19*sa_r)
						  + x*(gfpc.t11 + p_r*(gfpc.t14 + gfpc.t18*p_r) + sa_r*(gfpc.t16 + gfpc.t20*p_r
								  + gfpc.t22*sa_r)));
	/** Adjust for the effects of dissolved air */
	tf -= saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gtc.gsw_sso));
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

/**************************************************************************
==========================================================================
method: gsw_t_freezing_first_derivatives (sa, p, &
                            saturation_fraction, tfreezing_sa, tfreezing_p)
==========================================================================
   Calculates the first derivatives of the in-situ temperature at which
   seawater freezes with respect to Absolute Salinity SA and pressure P (in
   Pa).  These expressions come from differentiating the expression that
   defines the freezing temperature, namely the equality between the
   chemical potentials of water in seawater and in ice.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction = the saturation fraction of dissolved air in
                         seawater
   tfreezing_SA : the derivative of the in-situ freezing temperature
                  (ITS-90) with respect to Absolute Salinity at fixed
                  pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]
   tfreezing_P  : the derivative of the in-situ freezing temperature
                  (ITS-90) with respect to pressure (in Pa) at fixed
                  Absolute Salinity                                [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_t_freezing_first_derivatives(double sa,double p,double saturation_fraction,double *tfreezing_sa,double *tfreezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double rec_denom,tf,g_per_kg = 1000.0;
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	rec_denom = 1.0/
					(g_per_kg*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p)
					 + gsw_entropy_ice(tf,p));
	if (tfreezing_sa != NULL)
		*tfreezing_sa =
			gsw_dilution_coefficient_t_exact(sa,tf,p)*rec_denom
			+ saturation_fraction*(1e-3)/(2.0*gtc.gsw_sso);
	if (tfreezing_p != NULL)
		*tfreezing_p =
			-(gsw_specvol_t_exact(sa,tf,p) - sa*gsw_gibbs(1,0,1,sa,tf,p)
			  - gsw_specvol_ice(tf,p))*rec_denom;
}

/**************************************************************************
==========================================================================
method: gsw_t_freezing_first_derivatives_poly (sa, p, &
                            saturation_fraction, tfreezing_sa, tfreezing_p)
==========================================================================
   Calculates the first derivatives of the in-situ temperature at which
   seawater freezes with respect to Absolute Salinity SA and pressure P (in
   Pa).  These expressions come from differentiating the expression that
   defines the freezing temperature, namely the equality between the
   chemical potentials of water in seawater and in ice.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
                         seawater
   tfreezing_SA : the derivative of the in-situ freezing temperature
                  (ITS-90) with respect to Absolute Salinity at fixed
                  pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]
   tfreezing_P  : the derivative of the in-situ freezing temperature
                  (ITS-90) with respect to pressure (in Pa) at fixed
                  Absolute Salinity                                [ K/Pa ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_t_freezing_first_derivatives_poly(double sa,double p,double saturation_fraction,double *tfreezing_sa,double *tfreezing_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x, c = 1e-3/(2.0*gtc.gsw_sso);
	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;
	if (tfreezing_sa != NULL)
	{
		*tfreezing_sa =
			(gfpc.t1 + x*(1.5*gfpc.t2 + x*(2.0*gfpc.t3 + x*(2.5*gfpc.t4 + x*(3.0*gfpc.t5
													 + 3.5*gfpc.t6*x)))) + p_r*(gfpc.t10 + x*(1.5*gfpc.t11 + x*(2.0*gfpc.t13
															 + x*(2.5*gfpc.t16 + x*(3.0*gfpc.t19 + 3.5*gfpc.t22*x))))
															 + p_r*(gfpc.t12 + x*(1.5*gfpc.t14 + x*(2.0*gfpc.t17 + 2.5*gfpc.t20*x))
																	 + p_r*(gfpc.t15 + x*(1.5*gfpc.t18 + 2.0*gfpc.t21*x)))))*1e-2
			+ saturation_fraction*c;
	}
	if (tfreezing_p != NULL)
	{
		*tfreezing_p =
			(gfpc.t7 + sa_r*(gfpc.t10 + x*(gfpc.t11 + x*(gfpc.t13 + x*(gfpc.t16 + x*(gfpc.t19 + gfpc.t22*x)))))
			 + p_r*(2.0*gfpc.t8 + sa_r*(2.0*gfpc.t12 + x*(2.0*gfpc.t14 + x*(2.0*gfpc.t17
												 + 2.0*gfpc.t20*x))) + p_r*(3.0*gfpc.t9 + sa_r*(3.0*gfpc.t15 + x*(3.0*gfpc.t18
														 + 3.0*gfpc.t21*x)))))*1e-8;
	}
}

/**************************************************************************
==========================================================================
method: gsw_t_freezing_poly (sa, p, saturation_fraction)
==========================================================================
   Calculates the in-situ temperature at which seawater freezes from a
   computationally efficient polynomial.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
                         seawater
   t_freezing : in-situ temperature at which seawater freezes.    [ deg C ]
                (ITS-90)
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_t_freezing_poly(double sa,double p,double saturation_fraction)
{
	double ctf, return_value;
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	return_value = gsw_t_from_ct(sa,ctf,p);
	return (return_value);
}

/**************************************************************************
==========================================================================
method: gsw_t_freezing_poly (sa, p, saturation_fraction, polynomial)
==========================================================================
   Calculates the in-situ temperature at which seawater freezes from a
   computationally efficient polynomial.
   SA  :  Absolute Salinity                                        [ g/kg ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   saturation_fraction : the saturation fraction of dissolved air in
                         seawater
   t_freezing : in-situ temperature at which seawater freezes.    [ deg C ]
                (ITS-90)
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_t_freezing_poly(double sa,double p,double saturation_fraction,int polynomial)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */
	double p_r, sa_r, x, ctf, sfrac, return_value;
	int direct_poly=polynomial;
	if (direct_poly)
	{
		sfrac = saturation_fraction;
		ctf = gsw_ct_freezing_poly(sa,p,sfrac);
		return_value = gsw_t_from_ct(sa,ctf,p);
	}
	else
	{
		/** Alternative calculation ... */
		sa_r = sa*1e-2;
		x = sqrt(sa_r);
		p_r = p*1e-4;
		return_value = gfpc.t0
							+ sa_r*(gfpc.t1 + x*(gfpc.t2 + x*(gfpc.t3 + x*(gfpc.t4 + x*(gfpc.t5 + gfpc.t6*x)))))
							+ p_r*(gfpc.t7 + p_r*(gfpc.t8 + gfpc.t9*p_r))
							+ sa_r*p_r*(gfpc.t10 + p_r*(gfpc.t12 + p_r*(gfpc.t15 + gfpc.t21*sa_r))
											+ sa_r*(gfpc.t13 + gfpc.t17*p_r + gfpc.t19*sa_r)
											+ x*(gfpc.t11 + p_r*(gfpc.t14 + gfpc.t18*p_r) + sa_r*(gfpc.t16 + gfpc.t20*p_r
													+ gfpc.t22*sa_r)));
		/** Adjust for the effects of dissolved air */
		return_value = return_value -
							saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gtc.gsw_sso));
	}
	return (return_value);
}

/**************************************************************************
==========================================================================
method: gsw_t_from_ct(sa,ct,p)
==========================================================================
 Calculates in-situ temperature from Conservative Temperature of seawater
 sa      : Absolute Salinity                              [g/kg]
 ct      : Conservative Temperature                       [deg C]
 gsw_t_from_ct : in-situ temperature                      [deg C]
*************************************************************************/

double TeosBase::gsw_t_from_ct(double sa, double ct, double p)
{
	double pt0, p0=0.0;
	pt0 = gsw_pt_from_ct(sa,ct);
	return (gsw_pt_from_t(sa,pt0,p0,p));
}

/**************************************************************************
==========================================================================
method: gsw_util_indx(x,n,z,k)
==========================================================================
  Finds the index of the value in a monotonically increasing array
  x	 :  array of monotonically increasing values
  n     :  length of the array
  z     :  value to be indexed
  K      : index K - if X(K) <= Z < X(K+1), or
  N-1     		    - if Z = X(N)
*************************************************************************/

int TeosBase::gsw_util_indx(double *x, int n, double z)
{
	int k, ku, kl, km;
	if (z > x[0] && z < x[n-1])
	{
		kl	= 0;
		ku	= n-1;
		while (ku-kl > 1)
		{
			km	= (ku+kl)>>1;
			if (z > x[km])
				kl	= km;
			else
				ku	= km;
		}
		k	= kl;
		if (z == x[k+1])
			k++;
	}
	else if (z <= x[0])
		k	= 0;
	else
		k	= n-2;
	return (k);
}

/***********************************************************
   The default value for parameter "ASC" is "false"
   (see TeosBase.h) which means that "orderBySmallToLarge"
   is called by "sort" (dArray[0] set to smallest value).
   When parameter "ASC" is "true" "orderByLargeToSmall" is
   called by "sort" (dArray[0] set to largest value).
   -------------
   The exercise is to place the "dArray" and "iArray"
   into a list of "struct DI" (see TeosBase.h) objects.
   That list is then ordered by the values in the dArray
   according to the sequence setting ("ASC"). The sequence
   begins at dArray[0] and proceeds to the final value.
************************************************************/

void TeosBase::gsw_util_sort_dbl(double *dArray,int nx,int *iArray,bool ASC)
{
	/** instantiate a C++ list object for containing "struct DI" objects */
	list<DI> *pDI = new list<DI>;
	int i;
	/**
	   Transfer the values to the struct object "DI". Then
	   place those "DI" objects into the "pDI" list object.
	*/
	for (i = 0; i < nx; i++)
	{
		DI di;
		/** populate the struct DI */
		di.i = iArray[i];
		di.d = dArray[i];
		/** place DI object into list object */
		pDI->push_back(di);
	}
	/** call the selected sort function */
	if (ASC)
		pDI->sort(orderByLargeToSmall);
	else
		pDI->sort(orderBySmallToLarge);
	/**
	   Instantiate an iterator object for the list object
	   and set the iterator position to "begin".
	*/
	list<DI>::const_iterator itr = pDI->begin();
	/**
	   Transfer the sorted list of DI object
	   values back into double & int arrays.
	*/
	i = 0;
	while(itr != pDI->end())
	{
		iArray[i] = itr->i;
		dArray[i] = itr->d;
		++i;
		++itr;
	}
	/** clean up */
	if (pDI)
	{
		pDI->empty();
		delete pDI;
	}
}

/**************************************************************************
==========================================================================
method: gsw_util_xinterp1(x,y,n,x0)
==========================================================================
 Linearly interpolate a real array
 x      : y array (Must be monotonic)
 y      : y array
 n      : length of X and Y arrays
 x0     : value to be interpolated
 gsw_xinterp1 : Linearly interpolated value
*************************************************************************/

double TeosBase::gsw_util_xinterp1(double *x, double *y, int n, double x0)
{
	int k;
	double r;
	k = gsw_util_indx(x,n,x0);
	r = (x0-x[k])/(x[k+1]-x[k]);
	return (y[k] + r*(y[k+1]-y[k]));
}

/**************************************************************************
==========================================================================
method: gsw_z_from_p(p,lat,geo_strf_dyn_height,sea_surface_geopotential)
==========================================================================
 Calculates the height z from pressure p
 p      : sea pressure                                    [dbar]
 lat    : latitude                                        [deg]
 gsw_z_from_p : height                                    [m]
*************************************************************************/

double TeosBase::gsw_z_from_p(double p,double lat,double geo_strf_dyn_height,double sea_surface_geopotential)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double x, sin2, b, c, a;
	x = sin(lat*gtc.deg2rad);
	sin2 = x*x;
	b = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);
	a = -0.5*gtc.gamma*b;
	c = gsw_enthalpy_sso_0(p)
		 - (geo_strf_dyn_height + sea_surface_geopotential);
	return (-2.0*c/(b + sqrt(b*b - 4.0*a*c)));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_alpha_wrt_t_ice(double t, double p)
{
	return (gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(0,1,t,p));
}

/**************************************************************************
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
***************************************************************************/

double TeosBase::gsw_chem_potential_water_ice(double t,double p)
{
	return (gsw_gibbs_ice(0,0,t,p));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_enthalpy_diff(double sa,double ct,double p_shallow,double p_deep)
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

/**************************************************************************
==========================================================================
  method: gsw_helmholtz_energy_ice (t, p)
==========================================================================
   Calculates the Helmholtz energy of ice.
   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   Helmholtz_energy_ice  =  Helmholtz energy of ice      [ J/kg ]
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_helmholtz_energy_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (gsw_gibbs_ice(0,0,t,p) -
			  (gtc.db2pa*p + gtc.gsw_p0)*
			  gsw_gibbs_ice(0,1,t,p));
}

/**************************************************************************
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
*************************************************************************/

void TeosBase::gsw_ice_fraction_to_freeze_seawater(double sa,double ct,double p,double t_ih,double *sa_freeze,double *ct_freeze,double *w_ih)
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

/**************************************************************************
==========================================================================
  method: gsw_internal_energy_ice (t, p)
==========================================================================
   Calculates the specific internal energy of ice.
   t  :  in-situ temperature (ITS-90)              [ deg C ]
   p  =  sea pressure                     [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   internal_energy_ice  =  specific internal energy (u)         [J/kg]
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_internal_energy_ice(double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (gsw_gibbs_ice(0,0,t,p)
			  - (gtc.gsw_t0 + t)*gsw_gibbs_ice(1,0,t,p)
			  - (gtc.db2pa*p + gtc.gsw_p0)*gsw_gibbs_ice(0,1,t,p));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_kappa_const_t_ice(double t,double p)
{
	return (-gsw_gibbs_ice(0,2,t,p)/gsw_gibbs_ice(0,1,t,p));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_kappa_ice(double t,double p)
{
	double gi_tp, gi_tt;
	gi_tt = gsw_gibbs_ice(2,0,t,p);
	gi_tp = gsw_gibbs_ice(1,1,t,p);
	return ((gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))/
			  (gsw_gibbs_ice(0,1,t,p)*gi_tt));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_melting_ice_equilibrium_sa_ct_ratio(double sa,double p)
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_melting_ice_sa_ct_ratio(double sa,double ct,double p,double t_ih)
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_melting_seaice_equilibrium_sa_ct_ratio(double sa,double p)
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(double sa,double p)
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

/**************************************************************************
==========================================================================
elemental subroutine gsw_melting_seaice_into_seawater (sa, ct, p, &
          w_seaice, sa_seaice, t_seaice, sa_final, ct_final)
==========================================================================
   Calculates the Absolute Salinity and Conservative Temperature that
   results when a given mass of sea ice (or ice) melts and is mixed into a
   known mass of seawater (whose properties are (SA,CT,p)).
   If the ice contains no salt (e.g. if it is of glacial origin), then the
   input SA_seaice should be set to zero.
   Ice formed at the sea surface (sea ice) typically contains between 2 g/kg
   and 12 g/kg of salt (defined as the mass of salt divided by the mass of
   ice Ih plus brine) and this programme returns NaNs if the input
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
*************************************************************************/

void TeosBase::gsw_melting_seaice_into_seawater(double sa,double ct,double p,double w_seaice,double sa_seaice,double t_seaice,double *sa_final,double *ct_final)
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

/**************************************************************************
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
  ice Ih plus brine) and this programme returns NaNs if the input
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
*************************************************************************/

double TeosBase::gsw_melting_seaice_sa_ct_ratio_poly(double sa,double ct,double p,double sa_seaice,double t_seaice)
{
	double  ctf, delsa, h, h_brine, h_ih, sa_brine, ct_brine,
			  tf_sa_seaice, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;
	if (sa_seaice < 0.0 || sa_seaice > 15.0)
	{
		return (NAN);
	}
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	if (ct < ctf)
	{
		/**the seawater ct input is below the freezing temp*/
		return (cppGSW_INVALID_VALUE);
	}
	tf_sa_seaice = gsw_t_freezing_poly(sa_seaice,p,saturation_fraction) - 1e-6;
	if (t_seaice > tf_sa_seaice)
	{
		/**t_seaice exceeds the freezing sa*/
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_pot_enthalpy_from_pt_ice_poly(double pt0_ice)
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
	for (iteration = 1; iteration <= 5; iteration++)
	{
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_pressure_coefficient_ice(double t,double p)
{
	return (-gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(0,2,t,p));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_pressure_freezing_ct(double sa,double ct,double saturation_fraction)
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
	if (ct > ct_freezing_p0)
	{
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
	for (i_iter = 0; i_iter < number_of_iterations; i_iter++)
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_pt_from_pot_enthalpy_ice_poly_dh(double pot_enthalpy_ice)
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_pt_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice)
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_pt_from_t_ice(double t,double p,double p_ref)
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_rho_ice(double t,double p)
{
	return (1.0/gsw_gibbs_ice(0,1,t,p));
}

/**************************************************************************
==========================================================================
  method: gsw_sa_freezing_estimate (p, saturation_fraction,*ct,*t)
==========================================================================
  From an estimate of SA from a polynomial in CT and p
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_sa_freezing_estimate(double p,double saturation_fraction,double *ct,double *t)
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
		sa = gsw_max_d(-(*ct + 9e-4*p)/0.06, 0.0);
		ctx = *ct;
	}
	else if (t != NULL)
	{
		sa = gsw_max_d(-(*t + 9e-4*p)/0.06, 0.0);
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_sa_freezing_from_ct(double ct,double p,double saturation_fraction)
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
	sa = gsw_max_d(sa,0.0);
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_sa_freezing_from_t(double t,double p,double saturation_fraction)
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
	sa = gsw_max_d(sa,0.0);
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

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_sa_freezing_from_t_poly(double t,double p,double saturation_fraction)
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
	sa = gsw_max_d(sa,0.0);
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

/**************************************************************************
==========================================================================
  method: gsw_sa_p_inrange (sa, p)
==========================================================================
   Check for any values that are out of the TEOS-10 range ...
   SA    :  Absolute Salinity               [ g/kg ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
---------------------------------------------------------------------------
*************************************************************************/

int TeosBase::gsw_sa_p_inrange(double sa,double p)
{
	if (p > 10000.0 || sa > 120.0 ||
			(p + sa*71.428571428571402) > 13571.42857142857)
	{
		return (0);
	}
	return (1);
}

/**************************************************************************
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
*************************************************************************/

void TeosBase::gsw_seaice_fraction_to_freeze_seawater(double sa,double ct,double p,double sa_seaice,double t_seaice,double *sa_freeze,double *ct_freeze,double *w_seaice)
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
}

/**************************************************************************
==========================================================================
  method: gsw_sound_speed_ice (t, p)
==========================================================================
   Calculates the compression speed of sound in ice.
   t   :  in-situ temperature (ITS-90)             [ deg C ]
   p  :  sea pressure                    [ dbar ]
     ( i.e. absolute pressure - 10.1325 dbar )
   sound_speed_ice  :  compression speed of sound in ice       [ m/s ]
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_sound_speed_ice(double t,double p)
{
	double gi_tp, gi_tt;
	gi_tt = gsw_gibbs_ice(2,0,t,p);
	gi_tp = gsw_gibbs_ice(1,1,t,p);
	return (gsw_gibbs_ice(0,1,t,p) *
			  sqrt(gi_tt/(gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_t_from_pt0_ice(double pt0_ice,double p)
{
	double  p0 = 0.0;
	return (gsw_pt_from_t_ice(pt0_ice,p0,p));
}

/**************************************************************************
==========================================================================
method: gsw_adiabatic_lapse_rate_from_ct(sa,ct,p)
==========================================================================
  Calculates the adiabatic lapse rate from Conservative Temperature
  sa     : Absolute Salinity                                 [g/kg]
  ct     : Conservative Temperature                          [deg C]
  p      : sea pressure                                      [dbar]
  gsw_adiabatic_lapse_rate_from_ct : adiabatic lapse rate    [K/Pa]
*************************************************************************/

double TeosBase::gsw_adiabatic_lapse_rate_from_ct(double sa,double ct,double p)
{
	int n0=0, n1=1, n2=2;
	double pt0, pr0=0.0, t;
	pt0 = gsw_pt_from_ct(sa,ct);
	t = gsw_pt_from_t(sa,pt0,pr0,p);
	return (-gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n2,n0,sa,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_alpha_on_beta(sa,ct,p)
==========================================================================
   Calculates alpha divided by beta, where alpha is the thermal expansion
   coefficient and beta is the saline contraction coefficient of seawater
   from Absolute Salinity and Conservative Temperature.  This method uses
   the computationally-efficient expression for specific volume in terms of
   SA, CT and p (Roquet et al., 2014).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature                        [deg C]
  p      : sea pressure                                    [dbar]
  alpha_on_beta
         : thermal expansion coefficient with respect to   [kg g^-1 K^-1]
           Conservative Temperature divided by the saline
           contraction coefficient at constant Conservative
           Temperature
*************************************************************************/

double TeosBase::gsw_alpha_on_beta(double sa,double ct,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double xs, ys, z, v_ct_part, v_sa_part;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v_ct_part = gsw_gsvco_a(xs, ys, z);
	v_sa_part = gsw_gsvco_b(xs, ys, z);
	return (-(v_ct_part*xs)/(20.0*gtc.gsw_sfac*v_sa_part));
}

/**************************************************************************
==========================================================================
method: gsw_alpha_wrt_t_exact(sa,t,p)
==========================================================================
  Calculates thermal expansion coefficient of seawater with respect to
  in-situ temperature
  sa     : Absolute Salinity                               [g/kg]
  t      : insitu temperature                              [deg C]
  p      : sea pressure                                    [dbar]
  gsw_alpha_wrt_t_exact : thermal expansion coefficient    [1/K]
                          wrt (in-situ) temperature
*************************************************************************/

double TeosBase::gsw_alpha_wrt_t_exact(double sa,double t,double p)
{
	int n0=0, n1=1;
	return (gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_beta_const_t_exact(sa,t,p)
==========================================================================
  Calculates saline (haline) contraction coefficient of seawater at
  constant in-situ temperature.
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  beta_const_t_exact : haline contraction coefficient      [kg/g]
*************************************************************************/

double TeosBase::gsw_beta_const_t_exact(double sa,double t,double p)
{
	int n0=0, n1=1;
	return (-gsw_gibbs(n1,n0,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_c_from_sp(sp,t,p)
==========================================================================
   Calculates conductivity, C, from (SP,t,p) using PSS-78 in the range
   2 < SP < 42.  If the input Practical Salinity is less than 2 then a
   modified form of the Hill et al. (1986) fomula is used for Practical
   Salinity.  The modification of the Hill et al. (1986) expression is to
   ensure that it is exactly consistent with PSS-78 at SP = 2.
   The conductivity ratio returned by this method is consistent with the
   input value of Practical Salinity, SP, to 2x10^-14 psu over the full
   range of input parameters (from pure fresh water up to SP = 42 psu).
   This error of 2x10^-14 psu is machine precision at typical seawater
   salinities.  This accuracy is achieved by having four different
   polynomials for the starting value of Rtx (the square root of Rt) in
   four different ranges of SP, and by using one and a half iterations of
   a computationally efficient modified Newton-Raphson technique (McDougall
   and Wotherspoon, 2012) to find the root of the equation.
   Note that strictly speaking PSS-78 (Unesco, 1983) defines Practical
   Salinity in terms of the conductivity ratio, R, without actually
   specifying the value of C(35,15,0) (which we currently take to be
   42.9140 mS/cm).
  sp     : Practical Salinity                               [unitless]
  t      : in-situ temperature [ITS-90]                     [deg C]
  p      : sea pressure                                     [dbar]
  c      : conductivity                                     [ mS/cm ]
*************************************************************************/

double TeosBase::gsw_c_from_sp(double sp,double t,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SP_COEFFICIENTS use gspc */
	double  p0 = 4.577801212923119e-3,      p1 = 1.924049429136640e-1,
			  p2 = 2.183871685127932e-5,      p3 = -7.292156330457999e-3,
			  p4 = 1.568129536470258e-4,      p5 = -1.478995271680869e-6,
			  p6 = 9.086442524716395e-4,      p7 = -1.949560839540487e-5,
			  p8 = -3.223058111118377e-6,     p9 = 1.175871639741131e-7,
			  p10 = -7.522895856600089e-5,    p11 = -2.254458513439107e-6,
			  p12 = 6.179992190192848e-7,     p13 = 1.005054226996868e-8,
			  p14 = -1.923745566122602e-9,    p15 = 2.259550611212616e-6,
			  p16 = 1.631749165091437e-7,     p17 = -5.931857989915256e-9,
			  p18 = -4.693392029005252e-9,    p19 = 2.571854839274148e-10,
			  p20 = 4.198786822861038e-12,
			  q0 = 5.540896868127855e-5,      q1 = 2.015419291097848e-1,
			  q2 = -1.445310045430192e-5,     q3 = -1.567047628411722e-2,
			  q4 = 2.464756294660119e-4,      q5 = -2.575458304732166e-7,
			  q6 = 5.071449842454419e-3,      q7 = 9.081985795339206e-5,
			  q8 = -3.635420818812898e-6,     q9 = 2.249490528450555e-8,
			  q10 = -1.143810377431888e-3,    q11 = 2.066112484281530e-5,
			  q12 = 7.482907137737503e-7,     q13 = 4.019321577844724e-8,
			  q14 = -5.755568141370501e-10,   q15 = 1.120748754429459e-4,
			  q16 = -2.420274029674485e-6,    q17 = -4.774829347564670e-8,
			  q18 = -4.279037686797859e-9,    q19 = -2.045829202713288e-10,
			  q20 = 5.025109163112005e-12,
			  s0 = 3.432285006604888e-3,      s1 = 1.672940491817403e-1,
			  s2 = 2.640304401023995e-5,      s3 = 1.082267090441036e-1,
			  s4 = -6.296778883666940e-5,     s5 = -4.542775152303671e-7,
			  s6 = -1.859711038699727e-1,     s7 = 7.659006320303959e-4,
			  s8 = -4.794661268817618e-7,     s9 = 8.093368602891911e-9,
			  s10 = 1.001140606840692e-1,     s11 = -1.038712945546608e-3,
			  s12 = -6.227915160991074e-6,    s13 = 2.798564479737090e-8,
			  s14 = -1.343623657549961e-10,   s15 = 1.024345179842964e-2,
			  s16 = 4.981135430579384e-4,     s17 = 4.466087528793912e-6,
			  s18 = 1.960872795577774e-8,     s19 = -2.723159418888634e-10,
			  s20 = 1.122200786423241e-12,
			  u0 = 5.180529787390576e-3,      u1 = 1.052097167201052e-3,
			  u2 = 3.666193708310848e-5,      u3 = 7.112223828976632e0,
			  u4 = -3.631366777096209e-4,     u5 = -7.336295318742821e-7,
			  u6 = -1.576886793288888e+2,     u7 = -1.840239113483083e-3,
			  u8 = 8.624279120240952e-6,      u9 = 1.233529799729501e-8,
			  u10 = 1.826482800939545e+3,     u11 = 1.633903983457674e-1,
			  u12 = -9.201096427222349e-5,    u13 = -9.187900959754842e-8,
			  u14 = -1.442010369809705e-10,   u15 = -8.542357182595853e+3,
			  u16 = -1.408635241899082e0,     u17 = 1.660164829963661e-4,
			  u18 = 6.797409608973845e-7,     u19 = 3.345074990451475e-10,
			  u20 = 8.285687652694768e-13;
	double t68, ft68, x, rtx=0.0, dsp_drtx, sqrty, part1, part2, hill_ratio,
								sp_est, rtx_old, rt, aa, bb, cc, dd, ee, ra,r, rt_lc, rtxm,
								sp_hill_raw;
	t68 = t*1.00024e0;
	ft68 = (t68 - 15e0)/(1e0 + gspc.k*(t68 - 15e0));
	x = sqrt(sp);
	/**
	-------------------------------------------------------------------------
	   Finding the starting value of Rtx, the square root of Rt, using four
	   different polynomials of SP and t68.
	--------------------------------------------------------------------------
	*/
	if (sp >= 9.0)
	{
		rtx = p0
				+ x*(p1 + p4*t68 + x*(p3 + p7*t68 + x*(p6
											 + p11*t68 + x*(p10 + p16*t68 + x*p15))))
				+ t68*(p2+ t68*(p5 + x*x*(p12 + x*p17) + p8*x
									 + t68*(p9 + x*(p13 + x*p18)+ t68*(p14 + p19*x
											  + p20*t68))));
	}
	else if (sp >= 0.25 && sp < 9.0)
	{
		rtx = q0
				+ x*(q1 + q4*t68 + x*(q3 + q7*t68 + x*(q6
											 + q11*t68 + x*(q10 + q16*t68 + x*q15))))
				+ t68*(q2+ t68*(q5 + x*x*(q12 + x*q17) + q8*x
									 + t68*(q9 + x*(q13 + x*q18)+ t68*(q14 + q19*x
											  + q20*t68))));
	}
	else if (sp >= 0.003 && sp < 0.25)
	{
		rtx = s0
				+ x*(s1 + s4*t68 + x*(s3 + s7*t68 + x*(s6
											 + s11*t68 + x*(s10 + s16*t68 + x*s15))))
				+ t68*(s2+ t68*(s5 + x*x*(s12 + x*s17) + s8*x
									 + t68*(s9 + x*(s13 + x*s18)+ t68*(s14 + s19*x
											  + s20*t68))));
	}
	else if (sp < 0.003)
	{
		rtx = u0
				+ x*(u1 + u4*t68 + x*(u3 + u7*t68 + x*(u6
											 + u11*t68 + x*(u10 + u16*t68 + x*u15))))
				+ t68*(u2+ t68*(u5 + x*x*(u12 + x*u17) + u8*x
									 + t68*(u9 + x*(u13 + x*u18)+ t68*(u14 + u19*x
											  + u20*t68))));
	}
	/**
	   -------------------------------------------------------------------------
	   Finding the starting value of dSP_dRtx, the derivative of SP with respect
	   to Rtx.
	   -------------------------------------------------------------------------
	*/
	dsp_drtx =  gspc.a1
					+ (2e0*gspc.a2 + (3e0*gspc.a3
											+ (4e0*gspc.a4 + 5e0*gspc.a5*rtx)*rtx)*rtx)*rtx
					+ ft68*(gspc.b1 + (2e0*gspc.b2 + (3e0*gspc.b3 + (4e0*gspc.b4
											 + 5e0*gspc.b5*rtx)*rtx)*rtx)*rtx);
	if (sp < 2.0)
	{
		x = 400e0*(rtx*rtx);
		sqrty = 10.0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		dsp_drtx = dsp_drtx
					  + gspc.a0*800e0*rtx*(1.5e0 + 2e0*x)/(part1*part1)
					  + gspc.b0*ft68*(10e0 + sqrty*(20e0
											+ 30e0*sqrty))/(part2*part2);
		dsp_drtx = hill_ratio*dsp_drtx;
	}
	/**
	---------------------------------------------------------------------------
	   One iteration through the modified Newton-Raphson method (McDougall and
	   Wotherspoon, 2012) achieves an error in Practical Salinity of about
	   10^-12 for all combinations of the inputs.  One and a half iterations of
	   the modified Newton-Raphson method achevies a maximum error in terms of
	   Practical Salinity of better than 2x10^-14 everywhere.
	   We recommend one and a half iterations of the modified Newton-Raphson
	   method.
	   Begin the modified Newton-Raphson method.
	----------------------------------------------------------------------------
	*/
	sp_est = gspc.a0
				+ (gspc.a1 + (gspc.a2 + (gspc.a3
												 + (gspc.a4 + gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx
				+ ft68*(gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3
														+ (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);
	if (sp_est <  2.0)
	{
		x = 400e0*(rtx*rtx);
		sqrty = 10e0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		sp_hill_raw = sp_est - gspc.a0/part1 - gspc.b0*ft68/part2;
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		sp_est = hill_ratio*sp_hill_raw;
	}
	rtx_old = rtx;
	rtx = rtx_old - (sp_est - sp)/dsp_drtx;
	rtxm = 0.5e0*(rtx + rtx_old); /** This mean value of Rtx, Rtxm, is the
                  value of Rtx at which the derivative dSP_dRtx is evaluated. */
	dsp_drtx = gspc.a1
				  + (2e0*gspc.a2 + (3e0*gspc.a3 + (4e0*gspc.a4
										  + 5e0*gspc.a5*rtxm)*rtxm)*rtxm)*rtxm
				  + ft68*(gspc.b1 + (2e0*gspc.b2 + (3e0*gspc.b3
											+ (4e0*gspc.b4 + 5e0*gspc.b5*rtxm)*rtxm)*rtxm)*rtxm);
	if (sp_est <  2.0)
	{
		x = 400e0*(rtxm*rtxm);
		sqrty = 10e0*rtxm;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		dsp_drtx = dsp_drtx
					  + gspc.a0*800e0*rtxm*(1.5e0 + 2e0*x)/(part1*part1)
					  + gspc.b0*ft68*(10e0 + sqrty*(20e0
											+ 30e0*sqrty))/(part2*part2);
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		dsp_drtx = hill_ratio*dsp_drtx;
	}
	/**
	------------------------------------------------------------------------
	   The line below is where Rtx is updated at the end of the one full
	   iteration of the modified Newton-Raphson technique.
	-------------------------------------------------------------------------
	*/
	rtx = rtx_old - (sp_est - sp)/dsp_drtx;
	/**
	------------------------------------------------------------------------
	   Now we do another half iteration of the modified Newton-Raphson
	   technique, making a total of one and a half modified N-R iterations.
	-------------------------------------------------------------------------
	*/
	sp_est = gspc.a0
				+ (gspc.a1 + (gspc.a2 + (gspc.a3
												 + (gspc.a4 + gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx
				+ ft68*(gspc.b0 + (gspc.b1 + (gspc.b2 + (gspc.b3
														+ (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);
	if (sp_est <  2.0)
	{
		x = 400e0*(rtx*rtx);
		sqrty = 10e0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		sp_hill_raw = sp_est - gspc.a0/part1 - gspc.b0*ft68/part2;
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		sp_est = hill_ratio*sp_hill_raw;
	}
	rtx = rtx - (sp_est - sp)/dsp_drtx;
	/**
	----------------------------------------------------------------------------
	   Now go from Rtx to Rt and then to the conductivity ratio R at pressure p.
	----------------------------------------------------------------------------
	*/
	rt = rtx*rtx;
	aa = gspc.d3 + gspc.d4*t68;
	bb = 1e0 + t68*(gspc.d1 + gspc.d2*t68);
	cc = p*(gspc.e1 + p*(gspc.e2 + gspc.e3*p));
	/** rt_lc (i.e. rt_lower_case) corresponds to rt as defined in
	   the UNESCO 44 (1983) routines. */
	rt_lc = gspc.c0
			  + (gspc.c1 + (gspc.c2 + (gspc.c3
												+ gspc.c4*t68)*t68)*t68)*t68;
	dd = bb - aa*rt_lc*rt;
	ee = rt_lc*rt*aa*(bb + cc);
	ra = sqrt(dd*dd + 4e0*ee) - dd;
	r = 0.5e0*ra/aa;
	/**
	 The dimensionless conductivity ratio, R, is the conductivity input, C,
	   divided by the present estimate of C(SP=35, t_68=15, p=0) which is
	   42.9140 mS/cm (=4.29140 S/m^).
	*/
	return (gtc.gsw_c3515*r);
}

/**************************************************************************
==========================================================================
method: gsw_cp_t_exact(sa,t,p)
==========================================================================
  Calculates isobaric heat capacity of seawater
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  gsw_cp_t_exact : heat capacity                           [J/(kg K)]
*************************************************************************/

double TeosBase::gsw_cp_t_exact(double sa,double t,double p)
{
	int n0=0, n2=2;
	return (-(t+273.15e0)*gsw_gibbs(n0,n2,n0,sa,t,p));
}

/**************************************************************************
=========================================================================
method: gsw_ct_from_entropy (sa, entropy)
=========================================================================
   Calculates Conservative Temperature with entropy as an input variable.
   SA       =  Absolute Salinity                                   [ g/kg ]
   entropy  =  specific entropy                                   [ deg C ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_ct_from_entropy(double sa,double entropy)
{
	double pt;
	pt = gsw_pt_from_entropy(sa,entropy);
	return (gsw_ct_from_pt(sa,pt));
}

/**************************************************************************
==========================================================================
method: gsw_ct_first_derivatives (sa, pt, ct_sa, ct_pt)
==========================================================================
   Calculates the following two derivatives of Conservative Temperature
   (1) CT_SA, the derivative with respect to Absolute Salinity at
       constant potential temperature (with pr = 0 dbar), and
    2) CT_pt, the derivative with respect to potential temperature
       (the regular potential temperature which is referenced to 0 dbar)
       at constant Absolute Salinity.
   SA    :  Absolute Salinity                                        [ g/kg ]
   pt  =  potential temperature (ITS-90)                          [ deg C ]
          (whose reference pressure is 0 dbar)
   CT_SA  =  The derivative of Conservative Temperature with respect to
             Absolute Salinity at constant potential temperature
             (the regular potential temperature which has reference
             sea pressure of 0 dbar).
             The CT_SA output has units of:                     [ K/(g/kg)]
   CT_pt  =  The derivative of Conservative Temperature with respect to
             potential temperature (the regular one with pr = 0 dbar)
             at constant SA. CT_pt is dimensionless.           [ unitless ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_ct_first_derivatives(double sa,double pt,double *ct_sa,double *ct_pt)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double abs_pt, g_sa_mod, g_sa_t_mod, x, y_pt;
	abs_pt = gtc.gsw_t0 + pt ;
	if (ct_pt != NULL)
		*ct_pt = -(abs_pt*gsw_gibbs_pt0_pt0(sa,pt))/gtc.gsw_cp0;
	if (ct_sa == NULL)
		return;
	x = sqrt(gtc.gsw_sfac*sa);
	y_pt = 0.025*pt;
	g_sa_t_mod = 1187.3715515697959
					 + x*(-1480.222530425046
							+ x*(2175.341332000392 + x*(-980.14153344888
									+ 220.542973797483*x) + y_pt*(-548.4580073635929
											+ y_pt*(592.4012338275047 + y_pt*(-274.2361238716608
													  + 49.9394019139016*y_pt)))) + y_pt*(-258.3988055868252
															  + y_pt*(-90.2046337756875 + y_pt*10.50720794170734)))
					 + y_pt*(3520.125411988816 + y_pt*(-1351.605895580406
								+ y_pt*(731.4083582010072 + y_pt*(-216.60324087531103
										  + 25.56203650166196*y_pt))));
	g_sa_t_mod = 0.5*gtc.gsw_sfac*0.025*g_sa_t_mod;
	g_sa_mod = 8645.36753595126
				  + x*(-7296.43987145382
						 + x*(8103.20462414788 + y_pt*(2175.341332000392
								 + y_pt*(-274.2290036817964 + y_pt*(197.4670779425016
											+ y_pt*(-68.5590309679152 + 9.98788038278032*y_pt))))
								+ x*(-5458.34205214835 - 980.14153344888*y_pt
									  + x*(2247.60742726704 - 340.1237483177863*x
											 + 220.542973797483*y_pt))) + y_pt*(-1480.222530425046
													 + y_pt*(-129.1994027934126 + y_pt*(-30.0682112585625
															 + y_pt*(2.626801985426835 ))))) + y_pt*(1187.3715515697959
																	 + y_pt*(1760.062705994408 + y_pt*(-450.535298526802
																			 + y_pt*(182.8520895502518 + y_pt*(-43.3206481750622
																					 + 4.26033941694366*y_pt)))));
	g_sa_mod = 0.5*gtc.gsw_sfac*g_sa_mod;
	*ct_sa = (g_sa_mod - abs_pt*g_sa_t_mod)/gtc.gsw_cp0;
}

/**************************************************************************
===========================================================================
method: gsw_enthalpy_second_derivatives (sa, ct, p, h_sa_sa, &
                                                      h_sa_ct, h_ct_ct)
===========================================================================
   Calculates the following three second-order derivatives of specific
   enthalpy (h),using the computationally-efficient expression for
   specific volume in terms of SA, CT and p (Roquet et al., 2014).
    (1) h_SA_SA, second-order derivative with respect to Absolute Salinity
        at constant CT & p.
    (2) h_SA_CT, second-order derivative with respect to SA & CT at
        constant p.
    (3) h_CT_CT, second-order derivative with respect to CT at constant SA
        and p.
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   p     :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   h_SA_SA  :  The second derivative of specific enthalpy with respect to
               Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
   h_SA_CT  :  The second derivative of specific enthalpy with respect to
               SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
   h_CT_CT  :  The second derivative of specific enthalpy with respect to
               CT at constant SA and p.                      [ J/(kg K^2) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_enthalpy_second_derivatives(double sa,double ct,double p,double *h_sa_sa,double *h_sa_ct,double *h_ct_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double dynamic_h_ct_ct_part, dynamic_h_sa_ct_part, dynamic_h_sa_sa_part, xs, xs2, ys, z;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	if (h_sa_sa != NULL)
	{
		xs2 = pow(xs,2);
		dynamic_h_sa_sa_part =
			z*(-gsvco.h101 + xs2*(3.0*gsvco.h301 + xs*(8.0*gsvco.h401
										 + xs*(15.0*gsvco.h501 + 24.0*gsvco.h601*xs))) + ys*(- gsvco.h111
												 + xs2*(3.0*gsvco.h311 + xs*(8.0*gsvco.h411 + 15.0*gsvco.h511*xs)) + ys*(-gsvco.h121
														 + xs2*(3.0*gsvco.h321 + 8.0*gsvco.h421*xs) + ys*(-gsvco.h131 + 3.0*gsvco.h331*xs2
																 + ys*(-gsvco.h141 - gsvco.h151*ys)))) + z*(-gsvco.h102 + xs2*(3.0*gsvco.h302
																		 + xs*(8.0*gsvco.h402 + 15.0*gsvco.h502*xs)) + ys*(-gsvco.h112 + xs2*(3.0*gsvco.h312
																				 + 8.0*gsvco.h412*xs) + ys*(-gsvco.h122 + 3.0*gsvco.h322*xs2 + ys*(-gsvco.h132
																						 - gsvco.h142*ys ))) + z*(xs2*(8.0*gsvco.h403*xs + 3.0*gsvco.h313*ys)
																								 + z*(-gsvco.h103 + 3.0*gsvco.h303*xs2 + ys*(-gsvco.h113 + ys*(-gsvco.h123 - gsvco.h133*ys))
																										 + z*(-gsvco.h104 - gsvco.h114*ys - gsvco.h105*z)))));
		*h_sa_sa = 1e8*0.25*gtc.gsw_sfac*gtc.gsw_sfac*dynamic_h_sa_sa_part/pow(xs,3);
	}
	if (h_sa_ct != NULL)
	{
		dynamic_h_sa_ct_part =
			z*(gsvco.h111 + xs*(2.0*gsvco.h211 + xs*(3.0*gsvco.h311
									  + xs*(4.0*gsvco.h411 + 5.0*gsvco.h511*xs))) + ys*(2.0*gsvco.h121
											  + xs*(4.0*gsvco.h221 + xs*(6.0*gsvco.h321 + 8.0*gsvco.h421*xs))
											  + ys*(3.0*gsvco.h131 + xs*(6.0*gsvco.h231 + 9.0*gsvco.h331*xs)
													  + ys*(4.0*gsvco.h141 + 8.0*gsvco.h241*xs + 5.0*gsvco.h151*ys ))) + z*(gsvco.h112
															  + xs*(2.0*gsvco.h212 + xs*(3.0*gsvco.h312 + 4.0*gsvco.h412*xs))
															  + ys*(2.0*gsvco.h122 + xs*(4.0*gsvco.h222 + 6.0*gsvco.h322*xs)
																	  + ys*(3.0*gsvco.h132 + 6.0*gsvco.h232*xs + 4.0*gsvco.h142*ys)) + z*(gsvco.h113
																			  + xs*(2.0*gsvco.h213 + 3.0*gsvco.h313*xs) + ys*(2.0*gsvco.h123
																					  + 4.0*gsvco.h223*xs + 3.0*gsvco.h133*ys) + gsvco.h114*z)));
		*h_sa_ct = 1e8*0.025*0.5*gtc.gsw_sfac*dynamic_h_sa_ct_part/xs;
	}
	if (h_ct_ct != NULL)
	{
		dynamic_h_ct_ct_part =
			z*(2.0*gsvco.h021 + xs*(2.0*gsvco.h121 + xs*(2.0*gsvco.h221
											+ xs*(2.0*gsvco.h321 + 2.0*gsvco.h421*xs))) + ys*(6.0*gsvco.h031
													+ xs*(6.0*gsvco.h131 + xs*(6.0*gsvco.h231 + 6.0*gsvco.h331*xs))
													+ ys*(12.0*gsvco.h041 + xs*(12.0*gsvco.h141 + 12.0*gsvco.h241*xs)
															+ ys*(20.0*gsvco.h051 + 20.0*gsvco.h151*xs + 30.0*gsvco.h061*ys)))
				+ z*(2.0*gsvco.h022 + xs*(2.0*gsvco.h122 + xs*(2.0*gsvco.h222
												  + 2.0*gsvco.h322*xs)) + ys*(6.0*gsvco.h032 + xs*(6.0*gsvco.h132
														  + 6.0*gsvco.h232*xs) + ys*(12.0*gsvco.h042 + 12.0*gsvco.h142*xs
																  + 20.0*gsvco.h052*ys)) + z*(2.0*gsvco.h023 + xs*(2.0*gsvco.h123
																		  + 2.0*gsvco.h223*xs) + ys*(6.0*gsvco.h133*xs + 6.0*gsvco.h033
																				  + 12.0*gsvco.h043*ys) + 2.0*gsvco.h024*z)));
		*h_ct_ct = 1e8*6.25e-4*dynamic_h_ct_ct_part;
	}
}

/**************************************************************************
==========================================================================
method: gsw_entropy_from_pt (sa, pt)
==========================================================================
   Calculates specific entropy of seawater.
   SA    :  Absolute Salinity                                        [ g/kg ]
   pt  =  potential temperature (ITS-90)                          [ deg C ]
   entropy  =  specific entropy                                [ J/(kg*K) ]
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_entropy_from_pt(double sa,double pt)
{
	int n0 = 0, n1 = 1;
	double pr0 = 0.0;
	return (-gsw_gibbs(n0,n1,n0,sa,pt,pr0));
}

/**************************************************************************
==========================================================================
method: gsw_entropy_from_t(sa,t,p)
==========================================================================
  Calculates the specific entropy of seawater
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  gsw_entropy_from_t : specific entropy                    [J/(kg K)]
*************************************************************************/

double TeosBase::gsw_entropy_from_t(double sa,double t,double p)
{
	int n0=0, n1=1;
	return (-gsw_gibbs(n0,n1,n0,sa,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_latentheat_evap_ct(sa,ct)
==========================================================================
  Calculates latent heat, or enthalpy, of evaporation.
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature                        [deg C]
  latentheat_evaporation : latent heat of evaporation      [J/kg]
*************************************************************************/

double TeosBase::gsw_latentheat_evap_ct(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  c0  =  2.499065844825125e6, c1  = -1.544590633515099e-1,
			  c2  = -9.096800915831875e4, c3  =  1.665513670736000e2,
			  c4  =  4.589984751248335e1, c5  =  1.894281502222415e1,
			  c6  =  1.192559661490269e3, c7  = -6.631757848479068e3,
			  c8  = -1.104989199195898e2, c9  = -1.207006482532330e3,
			  c10 = -3.148710097513822e3, c11 =  7.437431482069087e2,
			  c12 =  2.519335841663499e3, c13 =  1.186568375570869e1,
			  c14 =  5.731307337366114e2, c15 =  1.213387273240204e3,
			  c16 =  1.062383995581363e3, c17 = -6.399956483223386e2,
			  c18 = -1.541083032068263e3, c19 =  8.460780175632090e1,
			  c20 = -3.233571307223379e2, c21 = -2.031538422351553e2,
			  c22 =  4.351585544019463e1, c23 = -8.062279018001309e2,
			  c24 =  7.510134932437941e2, c25 =  1.797443329095446e2,
			  c26 = -2.389853928747630e1, c27 =  1.021046205356775e2;
	double x, y;
	x = sqrt(gtc.gsw_sfac*sa);
	y = ct/40.0;
	return (c0 + x*(c1 + c4*y + x*(c3
											 + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))
													 + x*(c10 + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))
			  + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)
													+ y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y))))));
}

/**************************************************************************
==========================================================================
method: gsw_latentheat_evap_t(sa,t)
==========================================================================
  Calculates latent heat, or enthalpy, of evaporation.
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  gsw_latentheat_evap_t : latent heat of evaporation       [J/kg]
*************************************************************************/

double TeosBase::gsw_latentheat_evap_t(double sa,double t)
{
	double  ct = gsw_ct_from_pt(sa,t);
	return (gsw_latentheat_evap_ct(sa,ct));
}

/**************************************************************************
==========================================================================
method: gsw_pot_rho_t_exact(sa,t,p,p_ref)
==========================================================================
  Calculates the potential density of seawater
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  p_ref  : reference sea pressure                          [dbar]
  gsw_pot_rho_t_exact : potential density                  [kg/m^3]
*************************************************************************/

double TeosBase::gsw_pot_rho_t_exact(double sa,double t,double p,double p_ref)
{
	double pt = gsw_pt_from_t(sa,t,p,p_ref);
	return (gsw_rho_t_exact(sa,pt,p_ref));
}

/**************************************************************************
=========================================================================
method: gsw_pt_from_entropy (sa, entropy)
=========================================================================
   Calculates potential temperature with reference pressure p_ref = 0 dbar
   and with entropy as an input variable.
   SA       =  Absolute Salinity                                   [ g/kg ]
   entropy  =  specific entropy                                   [ deg C ]
   pt   =  potential temperature                                  [ deg C ]
           with reference sea pressure (p_ref) = 0 dbar.
   Note. The reference sea pressure of the output, pt, is zero dbar.
---------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_pt_from_entropy(double sa,double entropy)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int number_of_iterations;
	double c, dentropy, dentropy_dt, ent_sa, part1, part2, pt, ptm, pt_old;
	/**Find the initial value of pt*/
	part1 = 1.0 - sa/gtc.gsw_sso;
	part2 = 1.0 - 0.05*part1;
	ent_sa = (gtc.gsw_cp0/gtc.gsw_t0)*part1*(1.0 - 1.01*part1);
	c = (entropy - ent_sa)*(part2/gtc.gsw_cp0);
	pt = gtc.gsw_t0*(exp(c) - 1.0);
	dentropy_dt = gtc.gsw_cp0/((gtc.gsw_t0 + pt)*part2);
	for (number_of_iterations = 1; number_of_iterations <= 2;
			number_of_iterations++)
	{
		pt_old = pt;
		dentropy = gsw_entropy_from_pt(sa,pt_old) - entropy;
		pt = pt_old - dentropy/dentropy_dt;
		ptm = 0.5*(pt + pt_old);
		dentropy_dt = -gsw_gibbs_pt0_pt0(sa,ptm);
		pt = pt_old - dentropy/dentropy_dt;
	}
	/**
	  Maximum error of 2.2x10^-6 degrees C for one iteration.
	  Maximum error is 1.4x10^-14 degrees C for two iterations
	  (two iterations is the default, "for Number_of_iterations = 1:2").
	*/
	return (pt);
}

/**************************************************************************
==========================================================================
method: gsw_rho_first_derivatives_wrt_enthalpy(sa,ct,p,*rho_sa,*rho_h)
=========================================================================
   Calculates two first-order derivatives of specific volume (v).
   Note that this method uses the using the computationally-efficient
   expression for specific volume (Roquet et al., 2014).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   p     :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   rho_SA =  The first derivative of rho with respect to
               Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
   rho_h  =  The first derivative of rho with respect to
               SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_rho_first_derivatives_wrt_enthalpy(double sa,double ct,double p,double *rho_sa,double *rho_h)
{
	double rec_v2, v_h=0.0, v_sa;

	if ((rho_sa != NULL) && (rho_h != NULL))
	{
		gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,&v_h);
	}
	else if (rho_sa != NULL)
	{
		gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,NULL);
	}
	else if (rho_h != NULL)
	{
		gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,NULL,&v_h);
	}

	rec_v2 = pow(1.0/gsw_specvol(sa,ct,p), 2);

	if (rho_sa != NULL) *rho_sa = -v_sa*rec_v2;

	if (rho_h != NULL) *rho_h = -v_h*rec_v2;
}

/**************************************************************************
==========================================================================
method: gsw_rho_t_exact(sa,t,p)
==========================================================================
  Calculates in-situ density of seawater from Absolute Salinity and
  in-situ temperature.
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  gsw_rho_t_exact : in-situ density                        [kg/m^3]
*************************************************************************/

double TeosBase::gsw_rho_t_exact(double sa,double t,double p)
{
	int n0=0, n1=1;
	return (1.0/gsw_gibbs(n0,n0,n1,sa,t,p));
}

/**************************************************************************
==========================================================================
method: gsw_sa_from_rho(rho,ct,p)
==========================================================================
   Calculates the Absolute Salinity of a seawater sample, for given values
   of its density, Conservative Temperature and sea pressure (in dbar).
   rho :  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
    Note. This input has not had 1000 kg/m^3 subtracted from it.
      That is, it is density, not density anomaly.
   ct  :  Conservative Temperature (ITS-90)                      [ deg C ]
   p   :  sea pressure                                           [ dbar ]
   sa  :  Absolute Salinity                                      [g/kg]
*************************************************************************/

double TeosBase::gsw_sa_from_rho(double rho,double ct,double p)
{
	int no_iter;
	double sa, v_lab, v_0, v_50, v_sa, sa_old, delta_v, sa_mean;
	v_lab = 1.0/rho;
	v_0 = gsw_specvol(0.0,ct,p);
	v_50 = gsw_specvol(50.0,ct,p);
	sa = 50.0*(v_lab - v_0)/(v_50 - v_0);
	if (sa < 0.0 || sa > 50.0)
		return (cppGSW_INVALID_VALUE);
	v_sa = (v_50 - v_0)/50.0;
	for (no_iter=1; no_iter <= 2; no_iter++)
	{
		sa_old = sa;
		delta_v = gsw_specvol(sa_old,ct,p) - v_lab;
		sa = sa_old - delta_v/v_sa;
		sa_mean = 0.5*(sa + sa_old);
		gsw_specvol_first_derivatives(sa_mean,ct,p,&v_sa,NULL,NULL);
		sa = sa_old - delta_v/v_sa;
		if (sa < 0.0 || sa > 50.0)
			return (cppGSW_INVALID_VALUE);
	}
	return (sa);
}

/**************************************************************************
==========================================================================
method: gsw_sigma0(sa,ct)
==========================================================================
   Calculates potential density anomaly with reference pressure of 0 dbar,
   this being this particular potential density minus 1000 kg/m^3.  This
   method has inputs of Absolute Salinity and Conservative Temperature.
   This method uses the computationally-efficient 48-term expression for
   density in terms of SA, CT and p (IOC et al., 2010).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature                        [deg C]
  gsw_sigma0  : potential density anomaly with reference pressure of 0
                                                       (48 term equation)
*************************************************************************/

double TeosBase::gsw_sigma0(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double vp0, xs, ys;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	vp0 = gsvco.v000
			+ xs*(gsvco.v010 + xs*(gsvco.v020 + xs*(gsvco.v030
										  + xs*(gsvco.v040 + xs*(gsvco.v050 + gsvco.v060*xs)))))
			+ ys*(gsvco.v100 + xs*(gsvco.v110 + xs*(gsvco.v120
										  + xs*(gsvco.v130 + xs*(gsvco.v140 + gsvco.v150*xs))))
					+ ys*(gsvco.v200 + xs*(gsvco.v210 + xs*(gsvco.v220
												  + xs*(gsvco.v230 + gsvco.v240*xs))) + ys*(gsvco.v300
														  + xs*(gsvco.v310 + xs*(gsvco.v320 + gsvco.v330*xs))
														  + ys*(gsvco.v400 + xs*(gsvco.v410 + gsvco.v420*xs)
																  + ys*(gsvco.v500 + gsvco.v510*xs + gsvco.v600*ys)))));
	return (1.0/vp0 - 1000.0);
}

/**************************************************************************
==========================================================================
method: gsw_sound_speed(sa,ct,p)
==========================================================================
   Calculates the speed of sound in seawater.  This method has inputs of
   Absolute Salinity and Conservative Temperature.  This method uses the
   computationally-efficient expression for specific volume in terms of SA,
   CT and p (Roquet et al., 2014).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
  sound_speed  : speed of sound in seawater                [m/s]
*************************************************************************/

double TeosBase::gsw_sound_speed(double sa,double ct,double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double v, v_p, xs, ys, z;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	v = gsw_gsvco_v(xs, ys, z);
	v_p = gsw_gsvco_c(xs, ys, z);
	return (10000.0*sqrt(-v*v/v_p));
}

/**************************************************************************
==========================================================================
method: gsw_sound_speed_t_exact(sa,t,p)
==========================================================================
  Calculates the speed of sound in seawater
  sa     : Absolute Salinity                               [g/kg]
  t      : in-situ temperature                             [deg C]
  p      : sea pressure                                    [dbar]
  gsw_sound_speed_t_exact : sound speed                    [m/s]
*************************************************************************/

double TeosBase::gsw_sound_speed_t_exact(double sa,double t,double p)
{
	int n0=0, n1=1, n2=2;
	double g_tt, g_tp;
	g_tt = gsw_gibbs(n0,n2,n0,sa,t,p);
	g_tp = gsw_gibbs(n0,n1,n1,sa,t,p);
	return (gsw_gibbs(n0,n0,n1,sa,t,p) *
			  sqrt(g_tt/(g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p))));
}

/**************************************************************************
==========================================================================
method: gsw_sp_from_sk(sk)
==========================================================================
  Calculates Practical Salinity, SP, from SK
   SK    : Knudsen Salinity                        [parts per thousand, ppt]
  gsw_sp_from_sk  : Practical Salinity                              [unitless]
*************************************************************************/

double TeosBase::gsw_sp_from_sk(double sk)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double gsw_sp_from_sk_value = (sk - 0.03e0)*(gtc.gsw_soncl/1.805e0);
	if (gsw_sp_from_sk_value < 0e0)
		gsw_sp_from_sk_value = cppGSW_INVALID_VALUE;
	return (gsw_sp_from_sk_value);
}

/**************************************************************************
==========================================================================
method: gsw_sp_from_sr(sr)
==========================================================================
  Calculates Practical Salinity, sp, from Reference Salinity, sr.
  sr     : Reference Salinity                              [g/kg]
  gsw_sp_from_sr  : Practical Salinity                     [unitless]
*************************************************************************/

double TeosBase::gsw_sp_from_sr(double sr)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return(sr/gtc.gsw_ups);
}

/**************************************************************************
==========================================================================
method: gsw_sp_salinometer(rt,t)
==========================================================================
   Calculates Practical Salinity SP from a salinometer, primarily using the
   PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical Salinity
   is only valid in the range 2 < SP < 42.  If the PSS-78 algorithm
   produces a Practical Salinity that is less than 2 then the Practical
   Salinity is recalculated with a modified form of the Hill et al. (1986)
   formula.  The modification of the Hill et al. (1986) expression is to
   ensure that it is exactly consistent with PSS-78 at SP = 2.
   A laboratory salinometer has the ratio of conductivities, Rt, as an
   output, and the present method uses this conductivity ratio and the
   temperature t of the salinometer bath as the two input variables.
   rt  = C(SP,t_68,0)/C(SP=35,t_68,0)                          [ unitless ]
   t   = temperature of the bath of the salinometer,
         measured on the ITS-90 scale (ITS-90)                 [ deg C ]
   gsw_sp_salinometer = Practical Salinity on the PSS-78 scale [ unitless ]
*************************************************************************/

double TeosBase::gsw_sp_salinometer(double rt,double t)
{
	/** for GSW_SP_COEFFICIENTS use gspc */
	double t68, ft68, rtx, sp, hill_ratio, x, sqrty, part1, part2, sp_hill_raw;
	if (rt < 0)
	{
		return NAN;
	}
	t68 = t*1.00024;
	ft68 = (t68 - 15)/(1 + gspc.k*(t68 - 15));
	rtx = sqrt(rt);
	sp = gspc.a0
		  + (gspc.a1 + (gspc.a2 + (gspc.a3
											+ (gspc.a4 + gspc.a5 * rtx) * rtx) * rtx) * rtx) * rtx
		  + ft68 * (gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3
													 + (gspc.b4 + gspc.b5 * rtx) * rtx) * rtx) * rtx) * rtx);
	/**
	   The following section of the code is designed for SP < 2 based on the
	   Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
	   exactly equal to the PSS-78 algorithm at SP = 2.
	*/
	if (sp < 2)
	{
		hill_ratio = gsw_hill_ratio_at_sp2(t);
		x = 400e0*rt;
		sqrty = 10e0*rtx;
		part1 = 1e0 + x*(1.5e0 + x);
		part2 = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
		sp_hill_raw = sp - gspc.a0/part1 - gspc.b0*ft68/part2;
		sp = hill_ratio*sp_hill_raw;
	}
	return sp;
}

/**************************************************************************
==========================================================================
method: gsw_specvol_anom_standard(sa,ct,p)
==========================================================================
   Calculates specific volume anomaly of seawater.
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
  specvol_anom  :  specific volume anomaly of seawater
*************************************************************************/

double TeosBase::gsw_specvol_anom_standard(double sa,double ct,double p)
{
	return (gsw_specvol(sa,ct,p) - gsw_specvol_sso_0(p));
}

/**************************************************************************
==========================================================================
method: gsw_specvol_first_derivatives (sa,ct,p,*v_sa,*v_ct,*v_p)
=========================================================================
   Calculates three first-order derivatives of specific volume (v).
   Note that this method uses the computationally-efficient
   expression for specific volume (Roquet et al., 2014).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   p     :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   v_SA  =  The first derivative of specific volume with respect to
            Absolute Salinity at constant CT & p.       [ J/(kg (g/kg)^2) ]
   v_CT  =  The first derivative of specific volume with respect to
            CT at constant SA and p.                     [ J/(kg K(g/kg)) ]
   v_P   =  The first derivative of specific volume with respect to
            P at constant SA and CT.                         [ J/(kg K^2) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_specvol_first_derivatives(double sa,double ct,double p,double *v_sa,double *v_ct,double *v_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double  v_ct_part, v_p_part, v_sa_part, xs, ys, z;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	z = p*1e-4;
	if (v_sa != NULL)
	{
		v_sa_part = gsw_gsvco_b(xs, ys, z);
		*v_sa = 0.5*gtc.gsw_sfac*v_sa_part/xs;
	}
	if (v_ct != NULL)
	{
		v_ct_part =gsw_gsvco_a(xs, ys, z);
		*v_ct = 0.025*v_ct_part;
	}
	if (v_p != NULL)
	{
		v_p_part = gsw_gsvco_c(xs, ys, z);
		*v_p = 1e-8*v_p_part;
	}
}

/**************************************************************************
===========================================================================
method: gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,*v_sa,*v_h)
===========================================================================
   Calculates two first-order derivatives of specific volume (v).
   Note that this method uses the using the computationally-efficient
   expression for specific volume (Roquet et al., 2014).
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   v_SA  :  The first derivative of specific volume with respect to
               Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
   v_h  :  The first derivative of specific volume with respect to
               SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_specvol_first_derivatives_wrt_enthalpy(double sa,double ct,double p,double *v_sa,double *v_h)
{
	double h_ct=1.0, h_sa, rec_h_ct, vct_ct, vct_sa;
	if (v_sa != NULL)
	{
		gsw_specvol_first_derivatives(sa,ct,p,&vct_sa,&vct_ct,NULL);
		gsw_enthalpy_first_derivatives(sa,ct,p,&h_sa,&h_ct);
	}
	else if (v_h != NULL)
	{
		gsw_specvol_first_derivatives(sa,ct,p,NULL,&vct_ct,NULL);
		gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);
	}
	rec_h_ct = 1.0/h_ct;
	if (v_sa != NULL)
		*v_sa = vct_sa - (vct_ct*h_sa)*rec_h_ct;
	if (v_h != NULL)
		*v_h = vct_ct*rec_h_ct;
	return;
}

/**************************************************************************
==========================================================================
method: gsw_specvol_second_derivatives (sa, ct, p, v_sa_sa, &
                                   v_sa_ct, v_ct_ct, v_sa_p, v_ct_p, iflag)
=========================================================================
   Calculates five second-order derivatives of specific volume (v).
   Note that this method uses the computationally-efficient
   expression for specific volume (Roquet et al., 2014).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   p     :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   v_SA_SA  =  The second derivative of specific volume with respect to
               Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
   v_SA_CT  =  The second derivative of specific volume with respect to
               SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
   v_CT_CT  =  The second derivative of specific volume with respect to
               CT at constant SA and p.                      [ J/(kg K^2) ]
   v_SA_P  =  The second derivative of specific volume with respect to
               SA and P at constant CT.                  [ J/(kg K(g/kg)) ]
   v_CT_P  =  The second derivative of specific volume with respect to
               CT and P at constant SA.                  [ J/(kg K(g/kg)) ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_specvol_second_derivatives(double sa,double ct,double p,
		double *v_sa_sa,double *v_sa_ct,double *v_ct_ct,double *v_sa_p,double *v_ct_p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double v_ct_ct_part, v_ct_p_part, v_sa_ct_part, v_sa_p_part,
			 v_sa_sa_part, xs, xs2, ys, z;
	xs2 = gtc.gsw_sfac*sa + gtc.offset;
	xs = sqrt(xs2);
	ys = ct*0.025;
	z = p*1e-4;
	if (v_sa_sa != NULL)
	{
		v_sa_sa_part = (-gsvco.b000
							 + xs2*(gsvco.b200 + xs*(2.0*gsvco.b300 + xs*(3.0*gsvco.b400
															 + 4.0*gsvco.b500*xs))) + ys*(-gsvco.b010 + xs2*(gsvco.b210 + xs*(2.0*gsvco.b310
																	 + 3.0*gsvco.b410*xs)) + ys*(-gsvco.b020 + xs2*(gsvco.b220 + 2.0*gsvco.b320*xs)
																			 + ys*(-gsvco.b030 + gsvco.b230*xs2 + ys*(-gsvco.b040 - gsvco.b050*ys)))) + z*(-gsvco.b001
																					 + xs2*(gsvco.b201 + xs*(2.0*gsvco.b301 + 3.0*gsvco.b401*xs)) + ys*(-gsvco.b011
																							 + xs2*(gsvco.b211 + 2.0*gsvco.b311*xs) + ys*(-gsvco.b021 + gsvco.b221*xs2
																									 + ys*(-gsvco.b031 - gsvco.b041*ys))) + z*(-gsvco.b002 + xs2*(gsvco.b202 + 2.0*gsvco.b302*xs)
																											 + ys*(-gsvco.b012 + gsvco.b212*xs2 + ys*(-gsvco.b022 - gsvco.b032*ys)) + z*(-gsvco.b003
																													 - gsvco.b013*ys - gsvco.b004*z))))/xs2;

		*v_sa_sa = 0.25*gtc.gsw_sfac*gtc.gsw_sfac*v_sa_sa_part/xs;
	}
	if (v_sa_ct != NULL)
	{
		v_sa_ct_part = (gsvco.b010
							 + xs*(gsvco.b110 + xs*(gsvco.b210 + xs*(gsvco.b310 + gsvco.b410*xs)))
							 + ys*(2.0*(gsvco.b020 + xs*(gsvco.b120 + xs*(gsvco.b220 + gsvco.b320*xs)))
									 + ys*(3.0*(gsvco.b030 + xs*(gsvco.b130 + gsvco.b230*xs)) + ys*(4.0*(gsvco.b040
											 + gsvco.b140*xs) + 5.0*gsvco.b050*ys))) + z*(gsvco.b011 + xs*(gsvco.b111 + xs*(gsvco.b211
													 + gsvco.b311*xs)) + ys*(2.0*(gsvco.b021 + xs*(gsvco.b121 + gsvco.b221*xs))
															 + ys*(3.0*(gsvco.b031 + gsvco.b131*xs) + 4.0*gsvco.b041*ys)) + z*(gsvco.b012
																	 + xs*(gsvco.b112 + gsvco.b212*xs) + ys*(2.0*(gsvco.b022 + gsvco.b122*xs)
																			 + 3.0*gsvco.b032*ys) + gsvco.b013*z)))/xs;

		*v_sa_ct = 0.025*0.5*gtc.gsw_sfac*v_sa_ct_part;
	}
	if (v_ct_ct != NULL)
	{
		v_ct_ct_part = gsvco.a010
							+ xs*(gsvco.a110 + xs*(gsvco.a210 + xs*(gsvco.a310 + gsvco.a410*xs)))
							+ ys*(2.0*(gsvco.a020 + xs*(gsvco.a120 + xs*(gsvco.a220 + gsvco.a320*xs)))
									+ ys*(3.0*(gsvco.a030 + xs*(gsvco.a130 + gsvco.a230*xs)) + ys*(4.0*(gsvco.a040
											+ gsvco.a140*xs) + 5.0*gsvco.a050*ys))) + z*(gsvco.a011 + xs*(gsvco.a111 + xs*(gsvco.a211
													+ gsvco.a311*xs)) + ys*(2.0*(gsvco.a021 + xs*(gsvco.a121 + gsvco.a221*xs))
															+ ys*(3.0*(gsvco.a031 + gsvco.a131*xs) + 4.0*gsvco.a041*ys)) + z*(gsvco.a012
																	+ xs*(gsvco.a112 + gsvco.a212*xs) + ys*(2.0*(gsvco.a022 + gsvco.a122*xs)
																			+ 3.0*gsvco.a032*ys) + gsvco.a013*z));

		*v_ct_ct = 0.025*0.025*v_ct_ct_part;
	}
	if (v_sa_p != NULL)
	{
		v_sa_p_part = gsvco.b001
						  + xs*(gsvco.b101 + xs*(gsvco.b201 + xs*(gsvco.b301
														 + gsvco.b401*xs))) + ys*(gsvco.b011 + xs*(gsvco.b111 + xs*(gsvco.b211
																 + gsvco.b311*xs)) + ys*(gsvco.b021 + xs*(gsvco.b121 + gsvco.b221*xs)
																		 + ys*(gsvco.b031 + gsvco.b131*xs + gsvco.b041*ys))) + z*(2.0*(gsvco.b002 + xs*(gsvco.b102
																				 + xs*(gsvco.b202 + gsvco.b302*xs)) + ys*(gsvco.b012 + xs*(gsvco.b112
																						 + gsvco.b212*xs) + ys*(gsvco.b022
																								 + gsvco.b122*xs + gsvco.b032*ys))) + z*(3.0*(gsvco.b003 + gsvco.b103*xs + gsvco.b013*ys)
																										 + 4.0*gsvco.b004*z));

		*v_sa_p = 1e-8*0.5*gtc.gsw_sfac*v_sa_p_part;
	}
	if (v_ct_p != NULL)
	{
		v_ct_p_part = gsvco.a001
						  + xs*(gsvco.a101 + xs*(gsvco.a201 + xs*(gsvco.a301
														 + gsvco.a401*xs))) + ys*(gsvco.a011
																 + xs*(gsvco.a111 + xs*(gsvco.a211 + gsvco.a311*xs)) + ys*(gsvco.a021
																		 + xs*(gsvco.a121 + gsvco.a221*xs)
																		 + ys*(gsvco.a031 + gsvco.a131*xs + gsvco.a041*ys))) + z*(2.0*(gsvco.a002 + xs*(gsvco.a102
																				 + xs*(gsvco.a202 + gsvco.a302*xs)) + ys*(gsvco.a012 + xs*(gsvco.a112 + gsvco.a212*xs)
																						 + ys*(gsvco.a022 + gsvco.a122*xs + gsvco.a032*ys))) + z*(3.0*(gsvco.a003
																								 + gsvco.a103*xs + gsvco.a013*ys) + 4.0*gsvco.a004*z));

		*v_ct_p = 1e-8*0.025*v_ct_p_part;
	}
}

/**************************************************************************
==========================================================================
method: gsw_specvol_second_derivatives_wrt_enthalpy (sa, ct, &
                                          p, v_sa_sa, v_sa_h, v_h_h, iflag)
=========================================================================
   Calculates three first-order derivatives of specific volume (v) with
   respect to enthalpy. Note that this method uses the using the
   computationally-efficient expression for specific volume
   (Roquet et al., 2014).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   p     :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   v_SA_SA = The second-order derivative of specific volume with respect to
             Absolute Salinity at constant h & p.       [ J/(kg (g/kg)^2) ]
   v_SA_h  = The second-order derivative of specific volume with respect to
             SA and h at constant p.                     [ J/(kg K(g/kg)) ]
   v_h_h   = The second-order derivative with respect to h at
             constant SA & p.
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_specvol_second_derivatives_wrt_enthalpy(double sa,double ct,double p,
		double *v_sa_sa,double *v_sa_h,double *v_h_h)
{
	double h_ct, h_ct_ct, h_sa, h_sa_ct, h_sa_sa, rec_h_ct, v_h_h_part,
			  rec_h_ct2, v_ct, vct_ct_ct, vct_sa_ct, vct_sa_sa, v_sa_h_part;

	gsw_specvol_first_derivatives(sa,ct,p,NULL, &v_ct, NULL);

	if ((v_sa_sa != NULL) || (v_sa_h != NULL))
		gsw_enthalpy_first_derivatives(sa,ct,p,&h_sa,&h_ct);
	else
		gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);

	if (v_sa_sa != NULL)
		gsw_specvol_second_derivatives(sa,ct,p,&vct_sa_sa,&vct_sa_ct,
												 &vct_ct_ct, NULL, NULL);
	else if (v_sa_h != NULL)
		gsw_specvol_second_derivatives(sa,ct,p,NULL,&vct_sa_ct,&vct_ct_ct,
												 NULL, NULL);
	else
		gsw_specvol_second_derivatives(sa,ct,p,NULL,NULL,&vct_ct_ct,
												 NULL, NULL);

	if (v_sa_sa != NULL)
		gsw_enthalpy_second_derivatives(sa,ct,p,&h_sa_sa,&h_sa_ct,&h_ct_ct);
	else if (v_sa_h != NULL)
		gsw_enthalpy_second_derivatives(sa,ct,p,NULL,&h_sa_ct,&h_ct_ct);
	else
		gsw_enthalpy_second_derivatives(sa,ct,p,NULL,NULL,&h_ct_ct);

	rec_h_ct = 1.0/h_ct;
	rec_h_ct2 = rec_h_ct*rec_h_ct;

	v_h_h_part = (vct_ct_ct*h_ct - h_ct_ct*v_ct)*(rec_h_ct2*rec_h_ct);

	if (v_h_h != NULL) *v_h_h = v_h_h_part;

	if ((v_sa_sa != NULL) || (v_sa_h != NULL))
	{
		v_sa_h_part = (vct_sa_ct*h_ct - v_ct*h_sa_ct)*rec_h_ct2
						  - h_sa*v_h_h_part;

		if (v_sa_h != NULL) *v_sa_h = v_sa_h_part;

		if (v_sa_sa != NULL)
			*v_sa_sa = vct_sa_sa - (h_ct*(vct_sa_ct*h_sa
            - v_ct*h_sa_sa) + v_ct*h_sa*h_sa_ct)*rec_h_ct2
            - h_sa*v_sa_h_part;
	}
}

/**
==========================================================================
method: gsw_spiciness0(sa,ct)
==========================================================================
   Calculates spiciness from Absolute Salinity and Conservative
   Temperature at a pressure of 0 dbar, as described by McDougall and
   Krzysik (2015).  This routine is based on the computationally-efficient
   expression for specific volume in terms of SA, CT and p (Roquet et al.,
   2015).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                        [ deg C ]
   spiciness0  =  spiciness referenced to a pressure of 0 dbar,
                  i.e. the surface                                 [ kg/m^3 ]
*************************************************************************/

double TeosBase::gsw_spiciness0(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  s01 = -9.22982898371678e1,      s02 = -1.35727873628866e1,
			  s03 =  1.87353650994010e1,      s04 = -1.61360047373455e1,
			  s05 =  3.76112762286425e1,      s06 = -4.27086671461257e1,
			  s07 =  2.00820111041594e1,      s08 =  2.87969717584045e2,
			  s09 =  1.13747111959674e1,      s10 =  6.07377192990680e1,
			  s11 = -7.37514033570187e1,      s12 = -7.51171878953574e1,
			  s13 =  1.63310989721504e2,      s14 = -8.83222751638095e1,
			  s15 = -6.41725302237048e2,      s16 =  2.79732530789261e1,
			  s17 = -2.49466901993728e2,      s18 =  3.26691295035416e2,
			  s19 =  2.66389243708181e1,      s20 = -2.93170905757579e2,
			  s21 =  1.76053907144524e2,      s22 =  8.27634318120224e2,
			  s23 = -7.02156220126926e1,      s24 =  3.82973336590803e2,
			  s25 = -5.06206828083959e2,      s26 =  6.69626565169529e1,
			  s27 =  3.02851235050766e2,      s28 = -1.96345285604621e2,
			  s29 = -5.74040806713526e2,      s30 =  7.03285905478333e1,
			  s31 = -2.97870298879716e2,      s32 =  3.88340373735118e2,
			  s33 = -8.29188936089122e1,      s34 = -1.87602137195354e2,
			  s35 =  1.27096944425793e2,      s36 =  2.11671167892147e2,
			  s37 = -3.15140919876285e1,      s38 =  1.16458864953602e2,
			  s39 = -1.50029730802344e2,      s40 =  3.76293848660589e1,
			  s41 =  6.47247424373200e1,      s42 = -4.47159994408867e1,
			  s43 = -3.23533339449055e1,      s44 =  5.30648562097667,
			  s45 = -1.82051249177948e1,      s46 =  2.33184351090495e1,
			  s47 = -6.22909903460368,        s48 = -9.55975464301446,
			  s49 =  6.61877073960113;
	double xs, ys, spiciness0;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	spiciness0 = s01
					 +ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
					 +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
							+xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
								  +xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
										 +xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
												+xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
													  +xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))))))));
	return (spiciness0);
}

/**************************************************************************
==========================================================================
method: gsw_spiciness1(sa,ct)
==========================================================================
   Calculates spiciness from Absolute Salinity and Conservative
   Temperature at a pressure of 1000 dbar, as described by McDougall and
   Krzysik (2015).  This routine is based on the computationally-efficient
   expression for specific volume in terms of SA, CT and p (Roquet et al.,
   2015).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   spiciness1  =  spiciness referenced to a pressure of 1000 dbar [ kg/m^3 ]
*************************************************************************/

double TeosBase::gsw_spiciness1(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  s01 = -9.19874584868912e1,      s02 = -1.33517268529408e1,
			  s03 =  2.18352211648107e1,      s04 = -2.01491744114173e1,
			  s05 =  3.70004204355132e1,      s06 = -3.78831543226261e1,
			  s07 =  1.76337834294554e1,      s08 =  2.87838842773396e2,
			  s09 =  2.14531420554522e1,      s10 =  3.14679705198796e1,
			  s11 = -4.04398864750692e1,      s12 = -7.70796428950487e1,
			  s13 =  1.36783833820955e2,      s14 = -7.36834317044850e1,
			  s15 = -6.41753415180701e2,      s16 =  1.33701981685590,
			  s17 = -1.75289327948412e2,      s18 =  2.42666160657536e2,
			  s19 =  3.17062400799114e1,      s20 = -2.28131490440865e2,
			  s21 =  1.39564245068468e2,      s22 =  8.27747934506435e2,
			  s23 = -3.50901590694775e1,      s24 =  2.87473907262029e2,
			  s25 = -4.00227341144928e2,      s26 =  6.48307189919433e1,
			  s27 =  2.16433334701578e2,      s28 = -1.48273032774305e2,
			  s29 = -5.74545648799754e2,      s30 =  4.50446431127421e1,
			  s31 = -2.30714981343772e2,      s32 =  3.15958389253065e2,
			  s33 = -8.60635313930106e1,      s34 = -1.22978455069097e2,
			  s35 =  9.18287282626261e1,      s36 =  2.12120473062203e2,
			  s37 = -2.21528216973820e1,      s38 =  9.19013417923270e1,
			  s39 = -1.24400776026014e2,      s40 =  4.08512871163839e1,
			  s41 =  3.91127352213516e1,      s42 = -3.10508021853093e1,
			  s43 = -3.24790035899152e1,      s44 =  3.91029016556786,
			  s45 = -1.45362719385412e1,      s46 =  1.96136194246355e1,
			  s47 = -7.06035474689088,        s48 = -5.36884688614009,
			  s49 =  4.43247303092448;
	double xs, ys, spiciness1;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	spiciness1 = s01
					 +ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
					 +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
							+xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
								  +xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
										 +xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
												+xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
													  +xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))))))));
	return (spiciness1);
}

/**************************************************************************
==========================================================================
method: gsw_spiciness2(sa,ct)
==========================================================================
   Calculates spiciness from Absolute Salinity and Conservative
   Temperature at a pressure of 2000 dbar, as described by McDougall and
   Krzysik (2015).  This routine is based on the computationally-efficient
   expression for specific volume in terms of SA, CT and p (Roquet et al.,
   2015).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   spiciness2  =  spiciness referenced to a pressure of 2000 dbar [ kg/m^3 ]
*************************************************************************/

double TeosBase::gsw_spiciness2(double sa,double ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double  s01 = -9.17327320732265e1,      s02 = -1.31200235147912e1,
			  s03 =  2.49574345782503e1,      s04 = -2.41678075247398e1,
			  s05 =  3.61654631402053e1,      s06 = -3.22582164667710e1,
			  s07 =  1.45092623982509e1,      s08 =  2.87776645983195e2,
			  s09 =  3.13902307672447e1,      s10 =  1.69777467534459,
			  s11 = -5.69630115740438,        s12 = -7.97586359017987e1,
			  s13 =  1.07507460387751e2,      s14 = -5.58234404964787e1,
			  s15 = -6.41708068766557e2,      s16 = -2.53494801286161e1,
			  s17 = -9.86755437385364e1,      s18 =  1.52406930795842e2,
			  s19 =  4.23888258264105e1,      s20 = -1.60118811141438e2,
			  s21 =  9.67497898053989e1,      s22 =  8.27674355478637e2,
			  s23 =  5.27561234412133e-1,     s24 =  1.87440206992396e2,
			  s25 = -2.83295392345171e2,      s26 =  5.14485994597635e1,
			  s27 =  1.29975755062696e2,      s28 = -9.36526588377456e1,
			  s29 = -5.74911728972948e2,      s30 =  1.91175851862772e1,
			  s31 = -1.59347231968841e2,      s32 =  2.33884725744938e2,
			  s33 = -7.87744010546157e1,      s34 = -6.04757235443685e1,
			  s35 =  5.27869695599657e1,      s36 =  2.12517758478878e2,
			  s37 = -1.24351794740528e1,      s38 =  6.53904308937490e1,
			  s39 = -9.44804080763788e1,      s40 =  3.93874257887364e1,
			  s41 =  1.49425448888996e1,      s42 = -1.62350721656367e1,
			  s43 = -3.25936844276669e1,      s44 =  2.44035700301595,
			  s45 = -1.05079633683795e1,      s46 =  1.51515796259082e1,
			  s47 = -7.06609886460683,        s48 = -1.48043337052968,
			  s49 =  2.10066653978515;
	double xs, ys, spiciness2;
	xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys = ct*0.025;
	spiciness2 = s01
					 +ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
					 +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
							+xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
								  +xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
										 +xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
												+xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
													  +xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))))))));
	return (spiciness2);
}

/**************************************************************************
==========================================================================
method: gsw_util_indxPref(n,z)
==========================================================================
  calculates k in simulated array
  n     :  length of the array
  z     :  value to be indexed
  K      : index K - if X(K) <= Z < X(K+1), or
  N-1     		    - if Z = X(N)
*************************************************************************/

int TeosBase::gsw_util_indxPref(int n, double z)
{
	int k, ku, kl, km;
	if (z > pSaarData->get_p_ref(0) && z < pSaarData->get_p_ref(n-1))
	{
		kl  = 0;
		ku  = n-1;
		while (ku-kl > 1)
		{
			km  = (ku+kl)>>1;
			if (z > pSaarData->get_p_ref(km))
			{
				kl  = km;
			}
			else
			{
				ku  = km;
			}
		}
		k = kl;
		if (z == pSaarData->get_p_ref(k+1))
		{
			k++;
		}
	}
	else
	{
		if (z <= pSaarData->get_p_ref(0))
		{
			k = 0;
		}
		else
		{
			k = n-2;
		}
	}
	return k;
}

/**************************************
=======
method: gsw_gsvco_a(xs,ys,z)
=======
The code in this method was repeated
in various places in the C code.
In the C++ code it is replaced with
a call to this method.
**************************************/

double TeosBase::gsw_gsvco_a(double xs, double ys, double z)
{
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double rval = gsvco.a000
					  + xs*(gsvco.a100 + xs*(gsvco.a200 + xs*(gsvco.a300 + xs*(gsvco.a400 + gsvco.a500*xs))))
					  + ys*(gsvco.a010 + xs*(gsvco.a110 + xs*(gsvco.a210 + xs*(gsvco.a310 + gsvco.a410*xs)))
							  + ys*(gsvco.a020 + xs*(gsvco.a120 + xs*(gsvco.a220 + gsvco.a320*xs)) + ys*(gsvco.a030
									  + xs*(gsvco.a130 + gsvco.a230*xs) + ys*(gsvco.a040 + gsvco.a140*xs + gsvco.a050*ys ))))
					  + z*(gsvco.a001 + xs*(gsvco.a101 + xs*(gsvco.a201 + xs*(gsvco.a301 + gsvco.a401*xs)))
							 + ys*(gsvco.a011 + xs*(gsvco.a111 + xs*(gsvco.a211 + gsvco.a311*xs)) + ys*(gsvco.a021
									 + xs*(gsvco.a121 + gsvco.a221*xs) + ys*(gsvco.a031 + gsvco.a131*xs + gsvco.a041*ys)))
							 + z*(gsvco.a002 + xs*(gsvco.a102 + xs*(gsvco.a202 + gsvco.a302*xs)) + ys*(gsvco.a012
									 + xs*(gsvco.a112 + gsvco.a212*xs) + ys*(gsvco.a022 + gsvco.a122*xs + gsvco.a032*ys))
									+ z*(gsvco.a003 + gsvco.a103*xs + gsvco.a013*ys + gsvco.a004*z)));
	return rval;
}

/**************************************
=======
method: gsw_gsvco_b(xs,ys,z)
=======
The code in this method was repeated
in various places in the C code.
In the C++ code it is replaced with
a call to this method.
**************************************/

double TeosBase::gsw_gsvco_b(double xs, double ys, double z)
{
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double rval = gsvco.b000
					  + xs*(gsvco.b100 + xs*(gsvco.b200 + xs*(gsvco.b300
													 + xs*(gsvco.b400 + gsvco.b500*xs)))) + ys*(gsvco.b010
															 + xs*(gsvco.b110 + xs*(gsvco.b210 + xs*(gsvco.b310
																	 + gsvco.b410*xs))) + ys*(gsvco.b020 + xs*(gsvco.b120
																			 + xs*(gsvco.b220 + gsvco.b320*xs)) + ys*(gsvco.b030
																					 + xs*(gsvco.b130 + gsvco.b230*xs) + ys*(gsvco.b040
																							 + gsvco.b140*xs + gsvco.b050*ys)))) + z*(gsvco.b001
																									 + xs*(gsvco.b101 + xs*(gsvco.b201 + xs*(gsvco.b301
																											 + gsvco.b401*xs))) + ys*(gsvco.b011 + xs*(gsvco.b111
																													 + xs*(gsvco.b211 + gsvco.b311*xs)) + ys*(gsvco.b021
																															 + xs*(gsvco.b121 + gsvco.b221*xs) + ys*(gsvco.b031
																																	 + gsvco.b131*xs + gsvco.b041*ys))) + z*(gsvco.b002
																																			 + xs*(gsvco.b102 + xs*(gsvco.b202 + gsvco.b302*xs))
																																			 + ys*(gsvco.b012 + xs*(gsvco.b112 + gsvco.b212*xs)
																																					 + ys*(gsvco.b022 + gsvco.b122*xs + gsvco.b032*ys))
																																			 + z*(gsvco.b003 + gsvco.b103*xs + gsvco.b013*ys
																																					 + gsvco.b004*z)));
	return rval;
}

/**************************************
=======
method: gsw_gsvco_c(xs,ys,z)
=======
The code in this method was repeated
in various places in the C code.
In the C++ code it is replaced with
a call to this method.
**************************************/

double TeosBase::gsw_gsvco_c(double xs, double ys, double z)
{
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double rval = gsvco.c000
					  + xs*(gsvco.c100 + xs*(gsvco.c200 + xs*(gsvco.c300
													 + xs*(gsvco.c400 + gsvco.c500*xs)))) + ys*(gsvco.c010
															 + xs*(gsvco.c110 + xs*(gsvco.c210 + xs*(gsvco.c310
																	 + gsvco.c410*xs))) + ys*(gsvco.c020 + xs*(gsvco.c120
																			 + xs*(gsvco.c220 + gsvco.c320*xs)) + ys*(gsvco.c030
																					 + xs*(gsvco.c130 + gsvco.c230*xs) + ys*(gsvco.c040
																							 + gsvco.c140*xs + gsvco.c050*ys)))) + z*(gsvco.c001
																									 + xs*(gsvco.c101 + xs*(gsvco.c201 + xs*(gsvco.c301
																											 + gsvco.c401*xs))) + ys*(gsvco.c011 + xs*(gsvco.c111
																													 + xs*(gsvco.c211 + gsvco.c311*xs)) + ys*(gsvco.c021
																															 + xs*(gsvco.c121 + gsvco.c221*xs) + ys*(gsvco.c031
																																	 + gsvco.c131*xs + gsvco.c041*ys))) + z*(gsvco.c002
																																			 + xs*(gsvco.c102 + gsvco.c202*xs) + ys*(gsvco.c012
																																					 + gsvco.c112*xs + gsvco.c022*ys) + z*(gsvco.c003
																																							 + gsvco.c103*xs + gsvco.c013*ys + z*(gsvco.c004
																																									 + gsvco.c005*z))));
	return rval;
}

/**************************************
=======
method: gsw_gsvco_v(xs,ys,z)
=======
The code in this method was repeated
in various places in the C code.
In the C++ code it is replaced with
a call to this method.
**************************************/

double TeosBase::gsw_gsvco_v(double xs, double ys, double z)
{
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double rval = gsvco.v000
					  + xs*(gsvco.v010 + xs*(gsvco.v020 + xs*(gsvco.v030
													 + xs*(gsvco.v040 + xs*(gsvco.v050 + gsvco.v060*xs)))))
					  + ys*(gsvco.v100 + xs*(gsvco.v110 + xs*(gsvco.v120
													 + xs*(gsvco.v130 + xs*(gsvco.v140 + gsvco.v150*xs))))
							  + ys*(gsvco.v200 + xs*(gsvco.v210 + xs*(gsvco.v220
									  + xs*(gsvco.v230 + gsvco.v240*xs))) + ys*(gsvco.v300
											  + xs*(gsvco.v310 + xs*(gsvco.v320 + gsvco.v330*xs))
											  + ys*(gsvco.v400 + xs*(gsvco.v410 + gsvco.v420*xs)
													  + ys*(gsvco.v500 + gsvco.v510*xs + gsvco.v600*ys)))))
					  + z*(gsvco.v001 + xs*(gsvco.v011 + xs*(gsvco.v021
													+ xs*(gsvco.v031 + xs*(gsvco.v041 + gsvco.v051*xs))))
							 + ys*(gsvco.v101 + xs*(gsvco.v111 + xs*(gsvco.v121
									 + xs*(gsvco.v131 + gsvco.v141*xs))) + ys*(gsvco.v201
											 + xs*(gsvco.v211 + xs*(gsvco.v221 + gsvco.v231*xs))
											 + ys*(gsvco.v301 + xs*(gsvco.v311 + gsvco.v321*xs)
													 + ys*(gsvco.v401 + gsvco.v411*xs + gsvco.v501*ys))))
							 + z*(gsvco.v002 + xs*(gsvco.v012 + xs*(gsvco.v022
									 + xs*(gsvco.v032 + gsvco.v042*xs))) + ys*(gsvco.v102
											 + xs*(gsvco.v112 + xs*(gsvco.v122 + gsvco.v132*xs))
											 + ys*(gsvco.v202 + xs*(gsvco.v212 + gsvco.v222*xs)
													 + ys*(gsvco.v302 + gsvco.v312*xs + gsvco.v402*ys)))
									+ z*(gsvco.v003 + xs*(gsvco.v013 + gsvco.v023*xs)
										  + ys*(gsvco.v103 + gsvco.v113*xs + gsvco.v203*ys)
										  + z*(gsvco.v004 + gsvco.v014*xs + gsvco.v104*ys
												 + z*(gsvco.v005 + gsvco.v006*z)))));
	return rval;
}

/****************************************************************************
	method: gsw_max_d (double, double)
	previously an inline function
*****************************************************************************/

double TeosBase::gsw_max_d(double a,double b)
{
	double _max;

	if (a > b)
	{
		_max = a;
	}
	else
	{
		_max = b;
	}

	return _max;
}

/****************************************************************************
	method: gsw_min_d (double, double)
	previously an inline function
*****************************************************************************/

double TeosBase::gsw_min_d(double a,double b)
{
	double _min;

	if (a > b)
	{
		_min = b;
	}
	else
	{
		_min = a;
	}

	return _min;
}

/****************************************************************************
	method: gsw_max_u (unsigned, unsigned)
	previously an inline function
*****************************************************************************/

unsigned TeosBase::gsw_max_u(unsigned a,unsigned b)
{
	unsigned _max;

	if (a > b)
	{
		_max = a;
	}
	else
	{
		_max = b;
	}

	return _max;
}

/****************************************************************************
	method: gsw_min_u (unsigned, unsigned)
	previously an inline function
*****************************************************************************/

unsigned TeosBase::gsw_min_u(unsigned a,unsigned b)
{
	unsigned _min;

	if (a > b)
	{
		_min = b;
	}
	else
	{
		_min = a;
	}

	return _min;
}

/****************************************************************************
	method: gsw_max_i (int, int)
	previously an inline function
*****************************************************************************/

int TeosBase::gsw_max_i(int a,int b)
{
	int _max;

	if (a > b)
	{
		_max = a;
	}
	else
	{
		_max = b;
	}

	return _max;
}

/****************************************************************************
	method: gsw_min_i (int, int)
	previously an inline function
*****************************************************************************/

int TeosBase::gsw_min_i(int a,int b)
{
	int _min;

	if (a > b)
	{
		_min = b;
	}
	else
	{
		_min = a;
	}

	return _min;
}

/**************************************
=======
method: gsw_get_v_ct_part(xs,ys,z)
=======
This method supplements gsw_alpha ()
**************************************/

double TeosBase::gsw_get_v_ct_part(double xs, double ys, double z)
{
	double rval = gsw_gsvco_a(xs, ys, z);
	return rval;
}

/************************************************************
=======
method: gsw_get_pot_ent(x2,x,y)
=======
This method supplements gsw_ct_from_pt(double sa, double pt)
************************************************************/

double TeosBase::gsw_get_pot_ent(double x2, double x, double y)
{
	double rval = 61.01362420681071e0 + y*(168776.46138048015e0 +
														y*(-2735.2785605119625e0 + y*(2574.2164453821433e0 +
																y*(-1536.6644434977543e0 + y*(545.7340497931629e0 +
																		(-50.91091728474331e0 - 18.30489878927802e0*y)*y))))) +
					  x2*(268.5520265845071e0 + y*(-12019.028203559312e0 +
							  y*(3734.858026725145e0 + y*(-2046.7671145057618e0 +
									  y*(465.28655623826234e0 + (-0.6370820302376359e0 -
											  10.650848542359153e0*y)*y)))) +
							x*(937.2099110620707e0 + y*(588.1802812170108e0+
									y*(248.39476522971285e0 + (-3.871557904936333e0-
											2.6268019854268356e0*y)*y)) +
								x*(-1687.914374187449e0 + x*(246.9598888781377e0 +
										x*(123.59576582457964e0 - 48.5891069025409e0*x)) +
									y*(936.3206544460336e0 +
										y*(-942.7827304544439e0 + y*(369.4389437509002e0 +
												(-33.83664947895248e0 - 9.987880382780322e0*y)*y))))));
	return rval;
}

/************************************************************
=======
method: gsw_get_dyn_enth_prt(xs,ys,z)
=======
This method supplements method gsw_dynamic_enthalpy(sa,ct,p)
************************************************************/

double TeosBase::gsw_get_dyn_enth_prt(double xs, double ys, double z)
{
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double rval = z*(gsvco.h001 + xs*(gsvco.h101 + xs*(gsvco.h201 + xs*(gsvco.h301 + xs*(gsvco.h401
												 + xs*(gsvco.h501 + gsvco.h601*xs))))) + ys*(gsvco.h011 + xs*(gsvco.h111 + xs*(gsvco.h211 + xs*(gsvco.h311
														 + xs*(gsvco.h411 + gsvco.h511*xs)))) + ys*(gsvco.h021 + xs*(gsvco.h121 + xs*(gsvco.h221 + xs*(gsvco.h321
																 + gsvco.h421*xs))) + ys*(gsvco.h031 + xs*(gsvco.h131 + xs*(gsvco.h231 + gsvco.h331*xs)) + ys*(gsvco.h041
																		 + xs*(gsvco.h141 + gsvco.h241*xs) + ys*(gsvco.h051 + gsvco.h151*xs + gsvco.h061*ys))))) + z*(gsvco.h002
																				 + xs*(gsvco.h102 + xs*(gsvco.h202 + xs*(gsvco.h302 + xs*(gsvco.h402 + gsvco.h502*xs)))) + ys*(gsvco.h012
																						 + xs*(gsvco.h112 + xs*(gsvco.h212 + xs*(gsvco.h312 + gsvco.h412*xs))) + ys*(gsvco.h022 + xs*(gsvco.h122
																								 + xs*(gsvco.h222 + gsvco.h322*xs)) + ys*(gsvco.h032 + xs*(gsvco.h132 + gsvco.h232*xs) + ys*(gsvco.h042
																										 + gsvco.h142*xs + gsvco.h052*ys)))) + z*(gsvco.h003 + xs*(gsvco.h103 + xs*(gsvco.h203 + xs*(gsvco.h303
																												 + gsvco.h403*xs))) + ys*(gsvco.h013 + xs*(gsvco.h113 + xs*(gsvco.h213 + gsvco.h313*xs)) + ys*(gsvco.h023
																														 + xs*(gsvco.h123 + gsvco.h223*xs) + ys*(gsvco.h033 + gsvco.h133*xs + gsvco.h043*ys))) + z*(gsvco.h004
																																 + xs*(gsvco.h104 + gsvco.h204*xs) + ys*(gsvco.h014 + gsvco.h114*xs + gsvco.h024*ys) + z*(gsvco.h005
																																		 + gsvco.h105*xs + gsvco.h015*ys + z*(gsvco.h006 + gsvco.h007*z))))));
	return rval;
}

/***********************************************
=======
method: gsw_get_g03plusg08(z,y,x2,x)
=======
This method supplements gsw_entropy_part(sa,t,p)
************************************************/

double TeosBase::gsw_get_g03plusg08(double z, double y, double x2, double x)
{
	double rval;
	double g03 = z*(-270.983805184062 +
						 z*(776.153611613101 + z*(-196.51255088122 +
								 (28.9796526294175 - 2.13290083518327*z)*z))) +
					 y*(-24715.571866078 + z*(2910.0729080936 +
													  z*(-1513.116771538718 + z*(546.959324647056 +
															  z*(-111.1208127634436 + 8.68841343834394*z)))) +
						 y*(2210.2236124548363 + z*(-2017.52334943521 +
								 z*(1498.081172457456 + z*(-718.6359919632359 +
										 (146.4037555781616 - 4.9892131862671505*z)*z))) +
							 y*(-592.743745734632 + z*(1591.873781627888 +
									 z*(-1207.261522487504 + (608.785486935364 -
											 105.4993508931208*z)*z)) +
								 y*(290.12956292128547 + z*(-973.091553087975 +
										 z*(602.603274510125 + z*(-276.361526170076 +
												 32.40953340386105*z))) +
									 y*(-113.90630790850321 + y*(21.35571525415769 -
											 67.41756835751434*z) +
										 z*(381.06836198507096 + z*(-133.7383902842754 +
												 49.023632509086724*z)))))));
	double g08 = x2*(z*(729.116529735046 +
							  z*(-343.956902961561 + z*(124.687671116248 +
									  z*(-31.656964386073 + 7.04658803315449*z)))) +
						  x*( x*(y*(-137.1145018408982 + y*(148.10030845687618 +
										y*(-68.5590309679152 + 12.4848504784754*y))) -
									22.6683558512829*z) + z*(-175.292041186547 +
											(83.1923927801819 - 29.483064349429*z)*z) +
								y*(-86.1329351956084 + z*(766.116132004952 +
										z*(-108.3834525034224 + 51.2796974779828*z)) +
									y*(-30.0682112585625 - 1380.9597954037708*z +
										y*(3.50240264723578 + 938.26075044542*z)))) +
						  y*(1760.062705994408 + y*(-675.802947790203 +
								  y*(365.7041791005036 + y*(-108.30162043765552 +
										  12.78101825083098*y) +
									  z*(-1190.914967948748 + (298.904564555024 -
											  145.9491676006352*z)*z)) +
								  z*(2082.7344423998043 + z*(-614.668925894709 +
										  (340.685093521782 - 33.3848202979239*z)*z))) +
							  z*(-1721.528607567954 + z*(674.819060538734 +
									  z*(-356.629112415276 + (88.4080716616 -
											  15.84003094423364*z)*z)))));
	rval = (-(g03 + g08)*0.025);
	return rval;
}

/**************************************************************************
=======
method: gsw_g08_d(x2,z,y,x)
=======
This method supplements gsw_dilution_coefficient_t_exact
***************************************************************************/

double TeosBase::gsw_g08_d(double x2,double z,double y,double x)
{
	double g08_a = 2.0*(8103.20462414788 +
							  y*(2175.341332000392 +
								  y*(-274.2290036817964 +
									  y*(197.4670779425016 +
										  y*(-68.5590309679152 + 9.98788038278032*y))) -
								  90.6734234051316*z) + 1.5*
							  x*(-5458.34205214835 - 980.14153344888*y +
								  (4.0/3.0)*x*(2247.60742726704 -
													340.1237483177863*1.25*x + 220.542973797483*y) +
								  180.142097805543*z) +
							  z*(-219.1676534131548 +
								  (-16.32775915649044 - 120.7020447884644*z)*z));
	double g08_b = x2*g08_a +
						x*(-7296.43987145382 +
							z*(598.378809221703 +
								z*(-156.8822727844005 +
									(204.1334828179377 - 10.23755797323846*z)*z)) +
							y*(-1480.222530425046 +
								z*(-525.876123559641 +
									(249.57717834054571 - 88.449193048287*z)*z) +
								y*(-129.1994027934126 +
									z*(1149.174198007428 +
										z*(-162.5751787551336 + 76.9195462169742*z)) +
									y*(-30.0682112585625 - 1380.9597954037708*z +
										y*(2.626801985426835 + 703.695562834065*z))))) +
						11625.62913253464 + 1702.453469893412*y;
	return g08_b;
}

/**************************************************************************
=======
method: gsw_get_g03plusg08_h(z,y,x2,x)
=======
This method supplements method gsw_chem_potential_water_t_exact (sa, t, p)
*************************************************************************/

double TeosBase::gsw_get_g03plusg08_h(double z,double y,double x2,double x)
{
	double rval;
	double g03_g = 101.342743139674 + z*(100015.695367145 +
													 z*(-2544.5765420363 + z*(284.517778446287 +
															 z*(-33.3146754253611 + (4.20263108803084 - 0.546428511471039*z)*z)))) +
						y*(5.90578347909402 + z*(-270.983805184062 +
														 z*(776.153611613101 + z*(-196.51255088122 + (28.9796526294175 - 2.13290083518327*z)*z))) +
							y*(-12357.785933039 + z*(1455.0364540468 +
									z*(-756.558385769359 + z*(273.479662323528 +
											z*(-55.5604063817218 + 4.34420671917197*z)))) +
								y*(736.741204151612 + z*(-672.50778314507 +
										z*(499.360390819152 + z*(-239.545330654412 + (48.8012518593872 - 1.66307106208905*z)*z))) +
									y*(-148.185936433658 + z*(397.968445406972 +
											z*(-301.815380621876 + (152.196371733841 - 26.3748377232802*z)*z)) +
										y*(58.0259125842571 + z*(-194.618310617595 +
												z*(120.520654902025 + z*(-55.2723052340152 + 6.48190668077221*z))) +
											y*(-18.9843846514172 + y*(3.05081646487967 - 9.63108119393062*z) +
												z*(63.5113936641785 + z*(-22.2897317140459 + 8.17060541818112*z))))))));
	double g08_g = x2*(1416.27648484197 +
							 x*(-2432.14662381794 + x*(2025.80115603697 +
									 y*(543.835333000098 + y*(-68.5572509204491 +
											 y*(49.3667694856254 + y*(-17.1397577419788 + 2.49697009569508*y))) - 22.6683558512829*z) +
									 x*(-1091.66841042967 - 196.028306689776*y +
										 x*(374.60123787784 - 48.5891069025409*x + 36.7571622995805*y) + 36.0284195611086*z) +
									 z*(-54.7919133532887 + (-4.08193978912261 - 30.1755111971161*z)*z)) +
								 z*(199.459603073901 + z*(-52.2940909281335 + (68.0444942726459 - 3.41251932441282*z)*z)) +
								 y*(-493.407510141682 + z*(-175.292041186547 + (83.1923927801819 - 29.483064349429*z)*z) +
									 y*(-43.0664675978042 + z*(383.058066002476 +
											 z*(-54.1917262517112 + 25.6398487389914*z)) +
										 y*(-10.0227370861875 - 460.319931801257*z +
											 y*(0.875600661808945 + 234.565187611355*z))))) +
							 y*(168.072408311545));
	rval = g03_g + g08_g;
	return rval;
}

/*************************************************************************************
=======
method: gsw_gibbs_g03plusg08_a(double sa, double z, double y, double x2, double x)
=======
Supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g03plusg08_a(double sa, double z, double y, double x2, double x)
{
	double g03, g08, rval;
	g03 = 101.342743139674 +
			z*(100015.695367145 +
				z*(-2544.5765420363 +
					z*(284.517778446287 +
						z*(-33.3146754253611 + (4.20263108803084 -
														0.546428511471039*z)*z)))) +
			y*(5.90578347909402 +
				z*(-270.983805184062 +
					z*(776.153611613101 +
						z*(-196.51255088122 +
							(28.9796526294175 - 2.13290083518327*z)*z))) +
				y*(-12357.785933039 +
					z*(1455.0364540468 +
						z*(-756.558385769359 +
							z*(273.479662323528 +
								z*(-55.5604063817218 + 4.34420671917197*z)))) +
					y*(736.741204151612 +
						z*(-672.50778314507 +
							z*(499.360390819152 +
								z*(-239.545330654412 +
									(48.8012518593872 - 1.66307106208905*z)*z))) +
						y*(-148.185936433658 +
							z*(397.968445406972 +
								z*(-301.815380621876 + (152.196371733841 -
										26.3748377232802*z)*z)) +
							y*(58.0259125842571 +
								z*(-194.618310617595 +
									z*(120.520654902025 +
										z*(-55.2723052340152 +
											6.48190668077221*z))) +
								y*(-18.9843846514172 +
									y*(3.05081646487967 -
										9.63108119393062*z) +
									z*(63.5113936641785 +
										z*(-22.2897317140459 +
											8.17060541818112*z))))))));
	g08 = x2*(1416.27648484197 +
				 z*(-3310.49154044839 +
					 z*(384.794152978599 +
						 z*(-96.5324320107458 +
							 (15.8408172766824 - 2.62480156590992*z)*z))) +
				 x*(-2432.14662381794 +
					 x*(2025.80115603697 +
						 y*(543.835333000098 +
							 y*(-68.5572509204491 +
								 y*(49.3667694856254 +
									 y*(-17.1397577419788 +
										 2.49697009569508*y))) - 22.6683558512829*z) +
						 x*(-1091.66841042967 - 196.028306689776*y +
							 x*(374.60123787784 - 48.5891069025409*x +
								 36.7571622995805*y) + 36.0284195611086*z) +
						 z*(-54.7919133532887 + (-4.08193978912261 -
														 30.1755111971161*z)*z)) +
					 z*(199.459603073901 +
						 z*(-52.2940909281335 +
							 (68.0444942726459 - 3.41251932441282*z)*z)) +
					 y*(-493.407510141682 +
						 z*(-175.292041186547 +
							 (83.1923927801819 - 29.483064349429*z)*z) +
						 y*(-43.0664675978042 +
							 z*(383.058066002476 +
								 z*(-54.1917262517112 + 25.6398487389914*z)) +
							 y*(-10.0227370861875 - 460.319931801257*z +
								 y*(0.875600661808945 + 234.565187611355*z))))) +
				 y*(168.072408311545 +
					 z*(729.116529735046 +
						 z*(-343.956902961561 +
							 z*(124.687671116248 +
								 z*(-31.656964386073 + 7.04658803315449*z)))) +
					 y*(880.031352997204 +
						 y*(-225.267649263401 +
							 y*(91.4260447751259 +
								 y*(-21.6603240875311 +
									 2.13016970847183*y) +
								 z*(-297.728741987187 + (74.726141138756 -
										 36.4872919001588*z)*z)) +
							 z*(694.244814133268 +
								 z*(-204.889641964903 +
									 (113.561697840594 - 11.1282734326413*z)*z))) +
						 z*(-860.764303783977 +
							 z*(337.409530269367 +
								 z*(-178.314556207638 + (44.2040358308 -
										 7.92001547211682*z)*z))))));
	if (sa > 0.0)
	{
		g08 += (x2*(5812.81456626732 + 851.226734946706*y)*log(x));
	}
	rval = g03 + g08;
	return rval;
}

/*************************************************************************************
=======
method: gsw_gibbs_g03plusg08_b(sa,z,y,x2,x)
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g03plusg08_b(double sa,double z,double y,double x2,double x)
{
	double rval;
	double g03 = 5.90578347909402 +
					 z*(-270.983805184062 +
						 z*(776.153611613101 +
							 z*(-196.51255088122 +
								 (28.9796526294175 - 2.13290083518327*z)*z))) +
					 y*(-24715.571866078 +
						 z*(2910.0729080936 +
							 z*(-1513.116771538718 +
								 z*(546.959324647056 +
									 z*(-111.1208127634436 + 8.68841343834394*z)))) +
						 y*(2210.2236124548363 +
							 z*(-2017.52334943521 +
								 z*(1498.081172457456 +
									 z*(-718.6359919632359 +
										 (146.4037555781616 - 4.9892131862671505*z)*z))) +
							 y*(-592.743745734632 +
								 z*(1591.873781627888 +
									 z*(-1207.261522487504 + (608.785486935364 -
											 105.4993508931208*z)*z)) +
								 y*(290.12956292128547 +
									 z*(-973.091553087975 +
										 z*(602.603274510125 +
											 z*(-276.361526170076 + 32.40953340386105*z))) +
									 y*(-113.90630790850321 +
										 y*(21.35571525415769 - 67.41756835751434*z) +
										 z*(381.06836198507096 +
											 z*(-133.7383902842754 + 49.023632509086724*z)))))));
	double g08 = x2*(168.072408311545 +
						  z*(729.116529735046 +
							  z*(-343.956902961561 +
								  z*(124.687671116248 +
									  z*(-31.656964386073 + 7.04658803315449*z)))) +
						  x*(-493.407510141682 +
							  x*(543.835333000098 +
								  x*(-196.028306689776 + 36.7571622995805*x) +
								  y*(-137.1145018408982 +
									  y*(148.10030845687618 +
										  y*(-68.5590309679152 + 12.4848504784754*y))) -
								  22.6683558512829*z) +
							  z*(-175.292041186547 + (83.1923927801819 - 29.483064349429*z)*z) +
							  y*(-86.1329351956084 +
								  z*(766.116132004952 +
									  z*(-108.3834525034224 + 51.2796974779828*z)) +
								  y*(-30.0682112585625 - 1380.9597954037708*z +
									  y*(3.50240264723578 + 938.26075044542*z)))) +
						  y*(1760.062705994408 +
							  y*(-675.802947790203 +
								  y*(365.7041791005036 +
									  y*(-108.30162043765552 + 12.78101825083098*y) +
									  z*(-1190.914967948748 + (298.904564555024 -
											  145.9491676006352*z)*z)) +
								  z*(2082.7344423998043 +
									  z*(-614.668925894709 +
										  (340.685093521782 - 33.3848202979239*z)*z))) +
							  z*(-1721.528607567954 +
								  z*(674.819060538734 +
									  z*(-356.629112415276 + (88.4080716616 - 15.84003094423364*z)*z)))));
	if (sa > 0.0)
	{
		g08 += 851.226734946706*x2*log(x);
	}
	rval = (g03 + g08)*0.025;
	return rval;
}

/***************************************************************************************
=======
method: gsw_gibbs_g03plusg08_c(z,y,x2,x)
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g03plusg08_c(double z, double y, double x2, double x)
{
	double rval;
	double g03 = 100015.695367145 +
					 z*(-5089.1530840726 +
						 z*(853.5533353388611 +
							 z*(-133.2587017014444 +
								 (21.0131554401542 - 3.278571068826234*z)*z))) +
					 y*(-270.983805184062 +
						 z*(1552.307223226202 +
							 z*(-589.53765264366 + (115.91861051767 -
									 10.664504175916349*z)*z)) +
						 y*(1455.0364540468 +
							 z*(-1513.116771538718 +
								 z*(820.438986970584 +
									 z*(-222.2416255268872 + 21.72103359585985*z))) +
							 y*(-672.50778314507 +
								 z*(998.720781638304 +
									 z*(-718.6359919632359 + (195.2050074375488 -
											 8.31535531044525*z)*z)) +
								 y*(397.968445406972 +
									 z*(-603.630761243752 +
										 (456.589115201523 - 105.4993508931208*z)*z) +
									 y*(-194.618310617595 +
										 y*(63.5113936641785 - 9.63108119393062*y +
											 z*(-44.5794634280918 + 24.511816254543362*z)) +
										 z*(241.04130980405 +
											 z*(-165.8169157020456 + 25.92762672308884*z)))))));
	double g08 = x2*(-3310.49154044839 + z*(769.588305957198 +
														 z*(-289.5972960322374 + (63.3632691067296 -
																 13.1240078295496*z)*z)) +
						  x*(199.459603073901 + x*(-54.7919133532887 +
								  36.0284195611086*x - 22.6683558512829*y +
								  (-8.16387957824522 - 90.52653359134831*z)*z) +
							  z*(-104.588181856267 + (204.1334828179377 -
									  13.65007729765128*z)*z) +
							  y*(-175.292041186547 + (166.3847855603638 -
									  88.449193048287*z)*z +
								  y*(383.058066002476 + y*(-460.319931801257 +
										  234.565187611355*y) +
									  z*(-108.3834525034224 + 76.9195462169742*z)))) +
						  y*(729.116529735046 + z*(-687.913805923122 +
								  z*(374.063013348744 + z*(-126.627857544292 +
										  35.23294016577245*z))) +
							  y*(-860.764303783977 + y*(694.244814133268 +
									  y*(-297.728741987187 + (149.452282277512 -
											  109.46187570047641*z)*z) +
									  z*(-409.779283929806 + (340.685093521782 -
											  44.5130937305652*z)*z)) +
								  z*(674.819060538734 + z*(-534.943668622914 +
										  (176.8161433232 - 39.600077360584095*z)*z)))));
	rval = (g03 + g08)*1.0e-8;
	return rval;
}

/***************************************************************************************
=======
method: gsw_gibbs_g03plusg08_d(z,y,x2,x)
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g03plusg08_d(double z,double y,double x2,double x)
{
	double rval;
	double g03 = -24715.571866078 + z*(2910.0729080936 + z*
												  (-1513.116771538718 + z*(546.959324647056 +
														  z*(-111.1208127634436 + 8.68841343834394*z)))) +
					 y*(4420.4472249096725 + z*(-4035.04669887042 +
														 z*(2996.162344914912 + z*(-1437.2719839264719 +
																 (292.8075111563232 - 9.978426372534301*z)*z))) +
						 y*(-1778.231237203896 + z*(4775.621344883664 +
								 z*(-3621.784567462512 + (1826.356460806092 -
										 316.49805267936244*z)*z)) +
							 y*(1160.5182516851419 + z*(-3892.3662123519 +
									 z*(2410.4130980405 + z*(-1105.446104680304 +
											 129.6381336154442*z))) +
								 y*(-569.531539542516 + y*(128.13429152494615 -
										 404.50541014508605*z) +
									 z*(1905.341809925355 + z*(-668.691951421377 +
											 245.11816254543362*z))))));
	double g08 = x2*(1760.062705994408 + x*(-86.1329351956084 +
														 x*(-137.1145018408982 + y*(296.20061691375236 +
																 y*(-205.67709290374563 + 49.9394019139016*y))) +
														 z*(766.116132004952 + z*(-108.3834525034224 +
																 51.2796974779828*z)) +
														 y*(-60.136422517125 - 2761.9195908075417*z +
																 y*(10.50720794170734 + 2814.78225133626*z))) +
						  y*(-1351.605895580406 + y*(1097.1125373015109 +
								  y*(-433.20648175062206 + 63.905091254154904*y) +
								  z*(-3572.7449038462437 + (896.713693665072 -
										  437.84750280190565*z)*z)) +
							  z*(4165.4688847996085 + z*(-1229.337851789418 +
									  (681.370187043564 - 66.7696405958478*z)*z))) +
						  z*(-1721.528607567954 + z*(674.819060538734 +
								  z*(-356.629112415276 + (88.4080716616 -
										  15.84003094423364*z)*z))));
	rval = (g03 + g08)*0.000625;
	return rval;
}

/***************************************************************************************
=======
method: gsw_gibbs_g03plusg08_e(z,y,x2,x)
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g03plusg08_e(double z, double y, double x2, double x)
{
	double rval;
	double g03 = -270.983805184062 + z*(1552.307223226202 +
													z*(-589.53765264366 + (115.91861051767 -
															10.664504175916349*z)*z)) +
					 y*(2910.0729080936 + z*(-3026.233543077436 +
													 z*(1640.877973941168 + z*(-444.4832510537744 +
															 43.4420671917197*z))) +
						 y*(-2017.52334943521 + z*(2996.162344914912 +
								 z*(-2155.907975889708 + (585.6150223126464 -
										 24.946065931335752*z)*z)) +
							 y*(1591.873781627888 + z*(-2414.523044975008 +
									 (1826.356460806092 - 421.9974035724832*z)*z) +
								 y*(-973.091553087975 + z*(1205.20654902025 +
										 z*(-829.084578510228 + 129.6381336154442*z)) +
									 y*(381.06836198507096 - 67.41756835751434*y +
										 z*(-267.4767805685508 + 147.07089752726017*z))))));
	double g08 = x2*(729.116529735046 + z*(-687.913805923122 +
														z*(374.063013348744 + z*(-126.627857544292 +
																35.23294016577245*z))) +
						  x*(-175.292041186547 - 22.6683558512829*x +
							  (166.3847855603638 - 88.449193048287*z)*z +
							  y*(766.116132004952 + y*(-1380.9597954037708 +
									  938.26075044542*y) +
								  z*(-216.7669050068448 + 153.8390924339484*z))) +
						  y*(-1721.528607567954 + y*(2082.7344423998043 +
								  y*(-1190.914967948748 + (597.809129110048 -
										  437.84750280190565*z)*z) +
								  z*(-1229.337851789418 + (1022.055280565346 -
										  133.5392811916956*z)*z)) +
							  z*(1349.638121077468 + z*(-1069.887337245828 +
									  (353.6322866464 - 79.20015472116819*z)*z))));
	rval = (g03 + g08)*2.5e-10;
	return rval;
}

/***************************************************************************************
=======
method: gsw_gibbs_g03plusg08_f(z,y,x2,x);
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g03plusg08_f(double z, double y, double x2, double x)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double rval;
	double g08 = 2.0*(8103.20462414788 +
							y*(2175.341332000392 + y*(-274.2290036817964 +
									y*(197.4670779425016 + y*(-68.5590309679152 +
											9.98788038278032*y))) - 90.6734234051316*z) +
							1.5*x*(-5458.34205214835 - 980.14153344888*y +
									 (4.0/3.0)*x*(2247.60742726704 - 340.1237483177863*1.25*x +
													  220.542973797483*y) + 180.142097805543*z) +
							z*(-219.1676534131548 + (-16.32775915649044 -
									120.7020447884644*z)*z));
	if (x > 0.0)
	{
		g08 += (-7296.43987145382 + z*(598.378809221703 +
												 z*(-156.8822727844005 + (204.1334828179377 -
														 10.23755797323846*z)*z)) +
				  y*(-1480.222530425046 + z*(-525.876123559641 +
													  (249.57717834054571 - 88.449193048287*z)*z) +
					  y*(-129.1994027934126 + z*(1149.174198007428 +
														  z*(-162.5751787551336 + 76.9195462169742*z)) +
						  y*(-30.0682112585625 - 1380.9597954037708*z +
							  y*(2.626801985426835 + 703.695562834065*z)))))/x +
				 (11625.62913253464 + 1702.453469893412*y)/x2;
	}
	else
	{
		g08 = 0.0;
	}
	rval = 0.25*gtc.gsw_sfac*gtc.gsw_sfac*g08;
	return rval;
}

/***************************************************************************************
=======
method: gsw_gibbs_g03plusg08_g(z,y,x2,x);
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g03plusg08_g(double z, double y, double x2, double x)
{
	double rval;
	double g03 = -5089.1530840726 + z*(1707.1066706777221 +
												  z*(-399.7761051043332 + (84.0526217606168 -
														  16.39285534413117*z)*z)) +
					 y*(1552.307223226202 + z*(-1179.07530528732 +
														(347.75583155301 - 42.658016703665396*z)*z) +
						 y*(-1513.116771538718 + z*(1640.877973941168 +
								 z*(-666.7248765806615 + 86.8841343834394*z)) +
							 y*(998.720781638304 + z*(-1437.2719839264719 +
									 (585.6150223126464 - 33.261421241781*z)*z) +
								 y*(-603.630761243752 + (913.178230403046 -
										 316.49805267936244*z)*z +
									 y*(241.04130980405 + y*(-44.5794634280918 +
											 49.023632509086724*z) +
										 z*(-331.6338314040912 + 77.78288016926652*z))))));
	double g08 = x2*(769.588305957198 + z*(-579.1945920644748 +
														(190.08980732018878 - 52.4960313181984*z)*z) +
						  x*(-104.588181856267 + x*(-8.16387957824522 -
								  181.05306718269662*z) +
							  (408.2669656358754 - 40.95023189295384*z)*z +
							  y*(166.3847855603638 - 176.898386096574*z +
								  y*(-108.3834525034224 + 153.8390924339484*z))) +
						  y*(-687.913805923122 + z*(748.126026697488 +
								  z*(-379.883572632876 + 140.9317606630898*z)) +
							  y*(674.819060538734 + z*(-1069.887337245828 +
									  (530.4484299696 - 158.40030944233638*z)*z) +
								  y*(-409.779283929806 + y*(149.452282277512 -
										  218.92375140095282*z) +
									  (681.370187043564 - 133.5392811916956*z)*z))));
	rval = (g03 + g08)*1e-16;
	return rval;
}

/*************************************************************************************
=======
method: gsw_gibbs_g08_a(sa, z, y, x)
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g08_a(double sa, double z, double y, double x)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double rval;
	double g08 = 8645.36753595126 +
					 z*(-6620.98308089678 +
						 z*(769.588305957198 +
							 z*(-193.0648640214916 +
								 (31.6816345533648 - 5.24960313181984*z)*z))) +
					 x*(-7296.43987145382 +
						 x*(8103.20462414788 +
							 y*(2175.341332000392 +
								 y*(-274.2290036817964 +
									 y*(197.4670779425016 +
										 y*(-68.5590309679152 +
											 9.98788038278032*y))) - 90.6734234051316*z) +
							 x*(-5458.34205214835 - 980.14153344888*y +
								 x*(2247.60742726704 - 340.1237483177863*x +
									 220.542973797483*y) + 180.142097805543*z) +
							 z*(-219.1676534131548 + (-16.32775915649044 -
									 120.7020447884644*z)*z)) +
						 z*(598.378809221703 +
							 z*(-156.8822727844005 +
								 (204.1334828179377 - 10.23755797323846*z)*z)) +
						 y*(-1480.222530425046 +
							 z*(-525.876123559641 +
								 (249.57717834054571 - 88.449193048287*z)*z) +
							 y*(-129.1994027934126 +
								 z*(1149.174198007428 +
									 z*(-162.5751787551336 + 76.9195462169742*z)) +
								 y*(-30.0682112585625 - 1380.9597954037708*z +
									 y*(2.626801985426835 + 703.695562834065*z))))) +
					 y*(1187.3715515697959 +
						 z*(1458.233059470092 +
							 z*(-687.913805923122 +
								 z*(249.375342232496 +
									 z*(-63.313928772146 + 14.09317606630898*z)))) +
						 y*(1760.062705994408 +
							 y*(-450.535298526802 +
								 y*(182.8520895502518 +
									 y*(-43.3206481750622 + 4.26033941694366*y) +
									 z*(-595.457483974374 + (149.452282277512 -
											 72.9745838003176*z)*z)) +
								 z*(1388.489628266536 +
									 z*(-409.779283929806 +
										 (227.123395681188 - 22.2565468652826*z)*z))) +
							 z*(-1721.528607567954 +
								 z*(674.819060538734 +
									 z*(-356.629112415276 + (88.4080716616 -
											 15.84003094423364*z)*z)))));
	if (sa > 0.0)
	{
		g08 += (11625.62913253464 + 1702.453469893412*y)*log(x);
	}
	else
		g08 = 0.0;
	rval = 0.5*gtc.gsw_sfac*g08;
	return rval;
}

/***************************************************************************************
=======
method: gsw_gibbs_g08_b(z,y,x);
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g08_b(double z,double y,double x)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double rval;
	double g08 = -6620.98308089678 + z*(1539.176611914396 +
													z*(-579.1945920644748 + (126.7265382134592 -
															26.2480156590992*z)*z)) +
					 x*(598.378809221703 + x*(-219.1676534131548 +
													  180.142097805543*x - 90.6734234051316*y +
													  (-32.65551831298088 - 362.10613436539325*z)*z) +
						 z*(-313.764545568801 + (612.4004484538132 -
														 40.95023189295384*z)*z) +
						 y*(-525.876123559641 + (499.15435668109143 -
														 265.347579144861*z)*z +
							 y*(1149.174198007428 + y*(-1380.9597954037708 +
									 703.695562834065*y) +
								 z*(-325.1503575102672 + 230.7586386509226*z)))) +
					 y*(1458.233059470092 + z*(-1375.827611846244 +
														z*(748.126026697488 + z*(-253.255715088584 +
																70.4658803315449*z))) +
						 y*(-1721.528607567954 + y*(1388.489628266536 +
								 y*(-595.457483974374 + (298.904564555024 -
										 218.92375140095282*z)*z) +
								 z*(-819.558567859612 + (681.370187043564 -
										 89.0261874611304*z)*z)) +
							 z*(1349.638121077468 + z*(-1069.887337245828 +
									 (353.6322866464 - 79.20015472116819*z)*z))));
	rval = g08*gtc.gsw_sfac*0.5e-8;
	return rval;
}

/***************************************************************************************
=======
method: gsw_gibbs_g08_c(sa,z,y,x);
=======
This method supplements method gsw_gibbs(ns,nt,np,sa,t,p) for clarity.
****************************************************************************************/

double TeosBase::gsw_gibbs_g08_c(double sa, double z, double y, double x)
{
	/**  for GSW_TEOS10_CONSTANTS use gtc */
	double rval;
	double g08 = 1187.3715515697959 + z*(1458.233059470092 +
													 z*(-687.913805923122 + z*(249.375342232496 +
															 z*(-63.313928772146 + 14.09317606630898*z)))) +
					 x*(-1480.222530425046 + x*(2175.341332000392 +
														 x*(-980.14153344888 + 220.542973797483*x) +
														 y*(-548.4580073635929 + y*(592.4012338275047 +
																 y*(-274.2361238716608 + 49.9394019139016*y))) -
														 90.6734234051316*z) +
						 z*(-525.876123559641 + (249.57717834054571 -
														 88.449193048287*z)*z) +
						 y*(-258.3988055868252 + z*(2298.348396014856 +
								 z*(-325.1503575102672 + 153.8390924339484*z)) +
							 y*(-90.2046337756875 - 4142.8793862113125*z +
								 y*(10.50720794170734 + 2814.78225133626*z)))) +
					 y*(3520.125411988816 + y*(-1351.605895580406 +
														y*(731.4083582010072 + y*(-216.60324087531103 +
																25.56203650166196*y) +
																z*(-2381.829935897496 + (597.809129110048 -
																		291.8983352012704*z)*z)) +
														z*(4165.4688847996085 + z*(-1229.337851789418 +
																(681.370187043564 - 66.7696405958478*z)*z))) +
						 z*(-3443.057215135908 + z*(1349.638121077468 +
								 z*(-713.258224830552 + (176.8161433232 -
										 31.68006188846728*z)*z))));
	if (sa > 0.0)
	{
		g08 += 1702.453469893412*log(x);
	}
	rval = 0.5*gtc.gsw_sfac*0.025*g08;
	return rval;
}

/** MOVED HERE FROM TeosSea ******************************************/
/**************************************************************************
==========================================================================
method: gsw_turner_rsubrho(*sa,*ct,*p,nz,*tu,*rsubrho,*p_mid)
==========================================================================
   Calculates the Turner angle and the Rsubrho as a method of pressure
   down a vertical water column.  These quantities express the relative
   contributions of the vertical gradients of Conservative Temperature
   and Absolute Salinity to the vertical stability (the square of the
   Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at
   the mid pressure between the individual data points in the vertical.
   Note that in the double-diffusive literature, papers concerned with
   the "diffusive" form of double-diffusive convection often define the
   stability ratio as the reciprocal of what is defined here as the
   stability ratio.
  sa      : Absolute Salinity         (a profile (length nz))     [g/kg]
  ct      : Conservative Temperature  (a profile (length nz))     [deg C]
  p       : sea pressure              (a profile (length nz))     [dbar]
  nz      : number of bottles
  tu      : Turner angle, on the same (nz-1) grid as p_mid.
            Turner angle has units of:           [ degrees of rotation ]
  rsubrho : Stability Ratio, on the same (nz-1) grid as p_mid.
            Rsubrho is dimensionless.                       [ unitless ]
  p_mid   : Mid pressure between p grid  (length nz-1)           [dbar]
*************************************************************************/

void TeosBase::gsw_turner_rsubrho(double *sa,double *ct,double *p,int nz,double *tu,double *rsubrho,double *p_mid)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	int k;
	double dsa, sa_mid, dct, ct_mid, alpha_mid, beta_mid;
	if (nz < 2)
		return;
	for (k = 0; k < nz-1; k++)
	{
		dsa         = (sa[k] - sa[k+1]);
		sa_mid      = 0.5e0*(sa[k] + sa[k+1]);
		dct         = (ct[k] - ct[k+1]);
		ct_mid      = 0.5e0*(ct[k] + ct[k+1]);
		p_mid[k]    = 0.5e0*(p[k] + p[k+1]);
		gsw_specvol_alpha_beta(sa_mid,ct_mid,p_mid[k],NULL,&alpha_mid,
									  &beta_mid);
		tu[k] = gtc.rad2deg*atan2((alpha_mid*dct + beta_mid*dsa),
										  (alpha_mid*dct - beta_mid*dsa));
		if (dsa == 0.0)
			rsubrho[k] = cppGSW_INVALID_VALUE;
		else
			rsubrho[k] = (alpha_mid*dct)/(beta_mid*dsa);
	}
}

/**************************************************************************
==========================================================================
method: gsw_util_interp1q_int (x, iy, x_i)
==========================================================================
  Returns the value of the 1-D method iy (integer) at the points of column
  vector x_i using linear interpolation. The vector x specifies the
  coordinates of the underlying interval.
  result(y_i)
==========================================================================
*************************************************************************/

double *TeosBase::gsw_util_interp1q_int(int nx,double *x,int *iy,int nxi,double *x_i,double *y_i)
{
	char *in_rng=NULL;
	int *j=NULL, *k=NULL, *r=NULL, *jrev=NULL, *ki=NULL, imax_x, imin_x, i, n, m, ii;
	double *xi=NULL, *xxi=NULL, u, max_x, min_x;
	if (nx <= 0 || nxi <= 0)
		return (NULL);
	min_x = max_x = x[0];
	imin_x = imax_x = 0;
	for (i=0; i<nx; i++)
	{
		if (x[i] < min_x)
		{
			min_x = x[i];
			imin_x = i;
		}
		else if (x[i] > max_x)
		{
			max_x = x[i];
			imax_x = i;
		}
	}
	in_rng = new char[nxi];
	memset(in_rng, 0, nxi*sizeof (char));
	for (i=n=0; i<nxi; i++)
	{
		if (x_i[i] <= min_x)
		{
			y_i[i] = iy[imin_x];
		}
		else if (x_i[i] >= max_x)
		{
			y_i[i] = iy[imax_x];
		}
		else
		{
			in_rng[i] = 1;
			n++;
		}
	}
	if (n==0)
	{
		if (in_rng) delete []in_rng;
		return (y_i);
	}
	xi = new double[n];
	k = new int[3*n];
	ki = k+n;
	r = ki+n;
	m  = nx + n;
	xxi = new double[m];
	j = new int[2*m];
	jrev = j+m;
	ii = 0;
	for (i = 0; i<nxi; i++)
	{
		if (in_rng[i])
		{
			xi[ii] = x_i[i];
			ki[ii] = i;
			ii++;
		}
	}
	if (in_rng) delete []in_rng;
	/**
	 Note that the following operations on the index
	 vectors jrev and r depend on the sort utility
	 gsw_util_sort_dbl() consistently ordering the
	 sorting indexes either in ascending or descending
	 sequence for replicate values in the real vector.
	*/
	gsw_util_sort_dbl(xi, n, k);
	for (i = 0; i<nx; i++)
		xxi[i] = x[i];
	for (i = 0; i<n; i++)
		xxi[nx+i] = xi[k[i]];
	gsw_util_sort_dbl(xxi, nx+n, j);
	for (i = 0; i<nx+n; i++)
		jrev[j[i]] = i;
	for (i = 0; i<n; i++)
		r[k[i]] = jrev[nx+i] - i-1;
	for (i = 0; i<n; i++)
	{
		u = (xi[i]-x[r[i]])/(x[r[i]+1]-x[r[i]]);
		y_i[ki[i]] = iy[r[i]] + (iy[r[i]+1]-iy[r[i]])*u;
	}
	if (j) delete []j;
	if (xxi) delete []xxi;
	if (k) delete []k;
	if (xi) delete []xi;
	return (y_i);
}

/**************************************************************************
==========================================================================
method: gsw_util_linear_interp (x, y, x_i)
==========================================================================
  Returns the values of the methods y{ny} at the points of column
  vector x_i using linear interpolation. The vector x specifies the
  coordinates of the underlying interval, and the matrix y specifies
  the method values at each x coordinate. Note that y has dimensions
  nx x ny and y_i has dimensions nxi x ny.
  This method was adapted from Matlabs interp1q.
  result(y_i)
==========================================================================
*************************************************************************/

double *TeosBase::gsw_util_linear_interp(int nx,double *x,int ny,double *y,int nxi,double *x_i,double *y_i)
{
	char *in_rng=NULL;
	int *j=NULL, *k=NULL, *r=NULL, *jrev=NULL, *ki=NULL,
		  imax_x, imin_x, i, n, m, ii, jy, jy0, jyi0, r0;
	double *xi=NULL, *xxi=NULL, u, max_x, min_x;
	if (nx <= 0 || nxi <= 0 || ny <= 0)
		return (NULL);
	min_x = max_x = x[0];
	imin_x = imax_x = 0;
	for (i=0; i<nx; i++)
	{
		if (x[i] < min_x)
		{
			min_x = x[i];
			imin_x = i;
		}
		else if (x[i] > max_x)
		{
			max_x = x[i];
			imax_x = i;
		}
	}
	in_rng = new char[nxi];
	memset(in_rng, 0, nxi*sizeof (char));
	for (i=n=0; i<nxi; i++)
	{
		if (x_i[i] <= min_x)
		{
			for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
				y_i[jyi0+i] = y[jy0+imin_x];
		}
		else if (x_i[i] >= max_x)
		{
			for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
				y_i[jyi0+i] = y[jy0+imax_x];
		}
		else
		{
			in_rng[i] = 1;
			n++;
		}
	}
	if (n==0)
	{
		if (in_rng) delete []in_rng;
		return (y_i);
	}
	xi = new double[n];
	k = new int[3*n];
	ki = k+n;
	r = ki+n;
	m  = nx + n;
	xxi = new double[m];
	j = new int[2*m];
	jrev = j+m;
	ii = 0;
	for (i = 0; i<nxi; i++)
	{
		if (in_rng[i])
		{
			xi[ii] = x_i[i];
			ki[ii] = i;
			ii++;
		}
	}
	if (in_rng) delete []in_rng;
	/**
	   This algorithm mimics the Matlab interp1q method.
	   An explaination of this algorithm:
	   We have points we are interpolating from (x) and
	   points that we are interpolating to (xi).  We
	   sort the interpolating from points, concatenate
	   them with the interpolating to points and sort the result.
	   We then construct index r, the interpolation index in x for
	   each point in xi.
	   Note that the following operations on the index
	   vectors jrev and r depend on the sort utility
	   gsw_util_sort_dbl() consistently ordering the
	   sorting indexes either in ascending or descending
	   sequence for replicate values in the real vector.
	*/
	gsw_util_sort_dbl(xi, n, k);
	memmove(xxi, x, nx*sizeof (double));
	memmove(xxi+nx, xi, n*sizeof (double));
	gsw_util_sort_dbl(xxi, m, j);
	for (i = 0; i<m; i++)
		jrev[j[i]] = i;
	for (i = 0; i<n; i++)
		r[k[i]] = jrev[nx+i] - i - 1;
	/** this is now the interpolation index in x for a point in xi */
	for (jy=jy0=jyi0=0; jy < ny; jy++, jy0+=nx, jyi0+=nxi)
	{
		for (i = 0; i<n; i++)
		{
			u = (xi[i]-x[r[i]])/(x[r[i]+1]-x[r[i]]);
			r0 = jy0+r[i];
			y_i[jyi0+ki[i]] = y[r0] + (y[r0+1]-y[r0])*u;
		}
	}
	if (j) delete []j;
	if (xxi) delete []xxi;
	if (k) delete []k;
	if (xi) delete []xi;
	return (y_i);
}


/*************************************************************************
=================================================
method: gsw_util_pchip_interp(*x,*y,n,*xi,*yi,ni)
=================================================
   Piecewise-Hermite algorithm from
   https://en.wikipedia.org/wiki/Cubic_Hermite_spline
   Extrapolation to points outside the range is done by setting those
   points to the corresponding end values.
   The input x must be monotonically increasing; the interpolation points,
   xi, may be in any order, but the algorithm will be faster if they are
   monotonic, increasing or decreasing.
   Returns 0 on success, 1 if it fails because there are fewer than 2 points,
   2 if it fails because x is not increasing.
   Consistent with other GSW-C code at present, the memory allocations
   are assumed to succeed.
*************************************************************************/

int TeosBase::gsw_util_pchip_interp(double *x,double *y,int n,double *xi,double *yi,int ni)
{
	double *d=NULL;
	double t, tt, ttt, xx, dx;
	int i, j0, j1, err;
	double h00, h10, h01, h11;
	if (n<2)
	{
		return 1;
	}
	d = new double[n];
	err = gsw_pchip_derivs(x, y, n, d);
	if (err)
	{
		if (d) delete []d;
		return 2;
	}
	j0 = 0;
	for (i=0; i<ni; i++)
	{
		xx = xi[i];
		/** Linear search is appropriate and probably optimal for the
		   expected primary use case of interpolation to a finer grid.
		   It is inefficient but still methodal in the worst case of
		   randomly distributed xi.
		*/
		while (xx < x[j0] && j0 > 0)
		{
			j0--;
		}
		while (xx > x[j0+1] && j0 < n - 2)
		{
			j0++;
		}
		j1 = j0 + 1;
		if (xx >= x[j0] && xx <= x[j1])
		{
			dx = x[j1] - x[j0];
			t = (xx - x[j0]) / dx;
			tt = t * t;
			ttt = tt * t;
			/** Using intermediate variables for readability. */
			h00 = (2*ttt - 3*tt + 1);
			h10 =  (ttt - 2*tt + t);
			h01 = (-2*ttt + 3*tt);
			h11 = (ttt - tt);
			yi[i] = y[j0] * h00 + d[j0] * dx * h10 +
					  y[j1] * h01 + d[j1] * dx * h11;
		}
		else
		{
			/** extrapolate with constant end values */
			yi[i] = (xx < x[0]) ? y[0] : y[n-1];
		}
	}
	if (d) delete []d;
	return 0;
}

/*************************************************************************
=======
method: gsw_linear_interp_sa_ct_for_dh(*sa,*ct,*p,nz,*p_i,n_i,*sa_i,*ct_i)
=======
    Linearly interpolate to the grid made by define_grid_for_dh.
    We take advantage of what we know about the grids: they match
    at the end points, and both are monotonic.
**************************************************************************/

int TeosBase::gsw_linear_interp_sa_ct_for_dh(double *sa,double *ct,double *p,int nz,double *p_i,int n_i,double *sa_i,double *ct_i)
{
	int i, ii;
	double pfac;
	sa_i[0] = sa[0];
	sa_i[n_i-1] = sa[nz-1];
	ct_i[0] = ct[0];
	ct_i[n_i-1] = ct[nz-1];
	i = 1;
	for (ii=1; ii<n_i-1; ii++)
	{
		/** Find the second point of the pair in the original grid that
		   bracket the target.
		*/
		while (p[i] < p_i[ii])
		{
			i++;
			if (i == nz)
			{
				return -1;  /** error  */
			}
		}
		pfac = (p_i[ii] - p[i-1]) / (p[i] - p[i-1]);
		sa_i[ii] = sa[i-1] + pfac * (sa[i] - sa[i-1]);
		ct_i[ii] = ct[i-1] + pfac * (ct[i] - ct[i-1]);
	}
	return 0;
}

/**************************************************************************
==========================================================================
method: gsw_rr68_interp_sa_ct (sa, ct, p, p_i, sa_i, ct_i)
==========================================================================
   Interpolate Absolute Salinity and Conservative Temperature values to
   arbitrary pressures using the Reiniger and Ross (1968) interpolation
   scheme.
   Note that this interpolation scheme requires at least four observed
   bottles on the cast.
   SA   =  Absolute Salinity                                  [ g/kg ]
   CT   =  Conservative Temperature (ITS-90)                 [ deg C ]
   p    =  sea pressure                                       [ dbar ]
            ( i.e. absolute pressure - 10.1325 dbar )
   p_i  =  pressures to interpolate to.
   SA_i = interpolated SA values at pressures p_i.
   CT_i = interpolated CT values at pressures p_i.
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_rr68_interp_sa_ct(double *sa,double *ct,double *p,int mp,double *p_i,int mp_i,double *sa_i,double *ct_i)
{
	int i, j, nshallow, ncentral, ndeep, *ip=NULL, *ip_i=NULL,
													  *ip_ishallow=NULL, *ip_icentral=NULL, *ip_ideep=NULL;
	char *shallow=NULL, *central=NULL, *deep=NULL;
	double *ip_shallow=NULL, *ip_central=NULL, *ip_deep=NULL, *dp=NULL, *p_ii=NULL;
	if (mp < 4)
	{
		/** need at least four bottles to perform this interpolation */
		ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
		return;
	}
	dp = new double[mp];
	for (i=1; i<mp; i++)
	{
		if ((dp[i-1] = (p[i] - p[i-1])) <= 0.0)
		{
			if (dp) delete []dp;
			ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
			return;
		}
	}
	shallow = new char[3*mp_i];
	central = shallow+mp_i;
	deep = central+mp_i;
	nshallow=ncentral=ndeep=0;
	memset(shallow, 0, 3*mp_i*sizeof (char));
	for (i=0; i<mp_i; i++)
	{
		if (p_i[i] >= p[0] && p_i[i] <= p[1])
		{
			nshallow++;
			shallow[i] = 1;
		}
		if (p_i[i] >= p[1] && p_i[i] <= p[mp-2])
		{
			ncentral++;
			central[i] = 1;
		}
		if (p_i[i] >= p[mp-2] && p_i[i] <= p[mp-1])
		{
			ndeep++;
			deep[i] = 1;
		}
	}
	if ((nshallow == 0) || (ncentral == 0) || (ndeep == 0))
	{
		if (shallow) delete []shallow;
		if (dp) delete []dp;
		ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
		return;
	}
	ip = new int[mp+mp_i];
	ip_i = ip+mp;
	for (i=0; i<mp; i++)
		ip[i] = i;
	for (i=0; i<mp_i; i++)
		ip_i[i] = i;
	ip_ishallow = new int[nshallow+ncentral+ndeep];
	ip_icentral = ip_ishallow+nshallow;
	ip_ideep = ip_icentral+ncentral;
	ip_shallow = new double[2*(nshallow+ncentral+ndeep)];
	ip_central = ip_shallow+nshallow;
	ip_deep = ip_central+ncentral;
	p_ii = ip_deep+ndeep;
	/**
	  Calculate the 2 outer extrapolated values and the inner
	  interpolated values
	*/
	for (i=j=0; i<mp_i; i++)
	{
		if (central[i])
		{
			ip_icentral[j] = ip_i[i];
			j++;
		}
	}
	for (i=0; i<ncentral; i++)
		p_ii[i] = p_i[ip_icentral[i]];
	gsw_util_interp1q_int(mp,p,ip,ncentral,p_ii,ip_central);
	gsw_rr68_interp_section(0,sa,ct,p,mp,ncentral,ip_central,ip_icentral,
									p_i,sa_i,ct_i);
	for (i=j=0; i<mp_i; i++)
	{
		if (shallow[i])
		{
			ip_ishallow[j] = ip_i[i];
			j++;
		}
	}
	for (i=0; i<nshallow; i++)
		p_ii[i] = p_i[ip_ishallow[i]];
	gsw_util_interp1q_int(mp,p,ip,nshallow,p_ii,ip_shallow);
	gsw_rr68_interp_section(-1,sa,ct,p,mp,nshallow,ip_shallow,ip_ishallow,
									p_i,sa_i,ct_i);
	for (i=j=0; i<mp_i; i++)
	{
		if (deep[i])
		{
			ip_ideep[j] = ip_i[i];
			j++;
		}
	}
	for (i=0; i<ndeep; i++)
		p_ii[i] = p_i[ip_ideep[i]];
	gsw_util_interp1q_int(mp,p,ip,ndeep,p_ii,ip_deep);
	gsw_rr68_interp_section(1,sa,ct,p,mp,ndeep,ip_deep,ip_ideep,p_i,sa_i,ct_i);
	/**
	  Insert any observed bottles that are at the required interpolated
	  pressures
	*/
	for (i=0; i<mp_i; i++)
	{
		for (j=0; j<mp; j++)
		{
			if (p_i[i] == p[j])
			{
				sa_i[i] = sa[j];
				ct_i[i] = ct[j];
			}
		}
	}
	if (ip_shallow) delete []ip_shallow;
	if (ip_ishallow) delete []ip_ishallow;
	if (ip) delete []ip;
	if (shallow) delete []shallow;
	if (dp) delete []dp;
}

/**************************************************************************
==========================================================================
method: gsw_specvol_alpha_beta(sa,ct,p,*specvol,*alpha,*beta)
==========================================================================
   Calculates specific volume, the appropiate thermal expansion coefficient
   and the appropriate saline contraction coefficient of seawater from
   Absolute Salinity and Conservative Temperature.  This method uses the
   computationally-efficient expression for specific volume in terms of
   SA, CT and p (Roquet et al., 2014).
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   p     :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   specvol =  specific volume                                      [ m/kg ]
   alpha   =  thermal expansion coefficient                         [ 1/K ]
              with respect to Conservative Temperature
   beta    =  saline (i.e. haline) contraction                     [ kg/g ]
              coefficient at constant Conservative Temperature
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_specvol_alpha_beta(double sa,double ct,double p,double *specvol,double *alpha,double *beta)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */
	double  v, v_ct, v_sa_part, xs, ys, z;
	xs = sqrt(gtc.gsw_sfac * sa + gtc.offset);
	ys = ct * 0.025;
	z = p * gtc.rec_db2pa;
	v = gsw_gsvco_v(xs, ys, z);

	if (specvol != NULL)
		*specvol = v;

	if (alpha != NULL)
	{
		v_ct = gsw_gsvco_a(xs, ys, z);
		*alpha = 0.025*v_ct/v;
	}

	if (beta != NULL)
	{
		v_sa_part = gsw_gsvco_b(xs, ys, z);
		*beta = -v_sa_part*0.5*gtc.gsw_sfac/(v*xs);
	}
}

/**************************************************************************
==========================================================================
method: gsw_sigma1(sa,ct)
==========================================================================
   Calculates potential density anomaly with reference pressure of 1000 dbar,
   this being this particular potential density minus 1000 kg/m^3.  This
   method has inputs of Absolute Salinity and Conservative Temperature.
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature                        [deg C]
  sigma1 : potential density anomaly with reference pressure of 1000
*************************************************************************/

double TeosBase::gsw_sigma1(double sa,double ct)
{
	return (gsw_rho(sa,ct,1000.0) - 1000.0);
}

/**************************************************************************
==========================================================================
method: gsw_sigma2(sa,ct)
==========================================================================
   Calculates potential density anomaly with reference pressure of 2000 dbar,
   this being this particular potential density minus 1000 kg/m^3.  This
   method has inputs of Absolute Salinity and Conservative Temperature.
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature                        [deg C]
  sigma2 : potential density anomaly with reference pressure of 2000
**************************************************************************/

double TeosBase::gsw_sigma2(double sa,double ct)
{
	return (gsw_rho(sa,ct,2000.0) - 1000.0);
}

/**************************************************************************
==========================================================================
method: gsw_sigma3(sa,ct)
==========================================================================
   Calculates potential density anomaly with reference pressure of 3000 dbar,
   this being this particular potential density minus 1000 kg/m^3.  This
   method has inputs of Absolute Salinity and Conservative Temperature.
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature                        [deg C]
  sigma3 : potential density anomaly with reference pressure of 3000
*************************************************************************/

double TeosBase::gsw_sigma3(double sa,double ct)
{
	return (gsw_rho(sa,ct,3000.0) - 1000.0);
}

/**************************************************************************
==========================================================================
method: gsw_sigma4(sa,ct)
==========================================================================
   Calculates potential density anomaly with reference pressure of 4000 dbar,
   this being this particular potential density minus 1000 kg/m^3.  This
   method has inputs of Absolute Salinity and Conservative Temperature.
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature                        [deg C]
  sigma4  : potential density anomaly with reference pressure of 4000
*************************************************************************/

double TeosBase::gsw_sigma4(double sa,double ct)
{
	return (gsw_rho(sa,ct,4000.0) - 1000.0);
}

/**************************************************************************
==========================================================================
method: gsw_sp_from_sstar(sstar,p,lon,lat)
==========================================================================
  Calculates Practical Salinity, SP, from Preformed Salinity, Sstar.
  sstar  : Preformed Salinity                              [g/kg]
  p      : sea pressure                                    [dbar]
  lon   : longitude                                       [deg E]
  lat    : latitude                                        [deg N]
  gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]
*************************************************************************/

double TeosBase::gsw_sp_from_sstar(double sstar,double p,double lon,double lat)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double saar, sp_baltic;
	/**
	**  In the Baltic Sea, SA = Sstar.
	*/
	sp_baltic = gsw_sp_from_sa_baltic(sstar,lon,lat);
	if (sp_baltic < cppGSW_ERROR_LIMIT)
		return (sp_baltic);
	saar = gsw_saar(p,lon,lat);
	if (saar == cppGSW_INVALID_VALUE)
		return (saar);
	return ((sstar/gtc.gsw_ups)/(1.0 - 0.35e0*saar));
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_melting_ice_sa_ct_ratio_poly(double sa,double ct,double p,double t_ih)
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

/**************************************************************************
==========================================================================
method: gsw_cabbeling(sa,ct,p)
==========================================================================
   Calculates the cabbeling coefficient of seawater with respect to
   Conservative Temperature.  This method uses the computationally-
   efficient expression for specific volume in terms of SA, CT and p
   (Roquet et al., 2014).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
  cabbeling  : cabbeling coefficient with respect to       [1/K^2]
               Conservative Temperature.
*************************************************************************/
double TeosBase::gsw_cabbeling(double sa,double ct,double p)
{
	double alpha_ct, alpha_on_beta, alpha_sa, beta_sa, rho,
			 v_sa, v_ct, v_sa_sa, v_sa_ct, v_ct_ct;
	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct, NULL);
	gsw_specvol_second_derivatives(sa,ct,p,&v_sa_sa,&v_sa_ct,&v_ct_ct,
											 NULL, NULL);
	rho             = gsw_rho(sa,ct,p);
	alpha_ct        = rho*(v_ct_ct - rho*v_ct*v_ct);
	alpha_sa        = rho*(v_sa_ct - rho*v_sa*v_ct);
	beta_sa         = -rho*(v_sa_sa - rho*v_sa*v_sa);
	alpha_on_beta   = gsw_alpha_on_beta(sa,ct,p);
	return (alpha_ct +
			  alpha_on_beta*(2.0*alpha_sa - alpha_on_beta*beta_sa));
}
/**************************************************************************
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
  ice Ih plus brine) and this programme returns NaNs if the input
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
*************************************************************************/

double TeosBase::gsw_melting_seaice_sa_ct_ratio(double sa,double ct,double p,double sa_seaice,double t_seaice)
{
	double  ctf, delsa, h, h_brine, h_ih, sa_brine,
			  tf_sa_seaice, h_hat_sa, h_hat_ct;
	double  saturation_fraction = 0.0;
	if (sa_seaice < 0.0 || sa_seaice > 15.0)
	{
		return (cppGSW_INVALID_VALUE);
	}
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf)      /*the seawater ct input is below the freezing temp*/
	{
		return (cppGSW_INVALID_VALUE);
	}
	tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction) - 1e-6;
	if (t_seaice > tf_sa_seaice)     /*t_seaice exceeds the freezing sa*/
	{
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
	if (sa_brine > cppGSW_ERROR_LIMIT)
	{
		return (cppGSW_INVALID_VALUE);
	}
	h_brine = gsw_enthalpy_t_exact(sa_brine,t_seaice,p);
	delsa = sa - sa_seaice;
	return (h_hat_ct*delsa /
			  (h - h_ih - delsa*h_hat_sa - sa_seaice*(h_brine - h_ih)/sa_brine));
}

/**************************************************************************
==========================================================================
method: gsw_thermobaric(sa,ct,p)
==========================================================================
   Calculates the thermobaric coefficient of seawater with respect to
   Conservative Temperature.  This routine is based on the
   computationally-efficient expression for specific volume in terms of
   SA, CT and p (Roquet et al., 2014).
  sa     : Absolute Salinity                               [g/kg]
  ct     : Conservative Temperature (ITS-90)               [deg C]
  p      : sea pressure                                    [dbar]
  thermobaric  : thermobaric coefficient with              [1/(K Pa)]
                     respect to Conservative Temperature (48 term equation)
*************************************************************************/

double TeosBase::gsw_thermobaric(double sa,double ct,double p)
{
	double v_ct, v_ct_p, v_sa, v_sa_p;
	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,NULL);
	gsw_specvol_second_derivatives(sa,ct,p,NULL,NULL,NULL,&v_sa_p,&v_ct_p);
	return (gsw_rho(sa,ct,p)*(v_ct_p - (v_ct/v_sa)*v_sa_p));
}

/**************************************************************************
==========================================================================
  method: gsw_pt_from_pot_enthalpy_ice (pot_enthalpy_ice)
==========================================================================
   Calculates the potential temperature of ice from the potential enthalpy
   of ice.  The reference sea pressure of both the potential temperature
   and the potential enthalpy is zero dbar.
   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]
   pt0_ice  :  potential temperature of ice (ITS-90)         [ deg C ]
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_pt_from_pot_enthalpy_ice(double pot_enthalpy_ice)
{
	int iteration;
	double df_dt, f, mod_pot_enthalpy_ice, pt0_cold_ice, recip_df_dt,
			 pt0_cold_ice_old, pt0_ice, pt0_ice_old, ptm_cold_ice, ptm_ice;
	double  h00 = -6.320202333358860e5, /*gsw_enthalpy_ice(-gtc.gsw_t0,0)*/
			  p0 = 0.0;
	mod_pot_enthalpy_ice = gsw_max_d(pot_enthalpy_ice,h00);
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

/**************************************************************************
==========================================================================
  method: gsw_pt0_cold_ice_poly (pot_enthalpy_ice)
==========================================================================
   Calculates an initial estimate of pt0_ice when it is less than about
   -100 deg C.
   pot_enthalpy_ice  :  potential enthalpy of ice        [ J/kg ]
   pt0_cold_ice_poly  :  initial estimate of potential temperatur
          of very cold ice in dgress C (not K)     [ deg C ]
--------------------------------------------------------------------------
*************************************************************************/

double TeosBase::gsw_pt0_cold_ice_poly(double pot_enthalpy_ice)
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

/**************************************************************************
==========================================================================
  method: gsw_deltasa_from_sp(sp,p,lon,lat)
==========================================================================
  Calculates Absolute Salinity Anomaly, deltaSA, from Practical Salinity, SP.
  sp     : Practical Salinity               [unitless]
  p      : sea pressure                [dbar]
  lon    : longitude                   [deg E]
  lat    : latitude               [deg N]
  gsw_deltasa_from_sp : Absolute Salinty Anomaly      [g/kg]
*************************************************************************/

double TeosBase::gsw_deltasa_from_sp(double sp,double p,double lon,double lat)
{
	double  res;
	res = gsw_sa_from_sp(sp,p,lon,lat) - gsw_sr_from_sp(sp);
	if (res > cppGSW_ERROR_LIMIT)
	{
		res = cppGSW_INVALID_VALUE;
	}
	return (res);
}

/**************************************************************************
================================================================================================
method: gsw_rho_second_derivatives(sa,ct,p,*rho_sa_sa,*rho_sa_ct,*rho_ct_ct,*rho_sa_p,*rho_ct_p)
================================================================================================
   Calculates five second-order derivatives of rho. Note that this method
   uses the computationally-efficient expression for specific
   volume (Roquet et al., 2014).
   SA  :  Absolute Salinity                                        [ g/kg ]
   CT  :  Conservative Temperature (ITS-90)                       [ deg C ]
   p   :  sea pressure                                             [ dbar ]
          ( i.e. absolute pressure - 10.1325 dbar )
   rho_SA_SA : The second-order derivative of rho with respect to
               Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
   rho_SA_CT : The second-order derivative of rho with respect to
               SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
   rho_CT_CT : The second-order derivative of rho with respect to CT at
               constant SA & p
   rho_SA_P  : The second-order derivative with respect to SA & P at
               constant CT.
   rho_CT_P  : The second-order derivative with respect to CT & P at
               constant SA.
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_rho_second_derivatives(double sa,double ct,double p,
		double &rho_sa_sa,double &rho_sa_ct,double &rho_ct_ct,double &rho_sa_p,
		double &rho_ct_p)
{
	double rec_v, rec_v2, rec_v3, v_ct, v_ct_ct, v_ct_p, v_p, v_sa,
			 v_sa_ct, v_sa_p, v_sa_sa;

	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,&v_p);
	gsw_specvol_second_derivatives(sa,ct,p,&v_sa_sa,&v_sa_ct,&v_ct_ct,
											 &v_sa_p,&v_ct_p);

	rec_v = 1.0 / gsw_specvol(sa,ct,p);
	rec_v2 = pow(rec_v, 2);
	rec_v3 = rec_v2*rec_v;

	//if (rho_sa_sa != NULL)
	rho_sa_sa = -v_sa_sa*rec_v2 + 2.0*v_sa*v_sa*rec_v3;
	//if (rho_sa_ct != NULL)
	rho_sa_ct = -v_sa_ct*rec_v2 + 2.0*v_sa*v_ct*rec_v3;
	//if (rho_ct_ct != NULL)
	rho_ct_ct = -v_ct_ct*rec_v2 + 2.0*v_ct*v_ct*rec_v3;
	//if (rho_sa_p != NULL)
	rho_sa_p = -v_sa_p*rec_v2 + 2.0*v_sa*v_p*rec_v3;
	//if (rho_ct_p != NULL)
	rho_ct_p = -v_ct_p*rec_v2 + 2.0*v_ct*v_p*rec_v3;
}

/**************************************************************************
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
*************************************************************************/

double TeosBase::gsw_sa_freezing_from_ct_poly(double ct,double p,double saturation_fraction)
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
	sa = gsw_max_d(sa,0.0);
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

/**************************************************************************
==========================================================================
method: gsw_pt_first_derivatives (sa, ct, pt_sa, pt_ct)
=========================================================================
   Calculates the following two partial derivatives of potential temperature
   (the regular potential temperature whose reference sea pressure is 0 dbar)
   (1) pt_SA, the derivative with respect to Absolute Salinity at
        constant Conservative Temperature, and
   (2) pt_CT, the derivative with respect to Conservative Temperature at
        constant Absolute Salinity.
   SA    :  Absolute Salinity                                        [ g/kg ]
   CT    :  Conservative Temperature (ITS-90)                       [ deg C ]
   pt_SA =  The derivative of potential temperature with respect to
            Absolute Salinity at constant Conservative Temperature.
                                                                [ K/(g/kg)]
   pt_CT =  The derivative of potential temperature with respect to
            Conservative Temperature at constant Absolute Salinity.
            pt_CT is dimensionless.                            [ unitless ]
---------------------------------------------------------------------------
*************************************************************************/

void TeosBase::gsw_pt_first_derivatives(double sa,double ct,double *pt_sa,double *pt_ct)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	double abs_pt, ct_pt, ct_sa, pt, pr0 = 0.0;
	int n0=0, n1=1, n2=2;
	pt = gsw_pt_from_ct(sa,ct);
	abs_pt = (gtc.gsw_t0 + pt);
	ct_pt = -(abs_pt*gsw_gibbs(n0,n2,n0,sa,pt,pr0))/gtc.gsw_cp0;
	if (pt_sa != NULL)
	{
		ct_sa = (gsw_gibbs(n1,n0,n0,sa,pt,pr0) -
					abs_pt*gsw_gibbs(n1,n1,n0,sa,pt,pr0))/gtc.gsw_cp0;
		*pt_sa = -ct_sa/ct_pt;
	}
	if (pt_ct != NULL)
		*pt_ct = 1.0/ct_pt;
}
/** moved here from conversion of .m files to C++ */
/**
% gsw_t_from_rho_ice                        temperature from density of ice
% =========================================================================
%
% USAGE: gsw_t_from_rho_ice(rho_ice,p)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of ice, for given values of its
%  density and sea pressure (in dbar).
%
% INPUT:
%  rho_ice =  density of ice                                     [ kg/m^3 ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho_ice is MxN.
%
% OUTPUT:
%  t_ice  =  in-situ temperature of ice                           [ deg C ]
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.5 of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_t_from_rho_ice(double rho_ice, double p)
{
   double v_t_ice, t_ice, v_freezing, t_freezing, v_173_15,v_93_15,v_63_15;
   double t_ice_mean, delta_v_ice, t_ice_old, specvol_ice = 1.0 / rho_ice;

   v_173_15 = gsw_specvol_ice(-100.0 * 1.0, p);
   t_freezing = gsw_t_freezing(0.0, p, 0.0);
   v_freezing = gsw_specvol_ice(t_freezing, p);

   if (specvol_ice < v_173_15)
   {
      v_93_15 = gsw_specvol_ice(-180.0 * 1.0,p);
      v_63_15 = gsw_specvol_ice(-210.0 * 1.0,p);

      t_ice = -180.0 + 30.0 * (specvol_ice - v_93_15) / (v_93_15 - v_63_15);  // initial estimate.
      v_t_ice = (v_63_15 - v_93_15) / (-30.0); //initial estimate of v_t_ice, the t derivative of v

      t_ice = -100.0 + 80.0 * (specvol_ice - v_173_15) / (v_173_15 - v_93_15);  // initial estimate.
      v_t_ice = (v_93_15 - v_173_15) / (-80.0); //initial estimate of v_t_ice, the t derivative of v

      t_ice =  (100.0 - t_freezing * (specvol_ice -
                  v_freezing) / (v_freezing - v_173_15));  // initial estimate.

      v_t_ice = (v_freezing - v_173_15) /
                  (-100.0 - t_freezing); //initial estimate of v_t_ice, the t derivative of v
   }
   else
   {
      t_ice = (100.0 - t_freezing) * (specvol_ice -
                  v_freezing) / (v_freezing - v_173_15);  //% the initial estimate of t_ice

      v_t_ice = (v_freezing - v_173_15) /
                  (-100.0 - t_freezing); //%initial estimate of v_t_ice, the t derivative of v
   }

   /**
      --------------------------------------------------------------------------
            Begin the modified Newton-Raphson iterative procedure
      --------------------------------------------------------------------------
   */

   for (unsigned iter = 0; iter < 3; iter++) //Number_of_iterations = 1:3
   {
      t_ice_old = t_ice;
      delta_v_ice = gsw_specvol_ice(t_ice_old,p) - specvol_ice;
      t_ice = t_ice_old - delta_v_ice / v_t_ice ;

      // this is half way through the modified
      // N-R method (McDougall and Wotherspoon, 2012)

      t_ice_mean = 0.5 * (t_ice + t_ice_old);
      v_t_ice = gsw_gibbs_ice(1,1,t_ice_mean,p);
      t_ice = t_ice_old - delta_v_ice / v_t_ice;
   }

   /**
      After 3 iterations of this modified Newton-Raphson iteration,
      the error in t_ice is no larger than 2.3x10^-12 deg C, which
      is machine precision for this calculation.
   */
   return t_ice;
}
/**
% gsw_specvol_anom               specific volume anomaly (75-term equation)
%==========================================================================
%
% USAGE: gsw_specvol_anom(SA,CT,p,SA_ref,CT_ref)
%
% DESCRIPTION:
%  Calculates specific volume anomaly from Absolute Salinity, Conservative
%  Temperature and pressure. It uses the computationally-efficient
%  expression for specific volume as a function of SA, CT and p (Roquet
%  et al., 2015).  If the Absolute Salinity and Conservative Temperature
%  reference value (SA_ref and CT_ref) are not defined then a default
%  reference value of Absolute Salinity is SSO and the reference value of
%  Conservative Temperature is equal to 0 degress C.
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA     =  Absolute Salinity                                     [ g/kg ]
%  CT     =  Conservative Temperature (ITS-90)                    [ deg C ]
%  p      =  sea pressure                                          [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
% Optional:
%  SA_ref =  reference Absolute Salinity                           [ g/kg ]
%  CT_ref =  reference Conservative Temperature (ITS-90)          [ deg C ]
%
%  SA & CT  need to have the same dimensions.
%  SA_ref & CT_ref must be scalars and have dimensions 1x1.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  specvol_anom  =  specific volume anomaly                      [ m^3/kg ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_anom(double SA, double CT, double p, double SA_ref,
												 double CT_ref)
{
	double v010 = -1.5649734675e-5;
	double v011 =  1.8505765429e-5;
	double v012 = -1.1736386731e-6;
	double v013 = -3.6527006553e-7;
	double v014 =  3.1454099902e-7;
	double v020 =  2.7762106484e-5;
	double v021 = -1.1716606853e-5;
	double v022 =  2.1305028740e-6;
	double v023 =  2.8695905159e-7;
	double v030 = -1.6521159259e-5;
	double v031 =  7.9279656173e-6;
	double v032 = -4.6132540037e-7;
	double v040 =  6.9111322702e-6;
	double v041 = -3.4102187482e-6;
	double v042 = -6.3352916514e-8;
	double v050 = -8.0539615540e-7;
	double v051 =  5.0736766814e-7;
	double v060 =  2.0543094268e-7;
	double v100 = -3.1038981976e-4;
	double v101 =  2.4262468747e-5;
	double v102 = -5.8484432984e-7;
	double v103 =  3.6310188515e-7;
	double v104 = -1.1147125423e-7;
	double v110 =  3.5009599764e-5;
	double v111 = -9.5677088156e-6;
	double v112 = -5.5699154557e-6;
	double v113 = -2.7295696237e-7;
	double v120 = -3.7435842344e-5;
	double v121 = -2.3678308361e-7;
	double v122 =  3.9137387080e-7;
	double v130 =  2.4141479483e-5;
	double v131 = -3.4558773655e-6;
	double v132 =  7.7618888092e-9;
	double v140 = -8.7595873154e-6;
	double v141 =  1.2956717783e-6;
	double v150 = -3.3052758900e-7;
	double v200 =  6.6928067038e-4;
	double v201 = -3.4792460974e-5;
	double v202 = -4.8122251597e-6;
	double v203 =  1.6746303780e-8;
	double v210 = -4.3592678561e-5;
	double v211 =  1.1100834765e-5;
	double v212 =  5.4620748834e-6;
	double v220 =  3.5907822760e-5;
	double v221 =  2.9283346295e-6;
	double v222 = -6.5731104067e-7;
	double v230 = -1.4353633048e-5;
	double v231 =  3.1655306078e-7;
	double v240 =  4.3703680598e-6;
	double v300 = -8.5047933937e-4;
	double v301 =  3.7470777305e-5;
	double v302 =  4.9263106998e-6;
	double v310 =  3.4532461828e-5;
	double v311 = -9.8447117844e-6;
	double v312 = -1.3544185627e-6;
	double v320 = -1.8698584187e-5;
	double v321 = -4.8826139200e-7;
	double v330 =  2.2863324556e-6;
	double v400 =  5.8086069943e-4;
	double v401 = -1.7322218612e-5;
	double v402 = -1.7811974727e-6;
	double v410 = -1.1959409788e-5;
	double v411 =  2.5909225260e-6;
	double v420 =  3.8595339244e-6;
	double v500 = -2.1092370507e-4;
	double v501 =  3.0927427253e-6;
	double v510 =  1.3864594581e-6;
	double v600 =  3.1932457305e-5;
	double x2_ref,x2,xy_part_4_diff, xs,ys,xy_part_0, xy_part_1, xy_part_2, xy_part_3;
	double z, xs_ref, ys_ref, xy_part_0_ref, xy_part_1_ref, xy_part_2_ref, xy_part_3_ref;
	double sa = 0.0;

	if (SA > 0.0) sa = SA;

	x2 = gtc.gsw_sfac * sa;
	xs = sqrt(x2 + gtc.offset);
	ys = CT * 0.025;
	z = p * gtc.rec_db2pa;

	x2_ref = gtc.gsw_sfac * SA_ref;
	xs_ref = sqrt(x2_ref + gtc.offset);
	ys_ref = CT_ref * 0.025;

	xy_part_0 = xs *(v100
	               + xs *(v200 + xs *(v300 + xs *(v400 + xs *(v500
                  + v600 *xs))))) + ys *(v010 + xs *(v110 + xs *(v210
                  + xs *(v310 + xs *(v410 + v510 *xs)))) + ys *(v020
                  + xs *(v120 + xs *(v220 + xs *(v320 + v420 *xs)))
						+ ys *(v030 + xs *(v130 + xs *(v230 + v330 *xs))
						+ ys *(v040 + xs *(v140 + v240*xs) + ys *(v050
						+ v150 *xs + v060 *ys)))));

	xy_part_1 = xs *(v101 + xs *(v201 + xs *(v301 + xs *(v401 + v501 *xs))))
					   + ys *(v011 + xs *(v111 + xs *(v211 + xs *(v311
					   + v411 *xs))) + ys *(v021 + xs *(v121 + xs *(v221
					   + v321 *xs)) + ys *(v031 + xs *(v131 + v231 *xs)
						+ ys *(v041 + v141 *xs + v051 *ys))));

	xy_part_2 = xs *(v102 + xs *(v202 + xs *(v302 + v402 *xs)))
					+ ys *(v012 + xs *(v112 + xs *(v212 + v312 *xs))
					+ ys *(v022 + xs *(v122 + v222 *xs) + ys *(v032
					+ v132 *xs + v042 *ys)));

	xy_part_3 = xs * (v103 + v203 * xs) + ys * (v013 + v113 * xs + v023 * ys);

	xy_part_0_ref = xs_ref *(v100 + xs_ref *(v200 + xs_ref *(v300
	                  + xs_ref *(v400 + xs_ref *(v500 + v600 *xs_ref)))))
	                  + ys_ref * (v010 + xs_ref * (v110 + xs_ref *(v210
	                  + xs_ref *(v310 + xs_ref *(v410 + v510 *xs_ref))))
	                  + ys_ref *(v020 + xs_ref *(v120 + xs_ref *(v220
	                  + xs_ref *(v320 + v420 *xs_ref))) + ys_ref *(v030
	                  + xs_ref *(v130 + xs_ref *(v230 + v330 *xs_ref))
	                  + ys_ref *(v040 + xs_ref *(v140 + v240*xs_ref)
	                  + ys_ref *(v050 + v150 *xs_ref + v060 *ys_ref)))));

	xy_part_1_ref = xs_ref *(v101 + xs_ref *(v201 + xs_ref *(v301
	                  + xs_ref *(v401 + v501 *xs_ref)))) + ys_ref *(v011
	                  + xs_ref *(v111 + xs_ref *(v211
	                  + xs_ref *(v311 + v411 *xs_ref)))
	                  + ys_ref *(v021 + xs_ref *(v121
	                  + xs_ref *(v221 + v321 *xs_ref))
	                  + ys_ref *(v031 + xs_ref *(v131 + v231 *xs_ref)
                     + ys_ref *(v041 + v141 *xs_ref + v051 *ys_ref))));

	xy_part_2_ref = xs_ref *(v102 + xs_ref *(v202 + xs_ref *(v302
	                  + v402 *xs_ref))) + ys_ref * (v012 + xs_ref * (v112
	                  + xs_ref *(v212 + v312 *xs_ref)) + ys_ref *(v022
	                  + xs_ref *(v122 + v222 *xs_ref) + ys_ref *(v032
	                  + v132 *xs_ref + v042 *ys_ref)));

	xy_part_3_ref = xs_ref * (v103 + v203 * xs_ref) + ys_ref *
	                  (v013 + v113 *xs_ref + v023 *ys_ref);

	xy_part_4_diff = v104 * (xs - xs_ref) + v014 * (ys - ys_ref);

	return xy_part_0 - xy_part_0_ref + z *(xy_part_1 - xy_part_1_ref
               + z *(xy_part_2 - xy_part_2_ref + z *(xy_part_3 - xy_part_3_ref
                     + z *(xy_part_4_diff))));
}
/**
% gsw_internal_energy_first_derivatives       first derivatives of specific
%                             interal energy of seawater (75-term equation)
%==========================================================================
%
% USAGE: gsw_internal_energy_first_derivatives(SA,CT,p,&u_SA,
                                                &u_CT, &u_P)
%
% DESCRIPTION:
%  Calculates the first order derivates of specific internal energy of
%  seawater using the computationally-efficient expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  u_SA = The first derivative of internal energy with respect to
%           Absolute Salinity at constant CT & p.
%                                          [ (J/kg)(g/kg)^-1 ] i.e. [ J/g ]
%  u_CT = The first derivative of internal energy with respect to
%           Conservative Temperature at constant SA & p.    [ (J/kg) K^-1 ]
%  u_P = The first derivative of internal energy with respect to
%           pressure at constant SA & CT.                  [ (J/kg) Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_internal_energy_first_derivatives(double SA, double CT, double p,
		double &u_SA, double &u_CT,
		double &u_P)
{
	double v, P, h_SA, h_CT, v_SA, v_CT, v_P, sa = 0.0;

	if (SA > 0.0) sa = SA;

	P = (gtc.db2pa * p + gtc.gsw_p0);

	gsw_enthalpy_first_derivatives(sa, CT, p, &h_SA, &h_CT);
	v = gsw_specvol(sa, CT, p);
	gsw_specvol_first_derivatives(sa, CT, p, &v_SA, &v_CT, &v_P);

	u_SA = h_SA - (P * v_SA);

	u_CT = h_CT - (P * v_CT);

	u_P = v - (P * v_P);
}
/**
% gsw_z_from_depth                                    height, z, from depth
%==========================================================================
%
% USAGE: gsw_z_from_depth(depth)
%
% DESCRIPTION:
%  Calculates height, z, from depth.  Note that in general height is
%  negative in the ocean.
%
% INPUT:
%  depth  =  depth                                                    [ m ]
%
% OUTPUT:
%  z  =  height                                                       [ m ]
%
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_z_from_depth(double depth)
{
	return -depth;
}
/**
% gsw_beta_const_CT_t_exact                  saline contraction coefficient
%                                      at constant Conservative temperature
%==========================================================================
%
% USAGE: gsw_beta_const_CT_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the saline (i.e. haline) contraction coefficient of seawater
%  at constant Conservative Temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar)
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  beta_const_CT_t_exact  =  saline contraction coefficient        [ kg/g ]
%                            at constant Conservative Temperature
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.19.3) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_beta_const_CT_t_exact(double SA, double t, double p)
{
	double g_SA_mod, g_SA_T_mod, z, x2, x, y, y_pt, pt0, gp, factor, factora;

	pt0 = gsw_pt0_from_t(SA, t, p);
	x2 = gtc.gsw_sfac * SA;
	x = sqrt(x2);
	y = 0.025 * t;
	y_pt = 0.025 * pt0;
	z = gtc.rec_db2pa * p; //Note.The input pressure (p) is sea pressure in units of dbar.

	g_SA_T_mod = 1187.3715515697959 + z *(1458.233059470092 +
      z *(-687.913805923122 + z *(249.375342232496 + z *(-63.313928772146 + 14.09317606630898 *z)))) +
      x *(-1480.222530425046 + x *(2175.341332000392 + x *(-980.14153344888 + 220.542973797483 *x) +
      y *(-548.4580073635929 + y *(592.4012338275047 + y *(-274.2361238716608 + 49.9394019139016 *y))) -
      90.6734234051316 *z) + z *(-525.876123559641 + (249.57717834054571 - 88.449193048287 *z) *z) +
      y *(-258.3988055868252 + z *(2298.348396014856 + z *(-325.1503575102672 + 153.8390924339484 *z)) +
      y *(-90.2046337756875 - 4142.8793862113125 *z + y *(10.50720794170734 + 2814.78225133626 *z)))) +
      y *(3520.125411988816 + y *(-1351.605895580406 +
      y *(731.4083582010072 + y *(-216.60324087531103 + 25.56203650166196 *y) +
      z *(-2381.829935897496 + (597.809129110048 - 291.8983352012704 *z) *z)) +
      z *(4165.4688847996085 + z *(-1229.337851789418 + (681.370187043564 - 66.7696405958478 *z) *z))) +
      z *(-3443.057215135908 + z *(1349.638121077468 +
      z *(-713.258224830552 + (176.8161433232 - 31.68006188846728 *z) *z))));

	g_SA_T_mod = 0.5 * gtc.gsw_sfac * 0.025 * g_SA_T_mod;

	g_SA_mod = 8645.36753595126 +
            x *(-7296.43987145382 + x *(8103.20462414788 +
            y_pt *(2175.341332000392 + y_pt *(-274.2290036817964 +
				y_pt *(197.4670779425016 + y_pt *(-68.5590309679152 + 9.98788038278032 *y_pt)))) +
				x *(-5458.34205214835 - 980.14153344888 *y_pt +
				x *(2247.60742726704 - 340.1237483177863 *x + 220.542973797483 *y_pt))) +
				y_pt *(-1480.222530425046 +
				y_pt *(-129.1994027934126 +
				y_pt *(-30.0682112585625 + y_pt *(2.626801985426835 ))))) +
            y_pt *(1187.3715515697959 +
				y_pt *(1760.062705994408 + y_pt *(-450.535298526802 +
				y_pt *(182.8520895502518 + y_pt *(-43.3206481750622 + 4.26033941694366 *y_pt)))));

	g_SA_mod = 0.5 * gtc.gsw_sfac * g_SA_mod;

	gp = gsw_gibbs(0,0,1, SA, t, p);

	factora = g_SA_T_mod - g_SA_mod / (gtc.gsw_t0 + pt0);

	factor = factora / (gp * gsw_gibbs(0,2,0, SA, t, p));

	return gsw_gibbs(0,1,1, SA, t, p) * factor - gsw_gibbs(1,0,1, SA, t, p) / gp;
}

/**
% gsw_internal_energy_CT_exact         specific internal energy of seawater
%==========================================================================
%
% USAGE: gsw_internal_energy_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the specific internal energy of seawater from Absolute
%  Salinity, Conservative Temperature and pressure.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_internal_energy(SA,CT,p) which uses the computationally
%  efficient 75 term expression for specific volume in terms of SA, CT
%  and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  internal_energy_CT_exact  =  specific internal energy (u)       [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall                                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.11.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_internal_energy_CT_exact(double SA, double CT, double p)
{
   double t = gsw_t_from_ct(SA, CT, p);

   return gsw_internal_energy_t_exact(SA, t, p);
}

/**
% gsw_kappa_CT_exact                             isentropic compressibility
%==========================================================================
%
% USAGE: gsw_kappa_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the isentropic compressibility of seawater.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_kappa(SA,CT,p),
%  which uses the computationally efficient 75-term expression for density
%  in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  kappa_CT_exact  =  isentropic compressibility                   [ 1/Pa ]
%   Note. The output units are 1/Pa not 1/dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (2.16.1) and the row for kappa in Table P.1 of appendix P
%    of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_kappa_CT_exact(double SA, double CT, double p)
{
	double t = gsw_t_from_ct(SA, CT, p);

	double g_tt = gsw_gibbs(0,2,0, SA, t, p);
	double g_tp = gsw_gibbs(0,1,1, SA, t, p);

	/**
	   these two (g001,g002) are named according to
	   the 1st 3 params in the gsw_gibbs call
	*/
	double g001 = gsw_gibbs(0,0,1, SA, t, p);
	double g002 = gsw_gibbs(0,0,2, SA, t, p);

	return (g_tp * g_tp - g_tt * g002) / (g001 * g_tt);
}

/**
% gsw_enthalpy_first_derivatives_wrt_t_exact           first derivatives of
%                                                                  enthalpy
%==========================================================================
%
% USAGE: gsw_enthalpy_first_derivatives_wrt_t_exact(SA,t,p,
                                                      &h_SA_wrt_t,
                                                      &h_T_wrt_t,
                                                      &h_P_wrt_t)
%
% DESCRIPTION:
%  Calculates the following three derivatives of specific enthalpy, h.
%  These derivatives are done with respect to in-situ temperature t (in the
%  case of h_T_wrt_t) or at constant in-situ tempertature (in the cases of
%  h_SA_wrt_t and h_P_wrt_t).
%   (1) h_SA_wrt_t, the derivative with respect to Absolute Salinity at
%       constant t and p.
%   (2) h_T_wrt_t, derivative with respect to in-situ temperature t at
%       constant SA and p.
%   (3) h_P_wrt_t, derivative with respect to pressure P (in Pa) at constant
%       SA and t.  This output has the same dimensions as specific volume,
%       but is not equal to specific volume.
%
%  Note that this function uses the full Gibbs function.  This function
%  avoids the Nan that would exist in h_sub_SA at SA=0 if it were
%  evaluated in the straightforward way from the gibbs function (as in the
%  commented line 111 of the code below).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  h_SA_wrt_t =  The first derivative of specific enthalpy with respect to
%                Absolute Salinity at constant t and p.
%                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
%  h_T_wrt_t  =  The first derivative of specific enthalpy with respect to
%                in-situ temperature, t, at constant SA and p. [ J/(kg K) ]
%
%  h_P_wrt_t  =  The first derivative of specific enthalpy with respect to
%                pressure P (in Pa) at constant SA and t.        [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_enthalpy_first_derivatives_wrt_t_exact(double SA, double t, double p,
		double &h_SA_wrt_t,
		double &h_T_wrt_t,
		double &h_P_wrt_t)
{
   double gibbs_SA_T, g08_SA_T, gibbs_SA, g08_SA, x2, x, y, z, sa = 0.0;

   if (SA > 0.0) sa = SA;

   h_T_wrt_t = gsw_cp_t_exact(sa, t, p);
   h_P_wrt_t = gsw_specvol_t_exact(sa, t, p) - (gtc.gsw_t0 + t) * gsw_gibbs(0,1,1, sa, t, p);

   x2 = gtc.gsw_sfac * sa;
   x = sqrt(x2);
   y = t * 0.025;
   z = p * gtc.rec_db2pa; //1e-4; Note.The input pressure (p) is sea pressure in units of dbar.

   g08_SA = 8645.36753595126 + z *(-6620.98308089678 +
        z *(769.588305957198 + z *(-193.0648640214916 + (31.6816345533648 - 5.24960313181984 *z) *z))) +
        x *(-7296.43987145382 + x *(8103.20462414788 +
        y *(2175.341332000392 + y *(-274.2290036817964 +
        y *(197.4670779425016 + y *(-68.5590309679152 + 9.98788038278032 *y))) - 90.6734234051316 *z) +
        x *(-5458.34205214835 - 980.14153344888 *y +
        x *(2247.60742726704 - 340.1237483177863 *x + 220.542973797483 *y) + 180.142097805543 *z) +
        z *(-219.1676534131548 + (-16.32775915649044 - 120.7020447884644 *z) *z)) +
        z *(598.378809221703 + z *(-156.8822727844005 + (204.1334828179377 - 10.23755797323846 *z) *z)) +
        y *(-1480.222530425046 + z *(-525.876123559641 + (249.57717834054571 - 88.449193048287 *z) *z) +
        y *(-129.1994027934126 + z *(1149.174198007428 + z *(-162.5751787551336 + 76.9195462169742 *z)) +
        y *(-30.0682112585625 - 1380.9597954037708 *z + y *(2.626801985426835 + 703.695562834065 *z))))) +
        y *(1187.3715515697959 + z *(1458.233059470092 +
        z *(-687.913805923122 + z *(249.375342232496 + z *(-63.313928772146 + 14.09317606630898 *z)))) +
        y *(1760.062705994408 + y *(-450.535298526802 +
        y *(182.8520895502518 + y *(-43.3206481750622 + 4.26033941694366 *y) +
        z *(-595.457483974374 + (149.452282277512 - 72.9745838003176 *z) *z)) +
        z *(1388.489628266536 + z *(-409.779283929806 + (227.123395681188 - 22.2565468652826 *z) *z))) +
        z *(-1721.528607567954 + z *(674.819060538734 +
        z *(-356.629112415276 + (88.4080716616 - 15.84003094423364 *z) *z)))));

   gibbs_SA = 0.5 * gtc.gsw_sfac * g08_SA;

   g08_SA_T = 1187.3715515697959 + z *(1458.233059470092 +
        z *(-687.913805923122 + z *(249.375342232496 + z *(-63.313928772146 + 14.09317606630898 *z)))) +
        x *(-1480.222530425046 + x *(2175.341332000392 + x *(-980.14153344888 + 220.542973797483 *x) +
        y *(-548.4580073635929 + y *(592.4012338275047 + y *(-274.2361238716608 + 49.9394019139016 *y))) -
        90.6734234051316 *z) + z *(-525.876123559641 + (249.57717834054571 - 88.449193048287 *z) *z) +
        y *(-258.3988055868252 + z *(2298.348396014856 + z *(-325.1503575102672 + 153.8390924339484 *z)) +
        y *(-90.2046337756875 - 4142.8793862113125 *z + y *(10.50720794170734 + 2814.78225133626 *z)))) +
        y *(3520.125411988816 + y *(-1351.605895580406 +
        y *(731.4083582010072 + y *(-216.60324087531103 + 25.56203650166196 *y) +
        z *(-2381.829935897496 + (597.809129110048 - 291.8983352012704 *z) *z)) +
        z *(4165.4688847996085 + z *(-1229.337851789418 + (681.370187043564 - 66.7696405958478 *z) *z))) +
        z *(-3443.057215135908 + z *(1349.638121077468 +
        z *(-713.258224830552 + (176.8161433232 - 31.68006188846728 *z) *z))));

   gibbs_SA_T = 0.5 * gtc.gsw_sfac *0.025 * g08_SA_T;

   h_SA_wrt_t = gibbs_SA - (gtc.gsw_t0 + t) * gibbs_SA_T;
}
/**
% gsw_specvol_CT_exact                                      specific volume
%==========================================================================
%
% USAGE: gsw_specvol_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific volume from Absolute Salinity, Conservative
%  Temperature and pressure.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_specvol(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA and CT are MxN.
%
% OUTPUT:
%  specvol_CT_exact  =  specific volume                          [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.7.2) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_CT_exact(double SA, double CT, double p)
{
	double t = gsw_t_from_ct(SA, CT, p);

	return gsw_specvol_t_exact(SA, t, p);
}
/**
% gsw_t90_from_t48              ITS-90 temperature from IPTS-48 temperature
%==========================================================================
%
% USAGE: gsw_t90_from_t48(t48)
%
% DESCRIPTION:
%  Calculates IPTS-48 temperature to International Temperature Scale 1990
%  (ITS-90) temperature.  This conversion should be applied to all in-situ
%  data collected prior to 31/12/1967.
%
% INPUT:
%  t48  =  in-situ temperature  (IPTS-48)                         [ deg C ]
%
% OUTPUT:
%  t90  =  in-situ temperature  (ITS-90)                          [ deg C ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  International Temperature Scales of 1948, 1968 and 1990, an ICES
%   note, available from http://www.ices.dk/ocean/procedures/its.htm
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_t90_from_t48(double t48)
{
	return (t48 - 4.4E-6 * t48 * (100.0 - t48)) / 1.00024;
}
/**
% gsw_dynamic_enthalpy_t_exact                  dyamic enthalpy of seawater
%==========================================================================
%
% USAGE: gsw_dynamic_enthalpy_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the dynamic enthalpy of seawater from Absolute Salinity,
%  in situ temperature and pressure.  Dynamic enthalpy was defined by
%  Young (2010) as the difference between enthalpy and potential enthalpy.
%  Note that this function uses the full TEOS-10 Gibbs function (i.e. the
%  sum of the IAPWS-09 and IAPWS-08 Gibbs functions, see the TEOS-10
%  Manual, IOC et al. (2010)).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  dynamic_enthalpy_t_exact  =  dynamic enthalpy                   [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Young, W.R., 2010: Dynamic enthalpy, Conservative Temperature, and the
%   seawater Boussinesq approximation. Journal of Physical Oceanography,
%   40, 394-400.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_dynamic_enthalpy_t_exact(double SA, double t, double p)
{
   double CT = gsw_ct_from_t(SA, t, p);

   return gsw_enthalpy_t_exact(SA, t, p) - gtc.gsw_cp0 * CT;
}
/**
% gsw_pot_enthalpy_from_pt    potential enthalpy from potential temperature
%==========================================================================
%
% USAGE: gsw_pot_enthalpy_from_pt(SA,pt)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of seawater from potential
%  temperature (whose reference sea pressure is zero dbar).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]
%
%  SA & pt need to have the same dimensions.
%
% OUTPUT:
%  pot_enthalpy  =  potential enthalpy                             [ J/kg ]
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.2 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_pot_enthalpy_from_pt(double SA, double pt)
{
   double pot_enthalpy, x2, x, y, sa = 0.0;

   if (SA > 0.0) sa = SA;

   x2 = gtc.gsw_sfac * sa;
   x = sqrt(x2);
   y = pt * 0.025;                                // normalize for F03 and F08

   pot_enthalpy =  61.01362420681071 + y *(168776.46138048015 +
    y *(-2735.2785605119625 + y *(2574.2164453821433 +
    y *(-1536.6644434977543 + y *(545.7340497931629 +
    (-50.91091728474331 - 18.30489878927802 *y) *y))))) +
    x2 *(268.5520265845071 + y *(-12019.028203559312 +
    y *(3734.858026725145 + y *(-2046.7671145057618 +
    y *(465.28655623826234 + (-0.6370820302376359 -
    10.650848542359153 *y) *y)))) +
    x *(937.2099110620707 + y *(588.1802812170108 +
    y *(248.39476522971285 + (-3.871557904936333 -
    2.6268019854268356 *y) *y)) +
    x *(-1687.914374187449 + x *(246.9598888781377 +
    x *(123.59576582457964 - 48.5891069025409 *x)) +
    y *(936.3206544460336 +
    y *(-942.7827304544439 + y *(369.4389437509002 +
    (-33.83664947895248 - 9.987880382780322 *y) *y))))));

//--------------------------------------------------------------------------
// The above polynomial for pot_enthalpy is the full expression for
// potential entahlpy in terms of SA and pt, obtained from the Gibbs
// function as below.  The above polynomial has simply collected like powers
// of x and y so that it is computationally faster than calling the Gibbs
// function twice as is done in the commented code below.  When this code
// below is run, the results are identical to the present function
// gsw_pot_enthalpy_from_pt, to machine precision.
//
//  pr0 = zeros(size(SA));
//  pot_enthalpy = gsw_gibbs(0,0,0,SA,pt,pr0) -
//                       (gsw_T0 + pt) *gsw_gibbs(0,1,0,SA,pt,pr0);
//
//-----------------This is the end of the alternative code------------------

   return pot_enthalpy;
}

/**
% gsw_N2Osol_SP_pt                            solubility of N2O in seawater
%==========================================================================
%
% USAGE: gsw_N2Osol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the nitrous oxide, N2O, concentration expected at equilibrium
%  with air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
%  including saturated water vapor  This function uses the solubility
%  coefficients as listed in Weiss and Price (1980).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  pt  =  potential temperature (ITS-90) referenced               [ deg C ]
%         to one standard atmosphere (0 dbar).
%
%  SP & pt need to have the same dimensions.
%
% OUTPUT:
%  N2Osol = solubility of nitrous oxide                           [ mol/L ]
%
% AUTHOR:  Rich Pawlowicz, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Weiss, R.F., and B.A. Price, 1980: Nitrous oxide solubility in water and
%   seawater. Mar. Chem., 8, 347-359.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_N2Osol_SP_pt(double SP, double pt)
{
	double N2Osol_tmp, pt68, y_100;

	pt68 = pt * 1.00024;
	y_100 = (pt68 + gtc.gsw_t0) * 0.01;
	N2Osol_tmp = log(y_100);

	return exp((((22287.43 / (pt68 + gtc.gsw_t0) + -165.8806) +
					 92.0792 * N2Osol_tmp) + -1.48425 * (y_100 * y_100)) +
				      SP * (y_100 * (-0.0048472 * y_100 + 0.031619) + -0.056235)) /
                     (1.0 - exp(((24.4543 - 6745.09 / (pt68 + gtc.gsw_t0))
                     - 4.8489 * N2Osol_tmp) - 0.000544 * SP));
}
/**
% gsw_t_maxdensity_exact                     in-situ temperature of maximum
%                                                       density of seawater
% =========================================================================
%
% USAGE: gsw_t_maxdensity_exact(SA,p)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of maximum density of seawater.
%  This function returns the in-situ temperature at which the density
%  of seawater is a maximum, at given Absolute Salinity, SA, and sea
%  pressure, p (in dbar).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  t_maxdensity_exact  =  in-situ temperature at which            [ deg C ]
%                         the density of seawater is a maximum for
%                         given Absolute Salinity and pressure.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.42 of this TEOS-10 Manual.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2013: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_t_maxdensity_exact(double SA, double p)
{
	double t_old, gibbs_PT, t_mean;
	double dt = 0.001;
	double t = 3.978 - (0.22072 * SA);
	double gibbs_PTT = 1.1e-8;

	for (unsigned itnum = 0; itnum < 3; itnum++)
	{
		t_old = t;
		gibbs_PT = gsw_gibbs(0, 1, 1, SA, t_old, p);
		t = t_old - gibbs_PT / gibbs_PTT;
		t_mean = 0.5 * (t + t_old);

		gibbs_PTT = gsw_gibbs(0, 1, 1, SA, t_mean + dt, p);
		gibbs_PTT -= gsw_gibbs(0, 1, 1, SA, t_mean - dt, p);
		gibbs_PTT /= dt + dt;

		t = t_old - gibbs_PT / gibbs_PTT;
	}

	return t;
}

/**
% gsw_sigma2_CT_exact                        potential density anomaly with
%                                       reference sea pressure of 2000 dbar
%==========================================================================
%
% USAGE: gsw_sigma2_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 2000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma2(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma2_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 2000 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sigma2_CT_exact(double SA, double CT)
{
   double t = gsw_t_from_ct(SA, CT, 2000.0);

   return gsw_rho_t_exact(SA, t, 2000.0) - 1000.0;
}
/**
% gsw_Nsquared_lowerlimit             specified profile of minimum buoyancy
%                                                         frequency squared
%==========================================================================
%
% USAGE: gsw_Nsquared_lowerlimit(p,long,lat)
%
% DESCRIPTION:
%  Calculates the minimum Nsquared such that a cast is stable.
%
% INPUT:
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  long =  longitude in decimal degrees                      [ 0 +360 ]
%                                                     or  [ -180 +180 ]
%  lat  =  latitude in decimal degrees north                [ -90 +90 ]
%
%  lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN, where p is MxN.
%
% OUTPUT:
%  Nsquared_lowerlimit = Minimum Brunt-Vaisala Frequency squared  [ 1/s^2 ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Barker, P.M., and T.J. McDougall, 2017: Stabilizing hydrographic
%   profiles with minimal change to the water masses. J. Atmosph. Ocean.
%   Tech., pp. 1935 - 1945, http://dx.doi.org/10.1175/JTECH-D-16-0111.1
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.18.3) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Nsquared_lowerlimit(double p, double b_long, double lat)
{
	double unnamed_idx_0 = p;

	if (p < -1.5)
	{
		unnamed_idx_0 = rtNaN;
	}

	return (0.75 * exp(-unnamed_idx_0 / 1000.0) + 0.25) * 1.0E-7;
}
/**
% gsw_Krsol_SP_pt                              solubility of Kr in seawater
%==========================================================================
%
% USAGE: gsw_Krsol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the krypton, Kr, concentration expected at equilibrium with
%  air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
%  including saturated water vapor.  This function uses the solubility
%  coefficients derived from the data of Weiss (1971).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  pt  =  potential temperature (ITS-90) referenced               [ deg C ]
%         to one standard atmosphere (0 dbar).
%
%  SP & pt need to have the same dimensions.
%
% OUTPUT:
%  Krsol = solubility of krypton in micro-moles per kg          [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Weiss, R.F. and T.K. Kyser, 1978: Solubility of Krypton in Water and
%   Seawater. J. Chem. Thermodynamics, 23, 69-72.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Krsol_SP_pt(double SP, double pt)
{
	double pt68, y_100;

	pt68 = pt * 1.00024;
	y_100 = (pt68 + gtc.gsw_t0) * 0.01;

	return exp((((15358.170000000002 / (pt68 + gtc.gsw_t0) + -112.684) +
					 74.469 * log(y_100)) + -10.0189 * y_100) + SP *
					   (y_100 * (0.0011201 * y_100 + -0.001844) +
					      -0.011213)) * 44.7405273118549;
}
/**
% gsw_osmotic_pressure_t_exact                             osmotic pressure
%==========================================================================
%
% USAGE: gsw_osmotic_pressure_t_exact(SA,t,pw)
%
% DESCRIPTION:
%  Calculates the osmotic pressure of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  pw  =  sea pressure of the pure water side                      [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  pw may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  osmotic_pressure_t_exact  =  osmotic pressure of seawater       [ dbar ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.41 of this TEOS-10 Manual.
%
%  McDougall T.J. and S.J. Wotherspoon, 2013: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_osmotic_pressure_t_exact(double SA, double t, double pw)
{
	double p_mean, f, p_old, rec_g2kg = 1e3, df_dp, p, gibbs_pure_water, sa = 0.0;

	if (SA > 0.0) sa = SA;

	gibbs_pure_water = gsw_gibbs(0, 0, 0, 0.0, t, pw);

	p = pw + 235.4684;        //% Initial guess of p, in dbar

	df_dp = -gtc.db2pa * (gsw_gibbs(0, 0, 1, sa, t, p) -
	            sa * gsw_gibbs(1, 0, 1, sa, t, p)); //% Inital guess of df/dp

	for (unsigned iter = 0; iter < 2; iter++)
	{
		p_old = p;
		f = gibbs_pure_water - rec_g2kg * gsw_chem_potential_water_t_exact(sa, t, p_old);
		p = p_old - f / df_dp; //% this is half way through the modified N-R method
		p_mean = 0.5 * (p + p_old);

		df_dp = -gtc.db2pa * (gsw_gibbs(0, 0, 1, sa, t, p_mean) -
		            sa * gsw_gibbs(1, 0, 1, sa, t, p_mean));

		p = p_old - f / df_dp;
	}

	//% After two iterations though the modified Newton-Raphson technique
	//% (McDougall and Wotherspoon, 2013) the maximum error is 6x10^-12 dbar.

	return p - pw; //% osmotic pressure of seawater, in dbar.
}
/**
% gsw_alpha_CT_exact                          thermal expansion coefficient
%                                  with respect to Conservative Temperature
%==========================================================================
%
% USAGE: gsw_alpha_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of seawater with respect to
%  Conservative Temperature from Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_alpha(SA,CT,p)
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  alpha_CT_exact  =  thermal expansion coefficient                 [ 1/K ]
%                     with respect to Conservative Temperature
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.18.3) of this TEOS-10 manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_alpha_CT_exact(double SA, double CT, double p)
{
	double t = gsw_t_from_ct(SA, CT, p);

	return gsw_alpha_wrt_CT_t_exact(SA, t, p);
}
/**
% gsw_specvol_second_derivatives_CT_exact             second derivatives of
%                                                           specific volume
% =========================================================================
%
% USAGE: gsw_specvol_second_derivatives_CT_exact(SA,CT,p, *v_SA_SA,
                                                           *v_SA_CT,
                                                           *v_CT_CT,
                                                           *v_SA_P,
                                                           *v_CT_P)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of specific
%  volume (v),
%   (1) v_SA_SA, second-order derivative with respect to Absolute Salinity
%       at constant CT & p.
%   (2) v_SA_CT, second-order derivative with respect to SA & CT at
%       constant p.
%   (3) v_CT_CT, second-order derivative with respect to CT at constant SA
%       and p.
%   (4) v_SA_P, second-order derivative with respect to SA & P at
%       constant CT.
%   (5) v_CT_P, second-order derivative with respect to CT & P at
%       constant SA.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_specvol_second_derivatives(SA,CT,p) which uses the computationally
%  efficient 75 term expression for specific volume in terms of SA, CT
%  and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  v_SA_SA  =  The second derivative of specific volume with respect to
%              Absolute Salinity at constant CT & p.  [ (m^3/kg)(g/kg)^-2 ]
%  v_SA_CT  =  The second derivative of specific volume with respect to
%              SA and CT at constant p.           [ (m^3/kg)(g/kg)^-1 K^-1]
%  v_CT_CT  =  The second derivative of specific volume with respect to
%              CT at constant SA and p.                  [ (m^3/kg) K^-2) ]
%  v_SA_P  =  The second derivative of specific volume with respect to
%              SA and P at constant CT.                  [ (m^3/kg) Pa^-1 ]
%  v_CT_P  =  The second derivative of specific volume with respect to
%              CT and P at constant SA.             [ (m^3/kg) K^-1 Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_specvol_second_derivatives_CT_exact(double SA,
		double CT,
		double p,
		double &v_SA_SA,
		double &v_SA_CT,
		double &v_CT_CT,
		double &v_SA_P,
		double &v_CT_P)
{
	double rec_gTT_pt0, rec_gTT2, rec_gTT, gamma, rec_cp0, cp0_div_abs_pt0, pt0;
	double gTTP,gTTT,gSAT,gSATP,gSATT,gSAT_pt0,gSA_pt0,gSASA_pt0,gSASAP,gSASAT,gSAPP,gTPP;
	double t, rec_abs_pt0, pr0=0.0, part_a, part_b, part_c, sa = 0.0;

	if (SA > 0.0) sa = SA;

	rec_cp0 = 1.0 / gtc.gsw_cp0; // from Eqn. 3.3.3 of IOC et al. (2010).
	pt0 = gsw_pt_from_ct(sa, CT);
	rec_abs_pt0 = 1.0 / (gtc.gsw_t0 + pt0);
	cp0_div_abs_pt0 = gtc.gsw_cp0 * rec_abs_pt0;
	t = gsw_pt_from_t(sa, pt0, pr0, p);

	gamma = -gsw_gibbs(0, 1, 1, sa, t,p) / gsw_gibbs(0, 2, 0, sa, t,p);

	rec_gTT = 1.0 / gsw_gibbs(0, 2, 0, sa, t, p);
	rec_gTT2 = rec_gTT * rec_gTT;
	rec_gTT_pt0 = 1.0 / gsw_gibbs(0, 2, 0, sa, pt0, pr0);

	gTTP = gsw_gibbs(0,2,1, sa,t,p);
	gTTT = gsw_gibbs(0,3,0, sa,t,p);
	gSAT = gsw_gibbs(1,1,0, sa,t,p);
	gSATP = gsw_gibbs(1,1,1, sa,t,p);
	gSATT = gsw_gibbs(1,2,0, sa,t,p);
	gSAT_pt0 = gsw_gibbs(1,1,0, sa,pt0,pr0);
	gSA_pt0 = gsw_gibbs(1,0,0, sa,pt0,pr0);
	gSASA_pt0 = gsw_gibbs(2,0,0, sa,pt0,pr0);
	gSASAP = gsw_gibbs(2,0,1, sa,t,p);
	gSASAT = gsw_gibbs(2,1,0, sa,t,p);
	gSAPP = gsw_gibbs(1,0,2, sa,t,p);
	gTPP = gsw_gibbs(0,1,2, sa,t,p);

	part_a = (gTTP + gamma * gTTT) * rec_gTT2;
	part_b = (gSATP + gamma * gSATT) * rec_gTT;
	part_c = (gTPP + gamma * (2.0 * gTTP + gamma * gTTT)) * rec_gTT;

	v_CT_CT = (cp0_div_abs_pt0 * cp0_div_abs_pt0) * (gamma * rec_abs_pt0 * rec_gTT_pt0 + part_a);

	v_SA_CT = cp0_div_abs_pt0 *
				 (gamma * rec_abs_pt0 * gSAT_pt0 *
				 rec_gTT_pt0 - part_b + gSAT * part_a) -
				 gSA_pt0 * v_CT_CT * rec_cp0;

	v_SA_SA = gSASAP + gamma * gSASAT
				 - gamma * rec_abs_pt0 * gSASA_pt0
				 + gamma * rec_abs_pt0 * (gSAT_pt0 * gSAT_pt0) * rec_gTT_pt0
				 - 2.0 * gSAT * part_b
				 + (gSAT * gSAT) * part_a
				 - 2.0 * gSA_pt0 * rec_cp0 * v_SA_CT - (gSA_pt0 * gSA_pt0) *
				 (rec_cp0 * rec_cp0) * v_CT_CT;

	v_CT_P = -cp0_div_abs_pt0 * part_c;

	v_SA_P = gSAPP + gamma * (2.0 * gSATP + gamma *gSATT) +
				   part_c * (gSA_pt0 *rec_abs_pt0 - gSAT);
}
/**
% gsw_sigma4_CT_exact                        potential density anomaly with
%                                       reference sea pressure of 4000 dbar
%==========================================================================
%
% USAGE: gsw_sigma4_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 4000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma4(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma4_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 4000 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sigma4_CT_exact(double SA, double CT)
{
	double t = gsw_t_from_ct(SA, CT, 4000.0);

	return gsw_rho_t_exact(SA, t, 4000.0) - 1000.0;
}
/**
% gsw_specvol_anom_standard_t_exact                 specific volume anomaly
%==========================================================================
%
% USAGE: gsw_specvol_anom_standard_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates specific volume anomaly from Absolute Salinity, in-situ
%  temperature and pressure, using the full TEOS-10 Gibbs function.  The
%  reference value of Absolute Salinity is SSO and the temperature of the
%  reference value is equal to a Conservative Temperature of 0 degrees C.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions,
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  specvol_anom_t_exact  =  specific volume anomaly              [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_anom_standard_t_exact(double SA, double t, double p)
{
	double pr0 = 0.0, CT0 = pr0, t_zero, sa = 0.0;

	if (SA > 0.0) sa = SA;

	t_zero = gsw_t_from_ct(gtc.gsw_sso, CT0, p);

	return gsw_gibbs(0, 0, 1, sa, t, p) -
	         gsw_gibbs(0, 0, 1, gtc.gsw_sso, t_zero, p);
}

/**
% gsw_R_from_SP                                  conductivity ratio from SP
%==========================================================================
%
% USAGE: gsw_R_from_SP(SP,t,p)
%
% DESCRIPTION:
%  Calculates conductivity ratio from (SP,t,p) using PSS-78 in the range
%  2 < SP < 42.  If the input Practical Salinity is less than 2 then a
%  modified form of the Hill et al. (1986) fomula is used for Practical
%  Salinity.  The modification of the Hill et al. (1986) expression is to
%  ensure that it is exactly consistent with PSS-78 at SP = 2.
%
%  The conductivity ratio returned by this function is consistent with the
%  input value of Practical Salinity, SP, to 2x10^-14 psu over the full
%  range of input parameters (from pure fresh water up to SP = 42 psu).
%  This error of 2x10^-14 psu is machine precision at typical seawater
%  salinities.  This accuracy is achieved by having four different
%  polynomials for the starting value of Rtx (the square root of Rt) in
%  four different ranges of SP, and by using one and a half iterations of
%  a computationally efficient modified Newton-Raphson technique
%  (McDougall and Wotherspoon, 2013) to find the root of the equation.
%
%  Note that strictly speaking PSS-78 (Unesco, 1983) defines Practical
%  Salinity in terms of the conductivity ratio, R, without actually
%  specifying the value of C(35,15,0) (which we currently take to be
%  42.9140 mS cm^-1 (Culkin and Smith, 1980)).
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SP & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SP & t are MxN.
%
% OUTPUT:
%  R  =  conductivity ratio                                    [ unitless ]
%
% AUTHOR:
%  Trevor McDougall, Paul Barker and Rich Pawlowicz    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Culkin and Smith, 1980:  Determination of the Concentration of Potassium
%   Chloride Solution Having the Same Electrical Conductivity, at 15C and
%   Infinite Frequency, as Standard Seawater of Salinity 35.0000
%   (Chlorinity 19.37394), IEEE J. Oceanic Eng, 5, 22-23.
%
%  Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the
%   Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng.,
%   11, 109 - 112.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix E of this TEOS-10 Manual.
%
%  McDougall T.J., and S.J. Wotherspoon, 2013: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  Unesco, 1983: Algorithms for computation of fundamental properties of
%   seawater. Unesco Technical Papers in Marine Science, 44, 53 pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_R_from_SP(double SP, double t, double p)
{
	double p0 =   4.577801212923119e-3;
	double p1 =   1.924049429136640e-1;
	double p2 =   2.183871685127932e-5;
	double p3 =  -7.292156330457999e-3;
	double p4 =   1.568129536470258e-4;
	double p5 =  -1.478995271680869e-6;
	double p6 =   9.086442524716395e-4;
	double p7 =  -1.949560839540487e-5;
	double p8 =  -3.223058111118377e-6;
	double p9 =   1.175871639741131e-7;
	double p10 = -7.522895856600089e-5;
	double p11 = -2.254458513439107e-6;
	double p12 =  6.179992190192848e-7;
	double p13 =  1.005054226996868e-8;
	double p14 = -1.923745566122602e-9;
	double p15 =  2.259550611212616e-6;
	double p16 =  1.631749165091437e-7;
	double p17 = -5.931857989915256e-9;
	double p18 = -4.693392029005252e-9;
	double p19 =  2.571854839274148e-10;
	double p20 =  4.198786822861038e-12;

	double q0 =   5.540896868127855e-5;
	double q1 =   2.015419291097848e-1;
	double q2 =  -1.445310045430192e-5;
	double q3 =  -1.567047628411722e-2;
	double q4 =   2.464756294660119e-4;
	double q5 =  -2.575458304732166e-7;
	double q6 =   5.071449842454419e-3;
	double q7 =  -9.081985795339206e-5;
	double q8 =  -3.635420818812898e-6;
	double q9 =   2.249490528450555e-8;
	double q10 = -1.143810377431888e-3;
	double q11 =  2.066112484281530e-5;
	double q12 =  7.482907137737503e-7;
	double q13 =  4.019321577844724e-8;
	double q14 = -5.755568141370501e-10;
	double q15 =  1.120748754429459e-4;
	double q16 = -2.420274029674485e-6;
	double q17 = -4.774829347564670e-8;
	double q18 = -4.279037686797859e-9;
	double q19 = -2.045829202713288e-10;
	double q20 =  5.025109163112005e-12;

	double r0 =   3.432285006604888e-3;
	double r1 =   1.672940491817403e-1;
	double r2 =   2.640304401023995e-5;
	double r3 =   1.082267090441036e-1;
	double r4 =  -6.296778883666940e-5;
	double r5 =  -4.542775152303671e-7;
	double r6 =  -1.859711038699727e-1;
	double r7 =   7.659006320303959e-4;
	double r8 =  -4.794661268817618e-7;
	double r9 =   8.093368602891911e-9;
	double r10 =  1.001140606840692e-1;
	double r11 = -1.038712945546608e-3;
	double r12 = -6.227915160991074e-6;
	double r13 =  2.798564479737090e-8;
	double r14 = -1.343623657549961e-10;
	double r15 =  1.024345179842964e-2;
	double r16 =  4.981135430579384e-4;
	double r17 =  4.466087528793912e-6;
	double r18 =  1.960872795577774e-8;
	double r19 = -2.723159418888634e-10;
	double r20 =  1.122200786423241e-12;

	double u0 =    5.180529787390576e-3;
	double u1 =    1.052097167201052e-3;
	double u2 =    3.666193708310848e-5;
	double u3 =    7.112223828976632;
	double u4 =   -3.631366777096209e-4;
	double u5 =   -7.336295318742821e-7;
	double u6 =   -1.576886793288888e+2;
	double u7 =   -1.840239113483083e-3;
	double u8 =    8.624279120240952e-6;
	double u9 =    1.233529799729501e-8;
	double u10 =   1.826482800939545e+3;
	double u11 =   1.633903983457674e-1;
	double u12 =  -9.201096427222349e-5;
	double u13 =  -9.187900959754842e-8;
	double u14 =  -1.442010369809705e-10;
	double u15 =  -8.542357182595853e+3;
	double u16 =  -1.408635241899082;
	double u17 =   1.660164829963661e-4;
	double u18 =   6.797409608973845e-7;
	double u19 =   3.345074990451475e-10;
	double u20 =   8.285687652694768e-13;

   double Rtx_old, SP_Hill_raw, SP_est, Hill_ratio, part1, part2;
   double Ra, rt_lc, A, B, C, D, E, Rt, Rtxm, sqrty, x, Rtx=0.0, dSP_dRtx;
	double t68 = t * 1.00024;
	double ft68 = (t68 - 15.0) / (1.0 + gspc.k * (t68 - 15.0));

   x = sqrt(SP);

	//--------------------------------------------------------------------------
	// Finding the starting value of Rtx, the square root of Rt, using four
	// different polynomials of SP and t68.
	//--------------------------------------------------------------------------

   if (SP >= 9)
   {
	   Rtx =  p0 + x *(p1 + p4 * t68 + x * (p3 + p7 * t68
	               + x * (p6 + p11 * t68 + x * (p10 + p16 * t68
	               + x * p15)))) + t68 *(p2+ t68 *(p5
	               + x * x * (p12 + x *p17) + p8 * x + t68 * (p9
                  + x * (p13 + x *p18)+ t68 *(p14 + p19 * x + p20 * t68))));
   }

   if (SP >= 0.25 && SP < 9)
   {
	   Rtx =  q0 + x *(q1 + q4 * t68
	               + x * (q3 + q7 * t68 + x *(q6 + q11 * t68
                  + x * (q10 + q16*t68+ x *q15)))) + t68 * (q2 + t68 *(q5
                  + x * x * (q12 + x * q17) + q8 * x + t68 *(q9
                  + x * (q13 + x * q18) + t68 * (q14 + q19 * x + q20 * t68))));
   }

   if (SP >= 0.003 && SP < 0.25)
   {
	   Rtx =  r0 + x *(r1 + r4 * t68
	               + x * (r3 + r7*t68 + x *(r6 + r11*t68
	               + x * (r10 + r16*t68 + x * r15)))) + t68 * (r2 + t68 * (r5
	               + x * x * (r12 + x *r17) + r8 * x + t68 * (r9
                  + x *(r13 + x * r18) + t68 * (r14 + r19 *x + r20 * t68))));
   }

   if (SP < 0.003)
   {
	   Rtx =  u0 + x *(u1 + u4*t68
	               + x * (u3 + u7*t68
	               + x * (u6 + u11*t68
	               + x * (u10 + u16 * t68
	               + x * u15)))) + t68 * (u2 + t68 * (u5
	               + x * x * (u12 + x * u17) + u8 * x + t68 *(u9
	               + x * (u13 + x * u18) + t68 * (u14 + u19 * x + u20 * t68))));
   }

	//--------------------------------------------------------------------------
	// Finding the starting value of dSP_dRtx, the derivative of SP with respect
	// to Rtx.
	//--------------------------------------------------------------------------

   dSP_dRtx =  gspc.a1
               + (2.0 * gspc.a2
               + (3.0 * gspc.a3
               + (4.0 * gspc.a4
               + 5.0 * gspc.a5 * Rtx) * Rtx) * Rtx) * Rtx
               + ft68 * (gspc.b1
               + (2.0 * gspc.b2
               + (3.0 * gspc.b3
               + (4.0 * gspc.b4
               + 5.0 * gspc.b5 * Rtx) * Rtx) * Rtx) * Rtx);

   if (SP < 2.0)
   {
	   x = 400.0 * (Rtx * Rtx);
	   sqrty = 10.0 *Rtx;
	   part1 = 1.0 + x * (1.5 + x);
	   part2 = 1.0 + sqrty * (1.0 + sqrty * (1.0 + sqrty));
	   Hill_ratio = gsw_hill_ratio_at_sp2(t);

	   dSP_dRtx = dSP_dRtx + gspc.a0 * 800.0 * Rtx * (1.5 + 2.0 * x) / (part1 * part1)
                  + gspc.b0 * ft68 *
                  (10.0 + sqrty * (20.0 + 30.0 * sqrty)) / (part2 * part2);

	   dSP_dRtx = Hill_ratio * dSP_dRtx;
   }

	//--------------------------------------------------------------------------
	// One iteration through the modified Newton-Raphson method achieves an
	// error in Practical Salinity of about 10^-12 for all combinations of the
	// inputs.  One and a half iterations of the modified Newton-Raphson method
	// achevies a maximum error in terms of Practical Salinity of better than
	// 2x10^-14 everywhere.
	//
	// We recommend one and a half iterations of the modified Newton-Raphson
	// method.
	//
	// Begin the modified Newton-Raphson method.
	//--------------------------------------------------------------------------

   SP_est = gspc.a0 + (gspc.a1 + (gspc.a2 + (gspc.a3
            + (gspc.a4 + gspc.a5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx
			   + ft68 * (gspc.b0 + (gspc.b1 + (gspc.b2 + (gspc.b3
			   + (gspc.b4 + gspc.b5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx);

   if (SP_est < 2.0)
   {
	   x = 400.0 * (Rtx * Rtx);
	   sqrty = 10.0 * Rtx;
	   part1 = 1.0 + x * (1.5 + x) ;
	   part2 = 1.0 + sqrty * (1.0 + sqrty * (1.0 + sqrty));
	   SP_Hill_raw = SP_est - gspc.a0 / part1 - gspc.b0 * ft68 / part2;

	   Hill_ratio = gsw_hill_ratio_at_sp2(t);

	   SP_est = Hill_ratio * SP_Hill_raw;
   }

   Rtx_old = Rtx;
   Rtx = Rtx_old - (SP_est - SP) / dSP_dRtx;
   Rtxm = 0.5 * (Rtx + Rtx_old);
   // This mean value of Rtx, Rtxm, is the
	//	value of Rtx at which the derivative dSP_dRtx is evaluated.

   dSP_dRtx =  gspc.a1
               + (2.0 * gspc.a2
               + (3.0 * gspc.a3
               + (4.0 * gspc.a4
               + 5.0 * gspc.a5 * Rtxm) * Rtxm) * Rtxm) * Rtxm
               + ft68 * (gspc.b1
               + (2.0 * gspc.b2
               + (3.0 * gspc.b3
               + (4.0 * gspc.b4
               + 5.0 * gspc.b5 * Rtxm) * Rtxm) * Rtxm) * Rtxm);

   if (SP_est < 2.0)
   {
	   x = 400 * (Rtxm *Rtxm);
	   sqrty = 10.0 * Rtxm;
	   part1 = 1.0 + x * (1.5 + x) ;
	   part2 = 1.0 + sqrty * (1.0 + sqrty * (1.0 + sqrty));

	   dSP_dRtx = dSP_dRtx + gspc.a0 *800.0 *Rtxm *(1.5 + 2.0 * x) / (part1 * part1)
				      + gspc.b0 *ft68 *
				      (10.0 + sqrty * (20.0 + 30.0 * sqrty)) / (part2 * part2);

	   Hill_ratio = gsw_hill_ratio_at_sp2(t);

	   dSP_dRtx = Hill_ratio * dSP_dRtx;
   }

	//--------------------------------------------------------------------------
	// The line below is where Rtx is updated at the end of the one full
	// iteration of the modified Newton-Raphson technique.
	//--------------------------------------------------------------------------
   Rtx = Rtx_old - (SP_est - SP) /dSP_dRtx;
	//--------------------------------------------------------------------------
	// Now we do another half iteration of the modified Newton-Raphson
	// technique, making a total of one and a half modified N-R iterations.
	//--------------------------------------------------------------------------
   SP_est = gspc.a0
            + (gspc.a1 + (gspc.a2
            + (gspc.a3 + (gspc.a4
            + gspc.a5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx
			   + ft68 *(gspc.b0
			   + (gspc.b1 + (gspc.b2
			   + (gspc.b3 + (gspc.b4
			   + gspc.b5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx);

   if (SP_est < 2.0)
   {
	   x = 400.0 * (Rtx * Rtx);
	   sqrty = 10.0 * Rtx;
	   part1 = 1.0 + x * (1.5 + x) ;
	   part2 = 1.0 + sqrty * (1.0 + sqrty * (1.0 + sqrty));
	   SP_Hill_raw = SP_est - gspc.a0 / part1 - gspc.b0 * ft68 / part2;

	   Hill_ratio = gsw_hill_ratio_at_sp2(t);

	   SP_est = Hill_ratio * SP_Hill_raw;
   }

   Rtx = Rtx - (SP_est - SP) / dSP_dRtx;

	//--------------------------------------------------------------------------
	// The following lines of code are commented out, but when activated, return
	// the error, SP_error, in Rtx (in terms of psu).
	// SP_est = a0 + (a1 + (a2 + (a3 + (a4 + a5 *Rtx) *Rtx) *Rtx) *Rtx) *Rtx
	//     + ft68 *(b0 + (b1 + (b2+ (b3 + (b4 + b5 *Rtx) *Rtx) *Rtx) *Rtx) *Rtx);
	// if any(SP_est < 2)
	//     [I2] = find(SP_est < 2);
	//     x = 400 *(Rtx *Rtx);
	//     sqrty = 10 *Rtx;
	//     part1 = 1 + x *(1.5 + x) ;
	//     part2 = 1 + sqrty *(1 + sqrty *(1 + sqrty));
	//     SP_Hill_raw = SP_est - a0 /part1 - b0 *ft68 /part2;
	//     Hill_ratio = gsw_Hill_ratio_at_SP2(t);
	//     SP_est = Hill_ratio *SP_Hill_raw;
	// end
	//
	// SP_error = abs(SP - SP_est);
	//
	//--------------This is the end of the error testing------------------------


	//--------------------------------------------------------------------------
	// Now go from Rtx to Rt and then to the conductivity ratio R at pressure p.
	//--------------------------------------------------------------------------
   Rt = Rtx * Rtx;
   A  = gspc.d3 + gspc.d4 * t68;
   B  = 1.0 + gspc.d1 * t68 + gspc.d2 * (t68 * t68);
   C  = p * (gspc.e1 + gspc.e2 * p + gspc.e3 * (p * p));
	// rt_lc (i.e. rt_lower_case) corresponds to rt as defined in
	// the UNESCO 44 (1983) routines.
   rt_lc = gspc.c0
            + (gspc.c1
            + (gspc.c2
            + (gspc.c3
            + gspc.c4 * t68) * t68) * t68) * t68;

   D  = B - A * rt_lc * Rt;
   E  = rt_lc * Rt * A * (B + C);
   Ra = sqrt((D * D) + 4.0 * E) - D;

   return 0.5 * Ra / A;
}
/**
% gsw_specvol_alpha_beta_CT_exact        specific volume, thermal expansion
%                                          & saline contraction coefficient
%==========================================================================
%
% USAGE: gsw_specvol_alpha_beta_CT_exact(SA,CT,p, *specvol_CT_exact,
                                                   *alpha_CT_exact,
                                                   *beta_CT_exact)
%
% DESCRIPTION:
%  Calculates specific volume, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from
%  Absolute Salinity and Conservative Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_specvol_alpha_beta(SA,CT,p), which uses the computationally
%  efficient 75-term expression for density in terms of SA, CT and p
%  (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  specvol_CT_exact =  specific volume                             [ m/kg ]
%  alpha_CT_exact   =  thermal expansion coefficient                [ 1/K ]
%                      with respect to Conservative Temperature
%  beta_CT_exact    =  saline (i.e. haline) contraction            [ kg/g ]
%                      coefficient at constant Conservative Temperature
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_specvol_alpha_beta_CT_exact(double SA, double CT, double p,
		double &specvol_CT_exact,
		double &alpha_CT_exact,
		double &beta_CT_exact)
{
	double t = gsw_t_from_ct(SA, CT, p);

	specvol_CT_exact = gsw_specvol_t_exact(SA, t, p);
	alpha_CT_exact = gsw_alpha_wrt_CT_t_exact(SA, t, p);
	beta_CT_exact = gsw_beta_const_CT_t_exact(SA, t, p);
}
/**
% gsw_t_from_rho_exact                     in-situ temperature from density
% =========================================================================
%
% USAGE: gsw_t_from_rho_exact(rho,SA,p, *t, *t_multiple)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar).
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg m^-3 subtracted from it.
%     That is, it is |density|, not |density anomaly|.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & SA need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & SA are MxN.
%
% OUTPUT:
%  t  =  in-situ temperature  (ITS-90)                            [ deg C ]
%  t_multiple  =  in-situ temperature  (ITS-90)                   [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      temperatures for a single density.  This programme will output both
%      valid solutions.  To see this second solution the user must call the
%      programme with two outputs (i.e. [t,t_multiple]), if there is only
%      one possible solution and the programme has been called with two
%      outputs the second variable will be set to NaN.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_t_from_rho_exact(double rho, double SA, double p,
													double &t,
													double &t_multiple)
{
	double sa, P, rho_40;

	/** check range, if OOR set to NAN */
	if (SA < 0.0 || SA > 42.0)
	{
	   sa = NAN;
   }
   else
   {
      sa = SA;
   }

	if (p < -1.5 || p > 12000.0) P = NAN;

	rho_40 = gsw_rho_t_exact(sa, 40.0 * sa, p);
	if (rho - rho_40 < 0.0) rho_40 = NAN;

	if (sa == NAN || P == NAN || rho_40 == NAN)
	{
		t = t_multiple = NAN;

		return;
	}

	double t_max_rho, max_rho, v_t_t, rho_t_t, deriv_lin, discrim;
	double a, b, c, v_t, rho_t, t_old, t_upper, t_freezing;

	t_max_rho = gsw_t_maxdensity_exact(sa, p);
	max_rho = gsw_rho_t_exact(sa, t_max_rho, p);
	v_t_t = gsw_gibbs(0, 2, 1, sa, t_max_rho, p);
	rho_t_t = -v_t_t * rho * rho;   //This is approximate

	deriv_lin = rho_40 - max_rho / 40.0 - t_max_rho;
	discrim = (deriv_lin * deriv_lin) - 4.0 * 0.5 * rho_t_t * (max_rho - rho);

	if (discrim < 0.0) discrim = 0.0;

	t = t_max_rho + 2.0 * (max_rho - rho / -deriv_lin + sqrt(discrim));

	/**
	   Having found this initial value of t, begin the iterative solution
	   for the t_upper part, the solution warmer than the TMD
	*/
	for (unsigned iter = 0; iter < 3; iter++)
	{
		v_t = gsw_gibbs(0,1,1,sa,t,p);
		rho_t = -v_t * rho * rho;        // This is approximate

		v_t_t = gsw_gibbs(0,2,1,sa,t,p);
		rho_t_t = -v_t_t *rho *rho;      //This is approximate

		b = rho_t;
		a = 0.5 *rho_t_t;
		c = gsw_rho_t_exact(sa,t,p) - rho;
		discrim = b * b - 4.0 * a * c;
		if (discrim < 0.0) discrim = 0.0;
		t_old = t;
		t = t_old + (2.0 * c / -b + sqrt(discrim));
	} /** for */

	t_upper = t;

	/** Now start the t_multiple part, the solution cooler than the TMD */
	t = 2.0 * t_max_rho - t_upper;

	for (unsigned iter = 0; iter < 6; iter++)
	{
		v_t = gsw_gibbs(0,1,1,sa,t,p);
		rho_t = -v_t * rho * rho;           //This is approximate

		v_t_t = gsw_gibbs(0,2,1,sa,t,p);
		rho_t_t = -v_t_t * rho * rho;       //This is approximate

		b = rho_t;
		a = 0.5 * rho_t_t;
		c = gsw_rho_t_exact(sa,t,p) - rho;
		discrim = b * b - 4.0 * a * c;
		if (discrim < 0.0) discrim = 0.0;
		t_old = t;
		t = t_old + (2.0 * c / -b - sqrt(discrim)); //Note the sign change of the sqrt term
	} /** for */

	t_multiple = t;

	/**This assumes that the seawater is unsaturated with air */
	t_freezing = gsw_t_freezing(sa, p, 0.0);

	/** Set values outside the relevant bounds to NaNs */
	if (t_upper < t_freezing || t_upper < t_max_rho) t_upper = NAN;
	if (t_multiple > t_max_rho || t_multiple > t_max_rho) t_multiple = NAN;

	t = t_upper;
}
/**
% gsw_pot_enthalpy_from_specvol_ice                 potential enthalpy from
%                                                    specific volume of ice
%==========================================================================
%
% USAGE: gsw_pot_enthalpy_from_specvol_ice(specvol_ice,p)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from the specific volume
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar.
%
% INPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_pot_enthalpy_from_specvol_ice(double specvol_ice, double p)
{
	double rho_ice = 1.0 / specvol_ice;
	double t_ice = gsw_t_from_rho_ice(rho_ice,p);
	double pt0_ice = gsw_pt0_from_t_ice(t_ice,p);

	return gsw_pot_enthalpy_from_pt_ice(pt0_ice);
}
/**
% gsw_alpha_on_beta_CT_exact                                     alpha/beta
%==========================================================================
%
% USAGE: gsw_alpha_on_beta_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates alpha divided by beta, where alpha is the thermal expansion
%  coefficient and beta is the saline contraction coefficient of seawater
%  from Absolute Salinity and Conservative Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_alpha_on_beta(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  alpha_on_beta_CT_exact  =  thermal expansion coefficient with respect to
%                             Conservative Temperature divided by the
%                             saline contraction coefficient at constant
%                             Conservative Temperature     [ kg g^-1 K^-1 ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_alpha_on_beta_CT_exact(double SA, double CT, double p)
{
	double dummy, alpha_CT_exact, beta_CT_exact, sa = 0.0;

	if (SA > 0.0) sa = SA;

	gsw_specvol_alpha_beta_CT_exact(sa,CT,p,
											  dummy,
											  alpha_CT_exact,
											  beta_CT_exact);

	return alpha_CT_exact / beta_CT_exact;
}
/**
% gsw_f                                              Coriolis parameter (f)
%==========================================================================
%
% USAGE: gsw_f(lat)
%
% DESCRIPTION:
%  Calculates the Coriolis parameter (f) defined by:
%       f = 2*omega*sin(lat)
%  where,
%       omega = 7.292115e-5   (Groten, 2004)                  [ radians/s ]
%
% INPUT:
%  lat  =  latitude in decimal degrees North                [ -90 +90 ]
%
% OUTPUT:
%  f    =  Coriolis parameter                                 [ radians/s ]
%
% AUTHOR:
%  20th April 1993. Phil Morgan                        [ help@teos-10.org ]
%
% MODIFIED:
%  28th July, 2010 by Paul Barker
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCE:
%  Groten, E., 2004: Fundamental Parameters and Current (2004) Best
%   Estimates of the Parameters of Common Relevance to Astronomy, Geodesy,
%   and Geodynamics. Journal of Geodesy, 77, pp. 724-797.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_f(double lat)
{
	return 0.0001458423 * sin(lat * 0.017453292519943295);
}
/**
% gsw_rho_alpha_beta_indexed    in-situ density, thermal expansion & saline
%                                contraction coefficient (75-term equation)
%==========================================================================
%
% USAGE: gsw_rho_alpha_beta_indexed(SA,CT,p,double &rho,
                                double &alpha, double &beta)
%
% DESCRIPTION:
%  Calculates in-situ density, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from
%  Absolute Salinity and Conservative Temperature.  This function uses the
%  computationally-efficient expression for specific volume in terms of
%  SA, CT and p (Roquet et al., 2015).
%
%  Note that potential density (pot_rho) with respect to reference pressure
%  p_ref is obtained by calling this function with the pressure argument
%  being p_ref as in [pot_rho, ~, ~] = gsw_rho_alpha_beta(SA,CT,p_ref).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho    =  in-situ density                                       [ kg/m ]
%  alpha  =  thermal expansion coefficient                          [ 1/K ]
%            with respect to Conservative Temperature
%  beta   =  saline (i.e. haline) contraction                      [ kg/g ]
%            coefficient at constant Conservative Temperature
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_rho_alpha_beta_indexed(double SA, double CT, double p,
		double &rho,
		double &alpha,
		double &beta)
{
	double v_SA_part, v_CT, z, v, ys, xs, x2, sa = 0.0;

	if (SA > 0.0) sa = SA;

	x2 = gtc.gsw_sfac * sa;
	xs = sqrt(x2 + gtc.offset);
	ys = CT * 0.025;
	z = p * gtc.rec_db2pa;

	/* gsvco */
	v = gsvco.v000 + xs *(gsvco.v100 + xs *(gsvco.v200 + xs *(gsvco.v300 + xs *(gsvco.v400 + xs *(gsvco.v500
	 + gsvco.v600 *xs))))) + ys *(gsvco.v010 + xs *(gsvco.v110 + xs *(gsvco.v210 + xs *(gsvco.v310 + xs *(gsvco.v410
	 + gsvco.v510 *xs)))) + ys *(gsvco.v020 + xs *(gsvco.v120 + xs *(gsvco.v220 + xs *(gsvco.v320 + gsvco.v420 *xs)))
	 + ys *(gsvco.v030 + xs *(gsvco.v130 + xs *(gsvco.v230 + gsvco.v330 *xs)) + ys *(gsvco.v040 + xs *(gsvco.v140
	 + gsvco.v240 * xs) + ys * (gsvco.v050 + gsvco.v150 * xs + gsvco.v060 * ys))))) + z *(gsvco.v001 + xs *(gsvco.v101
	 + xs *(gsvco.v201 + xs *(gsvco.v301 + xs *(gsvco.v401 + gsvco.v501 *xs)))) + ys *(gsvco.v011 + xs *(gsvco.v111
	 + xs *(gsvco.v211 + xs *(gsvco.v311 + gsvco.v411 *xs))) + ys *(gsvco.v021 + xs *(gsvco.v121 + xs *(gsvco.v221
	 + gsvco.v321 *xs)) + ys *(gsvco.v031 + xs *(gsvco.v131 + gsvco.v231 *xs) + ys *(gsvco.v041 + gsvco.v141 *xs
	 + gsvco.v051 * ys)))) + z *(gsvco.v002 + xs *(gsvco.v102 + xs *(gsvco.v202 + xs *(gsvco.v302 + gsvco.v402 *xs)))
	 + ys *(gsvco.v012 + xs *(gsvco.v112 + xs *(gsvco.v212 + gsvco.v312 *xs)) + ys *(gsvco.v022 + xs *(gsvco.v122
	 + gsvco.v222 *xs) + ys *(gsvco.v032 + gsvco.v132 *xs + gsvco.v042 * ys))) + z *(gsvco.v003 + xs *(gsvco.v103
	 + gsvco.v203 *xs) + ys *(gsvco.v013 + gsvco.v113 *xs + gsvco.v023 *ys) + z *(gsvco.v004 + gsvco.v104 *xs
	 + gsvco.v014 *ys + z *(gsvco.v005 + gsvco.v006 *z)))));

	rho = 1.0 / v;

	v_CT = gsvco.a000 + xs *(gsvco.a100 + xs *(gsvco.a200 + xs *(gsvco.a300 + xs *(gsvco.a400 + gsvco.a500 *xs))))
		+ ys *(gsvco.a010 + xs *(gsvco.a110 + xs *(gsvco.a210 + xs *(gsvco.a310 + gsvco.a410 *xs)))
		+ ys *(gsvco.a020 + xs *(gsvco.a120 + xs *(gsvco.a220 + gsvco.a320 *xs)) + ys *(gsvco.a030
		+ xs *(gsvco.a130 + gsvco.a230 *xs) + ys *(gsvco.a040 + gsvco.a140 *xs + gsvco.a050 *ys ))))
		+ z *(gsvco.a001 + xs *(gsvco.a101 + xs *(gsvco.a201 + xs *(gsvco.a301 + gsvco.a401 *xs)))
		+ ys *(gsvco.a011 + xs *(gsvco.a111 + xs *(gsvco.a211 + gsvco.a311 *xs)) + ys *(gsvco.a021
		+ xs *(gsvco.a121 + gsvco.a221 *xs) + ys *(gsvco.a031 + gsvco.a131 *xs + gsvco.a041 *ys)))
		+ z *(gsvco.a002 + xs *(gsvco.a102 + xs *(gsvco.a202 + gsvco.a302 *xs)) + ys *(gsvco.a012
		+ xs *(gsvco.a112 + gsvco.a212 *xs) + ys *(gsvco.a022 + gsvco.a122 *xs + gsvco.a032 *ys))
		+ z *(gsvco.a003 + gsvco.a103 *xs + gsvco.a013 *ys + gsvco.a004 *z)));

	alpha = 0.025 * (v_CT / v);

	v_SA_part = gsvco.b000 + xs *(gsvco.b100 + xs *(gsvco.b200 + xs *(gsvco.b300
	   + xs * (gsvco.b400 + gsvco.b500 *xs)))) + ys *(gsvco.b010 + xs *(gsvco.b110
	   + xs * (gsvco.b210 + xs *(gsvco.b310 + gsvco.b410 *xs))) + ys *(gsvco.b020
	   + xs * (gsvco.b120 + xs *(gsvco.b220 + gsvco.b320 *xs)) + ys *(gsvco.b030
		+ xs * (gsvco.b130 + gsvco.b230 *xs) + ys *(gsvco.b040 + gsvco.b140 *xs + gsvco.b050 *ys))))
		+ z * (gsvco.b001 + xs *(gsvco.b101 + xs *(gsvco.b201 + xs *(gsvco.b301 + gsvco.b401 *xs)))
		+ ys * (gsvco.b011 + xs *(gsvco.b111 + xs *(gsvco.b211 + gsvco.b311 *xs)) + ys *(gsvco.b021
		+ xs * (gsvco.b121 + gsvco.b221 *xs) + ys *(gsvco.b031 + gsvco.b131 *xs + gsvco.b041 *ys)))
		+ z * (gsvco.b002 + xs *(gsvco.b102 + xs *(gsvco.b202 + gsvco.b302 *xs))+ ys *(gsvco.b012
		+ xs * (gsvco.b112 + gsvco.b212 *xs) + ys *(gsvco.b022 + gsvco.b122 *xs + gsvco.b032 *ys))
		+ z * (gsvco.b003 +  gsvco.b103 *xs + gsvco.b013 *ys + gsvco.b004 *z)));

	beta = -v_SA_part * (0.5 * (gtc.gsw_sfac / v * xs));
}
/**
% gsw_rho_first_derivatives_CT_exact       SA, CT and p partial derivatives
%                                                                of density
%==========================================================================
%
% USAGE: gsw_rho_first_derivatives_CT_exact(SA,CT,p,
                                             *rho_SA,
                                             *rho_CT,
                                             *rho_P)
%
% DESCRIPTION:
%  Calculates the three (3) partial derivatives of in-situ density with
%  respect to Absolute Salinity, Conservative Temperature and pressure.
%  Note that the pressure derivative is done with respect to pressure in
%  Pa, not dbar.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_rho_first_derivatives(SA,CT,p), which uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT and
%  p (Roquet et al., 2015)
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_SA  =  partial derivative of density           [ (kg/m^3)(g/kg)^-1 ]
%                 with respect to Absolute Salinity
%  rho_CT  =  partial derivative of density                  [ kg/(m^3 K) ]
%                 with respect to Conservative Temperature
%  rho_P   =  partial derivative of density                 [ kg/(m^3 Pa) ]
%                 with respect to pressure in Pa
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_rho_first_derivatives_CT_exact(double SA, double CT, double p,
		double &rho_SA, double &rho_CT,
		double &rho_P)
{
   double g_SA_mod_part, g_SA_T_mod, g_SA_T_mod_part;
   double x2, x, y, y_pt, z, pt0, t;
   double factor, factora, g_tp, g_psq_g_tt, g_p, g_SA_mod;
   double sa = 0.0;

   if (SA > 0.0) sa = SA;

   t = gsw_t_from_ct(sa,CT,p);

   pt0 = gsw_pt0_from_t(sa,t,p);

   x2 = gtc.gsw_sfac * sa;
   x = sqrt(x2);
   y = 0.025 * t;
   y_pt = 0.025 * pt0;

   z = gtc.rec_db2pa * p; //Note.The input pressure (p) is sea pressure in units of dbar.

   g_SA_T_mod_part = 1187.3715515697959 + z *(1458.233059470092 +
        z *(-687.913805923122 + z *(249.375342232496 + z *(-63.313928772146 + 14.09317606630898 *z)))) +
        x *(-1480.222530425046 + x *(2175.341332000392 + x *(-980.14153344888 + 220.542973797483 *x) +
        y *(-548.4580073635929 + y *(592.4012338275047 + y *(-274.2361238716608 + 49.9394019139016 *y))) -
        90.6734234051316 *z) + z *(-525.876123559641 + (249.57717834054571 - 88.449193048287 *z) *z) +
        y *(-258.3988055868252 + z *(2298.348396014856 + z *(-325.1503575102672 + 153.8390924339484 *z)) +
        y *(-90.2046337756875 - 4142.8793862113125 *z + y *(10.50720794170734 + 2814.78225133626 *z)))) +
        y *(3520.125411988816 + y *(-1351.605895580406 +
        y *(731.4083582010072 + y *(-216.60324087531103 + 25.56203650166196 *y) +
        z *(-2381.829935897496 + (597.809129110048 - 291.8983352012704 *z) *z)) +
        z *(4165.4688847996085 + z *(-1229.337851789418 + (681.370187043564 - 66.7696405958478 *z) *z))) +
        z *(-3443.057215135908 + z *(1349.638121077468 +
        z *(-713.258224830552 + (176.8161433232 - 31.68006188846728 *z) *z))));

   g_SA_T_mod = 0.5 * gtc.gsw_sfac * 0.025 * g_SA_T_mod_part;

   g_SA_mod_part = 8645.36753595126 +
        x * (-7296.43987145382 + x * (8103.20462414788 +
        y_pt * (2175.341332000392 + y_pt * (-274.2290036817964 +
        y_pt * (197.4670779425016 + y_pt * (-68.5590309679152 + 9.98788038278032 * y_pt)))) +
        x * (-5458.34205214835 - 980.14153344888 * y_pt +
        x * (2247.60742726704 - 340.1237483177863 * x + 220.542973797483 * y_pt))) +
        y_pt * (-1480.222530425046 +
        y_pt * (-129.1994027934126 +
        y_pt * (-30.0682112585625 + y_pt * (2.626801985426835 ))))) +
        y_pt * (1187.3715515697959 +
        y_pt * (1760.062705994408 + y_pt * (-450.535298526802 +
        y_pt * (182.8520895502518 + y_pt * (-43.3206481750622 + 4.26033941694366 * y_pt)))));

   g_SA_mod = 0.5 * gtc.gsw_sfac * g_SA_mod_part;

   g_p = gsw_gibbs(0, 0, 1, sa, t, p);

   g_psq_g_tt = g_p * g_p * gsw_gibbs(0, 2, 0, sa, t, p);

   g_tp = gsw_gibbs(0, 1, 1, sa, t, p);

   factora = g_SA_T_mod - g_SA_mod / (gtc.gsw_t0 + pt0);

   factor = factora / g_psq_g_tt;

   rho_SA = g_tp * factor - gsw_gibbs(1, 0, 1, sa, t, p) / (g_p * g_p);

   rho_CT = g_tp * gtc.gsw_cp0 / ((gtc.gsw_t0 + pt0) * g_psq_g_tt);

   rho_P = ((g_tp * g_tp) - gsw_gibbs(0, 2, 0, sa, t, p) *
               gsw_gibbs(0, 0, 2, sa, t, p)) / g_psq_g_tt;
}
/**
% gsw_Helmholtz_energy_t_exact                 Helmholtz energy of seawater
%==========================================================================
%
% USAGE: gsw_Helmholtz_energy_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the Helmholtz energy of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  Helmholtz_energy_t_exact  =  Helmholtz energy                   [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall                                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.13 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Helmholtz_energy_t_exact(double SA, double t, double p)
{
   return gsw_gibbs(0, 0, 0, SA, t, p) -
            (gtc.db2pa * p + gtc.gsw_p0) *
               gsw_gibbs(0, 0, 1, SA, t, p);
}
/**
% gsw_sigma1_CT_exact                        potential density anomaly with
%                                       reference sea pressure of 1000 dbar
%==========================================================================
%
% USAGE: gsw_sigma1_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 1000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma1(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma1_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 1000 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sigma1_CT_exact(double SA, double CT)
{
   double t = gsw_t_from_ct(SA, CT, 1000.0);

   return gsw_rho_t_exact(SA, t, 1000.0) - 1000.0;
}
/**
% gsw_specvol_anom_standard_CT_exact                specific volume anomaly
%==========================================================================
%
% USAGE: gsw_specvol_anom_standard_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific volume anomaly from Absolute Salinity, Conservative
%  Temperature and pressure.  The reference value of Absolute Salinity is
%  SSO and the reference value of Conservative Temperature is equal to
%  0 degress C.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_specvol_anom_standard(SA,CT,p), which uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT and
%  p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  specvol_anom_CT_exact  =  specific volume anomaly             [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_anom_standard_CT_exact(double SA, double CT, double p)
{
   double t = gsw_t_from_ct(SA,CT,p);

   return gsw_specvol_anom_standard_t_exact(SA, t, p);
}
/**
% gsw_N2sol_SP_pt                              solubility of N2 in seawater
%==========================================================================
%
% USAGE:
%  N2sol = gsw_N2sol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the nitrogen, N2, concentration expected at equilibrium with
%  air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
%  including saturated water vapor.  This function uses the solubility
%  coefficients as listed in Hamme and Emerson (2004).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  pt  =  potential temperature (ITS-90) referenced               [ deg C ]
%         to one standard atmosphere (0 dbar).
%
%  SP & pt need to have the same dimensions.
%
% OUTPUT:
%  N2sol = solubility of nitrogen in micro-moles per kg         [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Hamme, R., and S. Emerson, 2004: The solubility of neon, nitrogen and
%   argon in distilled water and seawater. Deep-Sea Research, 51,
%   1517-1528.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_N2sol_SP_pt(double SP, double pt)
{
	double y = log((298.15 - pt) / (pt + gtc.gsw_t0));

	return exp((y * (y * (4.69149 * y + 4.32531) + 2.92704) + 6.42931) +
				  SP * (y * (-0.0146775 * y + -0.00802566) + -0.00744129));
}
/**
% gsw_chem_potential_salt_t_exact                chemical potential of salt
%                                                               in seawater
%==========================================================================
%
% USAGE: gsw_chem_potential_salt_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the chemical potential of salt in seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  chem_potential_salt_t_exact  =  chemical potential of salt in seawater
%                                                                   [ J/g ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_chem_potential_salt_t_exact(double SA, double t, double p)
{
   return gsw_chem_potential_relative_t_exact(SA,t,p)
                               + gsw_chem_potential_water_t_exact(SA,t,p);
}
/**
% gsw_ntp_pt_vs_CT_ratio                    ratio of gradients of potential
%                             temperature and Conservative Temperature in a
%                            neutral tangent plane (in a locally-referenced
%                              potential density surface)(75-term equation)
% =========================================================================
%
% USAGE: gsw_ntp_pt_vs_CT_ratio(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the ratio of the two-dimensional gradient of potential
%  temperature versus that of Conservative Temperature, CT, along the
%  neutral tangent plane.  The potential temperature is the regular one
%  which has a reference sea pressure of 0 dbar.  Part of the calculation
%  uses the computationally-efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  ntp_pt_vs_CT_ratio  =  The ratio of the spatial gradient of
%                         potential temperature versus that of
%                         Conservative Temperature in the
%                         neutral tangent plane (ntp).         [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (A.14.5) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_ntp_pt_vs_CT_ratio(double SA, double CT, double p)
{
   double dummy, alpha, beta, pt_SA, pt_CT;

   gsw_specvol_alpha_beta(SA,CT,p, &dummy, &alpha, &beta);

/**
%--------------------------------------------------------------------------
% This function calculates the ntp_pt_vs_CT_ratio using the computationally
% efficient 75-term expression for specific volume in terms of SA, CT and
% p.  If one wanted to compute this with the full TEOS-10 Gibbs function
% expression for specific volume, the following lines of code will enable
% this.
%
%    t = gsw_t_from_CT(SA,CT,p);
%    beta = gsw_beta_const_CT_t_exact(SA,t,p);
%    alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p);
%
%--------- This is the end of the alternative code-------------------------
*/

   gsw_pt_first_derivatives(SA, CT, &pt_SA, &pt_CT);

   return pt_CT + pt_SA * (alpha / beta);
}
/**
% gsw_kappa_const_t_exact                        isothermal compressibility
%==========================================================================
%
% USAGE: gsw_kappa_const_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates isothermal compressibility of seawater.
%  Note. This is the compressibility of seawater AT CONSTANT IN-SITU
%    TEMPERATURE
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  kappa_const_t_exact  =  isothermal compressibility              [ 1/Pa ]
%   Note. The output units are 1/Pa not 1/dbar.
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%   See Eqn. (2.15.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_kappa_const_t_exact(double SA, double t, double p)
{
   return -gsw_gibbs(0, 0, 2, SA, t, p) / gsw_gibbs(0, 0, 1, SA, t, p);
}
/**
% gsw_specvol_from_pot_enthalpy_ice_poly               specific volume from
%                                          potential enthalpy of ice (poly)
%==========================================================================
%
% USAGE: gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_ice,p)
%
% DESCRIPTION:
%  Calculates the specific volume of ice from the potential enthalpy
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar.
%
%  Note that this is a comptationally efficient polynomial fit of the
%  programme gsw_specvol_from_pot_enthalpy_ice(pot_enthalpy_ice,p).
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice, double p)
{
	double v00 =  1.090843882189585e-3;
	double v01 = -2.993503247981132e-10;
	double v02 = -2.187697844100046e-15;
	double v03 = -6.597339965467078e-21;
	double v04 = -1.092035946035247e-26;
	double v05 = -9.164174131511919e-33;
	double v06 = -2.874852704079762e-39;
	double v07 = -1.558856898604863e-09;
	double v08 = -1.323880019968672e-15;
	double v09 = -2.685780994854041e-21;
	double v10 = -6.212785565052150e-27;
	double v11 = -4.797551830947803e-33;
	double v12 =  2.953039321948178e-15;
	double v13 = -2.403547886994993e-21;
	double v14 = -1.461069600352674e-20;

   return v00 + pot_enthalpy_ice
            * (v01 + pot_enthalpy_ice
            * (v02 + pot_enthalpy_ice
            * (v03 + pot_enthalpy_ice
            * (v04 + pot_enthalpy_ice
            * (v05 + v06 *pot_enthalpy_ice))))) + p
            * (v07 + pot_enthalpy_ice
            * (v08 + pot_enthalpy_ice
            * (v09 + pot_enthalpy_ice
            * (v10 + v11 *pot_enthalpy_ice))) + p
            * (v12 + v13 *pot_enthalpy_ice + v14 * p));
}
/**
% gsw_isopycnal_slope_ratio               ratio of the slopes of isopycnals
%                                      on the SA-CT diagram for p and p_ref
%                                                        (75-term equation)
% =========================================================================
%
% USAGE: gsw_isopycnal_slope_ratio(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates the ratio of alpha/beta at pressure, p, to that at reference
%  pressure, p_ref.  This function uses the computationally-efficient
%  75-term expression for specific volume in terms of SA, CT and p
%  (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  pr  =  reference pressure                                       [ dbar ]
%         ( i.e. absolute reference pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p and p_ref may have dimensions 1x1 or Mx1 or 1xN or MxN, where
%  SA and CT are MxN
%
% OUTPUT:
%  isopycnal_slope_ratio
%               =  The ratio of alpha/beta evaluated at        [ unitless ]
%                  pressure, p, to that at reference pressure, p_ref.
%
% AUTHOR:
%  Trevor McDougall, Paul Barker & David Jackett       [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (3.17.2) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_isopycnal_slope_ratio(double SA, double CT, double p, double p_ref)
{
   double *dummy1=new double, *alpha=new double, *beta=new double;
   double *dummy2=new double, *alpha_pref=new double, *beta_pref=new double;

   gsw_specvol_alpha_beta(SA, CT, p, dummy1, alpha, beta);
   gsw_specvol_alpha_beta(SA, CT, p_ref, dummy2, alpha_pref, beta_pref);

/**
%--------------------------------------------------------------------------
% This function calculates isopycnal_slope_ratio using the computationally
% efficient 75-term expression for specific volume as a function of SA, CT
% and p.  If one wanted to compute this with the full TEOS-10 Gibbs
% function expression for specific volume, the following lines of code will
% enable this.
%
%     t = gsw_pt_from_CT(SA,CT,p);
%     alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p);
%     beta = gsw_beta_const_CT_t_exact(SA,t,p);
%     tr = gsw_pt_from_t(SA,pt,p_ref0,p_ref);
%     alpha_pref = gsw_alpha_wrt_CT_t_exact(SA,tr,p_ref);
%     beta_pref = gsw_beta_const_CT_t_exact(SA,tr,p_ref);
%
%--------------This is the end of the alternative code---------------------
*/
   double v1 = *alpha, v2 = *beta_pref, v3 = *alpha_pref, v4 = *beta;
   double slope_ratio = v1 * v2 / v3 * v4;

   /** clean up */
   delete dummy1;
   delete alpha;
   delete beta;
   delete dummy2;
   delete alpha_pref;
   delete beta_pref;

   return slope_ratio;
}
/**
% gsw_chem_potential_relative_t_exact           relative chemical potential
%==========================================================================
%
% USAGE: gsw_chem_potential_relative_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the relative chemical potential of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  chem_potential_relative_t_exact  =  relative chemical potential  [ J/g ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_chem_potential_relative_t_exact(double SA, double t, double p)
{
   return gsw_gibbs(1, 0, 0, SA, t, p);
}
/**
% gsw_rho_CT_exact                                          in-situ density
%==========================================================================
%
% USAGE: gsw_rho_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates in-situ density from Absolute Salinity and Conservative
%  Temperature.
%
%  Note that potential density with respect to reference pressure, p_ref,
%  is obtained by calling this function with the pressure argument being
%  p_ref (i.e. "gsw_rho_CT_exact(SA,CT,p_ref)").
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_rho(SA,CT,p),
%  which uses the computationally efficient 75-term expression for density
%  in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_CT_exact  =  in-situ density                              [ kg/m^3 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.8.2) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_rho_CT_exact(double SA, double CT, double p)
{
	double t = gsw_t_from_ct(SA, CT, p);

	return gsw_rho_t_exact(SA, t, p);
}
/**
% gsw_sigma3_CT_exact                        potential density anomaly with
%                                       reference sea pressure of 3000 dbar
%==========================================================================
%
% USAGE: gsw_sigma3_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 3000
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma3(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma3_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 3000 dbar,
%                      that is, this potential density - 1000 kg m^-3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sigma3_CT_exact(double SA, double CT)
{
	double t = gsw_t_from_ct(SA, CT, 3000.0);

	return gsw_rho_t_exact(SA, t, 3000.0) - 1000.0;
}
/**
% gsw_rho_alpha_beta_CT_exact            in-situ density, thermal expansion
%                                  & saline contraction coefficient from CT
%==========================================================================
%
% USAGE: gsw_rho_alpha_beta_CT_exact(SA,CT,p,
                                       *rho_CT_exact,
                                       *alpha_CT_exact,
                                       *beta_CT_exact)
%
% DESCRIPTION:
%  Calculates in-situ density, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from
%  Absolute Salinity and Conservative Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_rho_alpha_beta(SA,CT,p), which uses the computationally-efficient
%  75-term expression for density in terms of SA, CT and p (Roquet et al.,
%  2015)
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_CT_exact    =  in-situ density                            [ kg/m^3 ]
%  alpha_CT_exact  =  thermal expansion coefficient                 [ 1/K ]
%                     with respect to Conservative Temperature
%  beta_CT_exact   =  saline contraction coefficient               [ kg/g ]
%                     at constant Conservative Temperature
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See sections (2.8), (2.18) and (2.19) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_rho_alpha_beta_CT_exact(double SA, double CT, double p,
		double &rho_CT_exact, double &alpha_CT_exact,
		double &beta_CT_exact)
{
	double t = gsw_t_from_ct(SA, CT, p);

	rho_CT_exact = gsw_rho_t_exact(SA, t, p);
	alpha_CT_exact = gsw_alpha_wrt_CT_t_exact(SA, t, p);
	beta_CT_exact = gsw_beta_const_CT_t_exact(SA, t, p);
}
/**
% gsw_rho_first_derivatives_wrt_enthalpy_CT_exact     second derivatives of
%                                  specific volume with respect to enthalpy
% =========================================================================
%
% USAGE: gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p, *rho_SA,
                                                         *rho_h)
%
% DESCRIPTION:
%  Calculates the following two first-order derivatives of rho with respect
%  to enthalpy h,
%   (1) rho_SA, first-order derivative with respect to Absolute Salinity
%       at constant h & p.
%   (2) rho_h, first-order derivative with respect to h at
%       constant SA and p.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_rho_first_derivatives_wrt_enthalpy(SA,CT,p) which uses the
%  computationally efficient 75-term expression for specific volume in
%  terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_SA  =  The second derivative of rho with respect to
%             Absolute Salinity at constant h & p.      [ J/(kg (g/kg)^2) ]
%  rho_h   =  The second derivative of rho with respect to
%             h at constant SA and p.                        [ J/(kg K^2) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(double SA,
		double CT,
		double p,
		double &rho_SA,
		double &rho_h)
{
	double rec_v, rec_v2, v_SA, v_h, sa = 0.0;

	if (SA > 0.0) sa = SA;

	rec_v = 1.0 / gsw_specvol_CT_exact(sa, CT, p);

	gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(sa, CT, p, v_SA, v_h);

	rec_v2 = rec_v * rec_v;

	rho_h = -v_h * rec_v2;

	rho_SA = -v_SA * rec_v2;
}
/**
% gsw_specvol_first_derivatives_CT_exact            first order derivatives
%                                                        of specific volume
% =========================================================================
%
% USAGE: gsw_specvol_first_derivatives_CT_exact(SA,CT,p, *v_SA,
                                                         *v_CT,
                                                         *v_P)
%
% DESCRIPTION:
%  Calculates the following three first-order derivatives of specific
%  volume (v),
%   (1) v_SA, first-order derivative with respect to Absolute Salinity
%       at constant CT & p.
%   (2) v_CT, first-order derivative with respect to SA & CT at
%       constant p.
%   (3) v_P, first-order derivative with respect to CT at constant SA
%       and p.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_specvol_first_derivatives(SA,CT,p) which uses the computationally
%  efficient 75 term expression for density in terms of SA, CT and p
%  (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  v_SA  =  The first derivative of specific volume with respect to
%           Absolute Salinity at constant CT & p.     [ (m^3/kg)(g/kg)^-1 ]
%  v_CT  =  The first derivative of specific volume with respect to
%           CT at constant SA and p.                         [ m^3/(K kg) ]
%  v_P   =  The first derivative of specific volume with respect to
%           P at constant SA and CT.                        [ m^3/(Pa kg) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_specvol_first_derivatives_CT_exact(double SA, double CT,
		double p,
		double &v_SA,
		double &v_CT,
		double &v_P)
{
	double rec_gTT, gSAP, gTP, gSAT,gSA_pt0, gPP;
	double t, rec_abs_pt0, pt0, pr0 = 0.0, sa = 0.0;

	if (SA > 0.0) sa = SA;

	pt0 = gsw_pt_from_ct(sa,CT);
	rec_abs_pt0 = 1.0 / (gtc.gsw_t0 + pt0);
	t = gsw_pt_from_t(sa,pt0,pr0,p);

	rec_gTT = 1.0 / gsw_gibbs(0,2,0, sa, t, p);

	gSAP = gsw_gibbs(1,0,1, sa, t, p);
	gTP = gsw_gibbs(0,1,1, sa, t, p);
	gSAT = gsw_gibbs(1,1,0, sa, t, p);
	gSA_pt0 = gsw_gibbs(1,0,0, sa, pt0, pr0);
	gPP = gsw_gibbs(0,0,2, sa, t, p);

	v_CT = -gtc.gsw_cp0 * gTP * rec_abs_pt0 * rec_gTT;

	v_SA = (gSAP - gTP) * (gSAT - rec_abs_pt0 * gSA_pt0) * rec_gTT;

	v_P = gPP - gTP * gTP *rec_gTT;
}
/**
% gsw_specvol_anom_CT_exact                         specific volume anomaly
%==========================================================================
%
% USAGE: gsw_specvol_anom_CT_exact(SA,CT,p,SA_ref,CT_ref)
%
% DESCRIPTION:
%  Calculates specific volume anomaly from Absolute Salinity, Conservative
%  Temperature and pressure.  If the Absolute Salinity and Conservative
%  Temperature reference value (SA_ref and CT_ref) are not defined then a
%  default reference value of Absolute Salinity is SSO and the reference
%  value of Conservative Temperature is equal to 0 degress C.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_specvol_anom(SA,CT,p,SA_ref,CT_ref), which uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT and
%  p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
% Optional:
%  SA_ref =  reference Absolute Salinity                           [ g/kg ]
%  CT_ref =  reference Conservative Temperature (ITS-90)          [ deg C ]
%
%  SA & CT  need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%  SA_ref & CT_ref muct be scalars and have dimensions 1x1.
%
% OUTPUT:
%  specvol_anom_CT_exact  =  specific volume anomaly             [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_anom_CT_exact(double SA, double CT, double p,
		double SA_ref, double CT_ref)
{
   double specvol_anom_CT_exact, t = gsw_t_from_ct(SA,CT,p);

   if (SA_ref == 0.0 && CT_ref == 0.0)
   {
      specvol_anom_CT_exact = gsw_specvol_anom_standard_t_exact(SA,t,p);
   }
   else
   {
      double t_ref_array, SA_ref_array = 1.0, CT_ref_array = 1.0;

      t_ref_array = gsw_t_from_ct(SA_ref_array, CT_ref_array, p);

      specvol_anom_CT_exact = gsw_gibbs(0,0,1, SA, t, p) -
                                 gsw_gibbs(0,0,1, SA_ref_array, t_ref_array, p);
   }

   return specvol_anom_CT_exact;
}
/**
% gsw_thermobaric_CT_exact                          thermobaric coefficient
%==========================================================================
%
% USAGE: gsw_thermobaric_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the thermobaric coefficient of seawater with respect to
%  Conservative Temperature.  This routine calculates the thermobaric
%  coefficient with the full TEOS-10 Gibbs function expression for density.
%  This function uses finite differences to calculate the temperature and
%  pressure derivatives.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_thermobaric(SA,CT,p)
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  thermobaric_CT_exact  =  thermobaric coefficient with       [ 1/(K Pa) ]
%                           respect to Conservative Temperature.
%  Note. The pressure derivative is taken with respect to
%    pressure in Pa not dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.8.2) and (P.2) of this TEOS-10 manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_thermobaric_CT_exact(double SA, double CT, double p)
{
	double v_SA_P, v_CT_P, v_SA, v_CT, dummy, dummy1, dummy2, dummy3, rho;
	double sa = 0.0;

	if (SA > 0.0) sa = SA;

	rho = gsw_rho_CT_exact(sa, CT, p);

	gsw_specvol_first_derivatives_CT_exact(sa, CT, p, v_SA, v_CT, dummy);

	gsw_specvol_second_derivatives_CT_exact(sa, CT, p, dummy1,
	                                          dummy2, dummy3, v_SA_P,
	                                          v_CT_P);

	return rho * (v_CT_P - (v_CT / v_SA) * v_SA_P);
}
/**
% gsw_t_from_pt0                         in-situ temperature from potential
%                                         temperature with a p_ref = 0 dbar
% =========================================================================
%
% USAGE: gsw_t_from_pt0(SA,pt0,p)
%
% DESCRIPTION:
%  Calculates in-situ temperature from potential temperature with a
%  reference pressure of 0 dbar.
%
%  It is also possible to calculate in-situ temperature from potential
%  temparature using the function gsw_pt_from_t.  In this case it would be
%  called with Absolute Salinity, SA, potential temperature referenced to
%  0 dbar, pt0, and the in-situ pressure (i.e. gsw_pt_from_t(SA,pt0,0,p) ).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt0 =  potential temperature with reference pressure = 0 dbar  [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & pt0 need to have the same dimensions. p may have dimensions 1x1 or
%  Mx1 or 1xN or MxN, where SA & pt0 are MxN.
%
% OUTPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.1 of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_t_from_pt0(double SA, double pt0, double p)
{
   return gsw_pt_from_t(SA, pt0, 0.0, p);
}
/*****************************************
   this is called from
   gsw_pot_enthalpy_from_specvol_ice_poly
   as a subfunction
*****************************************/
double TeosBase::gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(double pot_enthalpy_ice,
                                                         double p)
{
   double v01 = -2.993503247981132e-10;
   double v02 = -2.187697844100046e-15;
   double v03 = -6.597339965467078e-21;
   double v04 = -1.092035946035247e-26;
   double v05 = -9.164174131511919e-33;
   double v06 = -2.874852704079762e-39;
   double v08 = -1.323880019968672e-15;
   double v09 = -2.685780994854041e-21;
   double v10 = -6.212785565052150e-27;
   double v11 = -4.797551830947803e-33;
   double v13 = -2.403547886994993e-21;

   return v01 + pot_enthalpy_ice * (2.0 * v02 + pot_enthalpy_ice * (3.0 * v03
      + pot_enthalpy_ice * (4.0 * v04 + pot_enthalpy_ice * (5.0 * v05
         + 6.0 * v06 * pot_enthalpy_ice)))) + p * (v08 + pot_enthalpy_ice * (2.0 * v09
         + pot_enthalpy_ice * (3.0 * v10 + 4.0 * v11 * pot_enthalpy_ice)) + v13 * p);
}

/**
% gsw_pot_enthalpy_from_specvol_ice_poly            potential enthalpy from
%                                             specific volume of ice (poly)
%==========================================================================
%
% USAGE: gsw_pot_enthalpy_from_specvol_ice_poly(specvol_ice,p)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from the specific volume
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar.
%
%  Note that this is a comptationally efficient polynomial version of
%  specfic volume from potential enthalpy ice.
%
% INPUT:
%  specvol_ice =  specific volume                                [ m^3/kg ]
%  p           =  sea pressure                                     [ dbar ]
%              ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  pot_enthalpy_ice =  potential enthalpy of ice                   [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton|s method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_pot_enthalpy_from_specvol_ice_poly(double specvol_ice, double p)
{
	double h00 = -3.333548730778702e5;
	double h01 = -1.158533659785955e14;
	double h02 =  3.088110638980851e17;
	double h03 = -1.574635299523583e20;
	double h04 = -2.315779173662935e23;
	double h05 =  2.790387951606607e26;
	double h06 = -8.296646858006085e28;
	double h07 =  2.315408333461980e5;
	double h08 = -4.021450585146860e8;
	double h09 = -5.705978545141118e10;
	double h10 =  4.048693790458912e14;
	double h11 = -1.769028853661338e17;
	double h12 = -2.307016488735996e-2;
	double h13 =  2.135637957175861e1;
	double h14 =  8.593828063620491e-9;

	double pot_enthalpy_ice_old, delta_v_ice, pot_enthalpy_ice_mean;

   double pot_enthalpy_ice = h00 + specvol_ice *(h01 + specvol_ice *(h02 + specvol_ice *(h03
    + specvol_ice *(h04 + specvol_ice *(h05 + h06 *specvol_ice))))) + p *(h07
    + specvol_ice *(h08 + specvol_ice *(h09 + specvol_ice *(h10 + h11 *specvol_ice)))
    + p *(h12 + h13 *specvol_ice + h14 *p));

   double v_h_Ih = gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(pot_enthalpy_ice,p);
   // initial estimate of v_h_Ih, the pot_enthalpy_ice derivative of specvol_ice

   //--------------------------------------------------------------------------
   // Begin the modified Newton-Raphson iterative procedure
   //--------------------------------------------------------------------------

   pot_enthalpy_ice_old = pot_enthalpy_ice;
   delta_v_ice = gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old,p) - specvol_ice;
   pot_enthalpy_ice = pot_enthalpy_ice_old - delta_v_ice / v_h_Ih ;
   // this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)
   pot_enthalpy_ice_mean = 0.5*(pot_enthalpy_ice + pot_enthalpy_ice_old);
   v_h_Ih = gsw_specvol_from_pot_enthalpy_ice_first_derivatives_poly(pot_enthalpy_ice_mean,p);
   pot_enthalpy_ice = pot_enthalpy_ice_old - delta_v_ice / v_h_Ih;
   // this is the end of one loop.

   pot_enthalpy_ice_old = pot_enthalpy_ice;
   delta_v_ice = gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old,p) - specvol_ice;
   pot_enthalpy_ice = pot_enthalpy_ice_old - delta_v_ice /v_h_Ih ;
   // this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)

   // After 1.5 loops trhough the modified Newton-Raphson procedure the error
   // in pot_enthalpy_ice is 8 x 10^-9 J/kg which is machine precission for
   // this calculation.

   return pot_enthalpy_ice;
}

/**
% gsw_specvol_anom_standard_t_exact                 specific volume anomaly
%==========================================================================
% This function has changed name to gsw_specvol_anom_standard_t_exact
% Note that this is the last version of GSW (Version 3.05) that will
% contain the function "gsw_specvol_anom_t_exact".  Change your code to use
% "gsw_specvol_anom_standard_t_exact".
%==========================================================================
*/
double TeosBase::gsw_specvol_anom_t_exact(double SA, double t, double p)
{
   double t_zero = gsw_t_from_ct(gtc.gsw_sso,0.0, p);

   return gsw_gibbs(0,0,1,SA,t,p) -
            gsw_gibbs(0,0,1, gtc.gsw_sso, t_zero, p);

}
/**
% gsw_dynamic_enthalpy_CT_exact                 dyamic enthalpy of seawater
%==========================================================================
%
% USAGE: gsw_dynamic_enthalpy_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the dynamic enthalpy of seawater from Absolute Salinity and
%  Conservative Temperature and pressure.  Dynamic enthalpy is defined
%  as enthalpy minus potential enthalpy (Young, 2010).
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_dynamic_enthalpy(SA,CT,p), which uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT
%  and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  dynamic_enthalpy_CT_exact  =  dynamic enthalpy                  [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See apendix A.30 of this TEOS-10 Manual.
%
%  McDougall, T.J., 2003: Potential enthalpy: A conservative oceanic
%   variable for evaluating heat content and heat fluxes. Journal of
%   Physical Oceanography, 33, 945-963.
%    See Eqns. (18) and (22)
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  Young, W.R., 2010: Dynamic enthalpy, Conservative Temperature, and the
%   seawater Boussinesq approximation. Journal of Physical Oceanography,
%   40, 394-400.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_dynamic_enthalpy_CT_exact(double SA, double CT, double p)
{
   double t = gsw_t_from_ct(SA,CT,p);

   return gsw_enthalpy_t_exact(SA,t,p) - gtc.gsw_cp0 * CT;
}
/**
% gsw_cabbeling_CT_exact                              cabbeling coefficient
%==========================================================================
%
% USAGE: gsw_cabbeling_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the cabbeling coefficient of seawater with respect to
%  Conservative Temperature.  This routine calculates the cabbeling
%  coefficient with the full TEOS-10 Gibbs function expression for specific
%  volume.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_cabbeling(SA,CT,p)
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  cabbeling  =  cabbeling coefficient with respect to            [ 1/K^2 ]
%                Conservative Temperature.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.9.2) and (P.4) of this TEOS-10 manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_cabbeling_CT_exact(double SA, double CT, double p)
{
	double v_SA, v_CT, dummy1, rho, alpha_CT;
	double v_SA_SA, v_SA_CT, v_CT_CT, dummy2, dummy3;
	double alpha_SA, beta_SA, alpha_on_beta;

	gsw_specvol_first_derivatives_CT_exact(SA,CT,p, v_SA, v_CT, dummy1);

	gsw_specvol_second_derivatives_CT_exact(SA,CT,p,v_SA_SA, v_SA_CT, v_CT_CT, dummy2, dummy3);

	rho = gsw_rho_CT_exact(SA,CT,p);

	alpha_CT = rho * (v_CT_CT - rho * (v_CT * v_CT));

	alpha_SA = rho * (v_SA_CT - rho * v_SA * v_CT);

	beta_SA = -rho *(v_SA_SA - rho * (v_SA * v_SA));

	alpha_on_beta = gsw_alpha_on_beta_CT_exact(SA,CT,p);

	return alpha_CT + alpha_on_beta * (2.0 * alpha_SA - alpha_on_beta * beta_SA);
}
/**
% gsw_atomic_weight                 mole-weighted atomic weight of sea salt
%==========================================================================
%
% USAGE: gsw_atomic_weight
%
% DESCRIPTION:
%  This function returns the mole-weighted atomic weight of sea salt of
%  Reference Composition, which is 31.4038218 g/mol.  This has been
%  defined as part of the Reference-Composition Salinity Scale of 2008
%  (Millero et al., 2008).
%
% OUTPUT:
%  atomic_weight = mole-weighted atomic weight of sea salt of Reference
%                  Composition                                    [ g/mol ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Table D.4 of this TEOS-10 Manual.
%
%  Millero, F.J., R. Feistel, D.G. Wright, and T.J. McDougall, 2008:
%   The composition of Standard Seawater and the definition of the
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.
%     See Eqn. (5.3) of this paper.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_atomic_weight(void)
{
	return 31.4038218;
}
/**
% gsw_alpha_wrt_pt_t_exact                    thermal expansion coefficient
%                                     with respect to potential temperature
%==========================================================================
%
% USAGE: gsw_alpha_wrt_pt_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of seawater with respect to
%  potential temperature, with a reference pressure of zero.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  alpha_wrt_pt_t_exact  =  thermal expansion coefficient           [ 1/K ]
%                           with respect to potential temperature,
%                           with a reference pressure of zero dbar.
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.18.2) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_alpha_wrt_pt_t_exact(double SA, double t, double p)
{
   double pt0 = gsw_pt0_from_t(SA,t,p);

   double factor = gsw_gibbs(0,2,0,SA,pt0,0.0) / gsw_gibbs(0,2,0,SA,t,p);

   return factor * (gsw_gibbs(0,1,1,SA,t,p) / gsw_gibbs(0,0,1,SA,t,p));
}
/**
% gsw_CT_from_rho_exact               Conservative Temperature from density
% =========================================================================
%
% USAGE: gsw_CT_from_rho_exact(rho,SA,p, *CT, *CT_multiple)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar).
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_CT_from_rho(rho,SA,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it.
%     That is, it is |density|, not |density anomaly|.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & SA need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & SA are MxN.
%
% OUTPUT:
%  CT  =  Conservative Temperature                                [ deg C ]
%  CT_multiple  =  Conservative Temperature                       [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      temperatures for a single density.  This programme will output both
%      valid solutions.  To see this second solution the user must call the
%      programme with two outputs (i.e. [CT,CT_multiple]), if there is only
%      one possible solution and the programme has been called with two
%      outputs the second variable will be set to NaN.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_CT_from_rho_exact(double rho, double SA, double p,
													 double &CT,
													 double &CT_multiple)
{
	double t, t_multiple;

	gsw_t_from_rho_exact(rho, SA, p, t, t_multiple);

	CT = gsw_ct_from_t(SA, t, p);
	CT_multiple = gsw_ct_from_t(SA, t_multiple, p);
}
/**
% gsw_sigma0_pt0_exact                           potential density anomaly,
%                                 being potential density minus 1000 kg/m^3
%==========================================================================
%
% USAGE: gsw_sigma0_pt0_exact(SA,pt0)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference sea pressure of
%  zero (0) dbar.  The temperature input to this function is potential
%  temperature referenced to zero dbar.
%
% INPUT:
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  pt0  =  potential temperature with respect to a
%          reference sea pressure of 0 dbar (ITS-90)              [ deg C ]
%
%  SA & pt0 need to have the same dimensions.
%
% OUTPUT:
%  sigma0_pt0_exact  =  potential density anomaly with           [ kg/m^3 ]
%                       respect to a reference pressure of 0 dbar,
%                       that is, potential density minus 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.6.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sigma0_pt0_exact(double SA, double pt0)
{
	double sa = 0.0;

	if (SA > 0.0) sa = SA;

	return gsw_rho_t_exact(sa, pt0, 0.0) - 1000.0;
}
/**
% gsw_t_from_entropy                                    in-situ temperature
%                                                  as a function of entropy
% =========================================================================
%
% USAGE: gsw_t_from_entropy(SA,entropy,p)
%
% DESCRIPTION:
%  Calculates in-situ temperature with entropy as an input variable.
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  entropy  =  specific entropy                                [ J/(kg*K) ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & entropy need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & entropy are
%  MxN.
%
% OUTPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix  A.10 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_t_from_entropy(double SA, double entropy, double p)
{
	double p0 = 0.0, pt = gsw_pt_from_entropy(SA,entropy);

	/**
	   Note that pt is potential temperature
	   with a reference pressure of zero.
	*/

	return gsw_pt_from_t(SA, pt, p0, p);
}

/**
% gsw_internal_energy_first_derivatives_CT_exact       first derivatives of
%                                       specific interal energy of seawater
%==========================================================================
%
% USAGE: gsw_internal_energy_first_derivatives_CT_exact(SA,CT,p,
                                                         &u_SA,
                                                         &u_CT,
                                                         &u_P)
%
% DESCRIPTION:
%  Calculates the first order derivates of specific internal energy of
%  seawater.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_internal_energy_first_derivatives(SA,CT,p), which uses the
%  computationally efficient polynomial for specific volume in terms of
%  SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  u_SA = The first derivative of internal energy with respect to
%           Absolute Salinity at constant CT & p.
%                                          [ (J/kg)(g/kg)^-1 ] i.e. [ J/g ]
%  u_CT = The first derivative of internal energy with respect to
%           Conservative Temperature at constant SA & p.    [ (J/kg) K^-1 ]
%  u_P = The first derivative of internal energy with respect to
%           pressure at constant SA & CT.                  [ (J/kg) Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_internal_energy_first_derivatives_CT_exact(double SA,
		double CT,
		double p,
		double &u_SA,
		double &u_CT,
		double &u_P)
{
	double v_SA, v_CT, v_P, v, h_SA, h_CT, P, sa = 0.0;

	if (SA > 0.0) sa = SA;

	P = gtc.db2pa * p + gtc.gsw_p0;

	gsw_enthalpy_first_derivatives_ct_exact(sa, CT, p, &h_SA, &h_CT);
	v = gsw_specvol_CT_exact(sa,CT,p);
	gsw_specvol_first_derivatives_CT_exact(sa, CT, p, v_SA, v_CT, v_P);

	u_SA = h_SA - P * v_SA;

	u_CT = h_CT - P * v_CT;

	u_P = v - P * v_P;
}
/**
% gsw_alpha_wrt_CT_t_exact                    thermal expansion coefficient
%                                  with respect to Conservative Temperature
%==========================================================================
%
% USAGE: gsw_alpha_wrt_CT_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of seawater with respect to
%  Conservative Temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  alpha_wrt_CT_t_exact  =  thermal expansion coefficient           [ 1/K ]
%                           with respect to Conservative Temperature
%
% AUTHOR:
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.18.3) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_alpha_wrt_CT_t_exact(double SA, double t, double p)
{
	double pt0 = gsw_pt0_from_t(SA,t,p);
	double factor = -gtc.gsw_cp0 / ((gtc.gsw_t0 + pt0) * gsw_gibbs(0,2,0,SA,t,p));

	return factor * (gsw_gibbs(0,1,1,SA,t,p) / gsw_gibbs(0,0,1,SA,t,p));
}
/**
% gsw_internal_energy_second_derivatives              second derivatives of
%                                       specific interal energy of seawater
%                                                        (75-term equation)
%==========================================================================
%
% USAGE: gsw_internal_energy_second_derivatives(SA,CT,p,
                                                &u_SA_SA,
                                                &u_SA_CT,
                                                &u_CT_CT,
                                                &u_SA_P,
                                                &u_CT_P)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of
%  internal energy,
%  (1) u_SA_SA, second order derivative with respect to Absolute Salinity
%      at constant CT & p.
%  (2) u_SA_CT, second order derivative with respect to SA & CT at
%      constant p.
%  (3) u_CT_CT, second order derivative with respect to CT at constant
%      SA & p.
%  (4) u_SA_P, second-order derivative with respect to SA & P at
%      constant CT.
%  (5) u_CT_P, second-order derivative with respect to CT & P at
%      constant SA.
%
%  Note that this function uses the using the computationally-efficient
%  75-term expression for specific volume (Roquet et al., 2015).  There is
%  an alternative to calling this function, namely
%  gsw_internal_energy_second_derivatives_CT_exact(SA,CT,p) which uses the
%  full Gibbs function (IOC et al., 2010).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  u_SA_SA =  The second derivative of internal energy with respect to
%             Absolute Salinity at constant CT & p.     [ (J/kg)(g/kg)^-2 ]
%  u_SA_CT =  The second derivative of internal energy with respect to
%             SA & CT at constant p.                [ (J/kg)(g/kg)^-1 K^-1]
%  u_CT_CT =  The second derivative of internal energy with respect to
%             CT at constant SA and p.                      [ (J/kg) K^-2 ]
%  u_SA_P  =  The second derivative of internal energy with respect to
%             SA & P at constant CT.              [ (J/kg)(g/kg)^-1 Pa^-1 ]
%  u_CT_P  =  The second derivative of internal energy with respect to
%             CT & P at constant SA.                  [ (J/kg) K^-1 Pa^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_internal_energy_second_derivatives(double SA, double CT, double p,
		double &u_SA_SA, double &u_SA_CT,
		double &u_CT_CT, double &u_SA_P,
		double &u_CT_P)
{
	double *h_SA_SA=new double, *h_SA_CT=new double, *h_CT_CT=new double;
	double *v_SA_SA=new double, *v_SA_CT=new double, *v_CT_CT=new double;
	double *v_SA_P=new double, *v_CT_P=new double;
	double P, sa = 0.0;


	if (SA > 0.0) sa = SA;

	P = (gtc.db2pa * p + gtc.gsw_p0);

	gsw_enthalpy_second_derivatives(sa,CT,p, h_SA_SA, h_SA_CT, h_CT_CT);
	gsw_specvol_second_derivatives(sa,CT,p,v_SA_SA, v_SA_CT, v_CT_CT, v_SA_P, v_CT_P);

	u_SA_SA = *h_SA_SA - P * *v_SA_SA;

	u_SA_CT = *h_SA_CT - P * *v_SA_CT;

	u_CT_CT = *h_CT_CT - P * *v_CT_CT;

	u_SA_P = -P * *v_SA_P;

	u_CT_P = -P * *v_CT_P;

	/** clean up */
	delete h_SA_SA;
	delete h_SA_CT;
	delete h_CT_CT;
	delete v_SA_SA;
	delete v_SA_CT;
	delete v_CT_CT;
	delete v_SA_P;
	delete v_CT_P;
}
/**
   The function "entropy_t_exact"
   has changed to "gsw_entropy_from_t"
*/
double TeosBase::gsw_entropy_t_exact(double SA, double t, double p)
{
	return gsw_entropy_from_t(SA, t, p);
}


/**
% gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact     first derivatives
%                               of specific volume with respect to enthalpy
% =========================================================================
%
% USAGE: gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(SA,CT,p,
                                                               *v_SA,
                                                               *v_h)
%
% DESCRIPTION:
%  Calculates the following two first-order derivatives of specific
%  volume (v),
%   (1) v_SA, first-order derivative with respect to Absolute Salinity
%       at constant h & p.
%   (2) v_h, first-order derivative with respect to h at constant SA & p.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_specvol_first_derivatives_wrt_enthalpy(SA,CT,p) which uses the
%  computationally efficient 75 term expression for specific volume in
%  terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  v_SA =  The first derivative of specific volume with respect to
%          Absolute Salinity at constant h & p.         [ J/(kg (g/kg)^2) ]
%  v_h  =  The first derivative of specific volume with respect to
%          h at constant SA & p.                         [ J/(kg K(g/kg)) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(double SA,
		double CT,
		double p,
		double &v_SA,
		double &v_h)
{
	double rec_h_CT, vCT_SA, vCT_CT, dummy, h_SA, h_CT, sa = 0.0;

	 if (SA > 0.0) sa = SA;

	gsw_specvol_first_derivatives_CT_exact(sa, CT, p, vCT_SA, vCT_CT, dummy);
	gsw_enthalpy_first_derivatives_ct_exact(sa, CT, p, &h_SA, &h_CT);

	rec_h_CT = 1.0 / h_CT;

	v_SA = vCT_SA - (vCT_CT * h_SA) * rec_h_CT;

	v_h = vCT_CT * rec_h_CT;
}
/**
% gsw_SR_from_SP                 Reference Salinity from Practical Salinity
%==========================================================================
%
% USAGE: gsw_SR_from_SP(SP)
%
% DESCRIPTION:
%  Calculates Reference Salinity from Practical Salinity.
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%
% OUTPUT:
%  SR  =  Reference Salinity                                       [ g/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sr_from_sp(double sp)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	return (sp*gtc.gsw_ups);
}
/**
% gsw_p_from_Abs_Pressure                                      sea pressure
%==========================================================================
%
% USAGE: gsw_p_from_Abs_Pressure(Absolute_Pressure)
%
% DESCRIPTION:
%  Calculates sea pressure from Absolute Pressure.  Note that Absolute
%  Pressure is in Pa NOT dbar.
%
% INPUT:
%  Absolute_Pressure  =  Absolute Pressure                           [ Pa ]
%
% OUTPUT:
%  p  =  sea pressure                                              [ dbar ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.2.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_p_from_Abs_Pressure(double Absolute_Pressure)
{
	return (Absolute_Pressure - gtc.gsw_p0) * 0.0001;
}
/**
% gsw_SP_from_R                  Practical Salinity from conductivity ratio
%==========================================================================
%
% USAGE: gsw_SP_from_R(R,t,p)
%
% DESCRIPTION:
%  Calculates Practical Salinity, SP, from the conductivity ratio, R,
%  primarily using the PSS-78 algorithm.  Note that the PSS-78 algorithm
%  for Practical Salinity is only valid in the range 2 < SP < 42.  If the
%  PSS-78 algorithm produces a Practical Salinity that is less than 2 then
%  the Practical Salinity is recalculated with a modified form of the
%  Hill et al. (1986) formula.  The modification of the Hill et al. (1986)
%  expression are to ensure that it is exactly consistent with PSS-78
%  at SP = 2.
%
% INPUT:
%  R  =  conductivity ratio                                    [ unitless ]
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar )
%
%  t & p may have dimensions 1x1 or Mx1 or 1xN or MxN, where R is MxN.
%
% OUTPUT:
%  SP  =   Practical Salinity on the PSS-78 scale              [ unitless ]
%
% AUTHOR:
%  Paul Barker, Trevor McDougall and Rich Pawlowicz    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the
%   Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng.,
%   OE-11, 1, 109 - 112.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix E of this TEOS-10 Manual.
%
%  Unesco, 1983: Algorithms for computation of fundamental properties of
%   seawater. Unesco Technical Papers in Marine Science, 44, 53 pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_SP_from_R(double R, double t, double p)
{
	double SP, Rtx, Rt, Rp, rt_lc, ft68, t68 = t * 1.00024;

	ft68 = (t68 - 15.0) / (1.0 + gspc.k * (t68 - 15.0));

	//rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.
	rt_lc = gspc.c0 + (gspc.c1 + (gspc.c2
	         + (gspc.c3 + gspc.c4 * t68) * t68) * t68) * t68;

	Rp = 1.0 + (p * (gspc.e1 + p * (gspc.e2 + gspc.e3 * p))) / (1.0
		  + gspc.d1 * t68
		  + gspc.d2 * t68 * t68 + (gspc.d3
		  + gspc.d4 * t68) * R);

	Rt = R / (Rp * rt_lc);

	if (Rt < 0) Rt = NAN;

	Rtx = sqrt(Rt);

	SP = gspc.a0 + (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4
                  + gspc.a5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx
                  + ft68 * (gspc.b0 + (gspc.b1
                  + (gspc.b2 + (gspc.b3
                  + (gspc.b4 + gspc.b5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx);

	/**
	   The following section of the code is designed for SP < 2 based on the
	   Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
	   exactly equal to the PSS-78 algorithm at SP = 2.
	*/

	if (SP < 2.0)
	{
		double Hill_ratio,x,sqrty,part1,part2,SP_Hill_raw;

		Hill_ratio = gsw_hill_ratio_at_sp2(t);
		x = 400.0 * Rt;
		sqrty = 10.0 * Rtx;
		part1 = 1.0 + x * (1.5 + x);
		part2 = 1 + sqrty * (1.0 + sqrty * (1.0 + sqrty));
		SP_Hill_raw = SP - gspc.a0 / part1 - gspc.b0 * ft68 / part2;
		SP = Hill_ratio * SP_Hill_raw;
	}

	//This line ensures that SP is non-negative.
	if (SP < 0.0) SP = 0.0;

	return SP;
}
/**
% gsw_isochoric_heat_cap_t_exact        isochoric heat capacity of seawater
% =========================================================================
%
% USAGE: gsw_isochoric_heat_cap_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the isochoric heat capacity of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  isochoric_heat_cap_t_exact  =  isochoric heat capacity      [ J/(kg K) ]
%
% AUTHOR:
%  Trevor McDougall                                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.21 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_isochoric_heat_cap_t_exact(double SA, double t, double p)
{
	double g_tt, g_tp, g_pp;
	g_tt = gsw_gibbs(0,2,0,SA,t,p);
	g_tp = gsw_gibbs(0,1,1,SA,t,p);
	g_pp = gsw_gibbs(0,0,2,SA,t,p);

	return (-(t + gtc.gsw_t0) * (g_tt - g_tp * g_tp / g_pp));
}

/**
% gsw_beta_const_pt_t_exact                  saline contraction coefficient
%                                         at constant potential temperature
%==========================================================================
%
% USAGE: gsw_beta_const_pt_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the saline (i.e. haline) contraction coefficient of seawater
%  at constant potential temperature with a reference pressure of 0 dbar.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  beta_const_pt_t_exact  =  saline contraction coefficient        [ kg/g ]
%                            at constant potential temperature
%                            and with a reference pressure of 0 dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.19.2) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_beta_const_pt_t_exact(double SA, double t, double p)
{
	double factor,gp,g_SA_pt_mod,g_SA_T_mod,pt0,x2,x,y,y_pt,z;

	pt0 = gsw_pt0_from_t(SA,t,p);
	x2 = gtc.gsw_sfac * SA;
	x = sqrt(x2);
	y = 0.025 * t;
	y_pt = 0.025 * pt0;
	z = gtc.rec_db2pa * p; //Note.The input pressure (p) is sea pressure in units of dbar.

	g_SA_T_mod = 1187.3715515697959 + z *(1458.233059470092 +
      z *(-687.913805923122 + z *(249.375342232496 + z *(-63.313928772146 + 14.09317606630898 *z)))) +
      x *(-1480.222530425046 + x *(2175.341332000392 + x *(-980.14153344888 + 220.542973797483 *x) +
      y *(-548.4580073635929 + y *(592.4012338275047 + y *(-274.2361238716608 + 49.9394019139016 *y))) -
         90.6734234051316 *z) + z *(-525.876123559641 + (249.57717834054571 - 88.449193048287 *z) *z) +
      y *(-258.3988055868252 + z *(2298.348396014856 + z *(-325.1503575102672 + 153.8390924339484 *z)) +
      y *(-90.2046337756875 - 4142.8793862113125 *z + y *(10.50720794170734 + 2814.78225133626 *z)))) +
      y *(3520.125411988816 + y *(-1351.605895580406 +
      y *(731.4083582010072 + y *(-216.60324087531103 + 25.56203650166196 *y) +
      z *(-2381.829935897496 + (597.809129110048 - 291.8983352012704 *z) *z)) +
      z *(4165.4688847996085 + z *(-1229.337851789418 + (681.370187043564 - 66.7696405958478 *z) *z))) +
      z *(-3443.057215135908 + z *(1349.638121077468 +
      z *(-713.258224830552 + (176.8161433232 - 31.68006188846728 *z) *z))));

	g_SA_T_mod = 0.5 *gtc.gsw_sfac * 0.025 * g_SA_T_mod;

	g_SA_pt_mod = 1187.3715515697959 +
      x *(-1480.222530425046 + x *(2175.341332000392 + x *(-980.14153344888 + 220.542973797483 *x) +
      y_pt *(-548.4580073635929 + y_pt *(592.4012338275047 +
      y_pt *(-274.2361238716608 + 49.9394019139016 *y_pt)))) +
      y_pt *(-258.3988055868252 +
      y_pt *(-90.2046337756875 + y_pt *(10.50720794170734)))) +
      y_pt *(3520.125411988816 + y_pt *(-1351.605895580406 +
      y_pt *(731.4083582010072 + y_pt *(-216.60324087531103 + 25.56203650166196 *y_pt))));

	g_SA_pt_mod = 0.5 * gtc.gsw_sfac * 0.025 * g_SA_pt_mod;

	gp = gsw_gibbs(0,0,1,SA,t,p);

	factor = (g_SA_T_mod - g_SA_pt_mod) / (gp * gsw_gibbs(0,2,0,SA,t,p));

	return (gsw_gibbs(0,1,1,SA,t,p) * factor) - (gsw_gibbs(1,0,1,SA,t,p) / gp);
}

/***************************************
   returns true or false depending on
   equality or not of rtIsNaNF with val
***************************************/
bool TeosBase::gsw_rtIsNaNF(double val)
{
	return (val == rtNanF);
};

/**************************************************
   The function "gsw_adiabatic_lapse_rate_t_exact"
   has changed to "gsw_adiabatic_lapse_rate_from_t"
****************************************************/
double TeosBase::gsw_adiabatic_lapse_rate_t_exact(double SA, double t, double p)
{
	return gsw_adiabatic_lapse_rate_from_t(SA,t,p);
}
/**
% gsw_adiabatic_lapse_rate_from_t                      adiabatic lapse rate
%==========================================================================
%
% USAGE: gsw_adiabatic_lapse_rate_from_t(SA,t,p)
%
% DESCRIPTION:
%  Calculates the adiabatic lapse rate of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  adiabatic_lapse_rate  =  adiabatic lapse rate                   [ K/Pa ]
%    Note.  The output is in unit of degress Celsius per Pa,
%      (or equivilently K/Pa) not in units of K/dbar.
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.22.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_adiabatic_lapse_rate_from_t(double SA, double t, double p)
{
   return -gsw_gibbs(0,1,1,SA,t,p) / gsw_gibbs(0,2,0,SA,t,p);
}
/**
% gsw_osmotic_coefficient_t_exact                       osmotic coefficient
%==========================================================================
%
% USAGE: gsw_osmotic_coefficient_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the osmotic coefficient of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  osmotic_coefficient_t_exact  =  osmotic coefficient of seawater
%                                                              [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_osmotic_coefficient_t_exact(double SA, double t, double p)
{
   double sa = 0.0;
   double k = 3.777007343340624e-3; // k = M_S/R
   double part, x2, x, y, z, tl;
	double  oc01 =  7.231916621570606e1;
	double  oc02 =  1.059039593127674e1;
	double  oc03 = -3.025914794694813e1;
	double  oc04 =  5.040733670521486e1;
	double  oc05 = -4.074543321119333e1;
	double  oc06 =  1.864215613820487e1;
	double  oc07 = -3.022566485046178;
	double  oc08 = -6.138647522851840;
	double  oc09 =  1.353207379758663e1;
	double  oc10 = -7.316560781114737;
	double  oc11 =  1.829232499785750;
	double  oc12 = -5.358042980767074e-1;
	double  oc13 = -1.705887283375562;
	double  oc14 = -1.246962174707332e-1;
	double  oc15 =  1.228376913546017;
	double  oc16 =  1.089364009088042e-2;
	double  oc17 = -4.264828939262248e-1;
	double  oc18 =  6.213127679460041e-2;
	double  oc19 =  2.481543497315280;
	double  oc20 = -1.363368964861909;
	double  oc21 = -5.640491627443773e-1;
	double  oc22 =  1.344724779893754;
	double  oc23 = -2.180866793244492;
	double  oc24 =  4.765753255963401;
	double  oc25 = -5.726993916772165;
	double  oc26 =  2.918303792060746;
	double  oc27 = -6.506082399183509e-1;
	double  oc28 = -1.015695507663942e-1;
	double  oc29 =  1.035024326471108;
	double  oc30 = -6.742173543702397e-1;
	double  oc31 =  8.465642650849419e-1;
	double  oc32 = -7.508472135244717e-1;
	double  oc33 = -3.668086444057845e-1;
	double  oc34 =  3.189939162107803e-1;
	double  oc35 = -4.245629194309487e-2;

   if (SA > 0.0) sa = SA;

/******************************************************
   M_S = 0.0314038218;

   mole-weighted average atomic weight of the elements
   of Reference-Composition sea salt, in units of
   kg mol^-1.  Strictly speaking, the formula below
   applies only to seawater of Reference Composition.

   If molality and the osmotic coefficienit is
   required to an accuracy of better than 0.1% we
   suggest you contact the authors for further
   guidance.
*******************************************************/

   part = k *(1000.00 - sa) / (gtc.gsw_t0 + t);

   x2 = gtc.gsw_sfac * sa;
   x = sqrt(x2);
   y = t * 0.025;
   z = p * gtc.rec_db2pa; //Note that the input pressure (p) is sea pressure in units of dbar.

   tl = oc01 + oc02 * y
    + x *(oc03 + x *(oc04 + x *(oc05 + x *(oc06 + oc07*x)))
    + y *(oc08 + x *(oc09 + x *(oc10 + oc11*x))
    + y *(oc12 + oc13*x + y *(oc14 + oc15*x + y *(oc16 + x *(oc17 + oc18*y)))))
    + z *(oc19 + x *(oc20 + oc21*y + oc22*x) + y *(oc23 + y *(oc24 + y *(oc25 + oc26*y)))
    + z *(oc27 + oc28*x + y *(oc29 + oc30*y)
    + z *(oc31 + oc32*x + y *(oc33 + oc34*y) + oc35*z))));

   return tl * part;
}

/**
% gsw_Arsol_SP_pt                              solubility of Ar in seawater
%==========================================================================
%
% USAGE: gsw_Arsol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the argon, Ar, concentration expected at equilibrium with air
%  at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
%  saturated water vapor  This function uses the solubility coefficients
%  as listed in Hamme and Emerson (2004).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  pt  =  potential temperature (ITS-90) referenced               [ deg C ]
%         to one standard atmosphere (0 dbar).
%
%  SP & pt need to have the same dimensions.
%
% OUTPUT:
%  Arsol = solubility of argon                                  [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Hamme, R., and S. Emerson, 2004: The solubility of neon, nitrogen and
%   argon in distilled water and seawater. Deep-Sea Research, 51,
%   1517-1528.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Arsol_SP_pt(double SP, double pt)
{
	double y = log((298.15 - pt) / (pt + gtc.gsw_t0));

	return exp((y * (y * (4.90379 * y + 4.13116) + 3.17609) + 2.7915) +
				  SP * (y * (-0.0116888 * y + -0.0076667) + -0.00696233));
}

/**
% gsw_depth_from_z                                     depth from height, z
%==========================================================================
%
% USAGE: gsw_depth_from_z(z)
%
% DESCRIPTION:
%  Calculates depth from height, z.  Note that in general height is
%  negative in the ocean.
%
% INPUT:
%  z  =  height                                                       [ m ]
%
% OUTPUT:
%  depth  =  depth                                                    [ m ]
%
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_depth_from_z(double z)
{
	return -z;
}
/**
% gsw_molality_from_SA                                 molality of seawater
%==========================================================================
%
% USAGE: gsw_molality_from_SA(SA)
%
% DESCRIPTION:
%  Calculates the molality of seawater from Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%
% OUTPUT:
%  molality  =  molality of seawater                             [ mol/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_molality_from_SA(double SA)
{
	double unnamed_idx_0 = 0.0;

	if (SA > 0.0)
	{
		unnamed_idx_0 = SA;
	}

	return unnamed_idx_0 / (0.0314038218 * (1000.0 - unnamed_idx_0));
}
/**
% gsw_Abs_Pressure_from_p                                 Absolute Pressure
%==========================================================================
%
% USAGE: gsw_Abs_Pressure_from_p(p)
%
% DESCRIPTION:
%  Calculates Absolute Pressure from sea pressure.  Note that Absolute
%  Pressure is in Pa NOT dbar.
%
% INPUT:
%  p  =  sea pressure                                              [ dbar ]
%
% OUTPUT:
%  Absolute_Pressure  =  Absolute Pressure                           [ Pa ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.2.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Abs_Pressure_from_p(double p)
{
	return p * gtc.db2pa;
}
/**
% gsw_CT_maxdensity_exact               Conservative Temperature of maximum
%                                                       density of seawater
% =========================================================================
%
% USAGE: gsw_CT_maxdensity_exact(SA,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of maximum density of seawater.
%  This function returns the Conservative temperature at which the density
%  of seawater is a maximum, at given Absolute Salinity, SA, and sea
%  pressure, p (in dbar).
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_CT_maxdensity(SA,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA =  Absolute Salinity                                         [ g/kg ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar )
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  CT_maxdensity_exact  =  Conservative Temperature at which      [ deg C ]
%                          the density of seawater is a maximum for
%                          given Absolute Salinity and pressure.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.42 of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_CT_maxdensity_exact(double SA, double p)
{
	double t_maxdensity_exact = gsw_t_maxdensity_exact(SA,p);

	return gsw_ct_from_t(SA, t_maxdensity_exact, p);
}
/**
% gsw_internal_energy_t_exact          specific internal energy of seawater
%==========================================================================
%
% USAGE: gsw_internal_energy_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the specific internal energy of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  internal_energy_t_exact  =  specific internal energy (u)        [ J/kg ]
%
% AUTHOR:
%  Trevor McDougall                                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.11.1) of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_internal_energy_t_exact(double SA, double t, double p)
{
   return gsw_gibbs(0,0,0,SA,t,p) -
            (t + gtc.gsw_t0) * gsw_gibbs(0,1,0,SA,t,p) -
               (gtc.db2pa * p + gtc.gsw_p0 ) *
                  gsw_gibbs(0,0,1,SA,t,p);
}
/**
% gsw_sound_speed_CT_exact                                      sound speed
%==========================================================================
%
% USAGE: gsw_sound_speed_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the speed of sound in seawater from Absolute Salinity and
%  Conservative Temperature and pressure.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sound_speed(SA,CT,p)
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  sound_speed_CT_exact  =  speed of sound in seawater              [ m/s ]
%
% AUTHOR:
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.17.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sound_speed_CT_exact(double SA, double CT, double p)
{
	double t = gsw_t_from_ct(SA,CT,p);

	return gsw_sound_speed_t_exact(SA,t,p);
}
/**
% gsw_Hesol_SP_pt                              solubility of He in seawater
%==========================================================================
%
% USAGE: gsw_Hesol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the helium concentration expected at equilibrium with air at
%  an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
%  saturated water vapor.  This function uses the solubility coefficients
%  as listed in Weiss (1971).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  pt  =  potential temperature (ITS-90) referenced               [ deg C ]
%         to one standard atmosphere (0 dbar).
%
%  SP & pt need to have the same dimensions.
%
% OUTPUT:
%  Hesol = solubility of helium in micro-moles per kg           [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Dymond and Smith, 1980: The virial coefficients of pure gases and
%   mixtures. Clarendon Press, Oxford.
%
%  Weiss, R.F., 1971: Solubility of Helium and Neon in Water and Seawater.
%   J. Chem. and Engineer. Data, 16, 235-241.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Hesol_SP_pt(double SP, double pt)
{
	double pt68, y_100;

	pt68 = pt * 1.00024;
	y_100 = (pt68 + gtc.gsw_t0) * 0.01;

	return exp((((21634.42 / (pt68 + gtc.gsw_t0) + -167.2178) +
					 139.2032 * log(y_100)) + -22.6202 * y_100) +
				      SP * (y_100 * (-0.0034266 * y_100 + 0.023541)
				         + -0.044781)) * 44.558176715055367;
}
/**
% gsw_pt_second_derivatives     second derivatives of potential temperature
% =========================================================================
%
% USAGE: gsw_pt_second_derivatives(SA,CT, *pt_SA_SA, *pt_SA_CT, *pt_CT_CT)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of potential
%  temperature (the regular potential temperature which has a reference
%  sea pressure of 0 dbar),
%   (1) pt_SA_SA, the second derivative with respect to Absolute Salinity
%       at constant Conservative Temperature,
%   (2) pt_SA_CT, the derivative with respect to Conservative Temperature
%       and Absolute Salinity, and
%   (3) pt_CT_CT, the second derivative with respect to Conservative
%       Temperature at constant Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  pt_SA_SA  =  The second derivative of potential temperature (the
%               regular potential temperature which has reference sea
%               pressure of 0 dbar) with respect to Absolute Salinity
%               at constant Conservative Temperature.
%               pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
%  pt_SA_CT  =  The derivative of potential temperature with respect
%               to Absolute Salinity and Conservative Temperature.
%               pt_SA_CT has units of:                         [ 1/(g/kg) ]
%  pt_CT_CT  =  The second derivative of potential temperature (the
%               regular one with p_ref = 0 dbar) with respect to
%               Conservative Temperature at constant SA.
%               pt_CT_CT has units of:                              [ 1/K ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.9) and (A.12.10) of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_pt_second_derivatives(double SA, double CT, double &pt_SA_SA,
		double &pt_SA_CT, double &pt_CT_CT)
{
	double *pt_CT_l = new double;
	double *pt_CT_u = new double;
	double *pt_SA_l = new double;
	double *pt_SA_u = new double;
	double unnamed_idx_0;
	unnamed_idx_0 = rtNaN;

	if (SA >= 0.001)
	{
		unnamed_idx_0 = SA - 0.001;
	}

	if (SA < 0.001)
	{
		unnamed_idx_0 = 0.0;
	}

	gsw_pt_first_derivatives(unnamed_idx_0, CT, pt_SA_l, pt_CT_l);
	gsw_pt_first_derivatives(SA + 0.001, CT, pt_SA_u, pt_CT_u);
	pt_SA_SA = (pt_SA_u - pt_SA_l) / ((SA + 0.001) - unnamed_idx_0);
	gsw_pt_first_derivatives(SA, CT - 0.01, pt_SA_l, pt_CT_l);
	gsw_pt_first_derivatives(SA, CT + 0.01, pt_SA_u, pt_CT_u);
	unnamed_idx_0 = (CT + 0.01) - (CT - 0.01);
	pt_SA_CT = (pt_SA_u - pt_SA_l) / unnamed_idx_0;
	pt_CT_CT = (pt_CT_u - pt_CT_l) / unnamed_idx_0;

	/** clean up */
	delete pt_CT_l;
	delete pt_CT_u;
	delete pt_SA_l;
	delete pt_SA_u;
}
/**
% gsw_rho_first_derivatives                SA, CT and p partial derivatives
%                                             of density (75-term equation)
%==========================================================================
%
% USAGE: gsw_rho_first_derivatives(SA,CT,p, *rho_SA,  *rho_CT, *rho_P)
%
% DESCRIPTION:
%  Calculates the three (3) partial derivatives of in-situ density with
%  respect to Absolute Salinity, Conservative Temperature and pressure.
%  Note that the pressure derivative is done with respect to pressure in
%  Pa, not dbar.  This function uses the computationally-efficient
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_SA  =  partial derivative of density           [ (kg/m^3)(g/kg)^-1 ]
%                 with respect to Absolute Salinity
%  rho_CT  =  partial derivative of density                  [ kg/(m^3 K) ]
%                 with respect to Conservative Temperature
%  rho_P   =  partial derivative of density                 [ kg/(m^3 Pa) ]
%                 with respect to pressure in Pa
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_rho_first_derivatives(double SA, double CT, double p,
		double &rho_SA, double &rho_CT,
		double &rho_P)
{
   /**
      the following 'v' values differ from their
      counterparts in GSW_SPECVOL_COEFFICIENTS
   */
	double v010 = -1.5649734675e-5;
	double v011 =  1.8505765429e-5;
	double v012 = -1.1736386731e-6;
	double v013 = -3.6527006553e-7;
	double v014 =  3.1454099902e-7;
	double v020 =  2.7762106484e-5;
	double v021 = -1.1716606853e-5;
	double v022 =  2.1305028740e-6;
	double v023 =  2.8695905159e-7;
	double v030 = -1.6521159259e-5;
	double v031 =  7.9279656173e-6;
	double v032 = -4.6132540037e-7;
	double v040 =  6.9111322702e-6;
	double v041 = -3.4102187482e-6;
	double v042 = -6.3352916514e-8;
	double v050 = -8.0539615540e-7;
	double v051 =  5.0736766814e-7;
	double v060 =  2.0543094268e-7;
	double v100 = -3.1038981976e-4;
	double v101 =  2.4262468747e-5;
	double v102 = -5.8484432984e-7;
	double v103 =  3.6310188515e-7;
	double v104 = -1.1147125423e-7;
	double v110 =  3.5009599764e-5;
	double v111 = -9.5677088156e-6;
	double v112 = -5.5699154557e-6;
	double v113 = -2.7295696237e-7;
	double v120 = -3.7435842344e-5;
	double v121 = -2.3678308361e-7;
	double v122 =  3.9137387080e-7;
	double v130 =  2.4141479483e-5;
	double v131 = -3.4558773655e-6;
	double v132 =  7.7618888092e-9;
	double v140 = -8.7595873154e-6;
	double v141 =  1.2956717783e-6;
	double v150 = -3.3052758900e-7;
	double v200 =  6.6928067038e-4;
	double v201 = -3.4792460974e-5;
	double v202 = -4.8122251597e-6;
	double v203 =  1.6746303780e-8;
	double v210 = -4.3592678561e-5;
	double v211 =  1.1100834765e-5;
	double v212 =  5.4620748834e-6;
	double v220 =  3.5907822760e-5;
	double v221 =  2.9283346295e-6;
	double v222 = -6.5731104067e-7;
	double v230 = -1.4353633048e-5;
	double v231 =  3.1655306078e-7;
	double v240 =  4.3703680598e-6;
	double v300 = -8.5047933937e-4;
	double v301 =  3.7470777305e-5;
	double v302 =  4.9263106998e-6;
	double v310 =  3.4532461828e-5;
	double v311 = -9.8447117844e-6;
	double v312 = -1.3544185627e-6;
	double v320 = -1.8698584187e-5;
	double v321 = -4.8826139200e-7;
	double v330 =  2.2863324556e-6;
	double v400 =  5.8086069943e-4;
	double v401 = -1.7322218612e-5;
	double v402 = -1.7811974727e-6;
	double v410 = -1.1959409788e-5;
	double v411 =  2.5909225260e-6;
	double v420 =  3.8595339244e-6;
	double v500 = -2.1092370507e-4;
	double v501 =  3.0927427253e-6;
	double v510 =  1.3864594581e-6;
	double v600 =  3.1932457305e-5;
   double rho2, v_SA, v_CT, v_p, v;
   double x2,xs,ys,z,sa = 0.0;

   if (SA > 0.0) sa = SA;

   x2 = gtc.gsw_sfac * sa;
   xs = sqrt(x2 + gtc.offset);
   ys = CT * 0.025;
   z = p * gtc.rec_db2pa;

   v_SA = gsvco.b000 + xs *(gsvco.b100 + xs *(gsvco.b200 + xs *(gsvco.b300 + xs *(gsvco.b400 + gsvco.b500 *xs))))
    + ys *(gsvco.b010 + xs *(gsvco.b110 + xs *(gsvco.b210 + xs *(gsvco.b310 + gsvco.b410 *xs)))
    + ys *(gsvco.b020 + xs *(gsvco.b120 + xs *(gsvco.b220 + gsvco.b320 *xs)) + ys *(gsvco.b030
    + xs *(gsvco.b130 + gsvco.b230 *xs) + ys *(gsvco.b040 + gsvco.b140 *xs + gsvco.b050 *ys))))
    + z *(gsvco.b001 + xs *(gsvco.b101 + xs *(gsvco.b201 + xs *(gsvco.b301 + gsvco.b401 *xs)))
    + ys *(gsvco.b011 + xs *(gsvco.b111 + xs *(gsvco.b211 + gsvco.b311 *xs)) + ys *(gsvco.b021
    + xs *(gsvco.b121 + gsvco.b221 *xs) + ys *(gsvco.b031 + gsvco.b131 *xs + gsvco.b041 *ys)))
    + z *(gsvco.b002 + xs *(gsvco.b102 + xs *(gsvco.b202 + gsvco.b302 *xs))+ ys *(gsvco.b012
    + xs *(gsvco.b112 + gsvco.b212 *xs) + ys *(gsvco.b022 + gsvco.b122 *xs + gsvco.b032 *ys))
    + z *(gsvco.b003 +  gsvco.b103 *xs + gsvco.b013 *ys + gsvco.b004 *z)));

   v_CT = gsvco.a000 + xs *(gsvco.a100 + xs *(gsvco.a200 + xs *(gsvco.a300 + xs *(gsvco.a400 + gsvco.a500 *xs))))
    + ys *(gsvco.a010 + xs *(gsvco.a110 + xs *(gsvco.a210 + xs *(gsvco.a310 + gsvco.a410 *xs)))
    + ys *(gsvco.a020 + xs *(gsvco.a120 + xs *(gsvco.a220 + gsvco.a320 *xs)) + ys *(gsvco.a030
    + xs *(gsvco.a130 + gsvco.a230 *xs) + ys *(gsvco.a040 + gsvco.a140 *xs + gsvco.a050 *ys ))))
    + z *(gsvco.a001 + xs *(gsvco.a101 + xs *(gsvco.a201 + xs *(gsvco.a301 + gsvco.a401 *xs)))
    + ys *(gsvco.a011 + xs *(gsvco.a111 + xs *(gsvco.a211 + gsvco.a311 *xs)) + ys *(gsvco.a021
    + xs *(gsvco.a121 + gsvco.a221 *xs) + ys *(gsvco.a031 + gsvco.a131 *xs + gsvco.a041 *ys)))
    + z *(gsvco.a002 + xs *(gsvco.a102 + xs *(gsvco.a202 + gsvco.a302 *xs)) + ys *(gsvco.a012
    + xs *(gsvco.a112 + gsvco.a212 *xs) + ys *(gsvco.a022 + gsvco.a122 *xs + gsvco.a032 *ys))
    + z *(gsvco.a003 + gsvco.a103 *xs + gsvco.a013 *ys + gsvco.a004 *z))) ;

   v_p = gsvco.c000 + xs *(gsvco.c100 + xs *(gsvco.c200 + xs *(gsvco.c300 + xs *(gsvco.c400 + gsvco.c500 *xs))))
    + ys *(gsvco.c010 + xs *(gsvco.c110 + xs *(gsvco.c210 + xs *(gsvco.c310 + gsvco.c410 *xs))) + ys *(gsvco.c020
    + xs *(gsvco.c120 + xs *(gsvco.c220 + gsvco.c320 *xs)) + ys *(gsvco.c030 + xs *(gsvco.c130 + gsvco.c230 *xs)
    + ys *(gsvco.c040 + gsvco.c140 *xs + gsvco.c050 *ys)))) + z *(gsvco.c001 + xs *(gsvco.c101 + xs *(gsvco.c201
    + xs *(gsvco.c301 + gsvco.c401 *xs))) + ys *(gsvco.c011 + xs *(gsvco.c111 + xs *(gsvco.c211 + gsvco.c311 *xs))
    + ys *(gsvco.c021 + xs *(gsvco.c121 + gsvco.c221 *xs) + ys *(gsvco.c031 + gsvco.c131 *xs + gsvco.c041 *ys)))
    + z *(gsvco.c002 + xs *(gsvco.c102 + gsvco.c202 *xs) + ys *(gsvco.c012 + gsvco.c112 *xs + gsvco.c022 *ys)
    + z *(gsvco.c003 + gsvco.c103 *xs + gsvco.c013 *ys + z *(gsvco.c004 + gsvco.c005 *z))));

   v = gsvco.v000 + xs *(v100 + xs *(v200 + xs *(v300 + xs *(v400 + xs *(v500
    + v600 *xs))))) + ys *(v010 + xs *(v110 + xs *(v210 + xs *(v310 + xs *(v410
    + v510 *xs)))) + ys *(v020 + xs *(v120 + xs *(v220 + xs *(v320 +v420 *xs)))
    + ys *(v030 + xs *(v130 + xs *(v230 +v330 *xs)) + ys *(v040 + xs *(v140
    + v240*xs) + ys *(v050 +v150 *xs +v060 *ys))))) + z *(gsvco.v001 + xs *(v101
    + xs *(v201 + xs *(v301 + xs *(v401 +v501 *xs)))) + ys *(v011 + xs *(v111
    + xs *(v211 + xs *(v311 +v411 *xs))) + ys *(v021 + xs *(v121 + xs *(v221
    + v321 *xs)) + ys *(v031 + xs *(v131 +v231 *xs) + ys *(v041 +v141 *xs
    + v051 *ys)))) + z *(gsvco.v002 + xs *(v102 + xs *(v202 + xs *(v302 + v402 *xs)))
    + ys *(v012 + xs *(v112 + xs *(v212 +v312 *xs)) + ys *(v022 + xs *(v122
    + v222 *xs) + ys *(v032 +v132 *xs +v042 *ys))) + z *(gsvco.v003 + xs *(v103
    + v203 *xs) + ys *(v013 +v113 *xs +v023 *ys) + z *(gsvco.v004 +v104 *xs +v014 *ys
    + z *(gsvco.v005 + gsvco.v006 *z)))));

   rho2 = (1.0 / v) * (1.0 / v);

   rho_SA = -rho2 * 0.5 * gtc.gsw_sfac * v_SA / xs;

   rho_CT = -rho2 * 0.025 * v_CT;

   rho_P = -rho2 * gtc.rec_db2pa * gtc.rec_db2pa * v_p;
}
/**
% gsw_rho_second_derivatives_wrt_enthalpy                second derivatives
%                        of rho with respect to enthalpy (75-term equation)
% =========================================================================
%
% USAGE: gsw_rho_second_derivatives_wrt_enthalpy(SA,CT,p, *rho_SA_SA,
                                                            *rho_SA_h,
                                                            *rho_h_h)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of rho with
%  respect to enthalpy,
%   (1) rho_SA_SA, second-order derivative with respect to Absolute Salinity
%       at constant h & p.
%   (2) rho_SA_h, second-order derivative with respect to SA & h at
%       constant p.
%   (3) rho_h_h, second-order derivative with respect to h at
%       constant SA & p.
%
%  Note that this function uses the using the computationally-efficient
%  expression for specific volume (Roquet et al., 2015).  There is an
%  alternative to calling this function, namely
%  gsw_rho_second_derivatives_wrt_enthalpy_CT_exact(SA,CT,p) which uses
%  the full Gibbs function (IOC et al., 2010).
%
%  This 75-term equation has been fitted in a restricted range of parameter
%  space, and is most accurate inside the "oceanographic funnel" described
%  in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_SA_SA = The second-order derivative of rho with respect to
%              Absolute Salinity at constant h & p.     [ J/(kg (g/kg)^2) ]
%  rho_SA_h  = The second-order derivative of rho with respect to
%              SA and h at constant p.                   [ J/(kg K(g/kg)) ]
%  rho_h_h   = The second-order derivative of rho with respect to h at
%              constant SA & p
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_rho_second_derivatives_wrt_enthalpy(double SA, double CT,
		double p,
		double &rho_SA_SA,
		double &rho_SA_h,
		double &rho_h_h)
{
   double *v_SA_SA=new double, *v_SA_h=new double;
   double *v_h_h=new double, *v_SA=new double, *v_h=new double;
   double rec_v, rec_v2, rec_v3, sa = 0.0;

   if (SA > 0.0) sa = SA;

   rec_v = 1.0 / gsw_specvol(sa,CT,p);
   gsw_specvol_first_derivatives_wrt_enthalpy(sa,CT,p,v_SA, v_h);
   gsw_specvol_second_derivatives_wrt_enthalpy(sa,CT,p,v_SA_SA, v_SA_h, v_h_h);

   rec_v2 = rec_v * rec_v;
   rec_v3 = rec_v2 * rec_v;

   rho_h_h = -*v_h_h * rec_v2 + 2.0 * (*v_h * *v_h) * rec_v3;

   rho_SA_h = -*v_SA_h * rec_v2 + 2.0 * *v_SA * *v_h * rec_v3;

   rho_SA_SA = -*v_SA_SA * rec_v2 + 2.0 * (*v_SA * *v_SA) * rec_v3;

   /** clean up */
   delete v_SA_SA;
   delete v_SA_h;
   delete v_h_h;
   delete v_SA;
   delete v_h;
}

/**
% gsw_specvol_p_parts                        specific volume pressure terms
%                                                        (75-term equation)
%==========================================================================
%
% USAGE: gsw_specvol_p_parts(SA,CT, *vp0,  *vp1, *vp2,  *vp3,  *vp4,  *vp5,
												  *vp6)
%
% DESCRIPTION:
%  Calculates each of the pressure terms in specific volume from Absolute
%  Salinity and Conservative Temperature, using the computationally
%  efficient 75-term polynomial expression for specific volume (Roquet et
%  al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is available to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  vp0 = Absolute Salinity and Conservative Temperature terms for pressure
%         to the 0
%  vp1 = Absolute Salinity and Conservative Temperature terms for pressure
%         to the 1
%  vp2 = Absolute Salinity and Conservative Temperature terms for pressure
%         to the 2
%  vp3 = Absolute Salinity and Conservative Temperature terms for pressure
%         to the 3
%  vp4 = Absolute Salinity and Conservative Temperature terms for pressure
%         to the 4
%  vp5 = Absolute Salinity and Conservative Temperature terms for pressure
%         to the 5
%  vp6 = Absolute Salinity and Conservative Temperature terms for pressure
%         to the 6
%
% AUTHOR:
%  Fabien Roquet, Gurvan Madec, Trevor McDougall & Paul Barker
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmos. Ocean. Tech., 20,
%   730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
% ==========================================================================
*/
void TeosBase::gsw_specvol_p_parts(double SA, double CT, double &vp0, double &vp1,
												  double &vp2, double &vp3, double &vp4, double &vp5,
												  double &vp6)
{
   double ys,xs,x2,sa = 0.0;

   if (SA > 0.0) sa = SA;

   x2 = gtc.gsw_sfac * sa;
   xs = sqrt(x2 + gtc.offset);
   ys = CT * 0.025;

   vp0 = gsvco.v000 + xs *(gsvco.v100 + xs *(gsvco.v200 + xs *(gsvco.v300
      + xs *(gsvco.v400 + xs * (gsvco.v500 + xs * gsvco.v600)))))
      + ys *(gsvco.v010 + xs * (gsvco.v110 + xs * (gsvco.v210
      + xs *(gsvco.v310 + xs * (gsvco.v410 + xs * gsvco.v510))))
      + ys *(gsvco.v020 + xs * (gsvco.v120 + xs * (gsvco.v220
      + xs *(gsvco.v320 + xs * gsvco.v420)))
      + ys *(gsvco.v030 + xs *(gsvco.v130 + xs *(gsvco.v230
      + xs * gsvco.v330)) + ys *(gsvco.v040 + xs *(gsvco.v140 + xs * gsvco.v240)
      + ys *(gsvco.v050 + xs * gsvco.v150 + ys * gsvco.v060)))));

   vp1 =  gsvco.v001 + xs *(gsvco.v101 + xs *(gsvco.v201
      + xs * (gsvco.v301 + xs *(gsvco.v401 + xs * gsvco.v501))))
      + ys *(gsvco.v011 + xs *(gsvco.v111 + xs *(gsvco.v211
      + xs *(gsvco.v311 + xs * gsvco.v411)))
      + ys *(gsvco.v021 + xs *(gsvco.v121 + xs *(gsvco.v221 + xs * gsvco.v321))
      + ys *(gsvco.v031 + xs *(gsvco.v131 + xs * gsvco.v231)
      + ys *(gsvco.v041 + xs * gsvco.v141 + ys * gsvco.v051))));

   vp2 = gsvco.v002 + xs *(gsvco.v102 + xs *(gsvco.v202 + xs *(gsvco.v302 + xs * gsvco.v402)))
      + ys *(gsvco.v012 + xs *(gsvco.v112 + xs *(gsvco.v212 + xs * gsvco.v312))
      + ys *(gsvco.v022 + xs *(gsvco.v122 + xs * gsvco.v222)
      + ys *(gsvco.v032 + xs * gsvco.v132 + ys * gsvco.v042)));

   vp3 = gsvco.v003 + xs *(gsvco.v103 + xs * gsvco.v203)
      + ys *(gsvco.v013 + xs * gsvco.v113
      + ys * gsvco.v023);

   vp4 = gsvco.v004 + xs * gsvco.v104 + ys * gsvco.v014;

   vp5 = gsvco.v005;

   vp6 = gsvco.v006;

// Note that,
// specvol = vp0 + z.*(vp1 + z.*(vp2 + z.*(vp3
//               + z.*(vp4 + z.*(vp5 + z.*vp6)))));
}
/**
% gsw_t90_from_t68              ITS-90 temperature from IPTS-68 temperature
%==========================================================================
%
% USAGE: gsw_t90_from_t68(t68)
%
% DESCRIPTION:
%  Converts IPTS-68 temperature to International Temperature Scale 1990
%  (ITS-90) temperature.  This conversion should be applied to all in-situ
%  data collected between 1/1/1968 and 31/12/1989.
%
% INPUT:
%  t68  =  in-situ temperature  (IPTS-68)                         [ deg C ]
%
% OUTPUT:
%  t90  =  in-situ temperature  (ITS-90)                          [ deg C ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  International Temperature Scales of 1948, 1968 and 1990, an ICES
%   note, available from http://www.ices.dk/ocean/procedures/its.htm
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.1 and appendix A.1 of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_t90_from_t68(double t68)
{
	return t68 * 0.999760057586179;
}
/**
% gsw_weighted_nanmean                      weighted mean of non-NaN values
%==========================================================================
%
% USAGE: gsw_weighted_nanmean(data,weights, *weighted_mean, *weighted_weights)
%
% DESCRIPTION:
%  This function returns the weighted mean of the data, where the weights
%  are given by weights and a coresponding weighted weight for the weighted
%  mean.  If the weights are not supplied then the weights are considerd to
%  be all equal. This programme treats NaN|s as missing values.
%
% INPUT:
%  data    = data values
%
% Optional
%  weights = weights of the data
%
% OUTPUT:
%  weighted_mean     =  weighted mean of the data
%  weighted_weights  =  weighted weights of the weights
%
% AUTHOR:
%  Paul Barker and Trevor McDougall, based on nanmean from Matlab.
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_weighted_nanmean(double data, double weights,
													double &weighted_mean,
													double &weighted_weights)
{
	double unnamed_idx_0;
	double weighted_data;
	int nz;
	bool Inan;

	weighted_data = weights * data;
	Inan = gsw_rtIsNaN(weighted_data);
	unnamed_idx_0 = weights;

	if (Inan)
	{
		weighted_data = 0.0;
		unnamed_idx_0 = 0.0;
	}

	if (unnamed_idx_0 == 0.0)
	{
		unnamed_idx_0 = rtNaN;
	}

	weighted_mean = weighted_data / unnamed_idx_0;
	nz = !Inan;
	weighted_data = nz;

	if (nz == 0)
	{
		weighted_data = rtNaN;
	}

	weighted_weights = unnamed_idx_0 / weighted_data;
}
/**
% gsw_rho_second_derivatives_CT_exact              SA and CT partial second
%                                              order derivatives of density
%==========================================================================
%
% USAGE: gsw_rho_second_derivatives_CT_exact(SA,CT,p, *rho_SA_SA,
                                                      *rho_SA_CT,
                                                      *rho_CT_CT,
                                                      *rho_SA_P,
                                                      *rho_CT_P)
%
% DESCRIPTION:
%  Calculates the following five second-order derivatives of rho,
%   (1) rho_SA_SA, second-order derivative with respect to Absolute
%       Salinity at constant CT & p.
%   (2) rho_SA_CT, second-order derivative with respect to SA & CT at
%       constant p.
%   (3) rho_CT_CT, second-order derivative with respect to CT at
%       constant SA & p.
%   (4) rho_SA_P, second-order derivative with respect to SA & P at
%       constant CT.
%   (5) rho_CT_P, second-order derivative with respect to CT & P at
%       constant SA.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_rho_second_derivatives(SA,CT,p), which uses the computationally
%  efficient polynomial for specific volume in terms of SA, CT and p
%  (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_SA_SA = The second-order derivative of rho with respect to
%              Absolute Salinity at constant CT & p.  [ (kg/m^3)(g/kg)^-2 ]
%  rho_SA_CT = The second-order derivative of rho with respect to
%              SA and CT at constant p.          [ (kg/m^3)(g/kg)^-1 K^-1 ]
%  rho_CT_CT = The second-order derivative of rho with respect to CT at
%              constant SA & p                            [ (kg/m^3) K^-2 ]
%  rho_SA_P  = The second-order derivative with respect to SA & P at
%              constant CT.                    [ (kg/m^3) (g/kg)^-1 Pa^-1 ]
%  rho_CT_P  = The second-order derivative with respect to CT & P at
%              constant SA.                         [ (kg/m^3) K^-1 Pa^-1 ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
void TeosBase::gsw_rho_second_derivatives_CT_exact(double SA, double CT, double p,
		double *rho_SA_SA, double *rho_SA_CT,
		double *rho_CT_CT, double *rho_SA_P,
		double *rho_CT_P)
{
   double *v_SA_SA=new double, *v_SA_CT=new double, *v_CT_CT=new double;
   double *v_SA_P=new double, *v_CT_P=new double;
   double *v_SA=new double, *v_CT=new double, *v_P=new double;
   double rec_v, rec_v2, rec_v3, sa = 0.0;

   if (SA > 0.0) sa = SA;

   rec_v = 1.0 / gsw_specvol_CT_exact(sa,CT,p);
   gsw_specvol_first_derivatives_CT_exact(sa,CT,p,*v_SA, *v_CT, *v_P);
   gsw_specvol_second_derivatives_CT_exact(sa,CT,p,*v_SA_SA, *v_SA_CT, *v_CT_CT, *v_SA_P, *v_CT_P);

   rec_v2 = rec_v * rec_v;
   rec_v3 = rec_v2 * rec_v;

   *rho_CT_CT = -*v_CT_CT * rec_v2 + 2.0 * (*v_CT * *v_CT) * rec_v3;

   *rho_SA_CT = -*v_SA_CT * rec_v2 + 2.0 * *v_SA * *v_CT * rec_v3;

   *rho_SA_SA = -*v_SA_SA * rec_v2 + 2.0 * (*v_SA * *v_SA) * rec_v3;

   *rho_SA_P = -*v_SA_P * rec_v2 + 2.0 * *v_SA * *v_P * rec_v3;

   *rho_CT_P = -*v_CT_P * rec_v2 + 2.0 * *v_CT * *v_P * rec_v3;

   /** clean up */
   delete v_SA_SA;
   delete v_SA_CT;
   delete v_CT_CT;
   delete v_SA_P;
   delete v_CT_P;
   delete v_SA;
   delete v_CT;
   delete v_P;
}
/**
% gsw_specvol_diff         specific volume difference between two pressures
%                                                        (75-term equation)
%==========================================================================
%
% USAGE: gsw_specvol_diff(SA,CT,p_shallow,p_deep)
%
% DESCRIPTION:
%  Calculates the difference of the specific volume of seawater between
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  This
%  function uses the computationally-efficient expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).  The output
%  (specvol_diff) is the specific volume evaluated at (SA,CT,p_deep)
%  minus the specific volume at (SA,CT,p_shallow).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one|s data lies outside this "funnel".
%
% INPUT:
%  SA     =  Absolute Salinity                                     [ g/kg ]
%  CT     =  Conservative Temperature (ITS-90)                    [ deg C ]
%  p_shallow  =  upper sea pressure                                [ dbar ]
%                ( i.e. shallower absolute pressure - 10.1325 dbar )
%  p_deep     =  lower sea pressure                                [ dbar ]
%                ( i.e. deeper absolute pressure - 10.1325 dbar )
%
%  SA & CT  need to have the same dimensions.
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN,
%  where SA and CT are MxN.
%
% OUTPUT:
%  specvol_diff  =  difference of specific volume                [ m^3/kg ]
%                       (deep minus shallow)
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_diff(double SA, double CT, double p_shallow,
												 double p_deep)
{
	double v011 =  1.8505765429e-5;
	double v012 = -1.1736386731e-6;
	double v013 = -3.6527006553e-7;
	double v014 =  3.1454099902e-7;
	double v021 = -1.1716606853e-5;
	double v022 =  2.1305028740e-6;
	double v023 =  2.8695905159e-7;
	double v031 =  7.9279656173e-6;
	double v032 = -4.6132540037e-7;
	double v041 = -3.4102187482e-6;
	double v042 = -6.3352916514e-8;
	double v051 =  5.0736766814e-7;
	double v101 =  2.4262468747e-5;
	double v102 = -5.8484432984e-7;
	double v103 =  3.6310188515e-7;
	double v104 = -1.1147125423e-7;
	double v121 = -2.3678308361e-7;
	double v122 =  3.9137387080e-7;
	double v131 = -3.4558773655e-6;
	double v132 =  7.7618888092e-9;
	double v141 =  1.2956717783e-6;
	double v201 = -3.4792460974e-5;
	double v202 = -4.8122251597e-6;
	double v203 =  1.6746303780e-8;
	double v211 =  1.1100834765e-5;
	double v212 =  5.4620748834e-6;
	double v231 =  3.1655306078e-7;
	double v301 =  3.7470777305e-5;
	double v302 =  4.9263106998e-6;
	double v311 = -9.8447117844e-6;
	double v312 = -1.3544185627e-6;
	double v321 = -4.8826139200e-7;
	double v401 = -1.7322218612e-5;
	double v402 = -1.7811974727e-6;
	double v411 =  2.5909225260e-6;
	double v501 =  3.0927427253e-6;

	double dz,z_part_2,x2, xs, ys, z_shallow, z_shallow2, z_shallow3;
	double z_deep, z_deep2, z_deep_times_z_shallow,z_deep2_plus_z_shallow2;
   double xy_part_1, xy_part_2, xy_part_3, xy_part_4, z_part_3, z_part_4;
   double z_part_5, z_part_6;
   double sa = 0.0;

   if (SA > 0.0) sa = SA;

   x2 = gtc.gsw_sfac * sa;
   xs = sqrt(x2 + gtc.offset);
   ys = CT * 0.025;
   z_shallow = p_shallow * gtc.rec_db2pa;
   z_deep = p_deep * gtc.rec_db2pa;

   xy_part_1 = gsvco.v001 + xs *(v101 + xs *(v201 + xs *(v301 + xs *(v401 + v501 *xs))))
      + ys *(v011 + xs *(gsvco.v111 + xs *(v211 + xs *(v311 + v411 *xs))) + ys *(v021
      + xs *(v121 + xs *(gsvco.v221 + v321 *xs)) + ys *(v031 + xs *(v131 + v231 *xs)
      + ys *(v041 + v141 *xs + v051 *ys))));

   xy_part_2 = gsvco.v002 + xs *(v102 + xs *(v202 + xs *(v302 + v402 *xs)))
      + ys *(v012 + xs *(gsvco.v112 + xs *(v212 + v312 *xs)) + ys *(v022 + xs *(v122
      + gsvco.v222 *xs) + ys *(v032 + v132 *xs + v042 *ys)));

   xy_part_3 = gsvco.v003 + xs *(v103 + v203 *xs) + ys *(v013 + gsvco.v113 *xs + v023 *ys);

   xy_part_4 = gsvco.v004 + v104 *xs + v014 *ys;

   dz = z_deep - z_shallow;

   z_part_2 = z_deep + z_shallow;

   z_deep2 = z_deep * z_deep;
   z_shallow2 = z_shallow * z_shallow;
   z_deep_times_z_shallow = z_deep * z_shallow;
   z_deep2_plus_z_shallow2 = z_deep2 + z_shallow2;

   z_part_3 = z_deep2_plus_z_shallow2 + z_deep_times_z_shallow;

   z_part_4 = z_part_2 * z_deep2_plus_z_shallow2;

   z_shallow3 = z_shallow * z_shallow2;

   z_part_5 = z_deep2 * z_part_3 + z_shallow3 * z_part_2;

   z_part_6 = z_part_2 * z_part_3 * (z_deep2_plus_z_shallow2 - z_deep_times_z_shallow);

   return dz * (xy_part_1 + z_part_2 *xy_part_2 + z_part_3 * xy_part_3
            + z_part_4 * xy_part_4 + z_part_5 * gsvco.v005
            + z_part_6 * gsvco.v006);
}
/**
% gsw_sigma0_CT_exact                        potential density anomaly with
%                                         reference sea pressure of 0 dbar.
%==========================================================================
%
% USAGE: gsw_sigma0_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 0 dbar,
%  this being this particular potential density minus 1000 kg/m^3.  This
%  function has inputs of Absolute Salinity and Conservative Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_sigma0(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma0_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 0 dbar,
%                      that is, this potential density - 1000 kg/m^3.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_sigma0_CT_exact(double SA, double CT)
{
   double pt0 = gsw_pt_from_ct(SA, CT);

   return gsw_sigma0_pt0_exact(SA, pt0);
}


/**
% gsw_N2sol                                    solubility of N2 in seawater
%==========================================================================
%
% USAGE:
%  N2sol = gsw_N2sol(SA,CT,p,long,lat)
%
% DESCRIPTION:
%  Calculates the nitrogen, N2, concentration expected at equilibrium with
%  air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
%  including saturated water vapor.  This function uses the solubility
%  coefficients as listed in Hamme and Emerson (2004).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SA & CT need to have the same dimensions. p, lat and long may have
%  dimensions 1x1 or Mx1 or 1xN or MxN, where SA and CT are MxN.
%
% OUTPUT:
%  N2sol = solubility of nitrogen in micro-moles per kg         [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Hamme, R., and S. Emerson, 2004: The solubility of neon, nitrogen and
%   argon in distilled water and seawater. Deep-Sea Research, 51,
%   1517-1528.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_N2sol(double SA,double CT,double p,double lon,double lat)
{
   // The coefficents below are from Table 4 of Hamme and Emerson (2004)
   double a0 = 6.42931;
   double a1 = 2.92704;
   double a2 = 4.32531;
   double a3 = 4.69149;
   double b0 = -7.44129e-3;
   double b1 = -8.02566e-3;
   double b2 = -1.46775e-2;

   // Note that salinity argument is Practical Salinity, this is
   // beacuse the major ionic components of seawater related to Cl
   // are what affect the solubility of non-electrolytes in seawater.

   double y,pt,N2sol, SP = gsw_sp_from_sa(SA,p,lon,lat);

   pt = gsw_pt_from_ct(SA,CT); // pt is potential temperature referenced to the sea surface.

   y = log((298.15 - pt) / (gtc.gsw_t0 + pt));  // pt is the temperature in degress C
                                                // on the 1990 International Temperature Scale ITS-90.

   N2sol = exp(a0 + y * (a1 + y * (a2 + a3 * y)) + SP * (b0 + y * (b1 + b2 * y)));

   return N2sol;
}

/**
% gsw_N2Osol                                  solubility of N2O in seawater
%==========================================================================
%
% USAGE:
%  N2Osol = gsw_N2Osol(SA,CT,p,long,lat)
%
% DESCRIPTION:
%  Calculates the nitrous oxide, N2O, concentration expected at equilibrium
%  with air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
%  including saturated water vapor.  This function uses the solubility
%  coefficients as listed in Weiss and Price (1980).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SA & CT need to have the same dimensions. p, lat and long may have
%  dimensions 1x1 or Mx1 or 1xN or MxN, where SA and CT are MxN.
%
% OUTPUT:
%  N2Osol = solubility of nitrous oxide                           [ mol/L ]
%
% AUTHOR:  Rich Pawlowicz, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Weiss, R.F., and B.A. Price, 1980: Nitrous oxide solubility in water and
%   seawater. Mar. Chem., 8, 347-359.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_N2Osol(double SA,double CT,double p,double lon,double lat)
{
   double ph2odP,y_100,y,pt68,N2Osol,pt,SP = gsw_sp_from_sa(SA,p,lon,lat);

   // The coefficents below are from Table 2 of Weiss and Price (1980)
   double a0 = -165.8806;
   double a1 =  222.8743;
   double a2 =  92.0792;
   double a3 = -1.48425;

   double b1 = -0.056235;
   double b2 =  0.031619;
   double b3 = -0.0048472;

   double m0 = 24.4543;
   double m1 = 67.4509;
   double m2 = 4.8489;
   double m3 = 0.000544;

   // Note that salinity argument is Practical Salinity, this is
   // beacuse the major ionic components of seawater related to Cl
   // are what affect the solubility of non-electrolytes in seawater.

   pt = gsw_pt_from_ct(SA,CT); // pt is potential temperature referenced to the sea surface.

   pt68 = pt * 1.00024; // pt68 is the potential temperature in degress C on
                        // the 1968 International Practical Temperature Scale IPTS-68.
   y = pt68 + gtc.gsw_t0;
   y_100 = y * 1e-2;

   ph2odP = exp(m0 - m1  *100.0 / y - m2 * log(y_100) - m3 * SP); // Moist air correction at 1 atm.

   N2Osol = (exp(a0 + a1 * 100.0 / y + a2 * log(y_100) + a3 * (y_100 * y_100)
               + SP *(b1 + y_100 * (b2 + b3 * y_100)))) / (1.0 - ph2odP);

   return N2Osol;
}
/**
% gsw_Gibbs_energy_t_exact                         Gibbs energy of seawater
%==========================================================================
%
% USAGE:
%  Gibbs_energy_t_exact = gsw_Gibbs_energy_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the Gibbs energy of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  Gibbs_energy_t_exact = Gibbs energy                             [ J/kg ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.13 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Gibbs_energy_t_exact(double SA,double t,double p)
{
	return gsw_gibbs(0,0,0,SA,t,p);
}
/**
% gsw_Hesol                                    solubility of He in seawater
%==========================================================================
%
% USAGE:
%  Hesol = gsw_Hesol(SA,CT,p,long,lat)
%
% DESCRIPTION:
%  Calculates the helium concentration expected at equilibrium with air at
%  an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
%  saturated water vapor.  This function uses the solubility coefficients
%  as listed in Weiss (1971).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SA & CT need to have the same dimensions. p, lat and long may have
%  dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  Hesol = solubility of helium in micro-moles per kg           [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Dymond and Smith, 1980: The virial coefficients of pure gases and
%   mixtures. Clarendon Press, Oxford.
%
%  Weiss, R.F., 1971: Solubility of Helium and Neon in Water and Seawater.
%   J. Chem. and Engineer. Data, 16, 235-241.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Hesol(double SA,double CT,double p,double lon,double lat)
{
   double He_ml2umol, Hesol, y, y_100, pt68, pt, SP = gsw_sp_from_sa(SA,p,lon,lat);

   // Note that salinity argument is Practical Salinity, this is
   // beacuse the major ionic components of seawater related to Cl
   // are what affect the solubility of non-electrolytes in seawater.

   pt = gsw_pt_from_ct(SA,CT);// pt is potential temperature referenced to
                              // the sea surface.

   pt68 = pt * 1.00024; // pt68 is the potential temperature in degress C on
                        // the 1968 International Practical Temperature Scale IPTS-68.

   y = pt68 + gtc.gsw_t0;
   y_100 = y * 1e-2;

   // The coefficents below are from Table 3 of Weiss (1971)
   double a1 = -167.2178;
   double a2 =  216.3442;
   double a3 =  139.2032;
   double a4 = -22.6202;
   double b1 = -0.044781;
   double b2 =  0.023541;
   double b3 = -0.0034266;

   Hesol = exp(a1 + a2 * 100.0 / y + a3 * log(y_100) + a4 * y_100
            + SP * (b1 + y_100 * (b2 + b3 * y_100)));

   He_ml2umol = 4.455817671505537e1;// mL/kg to umol/kg for He (1/22.44257e-3)
                                    // Molar volume at STP (Dymond and Smith, 1980).
   Hesol = Hesol * He_ml2umol;
   return Hesol;
}

/**
% gsw_specvol_from_pot_enthalpy_ice                    specific volume from
%                                                 potential enthalpy of ice
%==========================================================================
%
% USAGE:
%  specvol_ice = gsw_specvol_from_pot_enthalpy_ice(pot_enthalpy_ice,p)
%
% DESCRIPTION:
%  Calculates the specific volume of ice from the potential enthalpy
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar.
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_specvol_from_pot_enthalpy_ice(double pot_enthalpy_ice, double p)
{
   double pt0_ice = gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice);

   double t_ice = gsw_t_from_pt0_ice(pt0_ice, p);

   return gsw_specvol_ice(t_ice, p);
}
/**
% gsw_Krsol                                    solubility of Kr in seawater
%==========================================================================
%
% USAGE:
%  Krsol = gsw_Krsol(SA,CT,p,long,lat)
%
% DESCRIPTION:
%  Calculates the krypton, Kr, concentration expected at equilibrium with
%  air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar)
%  including saturated water vapor.  This function uses the solubility
%  coefficients derived from the data of Weiss and Kyser (1978).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SA & CT need to have the same dimensions. p, lat and long may have
%  dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  Krsol = solubility of krypton in micro-moles per kg          [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Weiss, R.F., and T.K. Kyser, 1978: Solubility of Krypton in Water and
%   Seawater. J. Chem. Thermodynamics, 23, 69-72.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Krsol(double SA,double CT,double p,double lon,double lat)
{
   double Kr_ml2umol,Krsol,pt,pt68,y,y_100,SP = gsw_sp_from_sa(SA,p,lon,lat);

   // Note that salinity argument is Practical Salinity, this is
   // beacuse the major ionic components of seawater related to Cl
   // are what affect the solubility of non-electrolytes in seawater.

   pt = gsw_pt_from_ct(SA,CT);// pt is potential temperature referenced to
                              // the sea surface.

   pt68 = pt * 1.00024; // pt68 is the potential temperature in degress C on
                        // the 1968 International Practical Temperature Scale IPTS-68.

   y = pt68 + gtc.gsw_t0;
   y_100 = y * 1e-2;

   // Table 2 (Weiss and Kyser, 1978)
   double a1 = -112.6840;
   double a2 =  153.5817;
   double a3 =  74.4690;
   double a4 = -10.0189;
   double b1 = -0.011213;
   double b2 = -0.001844;
   double b3 =  0.0011201;

   Krsol = exp(a1 + a2*100 / y + a3 * log(y_100) + a4 * y_100
            + SP * (b1 + y_100 * (b2 + b3 * y_100)));

   Kr_ml2umol = 4.474052731185490e1;   // mL/kg to umol/kg for Kr (1/22.3511e-3)
                                       // Molar volume at STP (Dymond and Smith, 1980).

   Krsol = Krsol * Kr_ml2umol;

   return Krsol;
}

/**
% gsw_enthalpy_diff_CT_exact        difference of enthalpy at two pressures
%==========================================================================
%
% USAGE:
%  enthalpy_diff_CT_exact = gsw_enthalpy_diff_CT_exact(SA,CT,p_shallow,p_deep)
%
% DESCRIPTION:
%  Calculates the difference of the specific enthalpy of seawater between
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  The output
%  (enthalpy_diff_CT_exact) is the specific enthalpy evaluated at
%  (SA,CT,p_deep) minus the specific enthalpy at (SA,CT,p_shallow).
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely gsw_enthalpy_diff(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA        =  Absolute Salinity                                  [ g/kg ]
%  CT        =  Conservative Temperature (ITS-90)                 [ deg C ]
%  p_shallow =  upper sea pressure                                 [ dbar ]
%               ( i.e. shallower absolute pressure - 10.1325 dbar )
%  p_deep    =  lower sea pressure                                 [ dbar ]
%               ( i.e. deeper absolute pressure - 10.1325 dbar )
%
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN,
%  where SA and CT are MxN.
%
% OUTPUT:
%  enthalpy_diff_CT_exact = difference of specific enthalpy        [ J/kg ]
%                           (deep minus shallow)
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) of this TEOS-10 Manual.
%
%  McDougall, T.J., 2003: Potential enthalpy: A conservative oceanic
%   variable for evaluating heat content and heat fluxes. Journal of
%   Physical Oceanography, 33, 945-963.
%    See Eqns. (18) and (22)
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_enthalpy_diff_CT_exact(double SA,double CT,double p_shallow,double p_deep)
{
   double t_shallow = gsw_t_from_ct(SA, CT, p_shallow);
   double t_deep = gsw_t_from_ct(SA, CT, p_deep);

   return gsw_enthalpy_t_exact(SA, t_deep, p_deep) -
                              gsw_enthalpy_t_exact(SA, t_shallow, p_shallow);
}
/**
% gsw_Arsol                                    solubility of Ar in seawater
%==========================================================================
%
% USAGE:
%  Arsol = gsw_Arsol(SA,CT,p,long,lat)
%
% DESCRIPTION:
%  Calculates the argon, Ar, concentration expected at equilibrium with air
%  at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
%  saturated water vapor.  This function uses the solubility coefficients
%  as listed in Hamme and Emerson (2004).
%
%  Note that this algorithm has not been approved by IOC and is not work
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SA & CT need to have the same dimensions. p, lat and long may have
%  dimensions 1x1 or Mx1 or 1xN or MxN, where SA and CT are MxN.
%
% OUTPUT:
%  Arsol = solubility of argon                                  [ umol/kg ]
%
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Hamme, R., and S. Emerson, 2004: The solubility of neon, nitrogen and
%   argon in distilled water and seawater. Deep-Sea Research, 51,
%   1517-1528.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================
*/
double TeosBase::gsw_Arsol(double SA,double CT,double p,double lon,double lat)
{
   double y,pt,SP = gsw_sp_from_sa(SA,p,lon,lat);

   // Note that salinity argument is Practical Salinity, this is
   // because the major ionic components of seawater related to Cl
   // are what affect the solubility of non-electrolytes in seawater.

   pt = gsw_pt_from_ct(SA,CT);// pt is potential temperature referenced to
                              // the sea surface.

   y = log((298.15 - pt) / (273.15 + pt));// pt is the temperature in degress C
                                          // on the 1990 International Temperature Scale ITS-90.

   //%The coefficents below are from Table 4 of Hamme and Emerson (2004)
   double a0 =  2.79150;
   double a1 =  3.17609;
   double a2 =  4.13116;
   double a3 =  4.90379;
   double b0 = -6.96233e-3;
   double b1 = -7.66670e-3;
   double b2 = -1.16888e-2;

   double Arsol = exp(a0 + y * (a1 + y * (a2 + a3 * y)) + SP *(b0 + y * (b1 + b2 * y)));

	return Arsol;
}
/**************************************
   returns true or false depending on
   equality or not of rtIsNaN with val
**************************************/
bool TeosBase::gsw_rtIsNaN(double val)
{
	return val == rtNaN;
}

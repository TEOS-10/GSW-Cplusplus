#include "TeosBase.h"

/******************************************
  "TeosBase.cpp" Version 1.0
  by Randall Kent Whited
  rkwhited@gmail.com
  ------------------
  This is a modification of the TEOS-10
  file "gsw_oceanographic_toolbox.c"
  as adapted from the TEOS-10
  C version 3.05 http://www.teos-10.org
  -------------------------------------
  All copyrights and all license issues
  are the same as for the C version
*******************************************/

/*****************************************
  The C++ methods herein are simply
  implementations of C functions in the
  TEOS-10 "gsw_oceanographic_toolbox.c"
  source code file.
*****************************************/

/**************************************
  This base-class constructor is
  where the initialization of the
  "SaarDataHandler" object and
  all struct objects takes place.
  -------------------------------
  In methods that have calls into
  the struct data, each struct
  object is identified at the
  beginning of the method (example:
  "for GSW_TEOS10_CONSTANTS use gtc"
  ----------------------------------
  "use gtc" means use the object
  rather than the struct name, such
  as "gtc.gsw_sfac" or "gtc.offset"
  ----------------------------------
  "gtc" is the struct object, the data
  name follows the "."
--------------------------------------
  The "gsw_deltasa_atlas" C function
  was modified into two C++ methods
  "gsw_deltasa_atlas" and method
  "gsw_util_indxPref" to facilitate
  clear references & calls to the
  "p_ref" array pointer.
***************************************/

typedef struct {
        double  d;
        int     i;
} DI;

TeosBase::TeosBase()
: /** initialize the complex numbers */
   t1(3.68017112855051e-2, 5.10878114959572e-2),
   t2(3.37315741065416e-1, 3.35449415919309e-1),
   r1(4.47050716285388e1,  6.56876847463481e1),
   r20(-7.25974574329220e1, -7.81008427112870e1),
   r21(-5.57107698030123e-5, 4.64578634580806e-5),
   r22(2.34801409215913e-11, -2.85651142904972e-11)
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

    /** initialize GSW_FREEZING_POLY_COEFFICIENTS */
    gfpc.init_GSW_FREEZING_POLY_COEFFICIENTS();

    /** initialize GSW_GIBBS_ICE_COEFFICIENTS */
    gic.init_GSW_GIBBS_ICE_COEFFICIENTS();
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

/**
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
*/
void TeosBase::gsw_add_barrier(double *input_data, double lon, double lat,
                                  double long_grid, double lat_grid, double dlong_grid,
                                  double dlat_grid, double *output_data)
{
    /** for GSW_SAAR_DATA use gsaar */

	int	above_line[4];
	int	k, nmean, above_line0, kk;
	double	r, lats_line, data_mean;

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

	above_line[0]	= (lats_line <= lat_grid);
	above_line[3]	= (lats_line <= lat_grid+dlat_grid);

	k = gsw_util_indx(gsaar.longs_pan,6,long_grid+dlong_grid);
			/**the 2 & 3 lon/lat points */
	r = (long_grid+dlong_grid-gsaar.longs_pan[k])/
			(gsaar.longs_pan[k+1]-gsaar.longs_pan[k]);

	lats_line	= gsaar.lats_pan[k] + r*(gsaar.lats_pan[k+1]-gsaar.lats_pan[k]);

	above_line[1]	= (lats_line <= lat_grid);
	above_line[2]	= (lats_line <= lat_grid+dlat_grid);

	nmean = 0;
	data_mean	= 0.0;

	for (kk=0; kk<4; kk++)
	{
      if ((fabs(input_data[kk]) <= 100.0) && above_line0 == above_line[kk])
		{
            nmean	= nmean+1;
            data_mean	= data_mean+input_data[kk];
      }
	}
	if (nmean == 0)
	    data_mean	= 0.0;	/**error return*/
	else
	    data_mean	= data_mean/nmean;

	for (kk=0; kk<4; kk++)
	{
      if ((fabs(input_data[kk]) >= 1.0e10) || above_line0 != above_line[kk])
		{
         output_data[kk]	= data_mean;
      }
      else
         output_data[kk]	= input_data[kk];
	}
}

/**
==========================================================================
method: gsw_add_mean(data_in,data_out)
==========================================================================

 Replaces NaN's with non-nan mean of the 4 adjacent neighbours

 data_in   : data set of the 4 adjacent neighbours

 data_out : non-nan mean of the 4 adjacent neighbours     [unitless]
*/
void TeosBase::gsw_add_mean(double *data_in, double *data_out)
{
	int	k, nmean;
	double	data_mean;

	nmean = 0;
	data_mean	= 0.0;

	for (k=0; k<4; k++)
	{
	    if (fabs(data_in[k]) <= 100.0)
	    {
            nmean++;
            data_mean	= data_mean+data_in[k];
	    }
	}

	if (nmean == 0.0)
	    data_mean	= 0.0;    /**error return*/
	else
	    data_mean	= data_mean/nmean;

	for (k=0; k<4; k++)
	{
      if (fabs(data_in[k]) >= 100.0)
		{
         data_out[k]	= data_mean;
      }
      else
         data_out[k]	= data_in[k];
	}

	return;
}

/** this method supplements gsw_alpha () */
double TeosBase::gsw_get_v_ct_part(double xs, double ys, double z)
{
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

/**
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
*/
double TeosBase::gsw_alpha(double sa, double ct, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */

	double	xs, ys, z, v_ct_part, rval;

	xs	= sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys	= ct*0.025;
	z	= p*1e-4;

	v_ct_part = gsw_get_v_ct_part(xs, ys, z);

	rval = 0.025*v_ct_part/gsw_specvol(sa,ct,p);

	return rval;
}

/** supplements gsw_ct_from_pt */
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

/**
==========================================================================
method: gsw_ct_from_pt(sa,pt)
==========================================================================

 Calculates Conservative Temperature from potential temperature of seawater

 sa      : Absolute Salinity                              [g/kg]
 pt      : potential temperature with                     [deg C]
           reference pressure of 0 dbar

 gsw_ct_from_pt : Conservative Temperature                [deg C]
*/
double TeosBase::gsw_ct_from_pt(double sa, double pt)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	double	x2, x, y, pot_enthalpy, rval;

	x2 = gtc.gsw_sfac*sa;
	x = sqrt(x2);
	y = pt*0.025e0;	/*! normalize for F03 and F08 */

	pot_enthalpy =  gsw_get_pot_ent(x2, x, y);

   rval = pot_enthalpy/gtc.gsw_cp0;

	return rval;
}


/**
==========================================================================
method: gsw_ct_from_t(sa,t,p)
==========================================================================

 Calculates Conservative Temperature from in-situ temperature

 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]

 gsw_ct_from_t : Conservative Temperature                 [deg C]
*/
double TeosBase::gsw_ct_from_t(double sa, double t, double p)
{
	double	pt0;

	pt0	= gsw_pt0_from_t(sa,t,p);

	return (gsw_ct_from_pt(sa,pt0));
}

/** this method supplements gsw_dynamic_enthalpy */
double TeosBase::gsw_get_dyn_enth_prt(double xs, double ys, double z)
{
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

/**
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
*/
double TeosBase::gsw_dynamic_enthalpy(double sa, double ct, double p)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

	double	dynamic_enthalpy_part, xs, ys, z, rval;

	xs	= sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys	= ct*0.025;
	z	= p*1e-4;

	dynamic_enthalpy_part = gsw_get_dyn_enth_prt(xs, ys, z);

   rval = dynamic_enthalpy_part*gtc.db2pa*1e4;

	return rval;
}


/**
==========================================================================
elemental method: gsw_enthalpy_ct_exact (sa, ct, p)
==========================================================================

  Calculates specific enthalpy of seawater from Absolute Salinity and
  Conservative Temperature and pressure.

  Note that this method: uses the full Gibbs function.

  SA  =  Absolute Salinity                                        [ g/kg ]
  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
  p   =  sea pressure                                             [ dbar ]
         ( i.e. absolute pressure - 10.1325 dbar )

  enthalpy_CT_exact  =  specific enthalpy                         [ J/kg ]
--------------------------------------------------------------------------
*/
double TeosBase::gsw_enthalpy_ct_exact(double sa, double ct, double p)
{
	double	t;

	t = gsw_t_from_ct(sa,ct,p);

	return (gsw_enthalpy_t_exact(sa,t,p));
}

/**
==========================================================================
method: gsw_enthalpy_t_exact(sa,t,p)
==========================================================================

 Calculates the specific enthalpy of seawater

 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]

 gsw_enthalpy_t_exact : specific enthalpy                 [J/kg]
*/
double TeosBase::gsw_enthalpy_t_exact(double sa, double t, double p)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	int	n0=0, n1=1;

	return (gsw_gibbs(n0,n0,n0,sa,t,p) -
		(t+gtc.gsw_t0)*gsw_gibbs(n0,n1,n0,sa,t,p));
}
/** this method supplements gsw_entropy_part */
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

/**
==========================================================================
method: gsw_entropy_part(sa,t,p)
==========================================================================

 entropy minus the terms that are a function of only SA

 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]

 entropy_part : entropy part
*/
double TeosBase::gsw_entropy_part(double sa, double t, double p)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	double	x2, x, y, z, rval;

	x2	= gtc.gsw_sfac*sa;
	x	= sqrt(x2);
	y	= t*0.025;
	z	= p*1e-4;

   rval = gsw_get_g03plusg08(z, y, x2, x);

	return rval;
}

/**
   these next ten (protected) methods/member functions supplement gsw_gibbs clarity:

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
*/
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
double TeosBase::gsw_gibbs_g08_a(double sa, double z, double y, double x)
{
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
double TeosBase:: gsw_gibbs_g03plusg08_b(double sa, double z, double y, double x2, double x)
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
double TeosBase:: gsw_gibbs_g03plusg08_d(double z, double y, double x2, double x)
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
double TeosBase:: gsw_gibbs_g08_b(double z, double y, double x)
{
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

	rval = g08*gtc.gsw_sfac*0.5e-8;;

	return rval;
}
double TeosBase::gsw_gibbs_g08_c(double sa, double z, double y, double x)
{
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
double TeosBase::gsw_gibbs_g03plusg08_f(double z, double y, double x2, double x)
{
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
/**
==========================================================================
method: gsw_gibbs(ns,nt,np,sa,t,p)
==========================================================================

 seawater specific Gibbs free energy and derivatives up to order 2

 ns     : order of s derivative
 nt     : order of t derivative
 np     : order of p derivative
 sa     : Absolute Salinity                               [g/kg]
 t      : temperature                                     [deg C]
 p      : sea pressure                                    [dbar]
 								-1
 gsw_gibbs  : specific Gibbs energy or its derivative	   [J kg  ]
*/
double TeosBase::gsw_gibbs(int ns, int nt, int np, double sa, double t, double p)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	double	x2, x, y, z, return_value = cppGSW_INVALID_VALUE;

	x2	= gtc.gsw_sfac*sa;
	x	= sqrt(x2);
	y	= t*0.025;
	z	= p*1e-4;

	if (ns == 0  && nt == 0  && np == 0)
	{
      return_value	= gsw_gibbs_g03plusg08_a(sa, z, y, x2, x);
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

	} else if (ns == 0  && nt == 2  && np == 0) {

      return_value	= gsw_gibbs_g03plusg08_d(z, y, x2, x);

	} else if (ns == 1  && nt == 0  && np == 1) {

      return_value = gsw_gibbs_g08_b(z, y, x);

	} else if (ns == 0  && nt == 1  && np == 1) {

      return_value = gsw_gibbs_g03plusg08_e(z, y, x2, x);

	} else if (ns == 1  && nt == 1  && np == 0) {

      return_value = gsw_gibbs_g08_c(sa, z, y, x);

	} else if (ns == 2  && nt == 0  && np == 0) {

      return_value = gsw_gibbs_g03plusg08_f(z, y, x2, x);

	} else if (ns == 0  && nt == 0  && np == 2) {

      return_value = gsw_gibbs_g03plusg08_g(z, y, x2, x);
	}

	return return_value;
}

/**
==========================================================================
method: gsw_pt_from_ct(sa,ct)
==========================================================================

 potential temperature of seawater from conservative temperature

 sa     : Absolute Salinity                               [g/kg]
 ct     : Conservative Temperature                        [deg C]
 p      : sea pressure                                    [dbar]

 gsw_pt_from_ct : potential temperature with              [deg C]
                  reference pressure of  0 dbar
*/
double TeosBase::gsw_pt_from_ct(double sa, double ct)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	double	a5ct, b3ct, ct_factor, pt_num, pt_recden, ct_diff;
	double	pt, pt_old, ptm, dpt_dct, s1;
	double	a0	= -1.446013646344788e-2,
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
    **  Start the 1.5 iterations through the modified Newton-Raphson
    **  iterative method:.
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
/**
==========================================================================
method: gsw_pt_from_t(sa,t,p,p_ref)
==========================================================================

 Calculates potential temperature of seawater from in-situ temperature

 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]
 p_ref  : reference sea pressure                          [dbar]

 gsw_pt_from_t : potential temperature                    [deg C]
*/
double TeosBase::gsw_pt_from_t(double sa, double t, double p, double p_ref)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	int	n0=0, n2=2, no_iter;
	double	s1, pt, ptm, pt_old, dentropy, dentropy_dt,
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

	for (no_iter=1; no_iter <= 2; no_iter++) {
	    pt_old	= pt;
	    dentropy	= gsw_entropy_part(sa,pt_old,p_ref) - true_entropy_part;
	    pt = pt_old - dentropy/dentropy_dt;
	    ptm = 0.5*(pt + pt_old);
	    dentropy_dt	= -gsw_gibbs(n0,n2,n0,sa,ptm,p_ref);
	    pt = pt_old - dentropy/dentropy_dt;
	}

	return (pt);
}


/**
==========================================================================
method: gsw_pt0_from_t(sa,t,p)
==========================================================================

 Calculates potential temperature with reference pressure, p_ref = 0 dbar.

 sa     : Absolute Salinity                               [g/kg]
 t      : in-situ temperature                             [deg C]
 p      : sea pressure                                    [dbar]

 gsw_pt0_from_t : potential temperature, p_ref = 0        [deg C]
*/
double TeosBase::gsw_pt0_from_t(double sa, double t, double p)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	int	no_iter;
	double	pt0, pt0_old, dentropy, dentropy_dt;
	double	s1, true_entropy_part, pt0m;

	s1	= sa/gtc.gsw_ups;

	pt0	= t+p*( 8.65483913395442e-6  -
        	  s1 *  1.41636299744881e-6  -
		   p *  7.38286467135737e-9  +
		   t *(-8.38241357039698e-6  +
		  s1 *  2.83933368585534e-8  +
		   t *  1.77803965218656e-8  +
		   p *  1.71155619208233e-10));

	dentropy_dt	= gtc.gsw_cp0/((gtc.gsw_t0+pt0)*(1.0-0.05*(1.0-sa/gtc.gsw_sso)));

	true_entropy_part = gsw_entropy_part(sa,t,p);

	for (no_iter=1; no_iter <= 2; no_iter++) {
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

/**
==========================================================================
method: gsw_gibbs_pt0_pt0(sa,pt0)
==========================================================================

 gibbs_tt at (sa,pt,0)

 sa     : Absolute Salinity                            [g/kg]
 pt0    : potential temperature                        [deg C]

 gibbs_pt0_pt0 : gibbs_tt at (sa,pt,0)
*/
double TeosBase::gsw_gibbs_pt0_pt0(double sa, double pt0)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	double	x2, x, y, g03, g08;

	x2	= gtc.gsw_sfac*sa;
	x	= sqrt(x2);
	y	= pt0*0.025;

	g03	= -24715.571866078 +
		y*(4420.4472249096725 +
		y*(-1778.231237203896 +
		y*(1160.5182516851419 +
		y*(-569.531539542516 + y*128.13429152494615))));

	g08	= x2*(1760.062705994408 + x*(-86.1329351956084 +
		x*(-137.1145018408982 + y*(296.20061691375236 +
		y*(-205.67709290374563 + 49.9394019139016*y))) +
		y*(-60.136422517125 + y*10.50720794170734)) +
		y*(-1351.605895580406 + y*(1097.1125373015109 +
		y*(-433.20648175062206 + 63.905091254154904*y))));

	return ((g03 + g08)*0.000625);
}

/**
==========================================================================
method: gsw_t_from_ct(sa,ct,p)
==========================================================================

 Calculates in-situ temperature from Conservative Temperature of seawater

 sa      : Absolute Salinity                              [g/kg]
 ct      : Conservative Temperature                       [deg C]

 gsw_t_from_ct : in-situ temperature                      [deg C]
*/
double TeosBase::gsw_t_from_ct(double sa, double ct, double p)
{
	double	pt0, p0=0.0;

	pt0	= gsw_pt_from_ct(sa,ct);

	return (gsw_pt_from_t(sa,pt0,p0,p));
}
/**
==========================================================================
method: gsw_util_indxPref(n,z)
==========================================================================

  calculates k in simulated array

  n     :  length of the array
  z     :  value to be indexed

  K      : index K - if X(K) <= Z < X(K+1), or
  N-1     		    - if Z = X(N)

*/
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
/**
==========================================================================
method: gsw_util_indx(x,n,z,k)
==========================================================================

  Finds the index of the value in a monotonically increasing array

  x	 :  array of monotonically increasing values
  n     :  length of the array
  z     :  value to be indexed

  K      : index K - if X(K) <= Z < X(K+1), or
  N-1     		    - if Z = X(N)

*/
int TeosBase::gsw_util_indx(double *x, int n, double z)
{
	int	k, ku, kl, km;

	if (z > x[0] && z < x[n-1]) {
	    kl	= 0;
	    ku	= n-1;
	    while (ku-kl > 1) {
		km	= (ku+kl)>>1;
		if (z > x[km])
		    kl	= km;
		else
		    ku	= km;
	    }
	    k	= kl;
	    if (z == x[k+1])
		k++;
	} else if (z <= x[0])
	    k	= 0;
	else
	    k	= n-2;

	return (k);
}

/**
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

 specvol: specific volume                                 [m^3/kg]
*/
double TeosBase::gsw_specvol(double sa, double ct, double p)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	/** for GSW_SPECVOL_COEFFICIENTS use gsvco */

	double	xs, ys, z, value;

	xs	= sqrt(gtc.gsw_sfac*sa + gtc.offset);
	ys	= ct*0.025;
	z	= p*1e-4;

	value = gsvco.v000
    + xs*(gsvco.v010 + xs*(gsvco.v020 + xs*(gsvco.v030
    + xs*(gsvco.v040 + xs*(gsvco.v050
    + gsvco.v060*xs))))) + ys*(gsvco.v100 + xs*(gsvco.v110
    + xs*(gsvco.v120 + xs*(gsvco.v130 + xs*(gsvco.v140
    + gsvco.v150*xs)))) + ys*(gsvco.v200 + xs*(gsvco.v210
    + xs*(gsvco.v220 + xs*(gsvco.v230 + gsvco.v240*xs)))
    + ys*(gsvco.v300 + xs*(gsvco.v310 + xs*(gsvco.v320
    + gsvco.v330*xs)) + ys*(gsvco.v400 + xs*(gsvco.v410
    + gsvco.v420*xs) + ys*(gsvco.v500 + gsvco.v510*xs
    + gsvco.v600*ys))))) + z*(gsvco.v001 + xs*(gsvco.v011
    + xs*(gsvco.v021 + xs*(gsvco.v031 + xs*(gsvco.v041
    + gsvco.v051*xs)))) + ys*(gsvco.v101 + xs*(gsvco.v111
    + xs*(gsvco.v121 + xs*(gsvco.v131 + gsvco.v141*xs)))
    + ys*(gsvco.v201 + xs*(gsvco.v211 + xs*(gsvco.v221
    + gsvco.v231*xs)) + ys*(gsvco.v301 + xs*(gsvco.v311
    + gsvco.v321*xs) + ys*(gsvco.v401 + gsvco.v411*xs
    + gsvco.v501*ys)))) + z*(gsvco.v002 + xs*(gsvco.v012
    + xs*(gsvco.v022 + xs*(gsvco.v032 + gsvco.v042*xs)))
    + ys*(gsvco.v102 + xs*(gsvco.v112 + xs*(gsvco.v122
    + gsvco.v132*xs)) + ys*(gsvco.v202 + xs*(gsvco.v212
    + gsvco.v222*xs) + ys*(gsvco.v302 + gsvco.v312*xs
    + gsvco.v402*ys))) + z*(gsvco.v003 + xs*(gsvco.v013
    + gsvco.v023*xs) + ys*(gsvco.v103 + gsvco.v113*xs
    + gsvco.v203*ys) + z*(gsvco.v004 + gsvco.v014*xs
    + gsvco.v104*ys + z*(gsvco.v005 + gsvco.v006*z)))));

	return (value);
}

/**
==========================================================================
subroutine gsw_util_indx(x,n,z,k)
==========================================================================

!  Finds the index of the value in a monotonically increasing array
!
!  x     :  array of monotonically increasing values
!  n     :  length of the array
!  z     :  value to be indexed
!
!  K      : index K - if X(K) <= Z < X(K+1), or
!  N-1                      - if Z = X(N)
!
*/
/** missing 4 parameter version */



/**
--------------------------------------------------------------------------
 density and enthalpy, based on the 48-term expression for density
--------------------------------------------------------------------------

==========================================================================
method: gsw_rho(sa,ct,p)
==========================================================================

  Calculates in-situ density from Absolute Salinity and Conservative
  Temperature, using the computationally-efficient expression for
  specific volume in terms of SA, CT and p (Roquet et al., 2014).

 sa     : Absolute Salinity                               [g/kg]
 ct     : Conservative Temperature (ITS-90)               [deg C]
 p      : sea pressure                                    [dbar]
	   ( i.e. absolute pressure - 10.1325 dbar )

 rho    : in-situ density				   [kg/m]
*/
double TeosBase::gsw_rho(double sa, double ct, double p)
{
	return (1.0/gsw_specvol(sa,ct,p));
}
/**
==========================================================================
method: gsw_saar(p,lon,lat)
==========================================================================

 Calculates the Absolute Salinity Anomaly Ratio, SAAR.

 p      : sea pressure                                    [dbar]
 lon    : longitude                                       [deg E]
 lat    : latitude                                        [deg N]

 gsw_saar : Absolute Salinity Anomaly Ratio               [unitless]
*/
double TeosBase::gsw_saar(double p, double lon, double lat)
{
    /** for GSW_SAAR_DATA use gsaar */

        int     nx=gsw_nx, ny=gsw_ny, nz=gsw_nz;
        int     indx0, indy0, indz0, k, ndepth_index;
        double  saar[4], saar_old[4];
        double  sa_upper, sa_lower, dlong, dlat;
        double  r1, s1, t1, ndepth_max, return_value;

        return_value = cppGSW_INVALID_VALUE;

        if (isnan(lat) || isnan(lon) || isnan(p))
        {
            return (return_value);
        }

        if (lat  <  -86.0  ||  lat  >  90.0)
        {
            return (return_value);
        }

        if (lon  <  0.0)
        {
            lon += 360.0;
        }

        dlong   = pSaarData->get_longs_ref(1)-pSaarData->get_longs_ref(0);
        dlat    = pSaarData->get_lats_ref(1)-pSaarData->get_lats_ref(0);

        indx0   = floor(0 + (nx-1)*(lon-pSaarData->get_longs_ref(0))/
                        (pSaarData->get_longs_ref(nx-1)-pSaarData->get_longs_ref(0)));

        if (indx0 == nx-1)
        {
            indx0 = nx-2;
        }

        indy0 = floor(0 + (ny-1)*(lat-pSaarData->get_lats_ref(0))/
                        (pSaarData->get_lats_ref(ny-1)-pSaarData->get_lats_ref(0)));

        if(indy0 == ny-1)
        {
            indy0 = ny-2;
        }
/**
! Look for the maximum valid "ndepth_ref" value around our point.
! Note: invalid "ndepth_ref" values are NaNs (a hangover from the codes
! Matlab origins), but we have replaced the NaNs with a value of "9e90",
! hence we need an additional upper-limit check in the code below so they
! will not be recognised as valid values.
*/
        ndepth_max      = -1.0;

        for (k=0; k < 4; k++)
        {
            ndepth_index = indy0+gsaar.delj[k]+(indx0+gsaar.deli[k])*ny;

            if (pSaarData->get_ndepth_ref(ndepth_index) > 0.0 &&
                pSaarData->get_ndepth_ref(ndepth_index) < 1e90)
            {
                ndepth_max = gsw_max(ndepth_max, pSaarData->get_ndepth_ref(ndepth_index));
            }
        }
/**
! If we are a long way from the ocean then there will be no valid "ndepth_ref"
! values near the point (ie. surrounded by NaNs) - so just return SAAR = 0.0
*/
        if (ndepth_max == -1.0)
        {
            return (0.0);
        }

        if (p > pSaarData->get_p_ref((int)(ndepth_max)-1))
        {
            p = pSaarData->get_p_ref((int)(ndepth_max)-1);
        }

        indz0   = gsw_util_indxPref(nz,p);

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

        if (fabs(sa_lower) >=  cppGSW_ERROR_LIMIT)
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

/**
==========================================================================
method: gsw_sa_from_sp(sp,p,lon,lat)
==========================================================================

 Calculates Absolute Salinity, SA, from Practical Salinity, SP

 sp     : Practical Salinity                              [unitless]
 p      : sea pressure                                    [dbar]
 lon	 : longitude                                       [DEG E]
 lat    : latitude                                        [DEG N]

 gsw_sa_from_sp   : Absolute Salinity                     [g/kg]
*/
double TeosBase::gsw_sa_from_sp(double sp, double p, double lon, double lat)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	double	saar, gsw_sa_baltic;

	gsw_sa_baltic	= gsw_sa_from_sp_baltic(sp,lon,lat);

	if (gsw_sa_baltic < cppGSW_ERROR_LIMIT)
	    return (gsw_sa_baltic);

	saar = gsw_saar(p,lon,lat);
	if (saar == cppGSW_INVALID_VALUE)
	    return (saar);

	return (gtc.gsw_ups*sp*(1.e0 + saar));
}
/**
==========================================================================
method: gsw_sa_from_sp_baltic(sp,lon,lat)
==========================================================================

 For the Baltic Sea, calculates Absolute Salinity with a value
 computed analytically from Practical Salinity

 sp     : Practical Salinity                              [unitless]
 lon    : longitude                                       [deg E]
 lat    : latitude                                        [deg N]

 sa_from_sp_baltic : Absolute Salinity                    [g/kg]
*/
double TeosBase::gsw_sa_from_sp_baltic(double sp, double lon, double lat)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */
    /**  for GSW_BALTIC_DATA use gbalt */

	double	xx_left, xx_right, return_value;

	if (gbalt.xb_left[1] < lon  && lon < gbalt.xb_right[0]
            && gbalt.yb_left[0] < lat  && lat < gbalt.yb_left[2])
    {
	    xx_left	= gsw_util_xinterp1(gbalt.yb_left, gbalt.xb_left, 3, lat);

	    xx_right	= gsw_util_xinterp1(gbalt.yb_right, gbalt.xb_right, 2, lat);

	    if (xx_left <= lon  && lon <= xx_right)
            return_value	=((gtc.gsw_sso - 0.087)/35.0)*sp + 0.087;
	    else
            return_value	= cppGSW_INVALID_VALUE;
	} else
	    return_value	= cppGSW_INVALID_VALUE;

	return (return_value);
}

double TeosBase::gsw_sum(double *x, int n)
{
    int	i;
    double	val;

    for (val=0.0, i=0; i<n; val += x[i], i++);

    return (val);
}

/**
==========================================================================
method: gsw_util_xinterp1(x,y,n,x0)
==========================================================================

 Linearly interpolate a real array

 x      : y array (Must be monotonic)
 y      : y array
 n      : length of X and Y arrays
 x0     : value to be interpolated

 gsw_xinterp1 : Linearly interpolated value
*/
double TeosBase::gsw_util_xinterp1(double *x, double *y, int n, double x0)
{
	int	k;
	double	r;

	k	= gsw_util_indx(x,n,x0);
	r	= (x0-x[k])/(x[k+1]-x[k]);

	return (y[k] + r*(y[k+1]-y[k]));
}

/**
==========================================================================
method gsw_enthalpy_sso_0(p)
==========================================================================

!  This method calculates enthalpy at the Standard Ocean Salinity, SSO,
!  and at a Conservative Temperature of zero degrees C, as a method of
!  pressure, p, in dbar, using a streamlined version of the
!  computationally-efficient expression for specific volume, that is, a
!  streamlined version of the code "gsw_enthalpy(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
!
! enthalpy_sso_0 : enthalpy(sso,0,p)
*/
double TeosBase::gsw_enthalpy_sso_0(double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

    double  dynamic_enthalpy_sso_0_p, z;

    z = p*1.0e-4;

    dynamic_enthalpy_sso_0_p =
            z*( 9.726613854843870e-4 + z*(-2.252956605630465e-5
            + z*( 2.376909655387404e-6 + z*(-1.664294869986011e-7
            + z*(-5.988108894465758e-9 + z*(gsvco.h006 + gsvco.h007*z))))));

    return (dynamic_enthalpy_sso_0_p*gtc.db2pa*1.0e4);
}

/**
==========================================================================
method gsw_entropy_part_zerop(sa,pt0)
==========================================================================

! entropy part evaluated at the sea surface
!
! sa     : Absolute Salinity                               [g/kg]
! pt0    : insitu temperature                              [deg C]
!
! entropy_part_zerop : entropy part at the sea surface
*/
double TeosBase::gsw_entropy_part_zerop(double sa, double pt0)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  x2, x, y, g03, g08;

    x2      = gtc.gsw_sfac*sa;
    x       = sqrt(x2);
    y       = pt0*0.025;

    g03     = y*(-24715.571866078 + y*(2210.2236124548363 +
                y*(-592.743745734632 + y*(290.12956292128547 +
                y*(-113.90630790850321 + y*21.35571525415769)))));

    g08     = x2*(x*(x*(y*(-137.1145018408982 + y*(148.10030845687618 +
                y*(-68.5590309679152 + 12.4848504784754*y)))) +
                y*(-86.1329351956084 + y*(-30.0682112585625 +
                   y*3.50240264723578))) +
                y*(1760.062705994408 + y*(-675.802947790203 +
                y*(365.7041791005036 + y*(-108.30162043765552 +
                   12.78101825083098*y)))));

    return (-(g03 + g08)*0.025);
}

/**
==========================================================================
method: gsw_z_from_p(p,lat,geo_strf_dyn_height,sea_surface_geopotential)
==========================================================================

 Calculates the height z from pressure p

 p      : sea pressure                                    [dbar]
 lat    : latitude                                        [deg]

 gsw_z_from_p : height                                    [m]
*/
double TeosBase::gsw_z_from_p(double p, double lat, double geo_strf_dyn_height,
                              double sea_surface_geopotential)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

	double	x, sin2, b, c, a;

	x	= sin(lat*gtc.deg2rad);
	sin2	= x*x;
	b	= 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);
	a	= -0.5*gtc.gamma*b;
	c	= gsw_enthalpy_sso_0(p)
	      - (geo_strf_dyn_height + sea_surface_geopotential);

	return (-2.0*c/(b + sqrt(b*b - 4.0*a*c)));
}

/**
==========================================================================
method: gsw_p_from_z(z,lat, double geo_strf_dyn_height,
		               double sea_surface_geopotential)
==========================================================================

 Calculates the pressure p from height z

 z      : height                                          [m]
 lat    : latitude                                        [deg]

 gsw_p_from_z : pressure                                  [dbar]
*/
double TeosBase::gsw_p_from_z(double z, double lat,
                              double geo_strf_dyn_height,
                              double sea_surface_geopotential)
{
    /**  for GSW_TEOS10_CONSTANTS use gtc */

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
/**
==========================================================================
method: gsw_specvol_sso_0(p)
==========================================================================

  This method: calculates specific volume at the Standard Ocean Salinty,
  SSO, and at a Conservative Temperature of zero degrees C, as a function
  of pressure, p, in dbar, using a streamlined version of the CT version
  of specific volume, that is, a streamlined version of the code
  "gsw_specvol(SA,CT,p)".

 p      : sea pressure                                    [dbar]

 specvol_sso_0 : specvol(sso,0,p)                         [m  kg  ]
*/
double TeosBase::gsw_specvol_sso_0(double p)
{
    /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

	double	z, return_value;

	z = p*1.0e-4;

	return_value = 9.726613854843870e-04 + z*(-4.505913211160929e-05
		+ z*(7.130728965927127e-06 + z*(-6.657179479768312e-07
		+ z*(-2.994054447232880e-08 + z*(gsvco.v005 + gsvco.v006*z)))));

	return (return_value);

}

/**
==========================================================================
method: gsw_ct_maxdensity (sa, p)
==========================================================================
!
!  Calculates the Conservative Temperature of maximum density of seawater.
!  This method returns the Conservative temperature at which the density
!  of seawater is a maximum, at given Absolute Salinity, SA, and sea
!  pressure, p (in dbar).
!
!  SA =  Absolute Salinity                                         [ g/kg ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar )
!
!  CT_maxdensity  =  Conservative Temperature at which            [ deg C ]
!                    the density of seawater is a maximum for
!                    given Absolute Salinity and pressure.
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_ct_maxdensity(double sa, double p)
{
    int     number_of_iterations;
    double  alpha, ct, ct_mean, ct_old, dalpha_dct,
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
        ! After three iterations of this modified Newton-Raphson (McDougall and
        ! Wotherspoon, 2012) iteration, the error in CT_maxdensity is typically
        ! no larger than 1x10^-15 degress C.
        */

    return (ct);
}
/**
==========================================================================
method gsw_grav(lat,p)
==========================================================================

! Calculates acceleration due to gravity as a method of latitude and as
!  a method of pressure in the ocean.
!
! lat  =  latitude in decimal degress north                [ -90 ... +90 ]
! p    =  sea pressure                                     [ dbar ]
!
! grav : grav  =  gravitational acceleration               [ m s^-2 ]
*/
double TeosBase::gsw_grav(double lat, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  x, sin2, gs, z;

    x       = sin(lat*gtc.deg2rad);  /** convert to radians */
    sin2    = x*x;
    gs      = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);

    z       = gsw_z_from_p(p,lat,0.0,0.0);

    return (gs*(1.0 - gtc.gamma*z));    /** z is the height corresponding to p.
                                           Note. In the ocean z is negative. */
}
/**
--------------------------------------------------------------------------
! Practical Salinity (SP), PSS-78
!--------------------------------------------------------------------------

==========================================================================
method gsw_sp_from_c(c,t,p)
==========================================================================

!  Calculates Practical Salinity, SP, from conductivity, C, primarily using
!  the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical
!  Salinity is only valid in the range 2 < SP < 42.  If the PSS-78
!  algorithm produces a Practical Salinity that is less than 2 then the
!  Practical Salinity is recalculated with a modified form of the Hill et
!  al. (1986) formula.  The modification of the Hill et al. (1986)
!  expression is to ensure that it is exactly consistent with PSS-78
!  at SP = 2.  Note that the input values of conductivity need to be in
!  units of mS/cm (not S/m).
!
! c      : conductivity                                     [ mS/cm ]
! t      : in-situ temperature [ITS-90]                     [deg C]
! p      : sea pressure                                     [dbar]
!
! sp     : Practical Salinity                               [unitless]
*/
double TeosBase::gsw_sp_from_c(double c, double t, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_SP_COEFFICIENTS use gspc */

    double  sp, t68, ft68, r, rt_lc, rp, rt, rtx,
            hill_ratio, x, sqrty, part1, part2,
            sp_hill_raw;

    t68     = t*1.00024e0;
    ft68    = (t68 - 15e0)/(1e0 + gspc.k*(t68 - 15e0));

    /**
     ! The dimensionless conductivity ratio, R, is the conductivity input, C,
     ! divided by the present estimate of C(SP=35, t_68=15, p=0) which is
     ! 42.9140 mS/cm (=4.29140 S/m), (Culkin and Smith, 1980).
    */

    r = c/gtc.gsw_c3515;        /** 0.023302418791070513 = 1./42.9140 */

    /**rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.*/
    rt_lc   = gspc.c0 + (gspc.c1 + (gspc.c2 + (gspc.c3 + gspc.c4*t68)*t68)*t68)*t68;
    rp      = 1e0 + (p*(gspc.e1 + gspc.e2*p + gspc.e3*p*p))/(1e0 + gspc.d1*t68 + gspc.d2*t68*t68 +
                  (gspc.d3 + gspc.d4*t68)*r);
    rt      = r/(rp*rt_lc);

    if (rt < 0.0)
    {
        return (cppGSW_INVALID_VALUE);
    }

    rtx     = sqrt(rt);

    sp      = gspc.a0 + (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4 +
            gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx +
                  ft68*(gspc.b0 + (gspc.b1 + (gspc.b2 + (gspc.b3 +
                        (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);
    /**
     ! The following section of the code is designed for SP < 2 based on the
     ! Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
     ! exactly equal to the PSS-78 algorithm at SP = 2.
    */

    if (sp < 2)
    {
        hill_ratio  = gsw_hill_ratio_at_sp2(t);
        x           = 400e0*rt;
        sqrty       = 10e0*rtx;
        part1       = 1e0 + x*(1.5e0 + x);
        part2       = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
        sp_hill_raw = sp - gspc.a0/part1 - gspc.b0*ft68/part2;
        sp          = hill_ratio*sp_hill_raw;
    }

    /** This line ensures that SP is non-negative. */
    if (sp < 0.0)
    {
        sp  = cppGSW_INVALID_VALUE;
    }

    return (sp);
}
/**
==========================================================================
method  gsw_hill_ratio_at_sp2(t)
==========================================================================

!  Calculates the Hill ratio, which is the adjustment needed to apply for
!  Practical Salinities smaller than 2.  This ratio is defined at a
!  Practical Salinity = 2 and in-situ temperature, t using PSS-78. The Hill
!  ratio is the ratio of 2 to the output of the Hill et al. (1986) formula
!  for Practical Salinity at the conductivity ratio, Rt, at which Practical
!  Salinity on the PSS-78 scale is exactly 2.
!
!  t                 : in-situ temperature (ITS-90)              [deg C]
!  hill_ratio_at_sp2 : Hill ratio                                [dimensionless]
*/
double TeosBase::gsw_hill_ratio_at_sp2(double t)
{
    /** for GSW_SP_COEFFICIENTS use gspc */

    double  g0 = 2.641463563366498e-1, g1 = 2.007883247811176e-4,
            g2 = -4.107694432853053e-6, g3 = 8.401670882091225e-8,
            g4 = -1.711392021989210e-9, g5 = 3.374193893377380e-11,
            g6 = -5.923731174730784e-13, g7 = 8.057771569962299e-15,
            g8 = -7.054313817447962e-17, g9 = 2.859992717347235e-19,
            sp2 = 2.0;

    double  t68, ft68, rtx0, dsp_drtx, sp_est, rtx, rtxm, x, part1, part2;
    double  sqrty, sp_hill_raw_at_sp2;

    t68     = t*1.00024;
    ft68    = (t68 - 15.0)/(1.0 + gspc.k*(t68 - 15.0));

    /*!------------------------------------------------------------------------
    **! Find the initial estimates of Rtx (Rtx0) and of the derivative dSP_dRtx
    **! at SP = 2.
    **!------------------------------------------------------------------------
    */
    rtx0    = g0 + t68*(g1 + t68*(g2 + t68*(g3 + t68*(g4 + t68*(g5
                + t68*(g6 + t68*(g7 + t68*(g8 + t68*g9))))))));

    dsp_drtx= gspc.a1 + (2*gspc.a2 + (3*gspc.a3 +
        (4*gspc.a4 + 5*gspc.a5*rtx0)*rtx0)*rtx0)*rtx0 +
            ft68*(gspc.b1 + (2*gspc.b2 + (3*gspc.b3 +
                (4*gspc.b4 + 5*gspc.b5*rtx0)*rtx0)*rtx0)*rtx0);

    /*!-------------------------------------------------------------------------
    **! Begin a single modified Newton-Raphson iteration to find Rt at SP = 2.
    **!-------------------------------------------------------------------------
    */
    sp_est  = gspc.a0 + (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4 + gspc.a5*rtx0)*rtx0)*rtx0)*rtx0)*rtx0
                + ft68*(gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3 + (gspc.b4 + gspc.b5*rtx0)*rtx0)*rtx0)*
                  rtx0)*rtx0);

    rtx     = rtx0 - (sp_est - sp2)/dsp_drtx;

    rtxm    = 0.5*(rtx + rtx0);

    dsp_drtx= gspc.a1 + (2*gspc.a2 + (3*gspc.a3 + (4*gspc.a4 + 5*gspc.a5*rtxm)*rtxm)*rtxm)*rtxm
                + ft68*(gspc.b1 + (2*gspc.b2 + (3*gspc.b3 + (4*gspc.b4 + 5*gspc.b5*rtxm)*
                                                rtxm)*rtxm)*rtxm);
    rtx     = rtx0 - (sp_est - sp2)/dsp_drtx;
    /**
    **! This is the end of one full iteration of the modified Newton-Raphson
    **! iterative equation solver. The error in Rtx at this point is equivalent
    **! to an error in SP of 9e-16 psu.
    */

    x       = 400.0*rtx*rtx;
    sqrty   = 10.0*rtx;
    part1   = 1.0 + x*(1.5 + x);
    part2   = 1.0 + sqrty*(1.0 + sqrty*(1.0 + sqrty));

    sp_hill_raw_at_sp2 = sp2 - gspc.a0/part1 - gspc.b0*ft68/part2;

    return (2.0/sp_hill_raw_at_sp2);
}
/**
==========================================================================
elemental subroutine gsw_enthalpy_first_derivatives (sa, ct, p, h_sa, h_ct)
==========================================================================
!
!  Calculates the following two derivatives of specific enthalpy (h) of
!  seawater using the computationally-efficient expression for
!  specific volume in terms of SA, CT and p (Roquet et al., 2014).
!   (1) h_SA, the derivative with respect to Absolute Salinity at
!       constant CT and p, and
!   (2) h_CT, derivative with respect to CT at constant SA and p.
!  Note that h_P is specific volume (1/rho) it can be caclulated by calling
!  gsw_specvol(SA,CT,p).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_SA  =  The first derivative of specific enthalpy with respect to
!           Absolute Salinity at constant CT and p.
!                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
!  h_CT  =  The first derivative of specific enthalpy with respect to
!           CT at constant SA and p.                           [ J/(kg K) ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_enthalpy_first_derivatives(double sa, double ct, double p,
                                                double *h_sa, double *h_ct)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

    double  dynamic_h_ct_part, dynamic_h_sa_part, xs, ys, z;

    xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
    ys = ct*0.025;
    z = p*1e-4;

    if (h_sa != NULL)
    {
        dynamic_h_sa_part =
            z*(gsvco.h101 + xs*(2.0*gsvco.h201 + xs*(3.0*gsvco.h301
            + xs*(4.0*gsvco.h401 + xs*(5.0*gsvco.h501 + 6.0*gsvco.h601*xs)))) + ys*(gsvco.h111
            + xs*(2.0*gsvco.h211 + xs*(3.0*gsvco.h311 + xs*(4.0*gsvco.h411
            + 5.0*gsvco.h511*xs))) + ys*(gsvco.h121 + xs*(2.0*gsvco.h221 + xs*(3.0*gsvco.h321
            + 4.0*gsvco.h421*xs)) + ys*(gsvco.h131 + xs*(2.0*gsvco.h231 + 3.0*gsvco.h331*xs)
            + ys*(gsvco.h141 + 2.0*gsvco.h241*xs + gsvco.h151*ys)))) + z*(gsvco.h102
            + xs*(2.0*gsvco.h202 + xs*(3.0*gsvco.h302 + xs*(4.0*gsvco.h402
            + 5.0*gsvco.h502*xs))) + ys*(gsvco.h112 + xs*(2.0*gsvco.h212 + xs*(3.0*gsvco.h312
            + 4.0*gsvco.h412*xs)) + ys*(gsvco.h122 + xs*(2.0*gsvco.h222 + 3.0*gsvco.h322*xs)
            + ys*(gsvco.h132 + 2.0*gsvco.h232*xs + gsvco.h142*ys ))) + z*(gsvco.h103 + xs*(2.0*gsvco.h203
            + xs*(3.0*gsvco.h303 + 4.0*gsvco.h403*xs)) + ys*(gsvco.h113 + xs*(2.0*gsvco.h213
            + 3.0*gsvco.h313*xs) + ys*(gsvco.h123 + 2.0*gsvco.h223*xs + gsvco.h133*ys))
            + z*(gsvco.h104 + 2.0*gsvco.h204*xs + gsvco.h114*ys + gsvco.h105*z))));

      *h_sa = 1e8*0.5*gtc.gsw_sfac*dynamic_h_sa_part/xs;

    }

    if (h_ct != NULL)
    {
        dynamic_h_ct_part =
            z*(gsvco.h011 + xs*(gsvco.h111
            + xs*(gsvco.h211 + xs*(gsvco.h311 + xs*(gsvco.h411
            + gsvco.h511*xs)))) + ys*(2.0*(gsvco.h021 + xs*(gsvco.h121 + xs*(gsvco.h221 + xs*(gsvco.h321
            + gsvco.h421*xs)))) + ys*(3.0*(gsvco.h031 + xs*(gsvco.h131 + xs*(gsvco.h231 + gsvco.h331*xs)))
            + ys*(4.0*(gsvco.h041 + xs*(gsvco.h141 + gsvco.h241*xs)) + ys*(5.0*(gsvco.h051
            + gsvco.h151*xs) + 6.0*gsvco.h061*ys)))) + z*(gsvco.h012 + xs*(gsvco.h112 + xs*(gsvco.h212
            + xs*(gsvco.h312 + gsvco.h412*xs))) + ys*(2.0*(gsvco.h022 + xs*(gsvco.h122 + xs*(gsvco.h222
            + gsvco.h322*xs))) + ys*(3.0*(gsvco.h032 + xs*(gsvco.h132 + gsvco.h232*xs))
            + ys*(4.0*(gsvco.h042 + gsvco.h142*xs) + 5.0*gsvco.h052*ys))) + z*(gsvco.h013
            + xs*(gsvco.h113 + xs*(gsvco.h213 + gsvco.h313*xs)) + ys*(2.0*(gsvco.h023 + xs*(gsvco.h123
            + gsvco.h223*xs)) + ys*(3.0*(gsvco.h033 + gsvco.h133*xs) + 4.0*gsvco.h043*ys))
            + z*(gsvco.h014 + gsvco.h114*xs + 2.0*gsvco.h024*ys + gsvco.h015*z ))));

        *h_ct = gtc.gsw_cp0 + 1e8*0.025*dynamic_h_ct_part;
    }
}
/**
==========================================================================
method gsw_o2sol(sa, ct, p, lon, lat)
==========================================================================
!
! Calculates the oxygen concentration expected at equilibrium with air at
! an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
! saturated water vapor.  This method uses the solubility coefficients
! derived from the data of Benson and Krause (1984), as fitted by Garcia
! and Gordon (1992, 1993).
!
! Note that this algorithm has not been approved by IOC and is not work
! from SCOR/IAPSO Working Group 127. It is included in the GSW
! Oceanographic Toolbox as it seems to be oceanographic best practice.
!
! SA  :  Absolute Salinity of seawater                           [ g/kg ]
! CT  :  Conservative Temperature of seawater (ITS-90)           [ deg C ]
! p   :  sea pressure at which the melting occurs                [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
! lat : latitude                                                 [deg]
! lon : longitude                                                [deg]
!
! gsw_o2sol : olubility of oxygen in micro-moles per kg          [umol/kg]
*/
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
/**
==========================================================================
method gsw_o2sol_sp_pt(sp, pt)
==========================================================================
!
! Calculates the oxygen concentration expected at equilibrium with air at
! an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
! saturated water vapor.  This method uses the solubility coefficients
! derived from the data of Benson and Krause (1984), as fitted by Garcia
! and Gordon (1992, 1993).
!
! Note that this algorithm has not been approved by IOC and is not work
! from SCOR/IAPSO Working Group 127. It is included in the GSW
! Oceanographic Toolbox as it seems to be oceanographic best practice.
!
! SP  :  Practical Salinity  (PSS-78)                         [ unitless ]
! pt  :  potential temperature (ITS-90) referenced               [ dbar ]
!         to one standard atmosphere (0 dbar).
!
! gsw_o2sol_sp_pt : olubility of oxygen in micro-moles per kg     [umol/kg]
*/
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
/**
==========================================================================
method gsw_sp_from_sa_baltic(sa,lon,lat)
==========================================================================

! For the Baltic Sea, calculates Practical Salinity with a value
! computed analytically from Absolute Salinity
!
! sa     : Absolute Salinity                               [g/kg]
! lon    : longitude                                       [deg E]
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]
*/
double TeosBase::gsw_sp_from_sa_baltic(double sa, double lon, double lat)
{
      /** for GSW_TEOS10_CONSTANTS use gtc */
      /** for GSW_BALTIC_DATA use gbalt */

      double  xx_left, xx_right, return_value;

      if (gbalt.xb_left[1] < lon  && lon < gbalt.xb_right[0]  && gbalt.yb_left[0] < lat  &&
            lat < gbalt.yb_left[2])
      {
         xx_left     = gsw_util_xinterp1(gbalt.yb_left, gbalt.xb_left, 3, lat);

         xx_right    = gsw_util_xinterp1(gbalt.yb_right, gbalt.xb_right, 2, lat);

         if (xx_left <= lon  && lon <= xx_right)
            return_value    = (35.0/(gtc.gsw_sso - 0.087))*(sa - 0.087);
         else
            return_value   = cppGSW_INVALID_VALUE;
      } else
         return_value        = cppGSW_INVALID_VALUE;

      return (return_value);
}
/**
==========================================================================
method gsw_sp_from_sa(sa,p,lon,lat)
==========================================================================

! Calculates Practical salinity, sp, from Absolute salinity, sa
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! lon    : longitude                                       [DEG E]
! lat    : latitude                                        [DEG N]
!
! gsw_sp_from_sa      : Practical Salinity                 [unitless]
*/
double TeosBase::gsw_sp_from_sa(double sa, double p, double lon, double lat)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */

   double  saar, gsw_sp_baltic;

   gsw_sp_baltic   = gsw_sp_from_sa_baltic(sa,lon,lat);
   if (gsw_sp_baltic < cppGSW_ERROR_LIMIT)
      return (gsw_sp_baltic);

   saar = gsw_saar(p,lon,lat);
   if (saar == cppGSW_INVALID_VALUE)
      return (saar);

   return ((sa/gtc.gsw_ups)/(1e0 + saar));
}
/**
==========================================================================
method: gsw_t_deriv_chem_potential_water_t_exact (sa, t, p)
==========================================================================
!
!  Calculates the temperature derivative of the chemical potential of water
!  in seawater so that it is valid at exactly SA = 0.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  chem_potential_water_dt  =  temperature derivative of the chemical
!                           potential of water in seawater  [ J g^-1 K^-1 ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_t_deriv_chem_potential_water_t_exact(double sa, double t, double p)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */

   double  g03_t, g08_sa_t, x, x2, y, z, g08_t, kg2g = 1e-3;
/**
 Note. The kg2g, a factor of 1e-3, is needed to convert the output of this
! method into units of J/g. See section (2.9) of the TEOS-10 Manual.
*/

   x2 = gtc.gsw_sfac*sa;
   x = sqrt(x2);
   y = t*0.025;
   z = p*gtc.rec_db2pa;
   /** the input pressure (p) is sea pressure in units of dbar. */

   g03_t = 5.90578347909402 + z*(-270.983805184062 +
   z*(776.153611613101 + z*(-196.51255088122 + (28.9796526294175 -
   2.13290083518327*z)*z))) +
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

   g08_t = x2*(168.072408311545 +
        x*(-493.407510141682 + x*(543.835333000098 +
        x*(-196.028306689776 + 36.7571622995805*x) +
        y*(-137.1145018408982 + y*(148.10030845687618 +
        y*(-68.5590309679152 + 12.4848504784754*y))) -
        22.6683558512829*z) + z*(-175.292041186547 +
        (83.1923927801819 - 29.483064349429*z)*z) +
        y*(-86.1329351956084 + z*(766.116132004952 +
        z*(-108.3834525034224 + 51.2796974779828*z)) +
        y*(-30.0682112585625 - 1380.9597954037708*z +
        y*(3.50240264723578 + 938.26075044542*z)))));

        g08_sa_t = 1187.3715515697959 +
        x*(-1480.222530425046 + x*(2175.341332000392 +
        x*(-980.14153344888 + 220.542973797483*x) +
        y*(-548.4580073635929 + y*(592.4012338275047 +
        y*(-274.2361238716608 + 49.9394019139016*y))) -
        90.6734234051316*z) + z*(-525.876123559641 +
        (249.57717834054571 - 88.449193048287*z)*z) +
        y*(-258.3988055868252 + z*(2298.348396014856 +
        z*(-325.1503575102672 + 153.8390924339484*z)) +
        y*(-90.2046337756875 - 4142.8793862113125*z +
        y*(10.50720794170734 + 2814.78225133626*z))));

   return (kg2g*((g03_t + g08_t)*0.025 - 0.5*gtc.gsw_sfac*0.025*sa*g08_sa_t));
}
/**
==========================================================================
method gsw_specvol_t_exact(sa,t,p)
==========================================================================

! Calculates the specific volume of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! specvol_t_exact : specific volume                        [kg/m^3]
*/
double TeosBase::gsw_specvol_t_exact(double sa, double t, double p)
{
        int     n0=0, n1=1;

        return (gsw_gibbs(n0,n0,n1,sa,t,p));
}
/**
==========================================================================
method gsw_sr_from_sp(sp)
==========================================================================

! Calculates Reference Salinity, SR, from Practical Salinity, SP.
!
! sp     : Practical Salinity                              [unitless]
!
! gsw_sr_from_sp : Reference Salinity                      [g/kg]
*/
double TeosBase::gsw_sr_from_sp(double sp)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        return (sp*gtc.gsw_ups);
}
/**
==========================================================================
method gsw_enthalpy(sa,ct,p)
==========================================================================
!  Calculates specific enthalpy of seawater using the computationally-
!  efficient expression for specific volume in terms of SA, CT and p
!  (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
! enthalpy  :  specific enthalpy of seawater               [J/kg]
*/
double TeosBase::gsw_enthalpy(double sa, double ct, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    return (gtc.gsw_cp0*ct + gsw_dynamic_enthalpy(sa,ct,p));
}
/**
==========================================================================
method gsw_kappa(sa,ct,p)
==========================================================================

!  Calculates isentropic compressibility of seawater.  This method
!  has inputs of Absolute Salinity and Conservative Temperature.  This
!  method uses the computationally-efficient expression for
!  specific volume in terms of SA, CT and p (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!
! kappa  :  isentropic compressibility                     [1.0/Pa]
*/
double TeosBase::gsw_kappa(double sa, double ct, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

    double  v, v_p, xs, ys, z;

    xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
    ys      = ct*0.025;
    z       = p*1e-4;

    v   = gsvco.v000
        + xs*(gsvco.v010 + xs*(gsvco.v020 + xs*(gsvco.v030 + xs*(gsvco.v040 + xs*(gsvco.v050
        + gsvco.v060*xs))))) + ys*(gsvco.v100 + xs*(gsvco.v110 + xs*(gsvco.v120 + xs*(gsvco.v130 + xs*(gsvco.v140
        + gsvco.v150*xs)))) + ys*(gsvco.v200 + xs*(gsvco.v210 + xs*(gsvco.v220 + xs*(gsvco.v230 + gsvco.v240*xs)))
        + ys*(gsvco.v300 + xs*(gsvco.v310 + xs*(gsvco.v320 + gsvco.v330*xs)) + ys*(gsvco.v400 + xs*(gsvco.v410
        + gsvco.v420*xs) + ys*(gsvco.v500 + gsvco.v510*xs + gsvco.v600*ys))))) + z*(gsvco.v001 + xs*(gsvco.v011
        + xs*(gsvco.v021 + xs*(gsvco.v031 + xs*(gsvco.v041 + gsvco.v051*xs)))) + ys*(gsvco.v101 + xs*(gsvco.v111
        + xs*(gsvco.v121 + xs*(gsvco.v131 + gsvco.v141*xs))) + ys*(gsvco.v201 + xs*(gsvco.v211 + xs*(gsvco.v221
        + gsvco.v231*xs)) + ys*(gsvco.v301 + xs*(gsvco.v311 + gsvco.v321*xs) + ys*(gsvco.v401 + gsvco.v411*xs
        + gsvco.v501*ys)))) + z*(gsvco.v002 + xs*(gsvco.v012 + xs*(gsvco.v022 + xs*(gsvco.v032 + gsvco.v042*xs)))
        + ys*(gsvco.v102 + xs*(gsvco.v112 + xs*(gsvco.v122 + gsvco.v132*xs)) + ys*(gsvco.v202 + xs*(gsvco.v212
        + gsvco.v222*xs) + ys*(gsvco.v302 + gsvco.v312*xs + gsvco.v402*ys))) + z*(gsvco.v003 + xs*(gsvco.v013
        + gsvco.v023*xs) + ys*(gsvco.v103 + gsvco.v113*xs + gsvco.v203*ys) + z*(gsvco.v004 + gsvco.v014*xs
        + gsvco.v104*ys + z*(gsvco.v005 + gsvco.v006*z)))));


    v_p = gsvco.c000
        + xs*(gsvco.c100 + xs*(gsvco.c200 + xs*(gsvco.c300 + xs*(gsvco.c400 + gsvco.c500*xs))))
        + ys*(gsvco.c010 + xs*(gsvco.c110 + xs*(gsvco.c210 + xs*(gsvco.c310 + gsvco.c410*xs))) + ys*(gsvco.c020
        + xs*(gsvco.c120 + xs*(gsvco.c220 + gsvco.c320*xs)) + ys*(gsvco.c030 + xs*(gsvco.c130 + gsvco.c230*xs)
        + ys*(gsvco.c040 + gsvco.c140*xs + gsvco.c050*ys)))) + z*(gsvco.c001 + xs*(gsvco.c101 + xs*(gsvco.c201
        + xs*(gsvco.c301 + gsvco.c401*xs))) + ys*(gsvco.c011 + xs*(gsvco.c111 + xs*(gsvco.c211 + gsvco.c311*xs))
        + ys*(gsvco.c021 + xs*(gsvco.c121 + gsvco.c221*xs) + ys*(gsvco.c031 + gsvco.c131*xs + gsvco.c041*ys)))
        + z*(gsvco.c002 + xs*(gsvco.c102 + gsvco.c202*xs) + ys*(gsvco.c012 + gsvco.c112*xs + gsvco.c022*ys)
        + z*(gsvco.c003 + gsvco.c103*xs + gsvco.c013*ys + z*(gsvco.c004 + gsvco.c005*z))));

    return (-1e-8*v_p/v);
}
/**
==========================================================================
method gsw_fdelta(p,lon,lat)
==========================================================================

! Calculates fdelta.
!
! p      : sea pressure                                    [dbar]
! lon    : longitude                                       [deg E]
! lat    : latitude                                        [deg N]
!
! gsw_fdelta : Absolute Salinty Anomaly                    [unitless]
*/
double TeosBase::gsw_fdelta(double p, double lon, double lat)
{
    double  sa, saar;

    saar    = gsw_saar(p,lon,lat);

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
/**
==========================================================================
method: gsw_ct_from_enthalpy_exact (sa, h, p)
==========================================================================
!
!  Calculates the Conservative Temperature of seawater, given the Absolute
!  Salinity, specific enthalpy, h, and pressure p.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  h   =  specific enthalpy                                        [ J/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar )
!
!  CT  =  Conservative Temperature ( ITS-90)                      [ deg C ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_ct_from_enthalpy_exact(double sa, double h, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  ct, ct_freezing, ct_mean, ct_old, f, h_freezing,
            h_ct, h_40, ct_40 = 40.0;

    ct_freezing = gsw_t_freezing(sa,p,0.0);

    h_freezing = gsw_enthalpy_ct_exact(sa,ct_freezing,p);

    if (h < (h_freezing - gtc.gsw_cp0))
    {
        /**
            ! The input, seawater enthalpy h, is less than the enthalpy at the
            ! freezing temperature, i.e. the water is frozen.
        */
        return (cppGSW_INVALID_VALUE);
    }

    h_40 = gsw_enthalpy_ct_exact(sa,ct_40,p);
    if (h > h_40)
    {
        /**
            ! The input seawater enthalpy is greater than the enthalpy
            ! when CT is 40C
        */
        return (cppGSW_INVALID_VALUE);
    }

    /** First guess of ct */
    ct = ct_freezing + (ct_40 - ct_freezing)*
            (h - h_freezing)/(h_40 - h_freezing);

    gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,NULL,&h_ct);
    /**
        !------------------------------------------------------
        ! Begin the modified Newton-Raphson iterative procedure
        !------------------------------------------------------
    */
    ct_old = ct;
    f = gsw_enthalpy_ct_exact(sa,ct_old,p) - h;
    ct = ct_old - f/h_ct;
    ct_mean = 0.5*(ct + ct_old);
    gsw_enthalpy_first_derivatives_ct_exact(sa,ct_mean,p,NULL,&h_ct);
    ct = ct_old - f/h_ct;
    /**
        ! After 1 iteration of this modified Newton-Raphson iteration,
        ! the error in CT is no larger than 5x10^-14 degrees C, which
        ! is machine precision for this calculation.
    */
    return (ct);
}
/**
==========================================================================
method: gsw_ct_freezing (sa, p, saturation_fraction)
==========================================================================
!
!  Calculates the Conservative Temperature at which seawater freezes.  The
!  Conservative Temperature freezing point is calculated from the exact
!  in-situ freezing temperature which is found by a modified Newton-Raphson
!  iteration (McDougall and Wotherspoon, 2013) of the equality of the
!  chemical potentials of water in seawater and in ice.
!
!  An alternative GSW method, gsw_CT_freezing_poly, it is based on a
!  computationally-efficient polynomial, and is accurate to within -5e-4 K
!  and 6e-4 K, when compared with this method.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  CT_freezing = Conservative Temperature at freezing of seawater [ deg C ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_ct_freezing(double sa, double p, double saturation_fraction)
{
    double  t_freezing;

    t_freezing = gsw_t_freezing(sa,p,saturation_fraction);

    return (gsw_ct_from_t(sa,t_freezing,p));
}
/**
==========================================================================
elemental subroutine gsw_enthalpy_first_derivatives_ct_exact (sa, ct, p, &
                                                              h_sa, h_ct)
==========================================================================
!
!  Calculates the following two derivatives of specific enthalpy (h)
!   (1) h_SA, the derivative with respect to Absolute Salinity at
!       constant CT and p, and
!   (2) h_CT, derivative with respect to CT at constant SA and p.
!  Note that h_P is specific volume (1/rho) it can be calulated by calling
!  gsw_specvol_CT_exact(SA,CT,p). This method uses the full Gibbs method.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_SA  =  The first derivative of specific enthalpy with respect to
!           Absolute Salinity at constant CT and p.
!                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
!  h_CT  =  The first derivative of specific enthalpy with respect to
!           CT at constant SA and p.                           [ J/(kg K) ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_enthalpy_first_derivatives_ct_exact(double sa, double ct, double p,
                double *h_sa, double *h_ct)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  g_sa_mod_pt, g_sa_mod_t, pt0, t, temp_ratio, x, y, y_pt, z;

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
/**
==========================================================================
method gsw_t_freezing(sa,p,saturation_fraction)
==========================================================================

! Calculates the in-situ temperature at which seawater freezes
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
!         ( i.e. absolute pressure - 10.1325 dbar )
! saturation_fraction : the saturation fraction of dissolved air
!                       in seawater
!
! t_freezing : in-situ temperature at which seawater freezes.[deg C]
*/
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
 df_dt here is the initial value of the derivative of the method f whose
! zero (f = 0) we are finding (see Eqn. (3.33.2) of IOC et al (2010)).
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
 =========================================================================
method: gsw_gibbs_ice (nt, np, t, p)
! =========================================================================
!
!  Ice specific Gibbs energy and derivatives up to order 2.
!
!  nt  =  order of t derivative                      [ integers 0, 1 or 2 ]
!  np  =  order of p derivative                      [ integers 0, 1 or 2 ]
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!
!  gibbs_ice = Specific Gibbs energy of ice or its derivatives.
!            The Gibbs energy (when nt = np = 0) has units of:     [ J/kg ]
!            The temperature derivatives are output in units of:
!                                                      [ (J/kg) (K)^(-nt) ]
!            The pressure derivatives are output in units of:
!                                                     [ (J/kg) (Pa)^(-np) ]
!            The mixed derivatives are output in units of:
!                                           [ (J/kg) (K)^(-nt) (Pa)^(-np) ]
!  Note. The derivatives are taken with respect to pressure in Pa, not
!    withstanding that the pressure input into this routine is in dbar.
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_gibbs_ice (int nt, int np, double t, double p)
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

        } else
           return (cppGSW_INVALID_VALUE);
}

/** supplement to gsw_chem_potential_water_t_exact (sa, t, p) */
double TeosBase::gsw_get_g03plusg08_h(double z, double y, double x2, double x)
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
/**
==========================================================================
method: gsw_chem_potential_water_t_exact (sa, t, p)
==========================================================================
!
!  Calculates the chemical potential of water in seawater.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  chem_potential_water_t_exact  =  chemical potential of water in seawater
!                                                                   [ J/g ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_chem_potential_water_t_exact(double sa, double t, double p)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */

   double  g03g_plus_g08g, g_sa_part, x, x2, y, z, kg2g = 1e-3, rval;

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
/**
==========================================================================
method: gsw_ct_from_enthalpy (sa, h, p)
==========================================================================
!
!  Calculates the Conservative Temperature of seawater, given the Absolute
!  Salinity, specific enthalpy, h, and pressure p.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  h   =  specific enthalpy                                        [ J/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar )
!
!  CT  =  Conservative Temperature ( ITS-90)                      [ deg C ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_ct_from_enthalpy(double sa, double h, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  ct, ct_freezing, ct_mean, ct_old, f, h_freezing,
            h_ct, h_40, ct_40 = 40.0;

    ct_freezing = gsw_ct_freezing(sa,p,0.0);

    h_freezing = gsw_enthalpy(sa,ct_freezing,p);

    if (h < (h_freezing - gtc.gsw_cp0))
    {
        /**
        ! The input, seawater enthalpy h, is less than the enthalpy at the
        ! freezing temperature, i.e. the water is frozen.
        */
        return (cppGSW_INVALID_VALUE);
    }

    h_40 = gsw_enthalpy(sa,ct_40,p);

    if (h > h_40)
    {
        /**
        ! The input seawater enthalpy is greater than the enthalpy
        ! when CT is 40C
        */
        return (cppGSW_INVALID_VALUE);
    }

    /** first guess of ct */
    ct = ct_freezing + (ct_40 - ct_freezing)*
            (h - h_freezing)/(h_40 - h_freezing);

    gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);

    /**
        !------------------------------------------------------
        ! Begin the modified Newton-Raphson iterative procedure
        !------------------------------------------------------
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
        ! After 1.5d0 iterations of this modified Newton-Raphson iteration,
        ! the error in CT is no larger than 4x10^-13 degrees C, which
        ! is machine precision for this calculation.
    */
    return (ct);
}
/**
==========================================================================
method gsw_beta(sa,ct,p)
==========================================================================

!  Calculates the saline (i.e. haline) contraction coefficient of seawater
!  at constant Conservative Temperature using the computationally-efficient
!  expression for specific volume in terms of SA, CT and p
!  (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
! beta   : saline contraction coefficient of seawater      [kg/g]
!          at constant Conservative Temperature
*/
double TeosBase::gsw_beta(double sa, double ct, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

    double  xs, ys, z, v_sa_part, rval;

    xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
    ys      = ct*0.025;
    z       = p*1e-4;

    v_sa_part = gsvco.b000
         + xs*(gsvco.b100 + xs*(gsvco.b200 + xs*(gsvco.b300 + xs*(gsvco.b400 + gsvco.b500*xs))))
         + ys*(gsvco.b010 + xs*(gsvco.b110 + xs*(gsvco.b210 + xs*(gsvco.b310 + gsvco.b410*xs)))
         + ys*(gsvco.b020 + xs*(gsvco.b120 + xs*(gsvco.b220 + gsvco.b320*xs)) + ys*(gsvco.b030
         + xs*(gsvco.b130 + gsvco.b230*xs) + ys*(gsvco.b040 + gsvco.b140*xs + gsvco.b050*ys))))
         + z*(gsvco.b001 + xs*(gsvco.b101 + xs*(gsvco.b201 + xs*(gsvco.b301 + gsvco.b401*xs)))
         + ys*(gsvco.b011 + xs*(gsvco.b111 + xs*(gsvco.b211 + gsvco.b311*xs)) + ys*(gsvco.b021
         + xs*(gsvco.b121 + gsvco.b221*xs) + ys*(gsvco.b031 + gsvco.b131*xs + gsvco.b041*ys)))
         + z*(gsvco.b002 + xs*(gsvco.b102 + xs*(gsvco.b202 + gsvco.b302*xs))+ ys*(gsvco.b012
         + xs*(gsvco.b112 + gsvco.b212*xs) + ys*(gsvco.b022 + gsvco.b122*xs + gsvco.b032*ys))
         + z*(gsvco.b003 +  gsvco.b103*xs + gsvco.b013*ys + gsvco.b004*z)));

   rval = (-v_sa_part*0.5*gtc.gsw_sfac/(gsw_specvol(sa,ct,p)*xs));

   return rval;
}

/** this method supplements gsw_dilution_coefficient_t_exact */
double TeosBase::gsw_g08_d(double x2, double z, double y, double x)
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
/**
==========================================================================
method: gsw_dilution_coefficient_t_exact (sa, t, p)
==========================================================================
!
!  Calculates the dilution coefficient of seawater.  The dilution
!  coefficient of seawater is defined as the Absolute Salinity times the
!  second derivative of the Gibbs method with respect to Absolute
!  Salinity, that is, SA.*g_SA_SA.
!
!  SA =  Absolute Salinity                                         [ g/kg ]
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar )
!
!  dilution_coefficient_t_exact  =  dilution coefficient   [ (J/kg)(kg/g) ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_dilution_coefficient_t_exact(double sa, double t, double p)
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
! Note that this method avoids the singularity that occurs at SA = 0 if
! the straightforward expression for the dilution coefficient of seawater,
! SA*g_SA_SA is simply evaluated as SA.*gsw_gibbs(2,0,0,SA,t,p).
*/
}
/**
==========================================================================
method gsw_internal_energy(sa,ct,p)
==========================================================================

!  Calculates internal energy of seawater.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!
! internal_energy  :  internal_energy of seawater          [J/kg]
*/
double TeosBase::gsw_internal_energy(double sa, double ct, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    return (gsw_enthalpy(sa,ct,p) - (gtc.gsw_p0 + gtc.db2pa*p)
                *gsw_specvol(sa,ct,p));
}
/**
==========================================================================
method gsw_kappa_t_exact(sa,t,p)
==========================================================================

! isentropic compressibility of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_kappa_t_exact : isentropic compressibility           [1/Pa]
*/
double TeosBase::gsw_kappa_t_exact(double sa, double t, double p)
{
    int     n0=0, n1=1, n2=2;
    double  g_tt, g_tp;

    g_tt    = gsw_gibbs(n0,n2,n0,sa,t,p);
    g_tp    = gsw_gibbs(n0,n1,n1,sa,t,p);

    return ((g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p)) /
                (gsw_gibbs(n0,n0,n1,sa,t,p)*g_tt));
}
/**
=========================================================================
method: gsw_entropy_from_ct (sa, ct)
=========================================================================
!
!  Calculates specific entropy of seawater from Conservative Temperature.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  entropy  =  specific entropy                                   [ deg C ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_entropy_from_ct(double sa, double ct)
{
    double  pt0;

    pt0 = gsw_pt_from_ct(sa, ct);

    return (-gsw_gibbs(0,1,0,sa,pt0,0));
}
/**
==========================================================================
method: gsw_pot_enthalpy_from_pt_ice (pt0_ice)
==========================================================================
!
!  Calculates the potential enthalpy of ice from potential temperature of
!  ice (whose reference sea pressure is zero dbar).
!
!  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
!
!  gsw_pot_enthalpy_ice  =  potential enthalpy of ice              [ J/kg ]
!--------------------------------------------------------------------------
*/
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
/**
--------------------------------------------------------------------------
! isobaric melting enthalpy and isobaric evaporation enthalpy
!--------------------------------------------------------------------------

==========================================================================
method gsw_latentheat_melting(sa,p)
==========================================================================

! Calculates latent heat, or enthalpy, of melting.
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
!
! latentheat_melting : latent heat of melting              [kg/m^3]
*/
double TeosBase::gsw_latentheat_melting(double sa, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  tf = gsw_t_freezing(sa,p,0.0);

    return (1000.0*(gsw_chem_potential_water_t_exact(sa,tf,p)
           - (gtc.gsw_t0 + tf)*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p))
           - gsw_enthalpy_ice(tf,p));
}
/**
==========================================================================
method: gsw_enthalpy_ice (t, p)
==========================================================================
!
! Calculates the specific enthalpy of ice (h_Ih).
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar )
!
!  gsw_enthalpy_ice  :  specific enthalpy of ice                   [ J/kg ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_enthalpy_ice(double t, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_GIBBS_ICE_COEFFICIENTS use gic */

    double  tau, dzi, g0;
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
/**
==========================================================================
elemental subroutine gsw_entropy_first_derivatives (sa, ct, eta_sa, eta_ct)
! =========================================================================
!
!  Calculates the following two partial derivatives of specific entropy
!  (eta)
!   (1) eta_SA, the derivative with respect to Absolute Salinity at
!       constant Conservative Temperature, and
!   (2) eta_CT, the derivative with respect to Conservative Temperature at
!       constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  eta_SA =  The derivative of specific entropy with respect to
!            Absolute Salinity (in units of g kg^-1) at constant
!            Conservative Temperature.
!            eta_SA has units of:         [ J/(kg K(g/kg))]  or [ J/(g K) ]
!  eta_CT =  The derivative of specific entropy with respect to
!            Conservative Temperature at constant Absolute Salinity.
!            eta_CT has units of:                            [ J/(kg K^2) ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_entropy_first_derivatives(double sa, double ct, double *eta_sa,
        double *eta_ct)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  pt, pr0 = 0.0;
    int     n0=0, n1=1;

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
/**
==========================================================================
elemental subroutine gsw_ct_from_rho (rho, sa, p, ct, ct_multiple)
! =========================================================================
!
!  Calculates the Conservative Temperature of a seawater sample, for given
!  values of its density, Absolute Salinity and sea pressure (in dbar).
!
!  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
!   Note. This input has not had 1000 kg/m^3 subtracted from it.
!     That is, it is 'density', not 'density anomaly'.
!  SA   =  Absolute Salinity                                       [ g/kg ]
!  p    =  sea pressure                                            [ dbar ]
!          ( i.e. absolute pressure - 10.1325 dbar )
!
!  CT  =  Conservative Temperature  (ITS-90)                      [ deg C ]
!  CT_multiple  =  Conservative Temperature  (ITS-90)             [ deg C ]
!    Note that at low salinities, in brackish water, there are two possible
!      Conservative Temperatures for a single density.  This programme will
!      output both valid solutions.  To see this second solution the user
!      must call the programme with two outputs (i.e. [CT,CT_multiple]), if
!      there is only one possible solution and the programme has been
!      called with two outputs the second variable will be set to NaN.
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_ct_from_rho(double rho, double sa, double p,
                                double *ct, double *ct_multiple)
{
    int     number_of_iterations;
    double  a, alpha_freezing, alpha_mean, b, c, ct_a, ct_b, ct_diff,
                ct_freezing, ct_max_rho, ct_mean, ct_old,
                delta_ct, delta_v, factor, factorqa, factorqb,
                rho_40, rho_extreme, rho_freezing, rho_max, rho_mean,
                rho_old, sqrt_disc, top, v_ct, v_lab;

    /**
        ! alpha_limit is the positive value of the thermal expansion coefficient
        ! which is used at the freezing temperature to distinguish between
        ! salty and fresh water.
    */
    double  alpha_limit = 1e-5;

    /**
        ! rec_half_rho_TT is a constant representing the reciprocal of half the
        ! second derivative of density with respect to temperature near the
        ! temperature of maximum density.
    */
    double  rec_half_rho_tt = -110.0;

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

    /** Assumes that the seawater is always unsaturated with air*/
    ct_freezing = gsw_ct_freezing_poly(sa,p,0.0);

    gsw_rho_alpha_beta(sa,ct_freezing,p,&rho_freezing,&alpha_freezing,NULL);

    /** reset the extreme values*/
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
            /** Set the initial value of the quadratic solution roots.*/
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

            /** Set the initial value of the quadratic solution roots.*/
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
                ! After seven iterations of this quadratic iterative procedure,
                ! the error in rho is no larger than 4.6x10^-13 kg/m^3.
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

    /** Begin the modified Newton-Raphson iterative method*/

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
        ! After three iterations of this modified Newton-Raphson iteration,
        ! the error in rho is no larger than 1.6x10^-12 kg/m^3.
    */
    if (ct_multiple != NULL)
    {
        *ct_multiple = cppGSW_INVALID_VALUE;
    }

    return;
}
/**
==========================================================================
method: gsw_ct_freezing_poly (sa, p, saturation_fraction)
==========================================================================
!
!  Calculates the Conservative Temperature at which seawater freezes.
!  The error of this fit ranges between -5e-4 K and 6e-4 K when compared
!  with the Conservative Temperature calculated from the exact in-situ
!  freezing temperature which is found by a Newton-Raphson iteration of the
!  equality of the chemical potentials of water in seawater and in ice.
!  Note that the Conservative temperature freezing temperature can be found
!  by this exact method using the method gsw_CT_freezing.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  CT_freezing = Conservative Temperature at freezing of seawater [ deg C ]
!                That is, the freezing temperature expressed in
!                terms of Conservative Temperature (ITS-90).
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_ct_freezing_poly(double sa, double p, double saturation_fraction)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */

    double  p_r, sa_r, x, return_value;

    sa_r    = sa*1.0e-2;
    x       = sqrt(sa_r);
    p_r     = p*1.0e-4;

    return_value = gfpc.c0
        + sa_r*(gfpc.c1 + x*(gfpc.c2 + x*(gfpc.c3 + x*(gfpc.c4 + x*(gfpc.c5 + gfpc.c6*x)))))
        + p_r*(gfpc.c7 + p_r*(gfpc.c8 + gfpc.c9*p_r)) + sa_r*p_r*(gfpc.c10 + p_r*(gfpc.c12
        + p_r*(gfpc.c15 + gfpc.c21*sa_r)) + sa_r*(gfpc.c13 + gfpc.c17*p_r + gfpc.c19*sa_r)
        + x*(gfpc.c11 + p_r*(gfpc.c14 + gfpc.c18*p_r) + sa_r*(gfpc.c16 + gfpc.c20*p_r + gfpc.c22*sa_r)));

    /** Adjust for the effects of dissolved air */
    return_value = return_value - saturation_fraction*
        (1e-3)*(2.4 - gfpc.a*sa)*(1.0 + gfpc.b*(1.0 - sa/gtc.gsw_sso));

    return (return_value);
}
/**
==========================================================================
elemental subroutine gsw_rho_alpha_beta (sa, ct, p, rho, alpha, beta)
==========================================================================
!
!  Calculates in-situ density, the appropiate thermal expansion coefficient
!  and the appropriate saline contraction coefficient of seawater from
!  Absolute Salinity and Conservative Temperature.  This method uses the
!  computationally-efficient expression for specific volume in terms of
!  SA, CT and p (Roquet et al., 2014).
!
!  Note that potential density (pot_rho) with respect to reference pressure
!  p_ref is obtained by calling this method with the pressure argument
!  being p_ref as in [pot_rho, ~, ~] = gsw_rho_alpha_beta(SA,CT,p_ref).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho    =  in-situ density                                       [ kg/m ]
!  alpha  =  thermal expansion coefficient                          [ 1/K ]
!            with respect to Conservative Temperature
!  beta   =  saline (i.e. haline) contraction                      [ kg/g ]
!            coefficient at constant Conservative Temperature
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_rho_alpha_beta (double sa, double ct, double p, double *rho, double *alpha,
                        double *beta)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  v, v_ct_part, v_sa_part, xs, ys, z;

        xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys = ct*0.025;
        z = p*1e-4;

        v = gsvco.v000
        + xs*(gsvco.v010 + xs*(gsvco.v020 + xs*(gsvco.v030 + xs*(gsvco.v040 + xs*(gsvco.v050
        + gsvco.v060*xs))))) + ys*(gsvco.v100 + xs*(gsvco.v110 + xs*(gsvco.v120 + xs*(gsvco.v130 + xs*(gsvco.v140
        + gsvco.v150*xs)))) + ys*(gsvco.v200 + xs*(gsvco.v210 + xs*(gsvco.v220 + xs*(gsvco.v230 + gsvco.v240*xs)))
        + ys*(gsvco.v300 + xs*(gsvco.v310 + xs*(gsvco.v320 + gsvco.v330*xs)) + ys*(gsvco.v400 + xs*(gsvco.v410
        + gsvco.v420*xs) + ys*(gsvco.v500 + gsvco.v510*xs + gsvco.v600*ys))))) + z*(gsvco.v001 + xs*(gsvco.v011
        + xs*(gsvco.v021 + xs*(gsvco.v031 + xs*(gsvco.v041 + gsvco.v051*xs)))) + ys*(gsvco.v101 + xs*(gsvco.v111
        + xs*(gsvco.v121 + xs*(gsvco.v131 + gsvco.v141*xs))) + ys*(gsvco.v201 + xs*(gsvco.v211 + xs*(gsvco.v221
        + gsvco.v231*xs)) + ys*(gsvco.v301 + xs*(gsvco.v311 + gsvco.v321*xs) + ys*(gsvco.v401 + gsvco.v411*xs
        + gsvco.v501*ys)))) + z*(gsvco.v002 + xs*(gsvco.v012 + xs*(gsvco.v022 + xs*(gsvco.v032 + gsvco.v042*xs)))
        + ys*(gsvco.v102 + xs*(gsvco.v112 + xs*(gsvco.v122 + gsvco.v132*xs)) + ys*(gsvco.v202 + xs*(gsvco.v212
        + gsvco.v222*xs) + ys*(gsvco.v302 + gsvco.v312*xs + gsvco.v402*ys))) + z*(gsvco.v003 + xs*(gsvco.v013
        + gsvco.v023*xs) + ys*(gsvco.v103 + gsvco.v113*xs + gsvco.v203*ys) + z*(gsvco.v004 + gsvco.v014*xs + gsvco.v104*ys
        + z*(gsvco.v005 + gsvco.v006*z)))));

        if (rho != NULL)
            *rho = 1.0/v;

        if (alpha != NULL) {

            v_ct_part = gsvco.a000
                + xs*(gsvco.a100 + xs*(gsvco.a200 + xs*(gsvco.a300
                + xs*(gsvco.a400 + gsvco.a500*xs))))
                + ys*(gsvco.a010 + xs*(gsvco.a110 + xs*(gsvco.a210 + xs*(gsvco.a310 + gsvco.a410*xs)))
                + ys*(gsvco.a020 + xs*(gsvco.a120 + xs*(gsvco.a220 + gsvco.a320*xs)) + ys*(gsvco.a030
                + xs*(gsvco.a130 + gsvco.a230*xs) + ys*(gsvco.a040 + gsvco.a140*xs + gsvco.a050*ys ))))
                + z*(gsvco.a001 + xs*(gsvco.a101 + xs*(gsvco.a201 + xs*(gsvco.a301 + gsvco.a401*xs)))
                + ys*(gsvco.a011 + xs*(gsvco.a111 + xs*(gsvco.a211 + gsvco.a311*xs)) + ys*(gsvco.a021
                + xs*(gsvco.a121 + gsvco.a221*xs) + ys*(gsvco.a031 + gsvco.a131*xs + gsvco.a041*ys)))
                + z*(gsvco.a002 + xs*(gsvco.a102 + xs*(gsvco.a202 + gsvco.a302*xs)) + ys*(gsvco.a012
                + xs*(gsvco.a112 + gsvco.a212*xs) + ys*(gsvco.a022 + gsvco.a122*xs + gsvco.a032*ys))
                + z*(gsvco.a003 + gsvco.a103*xs + gsvco.a013*ys + gsvco.a004*z)));

            *alpha = 0.025*v_ct_part/v;

        }

        if (beta != NULL) {

            v_sa_part = gsvco.b000
                + xs*(gsvco.b100 + xs*(gsvco.b200 + xs*(gsvco.b300
                + xs*(gsvco.b400 + gsvco.b500*xs))))
                + ys*(gsvco.b010 + xs*(gsvco.b110 + xs*(gsvco.b210 + xs*(gsvco.b310 + gsvco.b410*xs)))
                + ys*(gsvco.b020 + xs*(gsvco.b120 + xs*(gsvco.b220 + gsvco.b320*xs)) + ys*(gsvco.b030
                + xs*(gsvco.b130 + gsvco.b230*xs) + ys*(gsvco.b040 + gsvco.b140*xs + gsvco.b050*ys))))
                + z*(gsvco.b001 + xs*(gsvco.b101 + xs*(gsvco.b201 + xs*(gsvco.b301 + gsvco.b401*xs)))
                + ys*(gsvco.b011 + xs*(gsvco.b111 + xs*(gsvco.b211 + gsvco.b311*xs)) + ys*(gsvco.b021
                + xs*(gsvco.b121 + gsvco.b221*xs) + ys*(gsvco.b031 + gsvco.b131*xs + gsvco.b041*ys)))
                + z*(gsvco.b002 + xs*(gsvco.b102 + xs*(gsvco.b202 + gsvco.b302*xs))+ ys*(gsvco.b012
                + xs*(gsvco.b112 + gsvco.b212*xs) + ys*(gsvco.b022 + gsvco.b122*xs + gsvco.b032*ys))
                + z*(gsvco.b003 + gsvco.b103*xs + gsvco.b013*ys + gsvco.b004*z)));

            *beta = -v_sa_part*0.5*gtc.gsw_sfac/(v*xs);

        }
}
/**
==========================================================================
elemental subroutine gsw_frazil_properties (sa_bulk, h_bulk, p, &
                                            sa_final, ct_final, w_ih_final)
==========================================================================
!
!  Calculates the mass fraction of ice (mass of ice divided by mass of ice
!  plus seawater), w_Ih_final, which results from given values of the bulk
!  Absolute Salinity, SA_bulk, bulk enthalpy, h_bulk, occuring at pressure
!  p.  The final values of Absolute Salinity, SA_final, and Conservative
!  Temperature, CT_final, of the interstitial seawater phase are also
!  returned.  This code assumes that there is no dissolved air in the
!  seawater (that is, saturation_fraction is assumed to be zero
!  throughout the code).
!
!  When the mass fraction w_Ih_final is calculated as being a positive
!  value, the seawater-ice mixture is at thermodynamic equlibrium.
!
!  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
!  is sufficiently large (i.e. sufficiently "warm") so that there is no ice
!  present in the final state.  In this case the final state consists of
!  only seawater rather than being an equlibrium mixture of seawater and
!  ice which occurs when w_Ih_final is positive.  Note that when
!  w_Ih_final = 0, the final seawater is not at the freezing temperature.
!
!  SA_bulk =  bulk Absolute Salinity of the seawater and ice mixture
!                                                                  [ g/kg ]
!  h_bulk  =  bulk enthalpy of the seawater and ice mixture        [ J/kg ]
!  p       =  sea pressure                                         [ dbar ]
!             ( i.e. absolute pressure - 10.1325 dbar )
!
!  SA_final    =  Absolute Salinity of the seawater in the final state,
!                 whether or not any ice is present.               [ g/kg ]
!  CT_final    =  Conservative Temperature of the seawater in the the final
!                 state, whether or not any ice is present.       [ deg C ]
!  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
!                 If this ice mass fraction is positive, the system is at
!                 thermodynamic equilibrium.  If this ice mass fraction is
!                 zero there is no ice in the final state which consists
!                 only of seawater which is warmer than the freezing
!                 temperature.                                   [unitless]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_frazil_properties(double sa_bulk, double h_bulk, double p,
                                    double *sa_final, double *ct_final,
                                    double *w_ih_final)
{
    int     number_of_iterations;
    double  cp_ih, ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly,
                func, func0, hf, h_hat_ct, h_hat_sa,
                h_ihf, sa, tf_sa, tf, w_ih_mean, w_ih_old, w_ih,
        /**
            ! Throughout this code seawater is taken to contain
            ! no dissolved air.
        */
    saturation_fraction = 0.0,
    num_f = 5.0e-2, num_f2 = 6.9e-7, num_p = 2.21;
    /**
        !---------------
        ! Finding func0
        !--------------
    */
    ctf = gsw_ct_freezing(sa_bulk,p,saturation_fraction);
    func0 = h_bulk - gsw_enthalpy_ct_exact(sa_bulk,ctf,p);
    /**
        !-----------------------------------------------------------------------
        ! When func0 is zero or positive we can immediately calculate the three
        ! outputs, as the bulk enthalpy, h_bulk, is too large to allow any ice
        ! at thermodynamic equilibrium. The result will be (warm) seawater with
        ! no frazil ice being present. The three outputs can be set and the rest
        ! of this code does not need to be performed.
        !-----------------------------------------------------------------------
    */
    if (func0 >= 0.0)
    {
        *sa_final = sa_bulk;
        *ct_final = gsw_ct_from_enthalpy_exact(sa_bulk,h_bulk,p);
        *w_ih_final = 0.0;

        return;
    }
    /**
        !-----------------------------------------------------------------------
        ! Begin to find the solution for those data points that have func0 < 0,
        ! implying that the output will be a positive ice mass fraction
        ! w_Ih_final.
        !
        ! Do a quasi-Newton step with a separate polynomial estimate of the
        ! derivative of func with respect to the ice mass fraction.  This
        ! section of the code delivers initial values of both w_Ih and SA to
        ! the rest of the more formal modified Newtons Method approach of
        ! McDougall and Wotherspoon (2014).
        !-----------------------------------------------------------------------
    */
    dfunc_dw_ih_mean_poly = 3.347814e+05
                                - num_f*func0*(1.0 + num_f2*func0) - num_p*p;

    w_ih = gsw_min(-func0/dfunc_dw_ih_mean_poly, 0.95);
    sa = sa_bulk/(1.0 - w_ih);

    if (sa < 0.0 || sa > 120.0)
    {
        *sa_final = cppGSW_INVALID_VALUE;
        *ct_final = *sa_final;
        *w_ih_final = *sa_final;

        return;
    }
    /**
        !-----------------------------------------------------------------------
        ! Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
        ! fed into the iterative Newton's Method.
        !-----------------------------------------------------------------------
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
        !-----------------------------------------------------------------------
        ! Enter the main McDougall-Wotherspoon (2014) modified Newton-Raphson
        | loop
        !-----------------------------------------------------------------------
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
            ! This will only trap cases that are smaller than zero by just
            ! machine precision
        */
        *sa_final = sa_bulk;
        *ct_final = gsw_ct_from_enthalpy_exact(*sa_final,h_bulk,p);
        *w_ih_final = 0.0;
    }
}
/**
==========================================================================
method: gsw_cp_ice (t, p)
==========================================================================
!
!  Calculates the isobaric heat capacity of seawater.
!
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!          ( i.e. absolute pressure - 10.1325 dbar )
!
!  gsw_cp_ice  =  heat capacity of ice                       [J kg^-1 K^-1]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_cp_ice(double t, double p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    return (-(t + gtc.gsw_t0)*gsw_gibbs_ice(2,0,t,p));
}
/**
==========================================================================
elemental subroutine gsw_ct_freezing_first_derivatives (sa, p, &
                          saturation_fraction, ctfreezing_sa, ctfreezing_p)
==========================================================================
!
!  Calculates the first derivatives of the Conservative Temperature at
!  which seawater freezes, with respect to Absolute Salinity SA and
!  pressure P (in Pa).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  CTfreezing_SA = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to Absolute Salinity at
!                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  CTfreezing_P  = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to pressure (in Pa) at
!                  fixed Absolute Salinity                         [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_ct_freezing_first_derivatives(double sa, double p,
                                                double saturation_fraction,
                                                double *ctfreezing_sa,
                                                double *ctfreezing_p)
{
    double  tf_sa, tf_p, ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t, tf;

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
/**
==========================================================================
elemental subroutine gsw_t_freezing_first_derivatives (sa, p, &
                            saturation_fraction, tfreezing_sa, tfreezing_p)
==========================================================================
!
!  Calculates the first derivatives of the in-situ temperature at which
!  seawater freezes with respect to Absolute Salinity SA and pressure P (in
!  Pa).  These expressions come from differentiating the expression that
!  defines the freezing temperature, namely the equality between the
!  chemical potentials of water in seawater and in ice.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  tfreezing_SA = the derivative of the in-situ freezing temperature
!                 (ITS-90) with respect to Absolute Salinity at fixed
!                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  tfreezing_P  = the derivative of the in-situ freezing temperature
!                 (ITS-90) with respect to pressure (in Pa) at fixed
!                 Absolute Salinity                                [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_t_freezing_first_derivatives(double sa, double p,
        double saturation_fraction, double *tfreezing_sa, double *tfreezing_p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        double  rec_denom, tf, g_per_kg = 1000.0;

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

        return;
}
/**
==========================================================================
elemental subroutine gsw_ct_first_derivatives_wrt_t_exact (sa, t, p, &
                                       ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
==========================================================================
!
!  Calculates the following three derivatives of Conservative Temperature.
!  These derivatives are done with respect to in-situ temperature t (in the
!  case of CT_T_wrt_t) or at constant in-situ tempertature (in the cases of
!  CT_SA_wrt_t and CT_P_wrt_t).
!   (1) CT_SA_wrt_t, the derivative of CT with respect to Absolute Salinity
!       at constant t and p, and
!   (2) CT_T_wrt_t, derivative of CT with respect to in-situ temperature t
!       at constant SA and p.
!   (3) CT_P_wrt_t, derivative of CT with respect to pressure P (in Pa) at
!       constant SA and t.
!
!  This method uses the full Gibbs method. Note that this method
!  avoids the NaN that would exist in CT_SA_wrt_t at SA = 0 if it were
!  evaluated in the straightforward way from the derivatives of the Gibbs
!  method method.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar)
!
!  CT_SA_wrt_t  =  The first derivative of Conservative Temperature with
!                  respect to Absolute Salinity at constant t and p.
!                                              [ K/(g/kg)]  i.e. [ K kg/g ]
!  CT_T_wrt_t  =  The first derivative of Conservative Temperature with
!                 respect to in-situ temperature, t, at constant SA and p.
!                                                              [ unitless ]
!  CT_P_wrt_t  =  The first derivative of Conservative Temperature with
!                 respect to pressure P (in Pa) at constant SA and t.
!                                                                  [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_ct_first_derivatives_wrt_t_exact(double sa, double t, double p,
                                                    double *ct_sa_wrt_t, double *ct_t_wrt_t,
                                                    double *ct_p_wrt_t)
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
/**
==========================================================================
method: gsw_entropy_ice (t, p)
==========================================================================
!
!  Calculates specific entropy of ice.
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  ice_entropy  =  specific entropy of ice                 [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_entropy_ice(double t, double p)
{
        return (-gsw_gibbs_ice(1,0,t,p));
}
/**
==========================================================================
method: gsw_specvol_ice (t, p)
==========================================================================
!
!  Calculates the specific volume of ice.
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  specvol_ice  =  specific volume                               [ m^3/kg ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_specvol_ice(double t, double p)
{
        return (gsw_gibbs_ice(0,1,t,p));
}
/**
==========================================================================
elemental subroutine gsw_frazil_properties_potential_poly (sa_bulk, &
                             h_pot_bulk, p, sa_final, ct_final, w_ih_final)
==========================================================================
!
!  Calculates the mass fraction of ice (mass of ice divided by mass of ice
!  plus seawater), w_Ih_final, which results from given values of the bulk
!  Absolute Salinity, SA_bulk, bulk potential enthalpy, h_pot_bulk,
!  occuring at pressure p.  The final equilibrium values of Absolute
!  Salinity, SA_final, and Conservative Temperature, CT_final, of the
!  interstitial seawater phase are also returned.  This code assumes that
!  there is no dissolved air in the seawater (that is, saturation_fraction
!  is assumed to be zero thoughout the code).
!
!  When the mass fraction w_Ih_final is calculated as being a positive
!  value, the seawater-ice mixture is at thermodynamic equlibrium.
!
!  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
!  is sufficiently large (i.e. sufficiently "warm") so that there is no ice
!  present in the final state.  In this case the final state consists of
!  only seawater rather than being an equlibrium mixture of seawater and
!  ice which occurs when w_Ih_final is positive.  Note that when
!  w_Ih_final = 0, the final seawater is not at the freezing temperature.
!
!  Note that this code uses the polynomial forms of CT_freezing and
!  pot_enthalpy_ice_freezing. This code is intended to be used in ocean
!  models where the model prognostic variables are SA_bulk and h_pot_bulk.
!
!  SA_bulk     =  bulk Absolute Salinity of the seawater and ice mixture
!                                                                  [ g/kg ]
!  h_pot_bulk  =  bulk potential enthalpy of the seawater and ice mixture
!                                                                  [ J/kg ]
!  p           =  sea pressure                                  [ dbar ]
!                  ( i.e. absolute pressure - 10.1325 dbar )
!
!  SA_final    =  Absolute Salinity of the seawater in the final state,
!                 whether or not any ice is present.               [ g/kg ]
!  CT_final    =  Conservative Temperature of the seawater in the the final
!                 state, whether or not any ice is present.       [ deg C ]
!  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
!                 If this ice mass fraction is positive, the system is at
!                 thermodynamic equilibrium.  If this ice mass fraction is
!                 zero there is no ice in the final state which consists
!                 only of seawater which is warmer than the freezing
!                 temperature.                                   [unitless]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_frazil_properties_potential_poly(double sa_bulk, double h_pot_bulk,
                                                    double p, double *sa_final,
                                                    double *ct_final, double *w_ih_final)
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
        !-----------------------------------------------------------------------
        ! Finding func0.  This is the value of the method, func, that would
        ! result in the output w_Ih_final being exactly zero.
        !-----------------------------------------------------------------------
    */
    func0 = h_pot_bulk - gtc.gsw_cp0
                *gsw_ct_freezing_poly(sa_bulk,p,saturation_fraction);
    /**
        !-----------------------------------------------------------------------
        ! Setting the three outputs for data points that have func0 non-negative
        !-----------------------------------------------------------------------
    */
    if (func0 >= 0.0)
    {
        /**
            ! When func0 is zero or positive then the final answer will contain
            ! no frazil ice; that is, it will be pure seawater that is warmer
            ! han the freezing temperature.  If func0 >= 0 we do not need to go
            ! through the modified Newton-Raphson procedure and we can simply
            ! write down the answer, as in the following 4 lines of code.
        */
        *sa_final = sa_bulk;
        *ct_final = h_pot_bulk/gtc.gsw_cp0;
        *w_ih_final = 0.0;

        return;
    }
    /**
        !-----------------------------------------------------------------------
        !Begin finding the solution for data points that have func0 < 0, so that
        !the output will have a positive ice mass fraction w_Ih_final.
        !-----------------------------------------------------------------------
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
            ! The ice mass fraction out of this code is restricted to be
            ! less than 0.9.
        */
        *sa_final = cppGSW_INVALID_VALUE;
        *ct_final = *sa_final;
        *w_ih_final = *sa_final;

        return;
    }
    /**
        ! The initial guess at the absolute salinity of the interstitial
        ! seawater
    */
    sa = sa_bulk/(1.0 - w_ih);
    /**
        !-----------------------------------------------------------------------
        ! Doing a Newton step with a separate polynomial estimate of the mean
        ! derivative dfunc_dw_Ih_mean_poly.
        !-----------------------------------------------------------------------
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
        !-----------------------------------------------------------------------
        ! Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
        ! fed into Newton's Method.
        !-----------------------------------------------------------------------
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
        *ct_final = gsw_ct_freezing_poly(sa,p,saturation_fraction);
        *w_ih_final = w_ih;
    }
}
/**
==========================================================================
method: gsw_pot_enthalpy_ice_freezing_poly (sa, p)
==========================================================================
!
!  Calculates the potential enthalpy of ice at which seawater freezes.
!  The error of this fit ranges between -2.5 and 1 J/kg with an rms of
!  1.07, between SA of 0 and 120 g/kg and p between 0 and 10,000 dbar (the
!  error in the fit is between -0.7 and 0.7 with an rms of
!  0.3, between SA of 0 and 120 g/kg and p between 0 and 5,000 dbar) when
!  compared with the potential enthalpy calculated from the exact in-situ
!  freezing temperature which is found by a Newton-Raphson iteration of the
!  equality of the chemical potentials of water in seawater and in ice.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  pot_enthalpy_ice_freezing = potential enthalpy of ice at freezing
!                              of seawater                         [ J/kg ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_pot_enthalpy_ice_freezing_poly(double sa, double p)
{
        double  p_r, sa_r, x,
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




/**
==========================================================================
elemental subroutine gsw_ct_freezing_first_derivatives_poly (sa, p, &
                          saturation_fraction, ctfreezing_sa, ctfreezing_p)
==========================================================================
!
!  Calculates the first derivatives of the Conservative Temperature at
!  which seawater freezes, with respect to Absolute Salinity SA and
!  pressure P (in Pa) of the comptationally efficient polynomial fit of the
!  freezing temperature (McDougall et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  CTfreezing_SA = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to Absolute Salinity at
!                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  CTfreezing_P  = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to pressure (in Pa) at
!                  fixed Absolute Salinity                         [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_ct_freezing_first_derivatives_poly(double sa, double p,
                                                        double saturation_fraction,
                                                        double *ctfreezing_sa,
                                                        double *ctfreezing_p)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */
    /** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */

    double  p_r, sa_r, x, d = -gfpc.a - gfpc.a*gfpc.b - 2.4*gfpc.b/gtc.gsw_sso,
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
            (gfpc.c7 + sa_r*(gfpc.c10 + x*(gfpc.c11 + x*(gfpc.c13 + x*(gfpc.c16 + x*(gfpc.c19 + gfpc.c22*x)))))
                + p_r*(2.0*gfpc.c8 + sa_r*(2.0*gfpc.c12 + x*(2.0*gfpc.c14 + x*(2.0*gfpc.c17
                + 2.0*gfpc.c20*x))) + p_r*(3.0*gfpc.c9 + sa_r*(3.0*gfpc.c15 + x*(3.0*gfpc.c18
                + 3.0*gfpc.c21*x)))))*1e-8;
    }
}
/**
==========================================================================
elemental subroutine gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(&
         sa, p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
==========================================================================
!
!  Calculates the first derivatives of the potential enthalpy of ice Ih at
!  which ice melts into seawater with Absolute Salinity SA and at pressure
!  p.  This code uses the comptationally efficient polynomial fit of the
!  freezing potential enthalpy of ice Ih (McDougall et al., 2015).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  pot_enthalpy_ice_freezing_SA = the derivative of the potential enthalpy
!                of ice at freezing (ITS-90) with respect to Absolute
!                salinity at fixed pressure  [ (J/kg)/(g/kg) ] i.e. [ J/g ]
!
!  pot_enthalpy_ice_freezing_P  = the derivative of the potential enthalpy
!                of ice at freezing (ITS-90) with respect to pressure
!                (in Pa) at fixed Absolute Salinity           [ (J/kg)/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(double sa, double p,
                                                                     double *pot_enthalpy_ice_freezing_sa,
                                                                     double *pot_enthalpy_ice_freezing_p)
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
/**
==========================================================================
elemental subroutine gsw_entropy_second_derivatives (sa, ct, eta_sa_sa, &
                                                     eta_sa_ct, eta_ct_ct)
! =========================================================================
!
!  Calculates the following three second-order partial derivatives of
!  specific entropy (eta)
!   (1) eta_SA_SA, the second derivative with respect to Absolute
!       Salinity at constant Conservative Temperature, and
!   (2) eta_SA_CT, the derivative with respect to Absolute Salinity and
!       Conservative Temperature.
!   (3) eta_CT_CT, the second derivative with respect to Conservative
!       Temperature at constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  eta_SA_SA =  The second derivative of specific entropy with respect
!               to Absolute Salinity (in units of g kg^-1) at constant
!               Conservative Temperature.
!               eta_SA_SA has units of:                 [ J/(kg K(g/kg)^2)]
!  eta_SA_CT =  The second derivative of specific entropy with respect
!               to Conservative Temperature at constant Absolute
!               Salinity. eta_SA_CT has units of:     [ J/(kg (g/kg) K^2) ]
!  eta_CT_CT =  The second derivative of specific entropy with respect
!               to Conservative Temperature at constant Absolute
!               Salinity.  eta_CT_CT has units of:           [ J/(kg K^3) ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_entropy_second_derivatives(double sa, double ct, double *eta_sa_sa,
                                                double *eta_sa_ct, double *eta_ct_ct)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  abs_pt, ct_pt, ct_sa, pt, ct_ct, pr0 = 0.0;
    int     n0=0, n1=1, n2=2;

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


/**
==========================================================================
elemental subroutine gsw_enthalpy_second_derivatives_ct_exact (sa, ct, p, &
                                                 h_sa_sa, h_sa_ct, h_ct_ct)
==========================================================================
!
!  Calculates three second-order derivatives of specific enthalpy (h).
!  Note that this method uses the full Gibbs method.
!
!  sa  =  Absolute Salinity                                        [ g/kg ]
!  ct  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_sa_sa  =  The second derivative of specific enthalpy with respect to
!              Absolute Salinity at constant ct & p.    [ J/(kg (g/kg)^2) ]
!  h_sa_ct  =  The second derivative of specific enthalpy with respect to
!              sa and ct at constant p.                  [ J/(kg K(g/kg)) ]
!  h_ct_ct  =  The second derivative of specific enthalpy with respect to
!              ct at constant sa and p.                      [ J/(kg K^2) ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_enthalpy_second_derivatives_ct_exact(double sa, double ct, double p,
                                                        double *h_sa_sa, double *h_sa_ct,
                                                        double *h_ct_ct)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    double  factor, gsa_pt0, gsat_pt0, gsat, part_b, pt0, h_ct_ct_val,
                rec_abs_pt0, rec_gtt_pt0, rec_gtt, t, temp_ratio,
                gsasa, gsasa_pt0, pr0 = 0.0, sa_small = 1e-100;

    int     n0=0, n1=1, n2=2;

    pt0 = gsw_pt_from_ct(sa,ct);
    rec_abs_pt0 = 1.0/(gtc.gsw_t0 + pt0);
    t = gsw_pt_from_t(sa,pt0,pr0,p);
    temp_ratio = (gtc.gsw_t0 + t)*rec_abs_pt0;

    rec_gtt_pt0 = 1.0/gsw_gibbs(n0,n2,n0,sa,pt0,pr0);
    rec_gtt = 1.0/gsw_gibbs(n0,n2,n0,sa,t,p);

    /**h_ct_ct is naturally well-behaved as sa approaches zero.*/
    h_ct_ct_val = gtc.gsw_cp0*gtc.gsw_cp0*
            (temp_ratio*rec_gtt_pt0 - rec_gtt)*(rec_abs_pt0*rec_abs_pt0);

    if (h_ct_ct != NULL)
    {
        *h_ct_ct = h_ct_ct_val;
    }

    if ((h_sa_sa == NULL) && (h_sa_ct == NULL))
    {
        return;
    }

    gsat_pt0 = gsw_gibbs(n1,n1,n0,sa,pt0,pr0);
    gsat = gsw_gibbs(n1,n1,n0,sa,t,p);
    gsa_pt0 = gsw_gibbs(n1,n0,n0,sa,pt0,pr0);

    part_b = (temp_ratio*gsat_pt0*rec_gtt_pt0 - gsat*rec_gtt)*rec_abs_pt0;
    factor = gsa_pt0/gtc.gsw_cp0;

    if (h_sa_sa != NULL)
    {
        gsasa = gsw_gibbs(n2,n0,n0,sa,t,p);
        gsasa_pt0 = gsw_gibbs(n2,n0,n0,sa,pt0,pr0);
        /**
          ! h_sa_sa has a singularity at sa = 0, and blows up as sa
          ! approaches zero.
        */
        *h_sa_sa = gsasa - temp_ratio*gsasa_pt0
                + temp_ratio*gsat_pt0*gsat_pt0*rec_gtt_pt0
                - gsat*gsat*rec_gtt
                - 2.0*gsa_pt0*part_b + (factor*factor)*h_ct_ct_val;

    }

    if (h_sa_ct == NULL)
    {
        return;
    }

    /**
        ! h_sa_ct should not blow up as sa approaches zero.
        ! The following lines of code ensure that the h_sa_ct output
        ! of this method does not blow up in this limit.
        ! That is, when sa < 1e-100 g/kg, we force the h_sa_ct
        ! output to be the same as if sa = 1e-100 g/kg.
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

    *h_sa_ct  = gtc.gsw_cp0*part_b - factor*h_ct_ct_val;
}
/**
==========================================================================
elemental subroutine gsw_frazil_ratios_adiabatic (sa, p, w_ih, &
                              dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
==========================================================================
!
!  Calculates the ratios of SA, CT and P changes when frazil ice forms or
!  melts in response to an adiabatic change in pressure of a mixture of
!  seawater and frazil ice crystals.
!
!  Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than
!  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero
!  whereas dCT/dSA would then be infinite.
!
!  Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
!  derivatives with the pressure measured in Pa not dbar.
!
!  SA  =  Absolute Salinity of seawater                            [ g/kg ]
!  p   =  sea pressure of seawater at which melting occurs         [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar )
!  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
!           sum of the masses of ice and seawater.  That is, the mass of
!           ice divided by the mass of the final mixed fluid.
!           w_Ih must be between 0 and 1.                      [ unitless ]
!
!  dSA_dCT_frazil =  the ratio of the changes in Absolute Salinity
!                    to that of Conservative Temperature       [ g/(kg K) ]
!  dSA_dP_frazil  =  the ratio of the changes in Absolute Salinity
!                    to that of pressure (in Pa)              [ g/(kg Pa) ]
!  dCT_dP_frazil  =  the ratio of the changes in Conservative Temperature
!                    to that of pressure (in Pa)                   [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_frazil_ratios_adiabatic (double sa, double p, double w_ih,
                                            double *dsa_dct_frazil, double *dsa_dp_frazil,
                                            double *dct_dp_frazil)
{
    double  bracket1, bracket2, cp_ih, gamma_ih, h, h_ih, part,
                rec_bracket3, tf, wcp, h_hat_sa, h_hat_ct, tf_sa, tf_p,
                ctf, ctf_sa, ctf_p;
    double  saturation_fraction = 0.0;

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
/**
==========================================================================
method: gsw_adiabatic_lapse_rate_ice (t, p)
==========================================================================
!
!  Calculates the adiabatic lapse rate of ice.
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!    Note.  The output is in unit of degress Celsius per Pa,
!      (or equivilently K/Pa) not in units of K/dbar.
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_adiabatic_lapse_rate_ice(double t, double p)
{
    return (-gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(2,0,t,p));
}
/**
==========================================================================
elemental subroutine gsw_frazil_properties_potential (sa_bulk, h_pot_bulk,&
                                         p, sa_final, ct_final, w_ih_final)
==========================================================================
!
!  Calculates the mass fraction of ice (mass of ice divided by mass of ice
!  plus seawater), w_Ih_final, which results from given values of the bulk
!  Absolute Salinity, SA_bulk, bulk potential enthalpy, h_pot_bulk,
!  occuring at pressure p.  The final equilibrium values of Absolute
!  Salinity, SA_final, and Conservative Temperature, CT_final, of the
!  interstitial seawater phase are also returned.  This code assumes that
!  there is no dissolved air in the seawater (that is, saturation_fraction
!  is assumed to be zero thoughout the code).
!
!  When the mass fraction w_Ih_final is calculated as being a positive
!  value, the seawater-ice mixture is at thermodynamic equlibrium.
!
!  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk,
!  is sufficiently large (i.e. sufficiently "warm") so that there is no ice
!  present in the final state.  In this case the final state consists of
!  only seawater rather than being an equlibrium mixture of seawater and
!  ice which occurs when w_Ih_final is positive.  Note that when
!  w_Ih_final = 0, the final seawater is not at the freezing temperature.
!
!  Note that this code uses the exact forms of CT_freezing and
!  pot_enthalpy_ice_freezing.
!
!  SA_bulk     =  bulk Absolute Salinity of the seawater and ice mixture
!                                                                  [ g/kg ]
!  h_pot_bulk  =  bulk potential enthalpy of the seawater and ice mixture
!                                                                  [ J/kg ]
!  p           =  sea pressure                                  [ dbar ]
!                  ( i.e. absolute pressure - 10.1325 dbar )
!
!  SA_final    =  Absolute Salinity of the seawater in the final state,
!                 whether or not any ice is present.               [ g/kg ]
!  CT_final    =  Conservative Temperature of the seawater in the the final
!                 state, whether or not any ice is present.       [ deg C ]
!  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
!                 If this ice mass fraction is positive, the system is at
!                 thermodynamic equilibrium.  If this ice mass fraction is
!                 zero there is no ice in the final state which consists
!                 only of seawater which is warmer than the freezing
!                 temperature.                                   [unitless]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_frazil_properties_potential(double sa_bulk, double h_pot_bulk, double p,
        double *sa_final, double *ct_final, double *w_ih_final)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    int     iterations, max_iterations;
    double  ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly,
            dpot_h_ihf_dsa, func, func0, h_pot_ihf, sa, w_ih_old, w_ih,
            x, xa, y, z;

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
        !-----------------------------------------------------------------------
        ! Finding func0.  This is the value of the method, func, that would
        ! result in the output w_Ih_final being exactly zero.
        !-----------------------------------------------------------------------
    */
    func0 = h_pot_bulk - gtc.gsw_cp0
                *gsw_ct_freezing(sa_bulk,p,saturation_fraction);
    /**
        !-----------------------------------------------------------------------
        ! Setting the three outputs for data points that have func0 non-negative
        !-----------------------------------------------------------------------
    */
    if (func0 >= 0.0)
    {
        /**
            ! When func0 is zero or positive then the final answer will contain
            ! no frazil ice; that is, it will be pure seawater that is warmer
            ! than the freezing temperature. If func0 >= 0 we do not need to go
            ! through the modified Newton-Raphson procedure and we can simply
            ! write down the answer, as in the following 4 lines of code.
        */
        *sa_final = sa_bulk;
        *ct_final = h_pot_bulk/gtc.gsw_cp0;
        *w_ih_final = 0.0;

        return;
    }
    /**
        !-----------------------------------------------------------------------
        !Begin finding the solution for data points that have func0 < 0, so that
        !the output will have a positive ice mass fraction w_Ih_final.
        !-----------------------------------------------------------------------
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
            ! The ice mass fraction out of this code is restricted to be
            !less than 0.9.
        */
        *sa_final = cppGSW_INVALID_VALUE;
        *ct_final = *sa_final;
        *w_ih_final = *sa_final;

        return;
    }

    /**
        !The initial guess at the absolute salinity of the interstitial seawater
    */
    sa = sa_bulk/(1.0 - w_ih);
    /**
        !-----------------------------------------------------------------------
        ! Doing a Newton step with a separate polynomial estimate of the mean
        ! derivative dfunc_dw_Ih_mean_poly.
        !-----------------------------------------------------------------------'
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
        !-----------------------------------------------------------------------
        ! Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
        ! fed into Newton's Method.
        !-----------------------------------------------------------------------
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
/**
==========================================================================
method: gsw_pot_enthalpy_ice_freezing (sa, p)
==========================================================================
!
!  Calculates the potential enthalpy of ice at which seawater freezes.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  pot_enthalpy_ice_freezing = potential enthalpy of ice at freezing
!                              of seawater                        [ deg C ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_pot_enthalpy_ice_freezing(double sa, double p)
{
        double  pt0_ice, t_freezing;

        t_freezing = gsw_t_freezing(sa,p,0.0) ;

        pt0_ice = gsw_pt0_from_t_ice(t_freezing,p);

        return (gsw_pot_enthalpy_from_pt_ice(pt0_ice));
}
/**
==========================================================================
elemental subroutine gsw_pot_enthalpy_ice_freezing_first_derivatives (sa, &
              p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
==========================================================================
!
!  Calculates the first derivatives of the potential enthalpy of ice at
!  which seawater freezes, with respect to Absolute Salinity SA and
!  pressure P (in Pa).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  pot_enthalpy_ice_freezing_SA = the derivative of the potential enthalpy
!                  of ice at freezing (ITS-90) with respect to Absolute
!                  salinity at fixed pressure  [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  pot_enthalpy_ice_freezing_P  = the derivative of the potential enthalpy
!                  of ice at freezing (ITS-90) with respect to pressure
!                  (in Pa) at fixed Absolute Salinity              [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa, double p,
                                                               double *pot_enthalpy_ice_freezing_sa,
                                                               double *pot_enthalpy_ice_freezing_p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        double  cp_ihf, pt_icef, ratio_temp, tf, tf_p, tf_sa;
        double  saturation_fraction = 0.0;

        tf = gsw_t_freezing(sa,p,saturation_fraction);
        pt_icef = gsw_pt0_from_t_ice(tf,p);
        ratio_temp = (gtc.gsw_t0 + pt_icef)/(gtc.gsw_t0 + tf);

        cp_ihf = gsw_cp_ice(tf,p);

        if ((pot_enthalpy_ice_freezing_sa != NULL) &&
            (pot_enthalpy_ice_freezing_p != NULL)) {
            gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
                        &tf_sa,&tf_p);
        } else if (pot_enthalpy_ice_freezing_sa != NULL) {
            gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
                                                  &tf_sa, NULL);
        } else if (pot_enthalpy_ice_freezing_p != NULL) {
            gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
                                                  NULL,&tf_p);
        }

        if (pot_enthalpy_ice_freezing_sa != NULL)
            *pot_enthalpy_ice_freezing_sa = ratio_temp*cp_ihf*tf_sa;

        if (pot_enthalpy_ice_freezing_p != NULL)
            *pot_enthalpy_ice_freezing_p = ratio_temp*cp_ihf*tf_p
                              - (gtc.gsw_t0 + pt_icef)*gsw_gibbs_ice(1,1,tf,p);
}
/**
 =========================================================================
method: gsw_pt0_from_t_ice (t, p)
! =========================================================================
!
!  Calculates potential temperature of ice Ih with a reference pressure of
!  0 dbar, from in-situ temperature, t.
!
!  t   =  in-situ temperature  (ITS-90)                           [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  pt0_ice  =  potential temperature of ice Ih with reference pressure of
!              zero dbar (ITS-90)                                 [ deg C ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_pt0_from_t_ice(double t, double p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        int     number_of_iterations;
        double  dentropy, dentropy_dt, pt0_ice,
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

        if (t < -45.0 && t > -273.0) {

            pt0_ice = t + p*(p1 + p*(p2 + p3*t) + t*(p4 + t*(p5 + p6*t)));

            if (pt0_ice < -gtc.gsw_t0) pt0_ice = -gtc.gsw_t0;
            /**
            ! we add 0.05d0 to the initial estimate of pt0_ice at
            ! temps less than -273 to ensure that it is never less than -273.15.
            */
            if (pt0_ice < -273.0) pt0_ice = pt0_ice + 0.05;

            dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);

            for (number_of_iterations = 1; number_of_iterations <= 3;
                number_of_iterations++) {
                pt0_ice_old = pt0_ice;
                dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
                pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
                ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
                dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
                pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
            }

        } else {

            pt0_ice = t + p*(s1 + t*(s2 + t*(s3 + t*s4)) + s5*p);
            dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);

            pt0_ice_old = pt0_ice;
            dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;

            pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
            ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
            dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
            pt0_ice = pt0_ice_old - dentropy/dentropy_dt;

        }

        if (pt0_ice < -273.0) {

            pt0_ice = t + p*(q1 + p*(q2 + q3*t) + t*(q4 + t*(q5 + q6*t)));
            /**
            ! add 0.01d0 to the initial estimate of pt_ice used in the
            ! derivative to ensure that it is never less than -273.15d0
            ! because the derivative approaches zero at absolute zero.
            */
            dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice+0.01);

            for (number_of_iterations = 1; number_of_iterations <= 3;
                number_of_iterations++) {
                pt0_ice_old = pt0_ice;
                dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
                pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
                ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
                /**
                ! add 0.01d0 to the estimate of ptm_ice for temperatures less
                | than -273 to ensure that they are never less than -273.15d0
                ! because the derivative approaches zero at absolute zero and
                ! the addition of 0.01d0 degrees c ensures that when we divide
                ! by the derivatve in the modified newton routine the method
                ! does not blow up.
                */
                ptm_ice = ptm_ice + 0.01;
                dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
                pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
            }

        }
        /**
        ! For temperatures less than -273.1 degsC the maximum error is less
        ! than 2x10^-7 degsC. For temperatures between -273.1 and 273 the
        ! maximum error is less than 8x10^-8 degsC, and for temperatures
        ! greater than -273 degsC the ! maximum error is 1.5x10^-12 degsC.
        ! These errors are over the whole ocean depths with p varying between
        ! 0 and 10,000 dbar, while the in-situ temperature varied independently
        ! between -273.15 and +2 degsC.
        */

        return (pt0_ice);
}
/**
! =========================================================================
! method: gsw_gibbs_ice_part_t (t, p)
! =========================================================================
!
!  part of the the first temperature derivative of Gibbs energy of ice
!  that is the outout is gibbs_ice(1,0,t,p) + S0
!
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!
!  gibbs_ice_part_t = part of temperature derivative       [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_gibbs_ice_part_t(double t, double p)
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
/**
 =========================================================================
    method: gsw_gibbs_ice_pt0_pt0 (pt0)
! =========================================================================
!
!  The second temperature derivative of Gibbs energy of ice at the
!  potential temperature with reference sea pressure of zero dbar.  That is
!  the output is gibbs_ice(2,0,pt0,0).
!
!  pt0  =  potential temperature with reference sea pressure of zero dbar
!                                                                 [ deg C ]
!
!  gsw_gibbs_ice_pt0_pt0 = temperature second derivative at pt0
!--------------------------------------------------------------------------
*/
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
/**
 =========================================================================
    method: gsw_gibbs_ice_pt0 (pt0)
! =========================================================================
!
!  Part of the the first temperature derivative of Gibbs energy of ice
!  that is the outout is "gibbs_ice(1,0,pt0,0) + s0"
!
!  pt0  =  potential temperature with reference sea pressure of zero dbar
!                                                                 [ deg C ]
!
!  gsw_gibbs_ice_pt0 = part of temperature derivative     [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------
*/
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

/**
==========================================================================
elemental subroutine gsw_frazil_ratios_adiabatic_poly (sa, p, w_ih, &
                              dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
==========================================================================
!
!  Calculates the ratios of SA, CT and P changes when frazil ice forms or
!  melts in response to an adiabatic change in pressure of a mixture of
!  seawater and frazil ice crystals.
!
!  Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than
!  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero
!  whereas dCT/dSA would then be infinite.
!
!  Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
!  derivatives with the pressure measured in Pa not dbar.
!
!  SA  =  Absolute Salinity of seawater                            [ g/kg ]
!  p   =  sea pressure of seawater at which melting occurs         [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar )
!  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
!           sum of the masses of ice and seawater.  That is, the mass of
!           ice divided by the mass of the final mixed fluid.
!           w_Ih must be between 0 and 1.                      [ unitless ]
!
!  dSA_dCT_frazil =  the ratio of the changes in Absolute Salinity
!                    to that of Conservative Temperature       [ g/(kg K) ]
!  dSA_dP_frazil  =  the ratio of the changes in Absolute Salinity
!                    to that of pressure (in Pa)              [ g/(kg Pa) ]
!  dCT_dP_frazil  =  the ratio of the changes in Conservative Temperature
!                    to that of pressure (in Pa)                   [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_frazil_ratios_adiabatic_poly(double sa, double p, double w_ih,
                                                double *dsa_dct_frazil, double *dsa_dp_frazil,
                                                double *dct_dp_frazil)
{
    double  bracket1, bracket2, cp_ih, gamma_ih, h, h_ih, part,
                rec_bracket3, tf, wcp, h_hat_sa, h_hat_ct, tf_sa, tf_p,
                ctf, ctf_sa, ctf_p;
    double  saturation_fraction = 0.0;

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
/**
==========================================================================
method: gsw_t_freezing_poly (sa, p, saturation_fraction)
==========================================================================
!
!  Calculates the in-situ temperature at which seawater freezes from a
!  computationally efficient polynomial.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
!               (ITS-90)
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_t_freezing_poly(double sa, double p, double saturation_fraction)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        double ctf, return_value;

    ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
    return_value = gsw_t_from_ct(sa,ctf,p);
        return (return_value);
}
/**
==========================================================================
elemental subroutine gsw_t_freezing_first_derivatives_poly (sa, p, &
                            saturation_fraction, tfreezing_sa, tfreezing_p)
==========================================================================
!
!  Calculates the first derivatives of the in-situ temperature at which
!  seawater freezes with respect to Absolute Salinity SA and pressure P (in
!  Pa).  These expressions come from differentiating the expression that
!  defines the freezing temperature, namely the equality between the
!  chemical potentials of water in seawater and in ice.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  tfreezing_SA = the derivative of the in-situ freezing temperature
!                 (ITS-90) with respect to Absolute Salinity at fixed
!                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  tfreezing_P  = the derivative of the in-situ freezing temperature
!                 (ITS-90) with respect to pressure (in Pa) at fixed
!                 Absolute Salinity                                [ K/Pa ]
!--------------------------------------------------------------------------
*/
void TeosBase::gsw_t_freezing_first_derivatives_poly(double sa, double p,
        double saturation_fraction, double *tfreezing_sa, double *tfreezing_p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */

        double  p_r, sa_r, x, c = 1e-3/(2.0*gtc.gsw_sso);

        sa_r = sa*1e-2;
        x = sqrt(sa_r);
        p_r = p*1e-4;

        if (tfreezing_sa != NULL)
            *tfreezing_sa =
            (gfpc.t1 + x*(1.5*gfpc.t2 + x*(2.0*gfpc.t3 + x*(2.5*gfpc.t4 + x*(3.0*gfpc.t5
                + 3.5*gfpc.t6*x)))) + p_r*(gfpc.t10 + x*(1.5*gfpc.t11 + x*(2.0*gfpc.t13
                + x*(2.5*gfpc.t16 + x*(3.0*gfpc.t19 + 3.5*gfpc.t22*x))))
                + p_r*(gfpc.t12 + x*(1.5*gfpc.t14 + x*(2.0*gfpc.t17 + 2.5*gfpc.t20*x))
                + p_r*(gfpc.t15 + x*(1.5*gfpc.t18 + 2.0*gfpc.t21*x)))))*1e-2
                + saturation_fraction*c;

        if (tfreezing_p != NULL)
            *tfreezing_p =
            (gfpc.t7 + sa_r*(gfpc.t10 + x*(gfpc.t11 + x*(gfpc.t13 + x*(gfpc.t16 + x*(gfpc.t19 + gfpc.t22*x)))))
                + p_r*(2.0*gfpc.t8 + sa_r*(2.0*gfpc.t12 + x*(2.0*gfpc.t14 + x*(2.0*gfpc.t17
                + 2.0*gfpc.t20*x))) + p_r*(3.0*gfpc.t9 + sa_r*(3.0*gfpc.t15 + x*(3.0*gfpc.t18
                + 3.0*gfpc.t21*x)))))*1e-8;

        return;
}
/**
!==========================================================================
elemental function gsw_t_freezing_poly (sa, p, saturation_fraction, polynomial)
!==========================================================================
!
!  Calculates the in-situ temperature at which seawater freezes from a
!  computationally efficient polynomial.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
!               (ITS-90)
!--------------------------------------------------------------------------
*/
double TeosBase::gsw_t_freezing_poly(double sa, double p, double saturation_fraction,
			int polynomial)
{
	/** for GSW_TEOS10_CONSTANTS use gtc */
	 /** for GSW_FREEZING_POLY_COEFFICIENTS use gfpc */

	double	p_r, sa_r, x, ctf, sfrac, return_value;
	int	direct_poly=polynomial;

	if (! direct_poly) {
	    sfrac = saturation_fraction;
	    ctf = gsw_ct_freezing_poly(sa,p,sfrac);
	    return_value = gsw_t_from_ct(sa,ctf,p);
	} else {
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
/**
==========================================================================
method gsw_sstar_from_sa(sa,p,lon,lat)
==========================================================================

! Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA.
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sa : Preformed Salinity                   [g/kg]
*/
double TeosBase::gsw_sstar_from_sa(double sa, double p, double lon, double lat)
{
        double  saar;

        saar    = gsw_saar(p,lon,lat);
    /**
        ! In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
        ! for saar in the Baltic.
    */
        if (saar == cppGSW_INVALID_VALUE)
        {
            return (saar);
        }

        return (sa*(1e0 - 0.35e0*saar)/(1e0 + saar));
}
/**
==========================================================================
method gsw_sa_from_sstar(sstar,p,lon,lat)
==========================================================================

! Calculates Absolute Salinity, SA, from Preformed Salinity, Sstar.
!
! Sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]
! lat    : latitude                                        [deg N]
!
! gsw_sa_from_sstar   : Absolute Salinity                  [g/kg]
*/
double TeosBase::gsw_sa_from_sstar(double sstar, double p, double lon, double lat)
{
        double  saar;

        saar    = gsw_saar(p,lon,lat);
        if (saar == cppGSW_INVALID_VALUE)
            return (saar);
    /**
    **! In the Baltic Sea, Sstar = SA, and note that gsw_saar returns zero
    **! for SAAR in the Baltic.
    */
        return (sstar*(1e0 + saar)/(1e0 - 0.35e0*saar));
}

/**
==========================================================================
method gsw_sstar_from_sp(sp,p,lon,lat)
==========================================================================

! Calculates Preformed Salinity, Sstar, from Practical Salinity, SP.
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon    : longitude                                       [deg E]
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sp  : Preformed Salinity                  [g/kg]
*/
double TeosBase::gsw_sstar_from_sp(double sp, double p, double lon, double lat)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        double  saar, sstar_baltic;

    /**
    !In the Baltic Sea, Sstar = SA.
    */
        sstar_baltic    = gsw_sa_from_sp_baltic(sp,lon,lat);

        if (sstar_baltic < cppGSW_ERROR_LIMIT)
        {
            return (sstar_baltic);
        }

        saar            = gsw_saar(p,lon,lat);

        if (saar == cppGSW_INVALID_VALUE)
        {
            return (saar);
        }

        return (gtc.gsw_ups*sp*(1 - 0.35e0*saar));
}
/**
==========================================================================
method: gsw_deltasa_atlas(p,lon,lat)
==========================================================================

 Calculates the Absolute Salinity Anomaly atlas value, delta_SA_atlas.

 p      : sea pressure                                    [dbar]
 lon    : longiture                                       [deg E]
 lat    : latitude                                        [deg N]

 deltasa_atlas : Absolute Salinity Anomaly atlas value    [g/kg]
*/
double TeosBase::gsw_deltasa_atlas(double p, double lon, double lat)
{
    /** for GSW_SAAR_DATA use gsaar */

	int	nx=gsw_nx, ny=gsw_ny, nz=gsw_nz;
	int	indx0, indy0, indz0, k, ndepth_index;

	double	dsar[4], dsar_old[4];
	double	dlong, dlat;
	double	return_value, sa_upper, sa_lower;
	double	r1, s1, t1, ndepth_max;

	return_value	= cppGSW_INVALID_VALUE;

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
        lon	+= 360.0;
    }

	dlong	= pSaarData->get_longs_ref(1)-pSaarData->get_longs_ref(0);
	dlat	= pSaarData->get_lats_ref(1)-pSaarData->get_lats_ref(0);

   /** these were set apart to resolve a C++ "problem"  calculating indx0 */
   double longsrefzero = pSaarData->get_longs_ref(0);
   double nxminusone = (double)nx-1.0;
   double longsrefnxminusone = pSaarData->get_longs_ref((unsigned)nxminusone);
   double lonminuslongsrefzero = lon-longsrefzero;

	indx0	= floor(nxminusone*(lonminuslongsrefzero/(longsrefnxminusone-longsrefzero)));

	if (indx0 == nx-1)
   {
      indx0	= nx-2;
   }

   /** these were set apart to resolve a C++ "problem" calculating indy0 */
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
         ndepth_max = gsw_max(ndepth_max, pSaarData->get_ndepth_ref(ndepth_index));
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
    Calculation of the derivatives is the key to the shape-preservation.
    There are other algorithms that could be used, but this appears to be
    the simplest, and adequate for our purposes.

    At minimal computational cost, we include here a check for increasing x.

    Returns 0 on success, 1 if x is not strictly increasing.
****************************************************************************/
int TeosBase::gsw_pchip_derivs(double *x, double *y, int n, double *d)
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
        /* change of sign, or either slope is zero */
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
double TeosBase::gsw_pchip_edge_case(double h0, double h1, double m0, double m1)
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
/***************************************************************************
    method: gsw_rr68_interp_section (sectnum, ip_sect, ip_isect, sa_i, ct_i)
****************************************************************************/
void TeosBase::gsw_rr68_interp_section(int sectnum, double *sa, double *ct, double *p, int mp,
        int nsect, double *ip_sect, int *ip_isect, double *p_i, double *sa_i,
        double *ct_i)
{
   int     i, *ip_1, *ip_2, *ip_3, *ip_4;
   double  m, *ct_12, *ct_13, *ct_23, *ct_34, ctp1,
                ctp2, *ct_ref, ctref_denom,
                ct_ref_minus_ctp1, ct_ref_minus_ctp2,
                ctref_num,
                gamma1_23, gamma1_24, gamma2_31,
                gamma2_34, gamma2_41, gamma3_12,
                gamma3_42, gamma4_12, gamma4_23,
                *sa_12, *sa_13, *sa_23, *sa_34, sap1,
                sap2, *sa_ref, saref_denom,
                sa_ref_minus_sap1, sa_ref_minus_sap2,
                saref_num, *p_ii;

   ip_1 = (int *) malloc(4*nsect*sizeof (int)); ip_2 = ip_1+nsect;
   ip_3 = ip_2+nsect; ip_4 = ip_3+nsect;

   ct_12 = (double *) malloc(12*nsect*sizeof (double));
   sa_12   = ct_12 +  1*nsect;
   sa_13   = ct_12 +  2*nsect;
   sa_23   = ct_12 +  3*nsect;
   sa_34   = ct_12 +  4*nsect;
   sa_ref  = ct_12 +  5*nsect;
   ct_13   = ct_12 +  6*nsect;
   ct_23   = ct_12 +  7*nsect;
   ct_34   = ct_12 +  8*nsect;
   ct_ref  = ct_12 +  9*nsect;
   p_ii    = ct_12 + 10*nsect;

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
   ! Construct the Reiniger & Ross reference curve equation.
   ! m = the power variable
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

         /** eqn (3a') */
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

      /** eqn (3b/3b') */
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

   free(ct_12); free(ip_1);
}
/**
==========================================================================
method: gsw_linear_interp_sa_ct (sa, ct, p, p_i, sa_i, ct_i)
==========================================================================
! This method interpolates the cast with respect to the interpolating
! variable p. This method finds the values of SA, CT at p_i on this cast.
!
! VERSION NUMBER: 3.05 (27th January 2015)
!
! This method was adapted from Matlab's interp1q.
==========================================================================
*/
void TeosBase::gsw_linear_interp_sa_ct(double *sa, double *ct, double *p, int np,
        double *p_i, int npi, double *sa_i, double *ct_i)
{
        char    *in_rng;
        int     *j, *k, *r, *jrev, *ki, imax_p, imin_p, i, n, m, ii;
        double  *xi, *xxi, u, max_p, min_p;

        min_p = max_p = p[0];
        imin_p = imax_p = 0;
        for (i=1; i<np; i++) {
            if (p[i] < min_p) {
                min_p = p[i];
                imin_p = i;
            } else if (p[i] > max_p) {
                max_p = p[i];
                imax_p = i;
            }
        }
        in_rng = (char *) malloc(npi*sizeof (char));
        memset(in_rng, 0, npi*sizeof (char));
        for (i=n=0; i<npi; i++) {
            if (p_i[i] <= min_p) {
                sa_i[i] = sa[imin_p];
                ct_i[i] = ct[imin_p];
            } else if (p_i[i] >= max_p) {
                sa_i[i] = sa[imax_p];
                ct_i[i] = ct[imax_p];
            } else {
                in_rng[i] = 1;
                n++;
            }
        }
        if (n==0)
            return;

        xi =(double *) malloc(n*sizeof (double));
        k  = (int *) malloc(3*n*sizeof (int)); ki = k+n; r = ki+n;
        m  = np + n;
        xxi = (double *) malloc(m*sizeof (double));
        j = (int *) malloc(2*m*sizeof (int)); jrev = j+m;

        ii = 0;
        for (i = 0; i<npi; i++) {
            if (in_rng[i]) {
                xi[ii] = p_i[i];
                ki[ii] = i;
                ii++;
            }
        }
        free(in_rng);
    /**
    **  Note that the following operations on the index
    **  vectors jrev and r depend on the sort utility
    **  gsw_util_sort_real() consistently ordering the
    **  sorting indexes either in ascending or descending
    **  sequence for replicate values in the real vector.
    */
        gsw_util_sort_real(xi, n, k);
        for (i = 0; i<np; i++)
            xxi[i] = p[i];
        for (i = 0; i<n; i++)
            xxi[np+i] = xi[k[i]];
        gsw_util_sort_real(xxi, m, j);

        for (i = 0; i<m; i++)
            jrev[j[i]] = i;
        for (i = 0; i<n; i++)
            r[k[i]] = jrev[np+i]-i-1;

        for (i = 0; i<n; i++) {
            u = (xi[i]-p[r[i]])/(p[r[i]+1]-p[r[i]]);
            sa_i[ki[i]] = sa[r[i]] + (sa[r[i]+1]-sa[r[i]])*u;
            ct_i[ki[i]] = ct[r[i]] + (ct[r[i]+1]-ct[r[i]])*u;
        }
        free(j); free(xxi); free(k); free(xi);
}


/** * Rank two items, by value if possible,
 * and by inverse index, if the values are
 * equal.
 * FIXME: decide if index method matches docs.
 */
int compareDI(const void *a, const void *b)
{
   DI *A = (DI*)a;
   DI *B = (DI*)b;

   if (A->d < B->d)
         return (-1);

   if (A->d > B->d)
         return (1);

   if (A->i < B->i)
         return (1);

   return (-1);
}
/**pure method gsw_util_sort_real (rarray) result(iarray)*/

/****
**  Sort the double array rarray into ascending value sequence
**  returning an index array of the sorted result.  This method
**  is thread-safe.
*/
void TeosBase::gsw_util_sort_real(double *rarray, int nx, int *iarray)
{
   int i;
   DI* di = (DI*)malloc(nx*sizeof(DI));

   for (i=0; i<nx; i++) {
      di[i].d = rarray[i];
      di[i].i = i;
   }

   qsort(di, nx, sizeof(DI), compareDI);

   for (i=0; i<nx; i++)
      iarray[i] = di[i].i;

   free(di);
}

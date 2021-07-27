#include "TeosSea.h"

/******************************************
  "TeosSea.cpp" Version 1.0
  by Randall Kent Whited
  rkwhited@gmail.com
  --------------------------------
  This is one modification of the TEOS-10
  file "gsw_oceanographic_toolbox.c"
  as adapted from the TEOS-10
  C version 3.05 http://www.teos-10.org
  ---------------------------------------
  All copyrights and all license issues
  are the same as for the C version
*******************************************/

/*****************************************
  Most C function calls in the original
  TEOS-10 "gsw_oceanographic_toolbox.c"
  and "gsw_saar.c" source code file that
  have been used are in the base-class
  "TeosBase.cpp" or in the descendant
  classes "TeosIce.cpp" or "TeosSea.cpp".
  ---------------------------------------
  Methods used by both "TeosIce" and
  "TeosSea" are implemented in the
  "TeosBase" class.
*****************************************/

/** descendant class constructor : base class constructor*/
TeosSea::TeosSea(): TeosBase()
{
}

/** descendant class destructor */
TeosSea::~TeosSea()
{
}



/**
==========================================================================
method gsw_adiabatic_lapse_rate_from_ct(sa,ct,p)
==========================================================================

! Calculates the adiabatic lapse rate from Conservative Temperature
!
! sa     : Absolute Salinity                                 [g/kg]
! ct     : Conservative Temperature                          [deg C]
! p      : sea pressure                                      [dbar]
!
! gsw_adiabatic_lapse_rate_from_ct : adiabatic lapse rate    [K/Pa]
*/
double TeosSea::gsw_adiabatic_lapse_rate_from_ct(double sa, double ct, double p)
{
        int     n0=0, n1=1, n2=2;
        double  pt0, pr0=0.0, t;

        pt0     = gsw_pt_from_ct(sa,ct);
        t       = gsw_pt_from_t(sa,pt0,pr0,p);

        return (-gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n2,n0,sa,t,p));

}
/**
==========================================================================
method gsw_cp_t_exact(sa,t,p)
==========================================================================

! Calculates isobaric heat capacity of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_cp_t_exact : heat capacity                           [J/(kg K)]
*/
double TeosSea::gsw_cp_t_exact(double sa, double t, double p)
{
        int     n0, n2;

        n0 = 0;
        n2 = 2;

        return (-(t+273.15e0)*gsw_gibbs(n0,n2,n0,sa,t,p));
}
/**
=========================================================================
method: gsw_ct_from_entropy (sa, entropy)
=========================================================================
!
!  Calculates Conservative Temperature with entropy as an input variable.
!
!  SA       =  Absolute Salinity                                   [ g/kg ]
!  entropy  =  specific entropy                                   [ deg C ]
!
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!--------------------------------------------------------------------------
*/
double TeosSea::gsw_ct_from_entropy(double sa, double entropy)
{
        double  pt;

        pt = gsw_pt_from_entropy(sa,entropy);
        return (gsw_ct_from_pt(sa,pt));
}
/**
=========================================================================
method: gsw_pt_from_entropy (sa, entropy)
=========================================================================
!
!  Calculates potential temperature with reference pressure p_ref = 0 dbar
!  and with entropy as an input variable.
!
!  SA       =  Absolute Salinity                                   [ g/kg ]
!  entropy  =  specific entropy                                   [ deg C ]
!
!  pt   =  potential temperature                                  [ deg C ]
!          with reference sea pressure (p_ref) = 0 dbar.
!  Note. The reference sea pressure of the output, pt, is zero dbar.
!--------------------------------------------------------------------------
*/
double TeosSea::gsw_pt_from_entropy(double sa, double entropy)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        int     number_of_iterations;
        double  c, dentropy, dentropy_dt, ent_sa, part1, part2, pt, ptm,
                pt_old;

        /**Find the initial value of pt*/
        part1 = 1.0 - sa/gtc.gsw_sso;
        part2 = 1.0 - 0.05*part1;
        ent_sa = (gtc.gsw_cp0/gtc.gsw_t0)*part1*(1.0 - 1.01*part1);
        c = (entropy - ent_sa)*(part2/gtc.gsw_cp0);
        pt = gtc.gsw_t0*(exp(c) - 1.0);
        dentropy_dt = gtc.gsw_cp0/((gtc.gsw_t0 + pt)*part2);

        for (number_of_iterations = 1; number_of_iterations <= 2;
            number_of_iterations++) {
            pt_old = pt;
            dentropy = gsw_entropy_from_pt(sa,pt_old) - entropy;
            pt = pt_old - dentropy/dentropy_dt;
            ptm = 0.5*(pt + pt_old);
            dentropy_dt = -gsw_gibbs_pt0_pt0(sa,ptm);
            pt = pt_old - dentropy/dentropy_dt;
        }
        /**
        ! Maximum error of 2.2x10^-6 degrees C for one iteration.
        ! Maximum error is 1.4x10^-14 degrees C for two iterations
        ! (two iterations is the default, "for Number_of_iterations = 1:2").
        */
        return (pt);
}
/**
==========================================================================
method: gsw_entropy_from_pt (sa, pt)
==========================================================================
!
!  Calculates specific entropy of seawater.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]
!
!  entropy  =  specific entropy                                [ J/(kg*K) ]
!--------------------------------------------------------------------------
*/
double TeosSea::gsw_entropy_from_pt(double sa, double pt)
{
        int     n0 = 0, n1 = 1;
        double  pr0 = 0.0;

        return (-gsw_gibbs(n0,n1,n0,sa,pt,pr0));
}
/**
==========================================================================
method gsw_entropy_from_t(sa,t,p)
==========================================================================

! Calculates the specific entropy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_entropy_from_t : specific entropy                    [J/(kg K)]
*/
double TeosSea::gsw_entropy_from_t(double sa, double t, double p)
{
        int     n0=0, n1=1;

        return (-gsw_gibbs(n0,n1,n0,sa,t,p));

}
/**
==========================================================================
elemental subroutine gsw_pt_first_derivatives (sa, ct, pt_sa, pt_ct)
! =========================================================================
!
!  Calculates the following two partial derivatives of potential temperature
!  (the regular potential temperature whose reference sea pressure is 0 dbar)
!  (1) pt_SA, the derivative with respect to Absolute Salinity at
!       constant Conservative Temperature, and
!  (2) pt_CT, the derivative with respect to Conservative Temperature at
!       constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  pt_SA =  The derivative of potential temperature with respect to
!           Absolute Salinity at constant Conservative Temperature.
!                                                               [ K/(g/kg)]
!  pt_CT =  The derivative of potential temperature with respect to
!           Conservative Temperature at constant Absolute Salinity.
!           pt_CT is dimensionless.                            [ unitless ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_pt_first_derivatives (double sa, double ct, double *pt_sa, double *pt_ct)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        double  abs_pt, ct_pt, ct_sa, pt, pr0 = 0.0;
        int     n0=0, n1=1, n2=2;

        pt = gsw_pt_from_ct(sa,ct);
        abs_pt = (gtc.gsw_t0 + pt);

        ct_pt = -(abs_pt*gsw_gibbs(n0,n2,n0,sa,pt,pr0))/gtc.gsw_cp0;

        if (pt_sa != NULL) {

            ct_sa = (gsw_gibbs(n1,n0,n0,sa,pt,pr0) -
                abs_pt*gsw_gibbs(n1,n1,n0,sa,pt,pr0))/gtc.gsw_cp0;

            *pt_sa = -ct_sa/ct_pt;

        }

        if (pt_ct != NULL)
            *pt_ct = 1.0/ct_pt;
}
/**
==========================================================================
method gsw_pt_from_ct(sa,ct)
==========================================================================

! potential temperature of seawater from conservative temperature
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_pt_from_ct : potential temperature with              [deg C]
!                  reference pressure of  0 dbar
*/
double TeosSea::gsw_pt_from_ct(double sa, double ct)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        double  a5ct, b3ct, ct_factor, pt_num, pt_recden, ct_diff;
        double  pt, pt_old, ptm, dpt_dct, s1;
        double  a0      = -1.446013646344788e-2,
                a1      = -3.305308995852924e-3,
                a2      =  1.062415929128982e-4,
                a3      =  9.477566673794488e-1,
                a4      =  2.166591947736613e-3,
                a5      =  3.828842955039902e-3,
                b0      =  1.000000000000000e0,
                b1      =  6.506097115635800e-4,
                b2      =  3.830289486850898e-3,
                b3      =  1.247811760368034e-6;

        s1      = sa/gtc.gsw_ups;

        a5ct    = a5*ct;
        b3ct    = b3*ct;

        ct_factor       = (a3 + a4*s1 + a5ct);
        pt_num          = a0 + s1*(a1 + a2*s1) + ct*ct_factor;
        pt_recden       = 1.0/(b0 + b1*s1 + ct*(b2 + b3ct));
        pt              = pt_num*pt_recden;

        dpt_dct = pt_recden*(ct_factor + a5ct - (b2 + b3ct + b3ct)*pt);

    /**
    **  Start the 1.5 iterations through the modified Newton-Raphson
    **  iterative method.
    */

        ct_diff = gsw_ct_from_pt(sa,pt) - ct;
        pt_old  = pt;
        pt      = pt_old - ct_diff*dpt_dct;
        ptm     = 0.5*(pt + pt_old);

        dpt_dct = -gtc.gsw_cp0/((ptm + gtc.gsw_t0)*gsw_gibbs_pt0_pt0(sa,ptm));

        pt      = pt_old - ct_diff*dpt_dct;
        ct_diff = gsw_ct_from_pt(sa,pt) - ct;
        pt_old  = pt;
        return (pt_old - ct_diff*dpt_dct);
}

/**
==========================================================================
elemental subroutine gsw_pt_second_derivatives (sa, ct, pt_sa_sa, &
                                                pt_sa_ct, pt_ct_ct)
! =========================================================================
!
!  Calculates the following three second-order derivatives of potential
!  temperature (the regular potential temperature which has a reference
!  sea pressure of 0 dbar),
!   (1) pt_SA_SA, the second derivative with respect to Absolute Salinity
!       at constant Conservative Temperature,
!   (2) pt_SA_CT, the derivative with respect to Conservative Temperature
!       and Absolute Salinity, and
!   (3) pt_CT_CT, the second derivative with respect to Conservative
!       Temperature at constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  pt_SA_SA  =  The second derivative of potential temperature (the
!               regular potential temperature which has reference sea
!               pressure of 0 dbar) with respect to Absolute Salinity
!               at constant Conservative Temperature.
!               pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
!  pt_SA_CT  =  The derivative of potential temperature with respect
!               to Absolute Salinity and Conservative Temperature.
!               pt_SA_CT has units of:                         [ 1/(g/kg) ]
!  pt_CT_CT  =  The second derivative of potential temperature (the
!               regular one with p_ref = 0 dbar) with respect to
!               Conservative Temperature at constant SA.
!               pt_CT_CT has units of:                              [ 1/K ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_pt_second_derivatives (double sa, double ct, double *pt_sa_sa,
        double *pt_sa_ct, double *pt_ct_ct)
{
        double  ct_l, ct_u, pt_ct_l, pt_ct_u, pt_sa_l, pt_sa_u, sa_l, sa_u,
                dct = 1e-2, dsa = 1e-3;

        if (pt_sa_sa != NULL) {

            if ((sa_l = sa - dsa) < 0.0)
                 sa_l = 0.0;
            sa_u = sa + dsa;

            gsw_pt_first_derivatives(sa_l,ct,&pt_sa_l,NULL);
            gsw_pt_first_derivatives(sa_u,ct,&pt_sa_u,NULL);

            *pt_sa_sa = (pt_sa_u - pt_sa_l)/(sa_u - sa_l);

        }

        if (pt_sa_ct != NULL || pt_ct_ct != NULL) {

            ct_l = ct - dct;
            ct_u = ct + dct;

            if ((pt_sa_ct != NULL) && (pt_ct_ct != NULL)) {

                gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,&pt_ct_l);
                gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,&pt_ct_u);

                *pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);
                *pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);

            } else if ((pt_sa_ct != NULL) && (pt_ct_ct == NULL)) {

                gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,NULL);
                gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,NULL);

                *pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);

            } else if ((pt_sa_ct == NULL) && (pt_ct_ct != NULL)) {

                gsw_pt_first_derivatives(sa,ct_l,NULL,&pt_ct_l);
                gsw_pt_first_derivatives(sa,ct_u,NULL,&pt_ct_u);

                *pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);
            }
        }
}
/**
==========================================================================
method gsw_pt_from_t(sa,t,p,p_ref)
==========================================================================

! Calculates potential temperature of seawater from in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
!
! gsw_pt_from_t : potential temperature                    [deg C]
*/
double TeosSea::gsw_pt_from_t(double sa, double t, double p, double p_ref)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        int     n0=0, n2=2, no_iter;
        double  s1, pt, ptm, pt_old, dentropy, dentropy_dt,
                true_entropy_part;

        s1      = sa/gtc.gsw_ups;
        pt      = t+(p-p_ref)*( 8.65483913395442e-6  -
                          s1 *  1.41636299744881e-6  -
                   (p+p_ref) *  7.38286467135737e-9  +
                          t  *(-8.38241357039698e-6  +
                          s1 *  2.83933368585534e-8  +
                          t  *  1.77803965218656e-8  +
                   (p+p_ref) *  1.71155619208233e-10));

        dentropy_dt     = gtc.gsw_cp0/((gtc.gsw_t0 + pt)*(1.0-0.05*(1.0 - sa/gtc.gsw_sso)));
        true_entropy_part       = gsw_entropy_part(sa,t,p);

        for (no_iter=1; no_iter <= 2; no_iter++) {
            pt_old      = pt;
            dentropy    = gsw_entropy_part(sa,pt_old,p_ref) - true_entropy_part;
            pt          = pt_old - dentropy/dentropy_dt;
            ptm         = 0.5*(pt + pt_old);
            dentropy_dt = -gsw_gibbs(n0,n2,n0,sa,ptm,p_ref);
            pt          = pt_old - dentropy/dentropy_dt;
        }
        return (pt);
}
/**
==========================================================================
elemental subroutine gsw_ct_second_derivatives (sa, pt, ct_sa_sa, ct_sa_pt, &
                                                ct_pt_pt)
==========================================================================
!
!  Calculates the following three, second-order derivatives of Conservative
!  Temperature
!   (1) CT_SA_SA, the second derivative with respect to Absolute Salinity
!       at constant potential temperature (with p_ref = 0 dbar),
!   (2) CT_SA_pt, the derivative with respect to potential temperature
!       (the regular potential temperature which is referenced to 0 dbar)
!       and Absolute Salinity, and
!   (3) CT_pt_pt, the second derivative with respect to potential
!       temperature (the regular potential temperature which is referenced
!       to 0 dbar) at constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]
!         (whose reference pressure is 0 dbar)
!
!  CT_SA_SA  =  The second derivative of Conservative Temperature with
!               respect to Absolute Salinity at constant potential
!               temperature (the regular potential temperature which
!               has reference sea pressure of 0 dbar).
!               CT_SA_SA has units of:                     [ K/((g/kg)^2) ]
!  CT_SA_pt  =  The derivative of Conservative Temperature with
!               respect to potential temperature (the regular one with
!               p_ref = 0 dbar) and Absolute Salinity.
!               CT_SA_pt has units of:                        [ 1/(g/kg) ]
!  CT_pt_pt  =  The second derivative of Conservative Temperature with
!               respect to potential temperature (the regular one with
!               p_ref = 0 dbar) at constant SA.
!               CT_pt_pt has units of:                              [ 1/K ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_ct_second_derivatives(double sa, double pt, double *ct_sa_sa,
        double *ct_sa_pt, double *ct_pt_pt)
{
        double  ct_pt_l, ct_pt_u, ct_sa_l, ct_sa_u, pt_l, pt_u, sa_l, sa_u,
                dsa = 1e-3, dpt = 1e-2;

        if ((ct_sa_sa != NULL)) {

            if ((sa_l = sa - dsa) < 0.0)
                sa_l = 0.0;
            sa_u = sa + dsa;

            gsw_ct_first_derivatives(sa_l,pt,&ct_sa_l,NULL);
            gsw_ct_first_derivatives(sa_u,pt,&ct_sa_u,NULL);

            *ct_sa_sa = (ct_sa_u - ct_sa_l)/(sa_u - sa_l);

        }

        if ((ct_sa_pt != NULL) || (ct_pt_pt != NULL)) {

            pt_l = pt - dpt;
            pt_u = pt + dpt;

            if ((ct_sa_pt != NULL) && (ct_pt_pt != NULL)) {

                gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,&ct_pt_l);
                gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,&ct_pt_u);

                *ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);
                *ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);

            } else if ((ct_sa_pt != NULL) && (ct_pt_pt == NULL)) {

                gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,NULL);
                gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,NULL);

                *ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);

            } else if ((ct_sa_pt == NULL) && (ct_pt_pt != NULL)) {

                gsw_ct_first_derivatives(sa,pt_l,NULL,&ct_pt_l);
                gsw_ct_first_derivatives(sa,pt_u,NULL,&ct_pt_u);

                *ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);

            }
        }
}
/**
==========================================================================
elemental subroutine gsw_ct_first_derivatives (sa, pt, ct_sa, ct_pt)
==========================================================================
!
!  Calculates the following two derivatives of Conservative Temperature
!  (1) CT_SA, the derivative with respect to Absolute Salinity at
!      constant potential temperature (with pr = 0 dbar), and
!   2) CT_pt, the derivative with respect to potential temperature
!      (the regular potential temperature which is referenced to 0 dbar)
!      at constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]
!         (whose reference pressure is 0 dbar)
!
!  CT_SA  =  The derivative of Conservative Temperature with respect to
!            Absolute Salinity at constant potential temperature
!            (the regular potential temperature which has reference
!            sea pressure of 0 dbar).
!            The CT_SA output has units of:                     [ K/(g/kg)]
!  CT_pt  =  The derivative of Conservative Temperature with respect to
!            potential temperature (the regular one with pr = 0 dbar)
!            at constant SA. CT_pt is dimensionless.           [ unitless ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_ct_first_derivatives(double sa, double pt, double *ct_sa, double *ct_pt)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        double  abs_pt, g_sa_mod, g_sa_t_mod, x, y_pt;

        abs_pt = gtc.gsw_t0 + pt ;

        if (ct_pt != NULL)
            *ct_pt = -(abs_pt*gsw_gibbs_pt0_pt0(sa,pt))/gtc.gsw_cp0;

        if (ct_sa == NULL)
            return;

        x = sqrt(gtc.gsw_sfac*sa);
        y_pt = 0.025*pt;

        g_sa_t_mod = 1187.3715515697959 + x*(-1480.222530425046
            + x*(2175.341332000392 + x*(-980.14153344888
            + 220.542973797483*x) + y_pt*(-548.4580073635929
            + y_pt*(592.4012338275047 + y_pt*(-274.2361238716608
            + 49.9394019139016*y_pt)))) + y_pt*(-258.3988055868252
            + y_pt*(-90.2046337756875 + y_pt*10.50720794170734)))
            + y_pt*(3520.125411988816  + y_pt*(-1351.605895580406
            + y_pt*(731.4083582010072  + y_pt*(-216.60324087531103
            + 25.56203650166196*y_pt))));
        g_sa_t_mod = 0.5*gtc.gsw_sfac*0.025*g_sa_t_mod;

        g_sa_mod = 8645.36753595126 + x*(-7296.43987145382
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
/**
==========================================================================
method gsw_latentheat_evap_t(sa,t)
==========================================================================
!
! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
!
! gsw_latentheat_evap_t : latent heat of evaporation       [J/kg]
*/
double TeosSea::gsw_latentheat_evap_t(double sa, double t)
{

        double  ct = gsw_ct_from_pt(sa,t);

        return (gsw_latentheat_evap_ct(sa,ct));
}
/**
==========================================================================
method gsw_latentheat_evap_ct(sa,ct)
==========================================================================

! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
!
! latentheat_evaporation : latent heat of evaporation      [J/kg]
*/
double TeosSea::gsw_latentheat_evap_ct(double sa, double ct)
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

        double  x, y;

        x       = sqrt(gtc.gsw_sfac*sa);
        y       = ct/40.0;

        return (c0 + x*(c1 + c4*y + x*(c3
                   + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))
                   + x*(c10 + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))
                   + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)
                   + y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y))))));
}

/**
    ** The following method is buggy; users are advised to use
    ** gsw_geo_strf_dyn_height_1 instead.

==========================================================================
pure method gsw_geo_strf_dyn_height (sa, ct, p, p_ref)
==========================================================================
!
!  Calculates dynamic height anomaly as the integral of specific volume
!  anomaly from the pressure p of the bottle to the reference pressure
!  p_ref.
!
!  Hence, geo_strf_dyn_height is the dynamic height anomaly with respect
!  to a given reference pressure.  This is the geostrophic streammethod
!  for the difference between the horizontal velocity at the pressure
!  concerned, p, and the horizontal velocity at p_ref.  Dynamic height
!  anomaly is the geostrophic streammethod in an isobaric surface.  The
!  reference values used for the specific volume anomaly are
!  SSO = 35.16504 g/kg and CT = 0 deg C.  This method calculates
!  specific volume anomaly using the computationally efficient
!  expression for specific volume of Roquet et al. (2015).
!
!  This method evaluates the pressure integral of specific volume using
!  SA and CT interpolated with respect to pressure using the method of
!  Reiniger and Ross (1968).  It uses a weighted mean of (i) values
!  obtained from linear interpolation of the two nearest data points, and
!  (ii) a linear extrapolation of the pairs of data above and below.  This
!  "curve fitting" method resembles the use of cubic splines.
!
!  SA    =  Absolute Salinity                                      [ g/kg ]
!  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
!  p     =  sea pressure                                           [ dbar ]
!           ( i.e. absolute pressure - 10.1325 dbar )
!  p_ref =  reference pressure                                     [ dbar ]
!           ( i.e. reference absolute pressure - 10.1325 dbar )
!
!  geo_strf_dyn_height  =  dynamic height anomaly               [ m^2/s^2 ]
!   Note. If p_ref exceeds the pressure of the deepest bottle on a
!     vertical profile, the dynamic height anomaly for each bottle
!     on the whole vertical profile is returned as NaN.
!--------------------------------------------------------------------------
*/
double  *TeosSea::gsw_geo_strf_dyn_height(double *sa, double *ct, double *p, double p_ref,
        int n_levels, double *dyn_height)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        int     m_levels = (n_levels <= 0) ? 1 : n_levels,
                p_cnt, top_pad, i, nz, ibottle, ipref, np_max, np, ibpr=0,
                *iidata;
        double  dp_min, dp_max, p_min, p_max, max_dp_i,
                *b, *b_av, *dp, *dp_i, *sa_i=NULL, *ct_i, *p_i=NULL,
                *geo_strf_dyn_height0;

/**
--------------------------------------------------------------------------
!  This max_dp_i is the limit we choose for the evaluation of specific
!  volume in the pressure integration.  That is, the vertical integration
!  of specific volume with respect to pressure is perfomed with the pressure
!  increment being no more than max_dp_i (the default value being 1 dbar).
!--------------------------------------------------------------------------
*/
        max_dp_i = 1.0;

        if ((nz = m_levels) <= 1)
            return (NULL);

        dp = (double *) malloc(nz*sizeof (double));
        dp_min = 11000.0;
        dp_max = -11000.0;
        for (i=0; i<nz-1; i++) {
            if ((dp[i] = p[i+1] - p[i]) < dp_min)
                dp_min = dp[i];
            if (dp[i] > dp_max)
                dp_max = dp[i];
        }

        if (dp_min <= 0.0) {
            /** pressure must be monotonic */
            free(dp);
            return (NULL);
        }
        p_min = p[0];
        p_max = p[nz-1];

        if (p_ref > p_max) {
            /**the reference pressure p_ref is deeper than all bottles*/
            free(dp);
            return (NULL);
        }

        /** Determine if there is a "bottle" at exactly p_ref */
        ipref = -1;
        for (ibottle = 0; ibottle < nz; ibottle++) {
            if (p[ibottle] == p_ref) {
                ipref = ibottle;
                break;
            }
        }
        if ((dp_max <= max_dp_i) && (p[0] == 0.0) && (ipref >= 0)) {
            /**
            !vertical resolution is good (bottle gap is no larger than max_dp_i)
            ! & the vertical profile begins at the surface (i.e. at p = 0 dbar)
            ! & the profile contains a "bottle" at exactly p_ref.
            */
            b = (double *) malloc(3*nz*sizeof (double));
            b_av = b+nz; geo_strf_dyn_height0 = b_av+nz;
            for (i=0; i<nz; i++) {
                b[i] = gsw_specvol_anom_standard(sa[i],ct[i],p[i]);
                if (i > 0)
                    b_av[i-1] = 0.5*(b[i] + b[i-1]);
            }
            /**
            ! "geo_strf_dyn_height0" is the dynamic height anomaly with respect
            ! to p_ref = 0 (the surface).
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
            free(b);
        } else {
        /**
        ! Test if there are vertical gaps between adjacent "bottles" which are
        ! greater than max_dp_i, and that there is a "bottle" exactly at the
        ! reference pressure.
        */
            iidata = (int *) malloc((nz+1)*sizeof (int));

            if ((dp_max <= max_dp_i) && (ipref >= 0)) {
            /**
            ! Vertical resolution is already good (no larger than max_dp_i), and
            ! there is a "bottle" at exactly p_ref.
            */
                sa_i = (double *) malloc(2*(nz+1)*sizeof (double));
                ct_i = sa_i+nz+1;
                p_i = (double *) malloc((nz+1)*sizeof (double));;

                if (p_min > 0.0) {
                /**
                ! resolution is fine and there is a bottle at p_ref, but
                ! there is not a bottle at p = 0. So add an extra bottle.
                */
                    for (i=0; i<nz; i++) {
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
                } else {
                /**
                ! resolution is fine, there is a bottle at p_ref, and
                ! there is a bottle at p = 0
                */
                    memmove(sa_i, sa, nz*sizeof (double));
                    memmove(ct_i, ct, nz*sizeof (double));
                    memmove(p_i, p, nz*sizeof (double));
                    ibpr = ipref;
                    for (i=0; i<nz; i++)
                        iidata[i] = i;
                    p_cnt = nz;
                }

            } else {
            /**
            ! interpolation is needed.
            */
                np_max = 2*rint(p[nz-1]/max_dp_i+0.5);
                p_i = (double *) malloc(np_max*sizeof (double));
                /** sa_i is allocated below, when its size is known */

                if (p_min > 0.0) {
                /**
                ! there is not a bottle at p = 0.
                */
                    if (p_ref < p_min) {
                    /**
                    ! p_ref is shallower than the minimum bottle pressure.
                    */
                        p_i[0] = 0.0;
                        gsw_p_sequence(p_i[0],p_ref,max_dp_i, p_i+1,&np);
                        ibpr = p_cnt = np;
                        p_cnt++;
                        gsw_p_sequence(p_ref,p_min,max_dp_i, p_i+p_cnt,&np);
                        p_cnt += np;
                        top_pad = p_cnt;
                    } else {
                    /**
                    ! p_ref is deeper than the minimum bottle pressure.
                    */
                        p_i[0] = 0.0;
                        p_i[1] = p_min;
                        top_pad = 2;
                        p_cnt = 2;
                    }
                } else {
                /**
                ! there is a bottle at p = 0.
                */
                    p_i[0] = p_min;
                    top_pad = 1;
                    p_cnt = 1;
                }

                for (ibottle=0; ibottle < nz-1; ibottle++) {

                    iidata[ibottle] = p_cnt-1;
                    if (p[ibottle] == p_ref) ibpr = p_cnt-1;

                    if (p[ibottle] < p_ref && p[ibottle+1] > p_ref) {
                    /**
                    ! ... reference pressure is spanned by bottle pairs -
                    ! need to include p_ref as an interpolated pressure.
                    */
                        gsw_p_sequence(p[ibottle],p_ref,max_dp_i, p_i+p_cnt,&np);
                        p_cnt += np;
                        ibpr = p_cnt-1;
                        gsw_p_sequence(p_ref,p[ibottle+1],max_dp_i,p_i+p_cnt,&np);
                        p_cnt += np;
                    } else {
                    /**
                    ! ... reference pressure is not spanned by bottle pairs.
                    */
                        gsw_p_sequence(p[ibottle],p[ibottle+1],max_dp_i,
                                p_i+p_cnt,&np);
                        p_cnt += np;
                    }

                }

                iidata[nz-1] = p_cnt-1;
                if (p[nz-1] == p_ref) ibpr = p_cnt-1;

                sa_i = (double *) malloc(2*p_cnt*sizeof (double));
                ct_i = sa_i+p_cnt;

                if (top_pad > 1) {
                    gsw_linear_interp_sa_ct(sa,ct,p,nz,
                        p_i,top_pad-1,sa_i,ct_i);
                }
                gsw_rr68_interp_sa_ct(sa,ct,p,nz,p_i+top_pad-1,p_cnt-top_pad+1,
                                      sa_i+top_pad-1,ct_i+top_pad-1);
            }

            b = (double *) malloc(4*p_cnt*sizeof (double));
            b_av = b+p_cnt; dp_i = b_av+p_cnt;
            geo_strf_dyn_height0 = dp_i+p_cnt;
            for (i=0; i<p_cnt; i++) {
                b[i] = gsw_specvol_anom_standard(sa_i[i],ct_i[i],p_i[i]);
                if (i > 0) {
                    dp_i[i-1] = p_i[i]-p_i[i-1];
                    b_av[i-1] = 0.5*(b[i] + b[i-1]);
                }
            }
            /**
            ! "geo_strf_dyn_height0" is the dynamic height anomaly with respect
            ! to p_ref = 0 (the surface).
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

            free(b);
            free(iidata);
            if (sa_i != NULL)
                free(sa_i);
            if (p_i != NULL)
                free(p_i);

        }
        free(dp);
        return (dyn_height);
}


void TeosSea::gsw_p_sequence(double p1, double p2, double max_dp_i, double *pseq, int *nps)
{
        double  dp, pstep;
        int             n, i;

        dp = p2 - p1;
        n = ceil(dp/max_dp_i);
        pstep = dp/n;

        if (nps != NULL) *nps = n;
        /**        ! Generate the sequence ensuring that the value of p2 is exact to
        ! avoid round-off issues, ie. don't do "pseq = p1+pstep*(i+1)".
        */
        for (i=0; i<n; i++)
            pseq[i] = p2-pstep*(n-1-i);

}
/**    ** This is a replacement for gsw_geo_strf_dyn_height, with a different
    ** signature and interpolation algorithms.

==========================================================================
int   (returns nonzero on error, 0 if OK)
gsw_geo_strf_dyn_height_1(double *sa, double *ct, double *p, double p_ref,
    int nz, double *dyn_height, double max_dp_i, int interp_method)
==========================================================================
!
!  Calculates dynamic height anomaly as the integral of specific volume
!  anomaly from the pressure p of the bottle to the reference pressure
!  p_ref.
!
!  Hence, geo_strf_dyn_height is the dynamic height anomaly with respect
!  to a given reference pressure.  This is the geostrophic streammethod
!  for the difference between the horizontal velocity at the pressure
!  concerned, p, and the horizontal velocity at p_ref.  Dynamic height
!  anomaly is the geostrophic streammethod in an isobaric surface.  The
!  reference values used for the specific volume anomaly are
!  SSO = 35.16504 g/kg and CT = 0 deg C.  This method calculates
!  specific volume anomaly using the computationally efficient
!  expression for specific volume of Roquet et al. (2015).
!
!  This method evaluates the pressure integral of specific volume using
!  SA and CT interpolated with respect to pressure. The interpolation method
!  may be chosen as linear or "PCHIP", piecewise cubic Hermite using a shape-
!  preserving algorithm for setting the derivatives.
!
!  SA    =  Absolute Salinity                                      [ g/kg ]
!  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
!  p     =  sea pressure  (increasing with index)                  [ dbar ]
!           ( i.e. absolute pressure - 10.1325 dbar )
!  nz    =  number of points in each array
!  p_ref =  reference pressure                                     [ dbar ]
!           ( i.e. reference absolute pressure - 10.1325 dbar )
!  geo_strf_dyn_height  =  dynamic height anomaly               [ m^2/s^2 ]
!  max_dp_i = maximum pressure difference between points for triggering
!              interpolation.
!  interp_method = 1 for linear, 2 for PCHIP
!
!   Note. If p_ref falls outside the range of a
!     vertical profile, the dynamic height anomaly for each bottle
!     on the whole vertical profile is returned as NaN.
!--------------------------------------------------------------------------
*/

/**    Make a new grid based on an original monotonic array, typically P, such that
    on the new monotonic grid, p_i:

        1) The first element is p[0], the last is p[nz-1].
        2) Approximate integer multiples of dp are included.
        3) All original p points are included.
        4) The value p_ref is included.

    Arguments:
        p : the original grid array
        p_ref : scalar reference pressure
        nz : size of p
        dp : target p interval for the output
        p_i : an array to hold the regridded pressures
        ni_max : size of p_i
        p_indices : an array of size nz
        p_ref_ind_ptr : pointer to an integer

    The method fills as much of p_i as it needs.  It returns with
    the array of indices (p_indices) of the original p grid points in
    p_i, and with a pointer to the index of p_ref within p_i.

    Its return argument is the number of elements used in p_i, or -1 on error.
*/
int TeosSea::gsw_refine_grid_for_dh(double *p, double p_ref, int nz,
    double dp,
    double *p_i, int ni_max,  /** size of p_i array; larger than needed */
    int *p_indices, int *p_ref_ind_ptr)
{
    int i, iuniform, iorig;
    double p_next;
    /** Don't add a new point if it is within p_tol of an original. */
    double p_tol = 0.001 * dp;

    p_i[0] = p[0];
    p_indices[0] = 0;
    *p_ref_ind_ptr = -1;  /** initialize to a flag value */
    if (p_ref <= p[0] + p_tol) {
        *p_ref_ind_ptr = 0;
    }
    for (i=1, iuniform=1, iorig=1; i<ni_max && iorig<nz; i++) {
        /** Candidate insertion based on uniform grid: */
        p_next = p[0] + dp * iuniform;

        /** See if we need to insert p_ref: */
        if (*p_ref_ind_ptr == -1 && p_ref <= p_next && p_ref <= p[iorig]) {
            p_i[i] = p_ref;
            *p_ref_ind_ptr = i;
            if (p_ref == p[iorig]) {
                p_indices[iorig] = i;
                iorig++;
            }
            if (p_ref > p_next - p_tol) {
                iuniform++;
            }
            continue;
        }

        /** We did not insert p_ref, so insert either p_next or p[iorig]. */
        if (p_next < p[iorig] - p_tol) {
            p_i[i] = p_next;
            iuniform++;
        }
        else {
            p_i[i] = p[iorig];
            p_indices[iorig] = i;
            /** Skip this p_next if it is close to the point we just added. */
            if (p_next < p[iorig] + p_tol) {
                iuniform++;
            }
            iorig++;
        }
    }

    if (i == ni_max) {
        return (-1);  /** error! */
    }
    return (i);  /** number of elements in p_i */
}

/******************************************************************
    Linearly interpolate to the grid made by define_grid_for_dh.
    We take advantage of what we know about the grids: they match
    at the end points, and both are monotonic.
*******************************************************************/
int TeosSea::gsw_linear_interp_SA_CT_for_dh(double *sa, double *ct, double *p, int nz,
                                        double *p_i, int n_i, double *sa_i, double *ct_i)
{
    int i, ii;
    double pfac;

    sa_i[0] = sa[0];
    sa_i[n_i-1] = sa[nz-1];
    ct_i[0] = ct[0];
    ct_i[n_i-1] = ct[nz-1];
    i = 1;
    for (ii=1; ii<n_i-1; ii++) {
        /** Find the second point of the pair in the original grid that
           bracket the target.
        */
        while (p[i] < p_i[ii]) {
            i++;
            if (i == nz) {
                return -1;  /** error! */
            }
        }
        pfac = (p_i[ii] - p[i-1]) / (p[i] - p[i-1]);
        sa_i[ii] = sa[i-1] + pfac * (sa[i] - sa[i-1]);
        ct_i[ii] = ct[i-1] + pfac * (ct[i] - ct[i-1]);
    }
    return 0;
}

/*************************************************
    this method returns nonzero on error, 0 if OK
    see: gsw_geo_strf_dyn_height_pc
**************************************************/
int  TeosSea::gsw_geo_strf_dyn_height_1(double *sa, double *ct, double *p, double p_ref,
    int nz, double *dyn_height, double max_dp_i, int interp_method)
{
    /** for GSW_TEOS10_CONSTANTS use gtc */

    int i, ipref, *p_indices, n_i, ni_max, err;
    double    dp_min, dp_max, p_min, p_max,
        *b, *b_av, *dp, *sa_i, *ct_i, *p_i,
        *dh_i;
    double dh_ref;

    if (nz < 2)
        return (1);

    dp = (double *)malloc((nz-1) * sizeof(double));
    dp_min = 11000.0;
    dp_max = -11000.0;

    for (i=0; i<nz-1; i++) {
        dp[i] = p[i+1] - p[i];
        if (dp[i] < dp_min) {
             dp_min = dp[i];
        }
        if (dp[i] > dp_max) {
            dp_max = dp[i];
        }
    }

    if (dp_min <= 0.0)
    {
        /** pressure must be monotonic */
        free(dp);
        return (2);
    }
    p_min = p[0];
    p_max = p[nz-1];

    if (p_ref > p_max || p_ref < p_min)
    {
        /** Reference pressure must be within the data range. */
        free(dp);
        return (3);
    }

    /** Determine if there is a sample at exactly p_ref */
    ipref = -1;
    for (i = 0; i < nz; i++) {
        if (p[i] == p_ref) {
            ipref = i;
            break;
        }
    }

    if ((dp_max <= max_dp_i) && (ipref >= 0))
    {
        /**
          vertical resolution is good (bottle gap is no larger than max_dp_i)
        ! & the profile contains a "bottle" at exactly p_ref.
         */
        b = (double *)malloc(nz*sizeof (double));
        b_av = (double *)malloc((nz-1) * sizeof(double));
        for (i=0; i<nz; i++) {
            b[i] = gsw_specvol_anom_standard(sa[i],ct[i],p[i]);
        }
        for (i=0; i<(nz-1); i++) {
            b_av[i] = 0.5*(b[i+1] + b[i]);
        }
        /**
          First calculate dynamic height relative to the
          first (shallowest) depth.
        */
        dyn_height[0] = 0.0;
        for (i=1; i<nz; i++) {
            dyn_height[i] = dyn_height[i-1] - b_av[i-1]*dp[i-1]*gtc.db2pa;
        }
        /**
          Then subtract out the value at the reference pressure.
        */
        dh_ref = dyn_height[ipref];
        for (i=0; i<nz; i++) {
            dyn_height[i] -= dh_ref;
        }

        /** release memory */
        free(b);
        free(b_av);
        free(dp);

        return (0);
    }

    /**
      If we got this far, then we need to interpolate: either or both of
      inserting a point for p_ref and subdividing the intervals to keep the max
      interval less than max_dp_i.
    */

    free(dp);  /** Need to recalculate, so free here and malloc when needed. */

    ni_max = nz + (int) ceil((p[nz-1] - p[0]) / max_dp_i) + 2;
    /**
      Maximum possible size of new grid: Original grid size plus
      the number of dp intervals plus 1 for the p_ref,
      plus 1 so that we can know we exited the loop before we hit it.
    */

    p_i = (double *) malloc(ni_max * sizeof(double));
    p_indices = (int *) malloc(nz * sizeof(int));

    n_i = gsw_refine_grid_for_dh(p, p_ref, nz, max_dp_i,
                             p_i, ni_max,
                             p_indices, &ipref);

    /** allocation reminder: if successful, this allocated p_i and p_indices. */
    if (n_i == -1) {
        free(p_i);
        free(p_indices);
        return (4);
    }

    ct_i = (double *)malloc(n_i * sizeof(double));
    sa_i = (double *)malloc(n_i * sizeof(double));

    if (interp_method == cppINTERP_METHOD_LINEAR)
    {
        err = gsw_linear_interp_SA_CT_for_dh(sa, ct, p, nz,
                                         p_i, n_i,
                                         sa_i, ct_i);
        if (err) err = 5;
    }
    else if (interp_method == cppINTERP_METHOD_PCHIP) {
        err = gsw_util_pchip_interp(p, sa, nz, p_i, sa_i, n_i);
        err = err || gsw_util_pchip_interp(p, ct, nz, p_i, ct_i, n_i);
        if (err) err = 6;
    }
    else {
        err = 7;
    }

    if (err) /** free allocated memory */
    {
        free(p_i);
        free(p_indices);
        free(ct_i);
        free(sa_i);
        return (err);
    }

    dh_i = (double *)malloc(n_i * sizeof(double));
    dp = (double *)malloc((n_i-1) * sizeof(double));
    b = (double *)malloc(n_i*sizeof (double));
    b_av = (double *)malloc((n_i-1) * sizeof(double));

    for (i=0; i<n_i; i++) {
        b[i] = gsw_specvol_anom_standard(sa_i[i], ct_i[i], p_i[i]);
    }
    /** release allocated memory */
    free(ct_i);
    free(sa_i);

    for (i=0; i<(n_i-1); i++) {
        b_av[i] = 0.5*(b[i+1] + b[i]);
        dp[i] = p_i[i+1] - p_i[i];
    }

    free(p_i); /** free allocated memory */

    /**
      First calculate dynamic height relative to the
      first (shallowest) depth.
    */
    dh_i[0] = 0.0;
    for (i=1; i<n_i; i++) {
        dh_i[i] = dh_i[i-1] - b_av[i-1]*dp[i-1]*gtc.db2pa;
    }
    /** free allocated memory */
    free(b);
    free(b_av);
    free(dp);

    dh_ref = dh_i[ipref];

    for (i=0; i<nz; i++) {
        dyn_height[i] = dh_i[p_indices[i]] - dh_ref;
    }
    /** free allocated memory */
    free(p_indices);
    free(dh_i);

    return 0;
}
/**
==========================================================================
method: gsw_rr68_interp_sa_ct (sa, ct, p, p_i, sa_i, ct_i)
==========================================================================
!
!  Interpolate Absolute Salinity and Conservative Temperature values to
!  arbitrary pressures using the Reiniger and Ross (1968) interpolation
!  scheme.
!  Note that this interpolation scheme requires at least four observed
!  bottles on the cast.
!
!  SA   =  Absolute Salinity                                  [ g/kg ]
!  CT   =  Conservative Temperature (ITS-90)                 [ deg C ]
!  p    =  sea pressure                                       [ dbar ]
!           ( i.e. absolute pressure - 10.1325 dbar )
!  p_i  =  pressures to interpolate to.
!
!  SA_i = interpolated SA values at pressures p_i.
!  CT_i = interpolated CT values at pressures p_i.
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_rr68_interp_sa_ct(double *sa, double *ct, double *p, int mp, double *p_i,
        int mp_i, double *sa_i, double *ct_i)
{
        int     i, j, nshallow, ncentral, ndeep,
                *ip, *ip_i, *ip_ishallow, *ip_icentral, *ip_ideep;
        char    *shallow, *central, *deep;
        double  *ip_shallow, *ip_central, *ip_deep, *dp, *p_ii;

        if (mp < 4) {
            /** need at least four bottles to perform this interpolation */
            ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
            return;
        }

        dp = (double *) malloc(mp*sizeof (double));
        for (i=1; i<mp; i++) {
            if ((dp[i-1] = (p[i] - p[i-1])) <= 0.0) {
                free(dp);
                ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
                return;
            }
        }

        shallow = (char *) malloc(3*mp_i*sizeof (char));
        central = shallow+mp_i; deep = central+mp_i;
        nshallow=ncentral=ndeep=0;
        memset(shallow, 0, 3*mp_i*sizeof (char));
        for (i=0; i<mp_i; i++) {
            if (p_i[i] >= p[0] && p_i[i] <= p[1]) {
                nshallow++;
                shallow[i] = 1;
            }
            if (p_i[i] >= p[1] && p_i[i] <= p[mp-2]) {
                ncentral++;
                central[i] = 1;
            }
            if (p_i[i] >= p[mp-2] && p_i[i] <= p[mp-1]) {
                ndeep++;
                deep[i] = 1;
            }
        }

        if ((nshallow == 0) || (ncentral == 0) || (ndeep == 0)) {
            free(shallow); free(dp);
            ct_i[0] = sa_i[0] = cppGSW_INVALID_VALUE;
            return;
        }

        ip = (int *) malloc((mp+mp_i)*sizeof (int)); ip_i = ip+mp;
        for (i=0; i<mp; i++)
            ip[i] = i;
        for (i=0; i<mp_i; i++)
            ip_i[i] = i;

        ip_ishallow = (int *) malloc((nshallow+ncentral+ndeep)*sizeof (int));
        ip_icentral = ip_ishallow+nshallow; ip_ideep = ip_icentral+ncentral;
        ip_shallow = (double *) malloc(2*(nshallow+ncentral+ndeep)*sizeof (double));
        ip_central = ip_shallow+nshallow; ip_deep = ip_central+ncentral;
        p_ii = ip_deep+ndeep;
        /**
        ! Calculate the 2 outer extrapolated values and the inner
        ! interpolated values
        */
        for (i=j=0; i<mp_i; i++) {
            if (central[i]) {
                ip_icentral[j] = ip_i[i];
                j++;
            }
        }
        for (i=0; i<ncentral; i++)
            p_ii[i] = p_i[ip_icentral[i]];
        gsw_util_interp1q_int(mp,p,ip,ncentral,p_ii,ip_central);
        gsw_rr68_interp_section(0,sa,ct,p,mp,ncentral,ip_central,ip_icentral,
                                p_i,sa_i,ct_i);

        for (i=j=0; i<mp_i; i++) {
            if (shallow[i]) {
                ip_ishallow[j] = ip_i[i];
                j++;
            }
        }
        for (i=0; i<nshallow; i++)
            p_ii[i] = p_i[ip_ishallow[i]];
        gsw_util_interp1q_int(mp,p,ip,nshallow,p_ii,ip_shallow);
        gsw_rr68_interp_section(-1,sa,ct,p,mp,nshallow,ip_shallow,ip_ishallow,
                                p_i,sa_i,ct_i);

        for (i=j=0; i<mp_i; i++) {
            if (deep[i]) {
                ip_ideep[j] = ip_i[i];
                j++;
            }
        }
        for (i=0; i<ndeep; i++)
            p_ii[i] = p_i[ip_ideep[i]];
        gsw_util_interp1q_int(mp,p,ip,ndeep,p_ii,ip_deep);
        gsw_rr68_interp_section(1,sa,ct,p,mp,ndeep,ip_deep,ip_ideep,p_i,sa_i,ct_i);

        /**
        ! Insert any observed bottles that are at the required interpolated
        ! pressures
        */
        for (i=0; i<mp_i; i++) {
            for (j=0; j<mp; j++) {
                if (p_i[i] == p[j]) {
                    sa_i[i] = sa[j];
                    ct_i[i] = ct[j];
                }
            }
        }
        free(ip_shallow); free(ip_ishallow); free(ip); free(shallow); free(dp);
}

/**
==========================================================================
pure method gsw_util_interp1q_int (x, iy, x_i) result(y_i)
==========================================================================
! Returns the value of the 1-D method iy (integer) at the points of column
! vector x_i using linear interpolation. The vector x specifies the
! coordinates of the underlying interval.
==========================================================================
*/
double *TeosSea::gsw_util_interp1q_int(int nx, double *x, int *iy, int nxi, double *x_i,
        double *y_i)
{
        char    *in_rng;
        int     *j, *k, *r, *jrev, *ki, imax_x, imin_x, i, n, m, ii;
        double  *xi, *xxi, u, max_x, min_x;

        if (nx <= 0 || nxi <= 0)
            return (NULL);

        min_x = max_x = x[0];
        imin_x = imax_x = 0;
        for (i=0; i<nx; i++) {
            if (x[i] < min_x) {
                min_x = x[i];
                imin_x = i;
            } else if (x[i] > max_x) {
                max_x = x[i];
                imax_x = i;
            }
        }
        in_rng = (char *) malloc(nxi*sizeof (char));
        memset(in_rng, 0, nxi*sizeof (char));

        for (i=n=0; i<nxi; i++) {
            if (x_i[i] <= min_x) {
                y_i[i] = iy[imin_x];
            } else if (x_i[i] >= max_x) {
                y_i[i] = iy[imax_x];
            } else {
                in_rng[i] = 1;
                n++;
            }
        }
        if (n==0)
            return (y_i);

        xi = (double *) malloc(n*sizeof (double));
        k  = (int *) malloc(3*n*sizeof (int)); ki = k+n; r = ki+n;
        m  = nx + n;
        xxi = (double *) malloc(m*sizeof (double));
        j = (int *) malloc(2*m*sizeof (int)); jrev = j+m;

        ii = 0;
        for (i = 0; i<nxi; i++) {
            if (in_rng[i]) {
                xi[ii] = x_i[i];
                ki[ii] = i;
                ii++;
            }
        }

        free(in_rng);

    /**    **  Note that the following operations on the index
    **  vectors jrev and r depend on the sort utility
    **  gsw_util_sort_real() consistently ordering the
    **  sorting indexes either in ascending or descending
    **  sequence for replicate values in the real vector.
    */
        gsw_util_sort_real(xi, n, k);
        for (i = 0; i<nx; i++)
            xxi[i] = x[i];
        for (i = 0; i<n; i++)
            xxi[nx+i] = xi[k[i]];
        gsw_util_sort_real(xxi, nx+n, j);

        for (i = 0; i<nx+n; i++)
            jrev[j[i]] = i;
        for (i = 0; i<n; i++)
            r[k[i]] = jrev[nx+i] - i-1;

        for (i = 0; i<n; i++)
        {
            u = (xi[i]-x[r[i]])/(x[r[i]+1]-x[r[i]]);
            y_i[ki[i]] = iy[r[i]] + (iy[r[i]+1]-iy[r[i]])*u;
        }

        free(j); free(xxi); free(k); free(xi);
        return (y_i);
}
/**
==========================================================================
pure method gsw_util_linear_interp (x, y, x_i) result(y_i)
==========================================================================
! Returns the values of the methods y{ny} at the points of column
! vector x_i using linear interpolation. The vector x specifies the
! coordinates of the underlying interval, and the matrix y specifies
| the method values at each x coordinate. Note that y has dimensions
| nx x ny and y_i has dimensions nxi x ny.
! This method was adapted from Matlab's interp1q.
==========================================================================
*/
double *TeosSea::gsw_util_linear_interp(int nx, double *x, int ny, double *y, int nxi,
        double *x_i, double *y_i)
{
        char    *in_rng;
        int     *j, *k, *r, *jrev, *ki, imax_x, imin_x, i, n, m, ii, jy,
                jy0, jyi0, r0;
        double  *xi, *xxi, u, max_x, min_x;

        if (nx <= 0 || nxi <= 0 || ny <= 0)
            return (NULL);

        min_x = max_x = x[0];
        imin_x = imax_x = 0;
        for (i=0; i<nx; i++) {
            if (x[i] < min_x) {
                min_x = x[i];
                imin_x = i;
            } else if (x[i] > max_x) {
                max_x = x[i];
                imax_x = i;
            }
        }
        in_rng = (char *) malloc(nxi*sizeof (char));
        memset(in_rng, 0, nxi*sizeof (char));

        for (i=n=0; i<nxi; i++) {
            if (x_i[i] <= min_x) {
                for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
                    y_i[jyi0+i] = y[jy0+imin_x];
            } else if (x_i[i] >= max_x) {
                for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
                    y_i[jyi0+i] = y[jy0+imax_x];
            } else {
                in_rng[i] = 1;
                n++;
            }
        }
        if (n==0)
            return (y_i);
        xi = (double *) malloc(n*sizeof (double));
        k  = (int *) malloc(3*n*sizeof (int)); ki = k+n; r = ki+n;
        m  = nx + n;
        xxi = (double *) malloc(m*sizeof (double));
        j = (int *) malloc(2*m*sizeof (int)); jrev = j+m;

        ii = 0;
        for (i = 0; i<nxi; i++) {
            if (in_rng[i]) {
                xi[ii] = x_i[i];
                ki[ii] = i;
                ii++;
            }
        }
        free(in_rng);

    /**
    **  This algorithm mimics the Matlab interp1q method.
    **
    **  An explaination of this algorithm:
    **  We have points we are interpolating from (x) and
    **  points that we are interpolating to (xi).  We
    **  sort the interpolating from points, concatenate
    **  them with the interpolating to points and sort the result.
    **  We then construct index r, the interpolation index in x for
    **  each point in xi.
    **
    **  Note that the following operations on the index
    **  vectors jrev and r depend on the sort utility
    **  gsw_util_sort_real() consistently ordering the
    **  sorting indexes either in ascending or descending
    **  sequence for replicate values in the real vector.
    */
        gsw_util_sort_real(xi, n, k);
        memmove(xxi, x, nx*sizeof (double));
        memmove(xxi+nx, xi, n*sizeof (double));
        gsw_util_sort_real(xxi, m, j);

        for (i = 0; i<m; i++)
            jrev[j[i]] = i;
        for (i = 0; i<n; i++)
            r[k[i]] = jrev[nx+i] - i - 1;
            /* this is now the interpolation index in x for a point in xi */

        for (jy=jy0=jyi0=0; jy < ny; jy++, jy0+=nx, jyi0+=nxi)
        {
            for (i = 0; i<n; i++) {
                u = (xi[i]-x[r[i]])/(x[r[i]+1]-x[r[i]]);
                r0 = jy0+r[i];
                y_i[jyi0+ki[i]] = y[r0] + (y[r0+1]-y[r0])*u;
            }
        }

        free(j); free(xxi); free(k); free(xi);

        return (y_i);
}
/**
==========================================================================
method gsw_specvol_anom_standard(sa,ct,p)
==========================================================================
!
!  Calculates specific volume anomaly of seawater.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!
! specvol_anom  :  specific volume anomaly of seawater
*/
double TeosSea::gsw_specvol_anom_standard(double sa, double ct, double p)
{
        return (gsw_specvol(sa,ct,p) - gsw_specvol_sso_0(p));
}
/*************************************************************************
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
*/
int TeosSea::gsw_util_pchip_interp(double *x, double *y, int n,
                          double *xi, double *yi, int ni)
{
    double *d;
    double t, tt, ttt, xx, dx;
    int i, j0, j1, err;
    double h00, h10, h01, h11;

    if (n<2)
    {
        return 1;
    }
    d = (double *)calloc(n, sizeof(double));
    err = gsw_pchip_derivs(x, y, n, d);
    if (err)
    {
                free(d);
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
    free(d);
    return 0;
}

/**
==========================================================================
subroutine gsw_turner_rsubrho(sa,ct,p,nz,tu,rsubrho,p_mid)
==========================================================================

!  Calculates the Turner angle and the Rsubrho as a method of pressure
!  down a vertical water column.  These quantities express the relative
!  contributions of the vertical gradients of Conservative Temperature
!  and Absolute Salinity to the vertical stability (the square of the
!  Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at
!  the mid pressure between the individual data points in the vertical.
!
!  Note that in the double-diffusive literature, papers concerned with
!  the "diffusive" form of double-diffusive convection often define the
!  stability ratio as the reciprocal of what is defined here as the
!  stability ratio.

! sa      : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct      : Conservative Temperature  (a profile (length nz))     [deg C]
! p       : sea pressure              (a profile (length nz))     [dbar]
! nz      : number of bottles
! tu      : Turner angle, on the same (nz-1) grid as p_mid.
!           Turner angle has units of:           [ degrees of rotation ]
! rsubrho : Stability Ratio, on the same (nz-1) grid as p_mid.
!           Rsubrho is dimensionless.                       [ unitless ]
! p_mid   : Mid pressure between p grid  (length nz-1)           [dbar]
*/
void TeosSea::gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz,
        double *tu, double *rsubrho, double *p_mid)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        int     k;
        double  dsa, sa_mid, dct, ct_mid, alpha_mid, beta_mid;

        if (nz < 2)
            return;

        for (k = 0; k < nz-1; k++) {
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

/**
==========================================================================
elemental subroutine gsw_specvol_alpha_beta (sa, ct, p, specvol, alpha, &
                                             beta)
==========================================================================
!
!  Calculates specific volume, the appropiate thermal expansion coefficient
!  and the appropriate saline contraction coefficient of seawater from
!  Absolute Salinity and Conservative Temperature.  This method uses the
!  computationally-efficient expression for specific volume in terms of
!  SA, CT and p (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  specvol =  specific volume                                      [ m/kg ]
!  alpha   =  thermal expansion coefficient                         [ 1/K ]
!             with respect to Conservative Temperature
!  beta    =  saline (i.e. haline) contraction                     [ kg/g ]
!             coefficient at constant Conservative Temperature
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_specvol_alpha_beta(double sa, double ct, double p, double *specvol,
                        double *alpha, double *beta)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  v, v_ct, v_sa_part, xs, ys, z;

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

        if (specvol != NULL)
            *specvol = v;

        if (alpha != NULL) {

            v_ct = gsvco.a000
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

            *alpha = 0.025*v_ct/v;
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
elemental subroutine gsw_specvol_first_derivatives (sa, ct, p, v_sa, v_ct, &
                                                    v_p, iflag)
! =========================================================================
!
!  Calculates three first-order derivatives of specific volume (v).
!  Note that this method uses the computationally-efficient
!  expression for specific volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  v_SA  =  The first derivative of specific volume with respect to
!           Absolute Salinity at constant CT & p.       [ J/(kg (g/kg)^2) ]
!  v_CT  =  The first derivative of specific volume with respect to
!           CT at constant SA and p.                     [ J/(kg K(g/kg)) ]
!  v_P   =  The first derivative of specific volume with respect to
!           P at constant SA and CT.                         [ J/(kg K^2) ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_specvol_first_derivatives (double sa, double ct, double p,
                                double *v_sa, double *v_ct, double *v_p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  v_ct_part, v_p_part, v_sa_part, xs, ys, z;

        xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys = ct*0.025;
        z = p*1e-4;

        if (v_sa != NULL) {
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

            *v_sa = 0.5*gtc.gsw_sfac*v_sa_part/xs;
        }

        if (v_ct != NULL) {

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

            *v_ct = 0.025*v_ct_part;

        }

        if (v_p != NULL) {

        v_p_part = gsvco.c000
        + xs*(gsvco.c100 + xs*(gsvco.c200 + xs*(gsvco.c300
        + xs*(gsvco.c400 + gsvco.c500*xs))))
        + ys*(gsvco.c010 + xs*(gsvco.c110 + xs*(gsvco.c210 + xs*(gsvco.c310 + gsvco.c410*xs))) + ys*(gsvco.c020
        + xs*(gsvco.c120 + xs*(gsvco.c220 + gsvco.c320*xs)) + ys*(gsvco.c030 + xs*(gsvco.c130 + gsvco.c230*xs)
        + ys*(gsvco.c040 + gsvco.c140*xs + gsvco.c050*ys)))) + z*(gsvco.c001 + xs*(gsvco.c101 + xs*(gsvco.c201
        + xs*(gsvco.c301 + gsvco.c401*xs))) + ys*(gsvco.c011 + xs*(gsvco.c111 + xs*(gsvco.c211 + gsvco.c311*xs))
        + ys*(gsvco.c021 + xs*(gsvco.c121 + gsvco.c221*xs) + ys*(gsvco.c031 + gsvco.c131*xs + gsvco.c041*ys)))
        + z*(gsvco.c002 + xs*(gsvco.c102 + gsvco.c202*xs) + ys*(gsvco.c012 + gsvco.c112*xs + gsvco.c022*ys)
        + z*(gsvco.c003 + gsvco.c103*xs + gsvco.c013*ys + z*(gsvco.c004 + gsvco.c005*z))));

            *v_p = 1e-8*v_p_part;

        }
}
/**
==========================================================================
method gsw_cabbeling(sa,ct,p)
==========================================================================

!  Calculates the cabbeling coefficient of seawater with respect to
!  Conservative Temperature.  This method uses the computationally-
!  efficient expression for specific volume in terms of SA, CT and p
!  (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!
! cabbeling  : cabbeling coefficient with respect to       [1/K^2]
!              Conservative Temperature.
*/
double TeosSea::gsw_cabbeling(double sa, double ct, double p)
{
        double  alpha_ct, alpha_on_beta, alpha_sa, beta_sa, rho,
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
/**
==========================================================================
elemental subroutine gsw_rho_second_derivatives_wrt_enthalpy (sa, ct, p, &
                                              rho_sa_sa, rho_sa_h, rho_h_h)
! =========================================================================
!
!  Calculates three second-order derivatives of rho with respect to enthalpy.
!  Note that this method uses the using the computationally-efficient
!  expression for specific volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_SA_SA = The second-order derivative of rho with respect to
!              Absolute Salinity at constant h & p.     [ J/(kg (g/kg)^2) ]
!  rho_SA_h  = The second-order derivative of rho with respect to
!              SA and h at constant p.                   [ J/(kg K(g/kg)) ]
!  rho_h_h   = The second-order derivative of rho with respect to h at
!              constant SA & p
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_rho_second_derivatives_wrt_enthalpy(double sa, double ct, double p,
        double *rho_sa_sa, double *rho_sa_h, double *rho_h_h)
{
        double  rec_v, rec_v2, rec_v3, v_h, v_h_h, v_sa, v_sa_h, v_sa_sa,
                *pv_sa, *pv_h, *pv_sa_sa, *pv_sa_h, *pv_h_h;

        pv_sa   = ((rho_sa_sa != NULL) || (rho_sa_h != NULL)) ? &v_sa : NULL;
        pv_h    = ((rho_sa_h != NULL) || (rho_h_h != NULL)) ?  &v_h : NULL;

        gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,pv_sa,pv_h);

        pv_sa_sa = ((rho_sa_sa != NULL)) ? &v_sa_sa : NULL;
        pv_sa_h  = ((rho_sa_h != NULL)) ? &v_sa_h : NULL;
        pv_h_h   = ((rho_h_h != NULL)) ? &v_h_h : NULL;

        gsw_specvol_second_derivatives_wrt_enthalpy(sa,ct,p,
                                                pv_sa_sa,pv_sa_h,pv_h_h);

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
/**
==========================================================================
elemental subroutine gsw_rho_second_derivatives (sa, ct, p, rho_sa_sa, &
                                  rho_sa_ct, rho_ct_ct, rho_sa_p, rho_ct_p)
==========================================================================
!
!  Calculates five second-order derivatives of rho. Note that this method
!  uses the computationally-efficient expression for specific
!  volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_SA_SA = The second-order derivative of rho with respect to
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  rho_SA_CT = The second-order derivative of rho with respect to
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!  rho_CT_CT = The second-order derivative of rho with respect to CT at
!              constant SA & p
!  rho_SA_P  = The second-order derivative with respect to SA & P at
!              constant CT.
!  rho_CT_P  = The second-order derivative with respect to CT & P at
!              constant SA.
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_rho_second_derivatives(double sa, double ct, double p,
                                            double *rho_sa_sa, double *rho_sa_ct,
                                            double *rho_ct_ct, double *rho_sa_p,
                                            double *rho_ct_p)
{
        double  rec_v, rec_v2, rec_v3, v_ct, v_ct_ct, v_ct_p, v_p, v_sa,
                v_sa_ct, v_sa_p, v_sa_sa;

        gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,&v_p);
        gsw_specvol_second_derivatives(sa,ct,p,&v_sa_sa,&v_sa_ct,&v_ct_ct,
                                    &v_sa_p,&v_ct_p);

        rec_v = 1.0/gsw_specvol(sa,ct,p);
        rec_v2 = pow(rec_v, 2);
        rec_v3 = rec_v2*rec_v;

        if (rho_sa_sa != NULL)
            *rho_sa_sa = -v_sa_sa*rec_v2 + 2.0*v_sa*v_sa*rec_v3;

        if (rho_sa_ct != NULL)
            *rho_sa_ct = -v_sa_ct*rec_v2 + 2.0*v_sa*v_ct*rec_v3;

        if (rho_ct_ct != NULL)
            *rho_ct_ct = -v_ct_ct*rec_v2 + 2.0*v_ct*v_ct*rec_v3;

        if (rho_sa_p != NULL)
            *rho_sa_p = -v_sa_p*rec_v2 + 2.0*v_sa*v_p*rec_v3;

        if (rho_ct_p != NULL)
            *rho_ct_p = -v_ct_p*rec_v2 + 2.0*v_ct*v_p*rec_v3;
}
/**
==========================================================================
method gsw_sa_from_rho(rho,ct,p)
==========================================================================

!  Calculates the Absolute Salinity of a seawater sample, for given values
!  of its density, Conservative Temperature and sea pressure (in dbar).
!
!  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
!   Note. This input has not had 1000 kg/m^3 subtracted from it.
!     That is, it is 'density', not 'density anomaly'.
!  ct  =  Conservative Temperature (ITS-90)                      [ deg C ]
!  p   =  sea pressure                                           [ dbar ]
!
!  sa  =  Absolute Salinity                                      [g/kg]
*/
double TeosSea::gsw_sa_from_rho(double rho, double ct, double p)
{
        int     no_iter;

        double  sa, v_lab, v_0, v_50, v_sa, sa_old, delta_v, sa_mean;

        v_lab   = 1.0/rho;
        v_0     = gsw_specvol(0.0,ct,p);
        v_50    = gsw_specvol(50.0,ct,p);

        sa      = 50.0*(v_lab - v_0)/(v_50 - v_0);
        if (sa < 0.0 || sa > 50.0)
            return (cppGSW_INVALID_VALUE);

        v_sa    = (v_50 - v_0)/50.0;

        for (no_iter=1; no_iter <= 2; no_iter++) {
            sa_old      = sa;
            delta_v     = gsw_specvol(sa_old,ct,p) - v_lab;
            sa          = sa_old - delta_v/v_sa;
            sa_mean     = 0.5*(sa + sa_old);
            gsw_specvol_first_derivatives(sa_mean,ct,p,&v_sa,NULL,NULL);
            sa          = sa_old - delta_v/v_sa;
            if (sa < 0.0 || sa > 50.0)
                return (cppGSW_INVALID_VALUE);
        }
        return (sa);
}
/**
==========================================================================
elemental subroutine gsw_specvol_first_derivatives_wrt_enthalpy (sa, ct, &
                                                       p, v_sa, v_h, iflag)
! =========================================================================
!
!  Calculates two first-order derivatives of specific volume (v).
!  Note that this method uses the using the computationally-efficient
!  expression for specific volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  v_SA  =  The first derivative of specific volume with respect to
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  v_h  =  The first derivative of specific volume with respect to
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_specvol_first_derivatives_wrt_enthalpy(double sa, double ct, double p,
        double *v_sa, double *v_h)
{
        double  h_ct=1.0, h_sa, rec_h_ct, vct_ct, vct_sa;

        if (v_sa != NULL) {

            gsw_specvol_first_derivatives(sa,ct,p,&vct_sa,&vct_ct,NULL);
            gsw_enthalpy_first_derivatives(sa,ct,p,&h_sa,&h_ct);

        } else if (v_h != NULL) {

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
/**
==========================================================================
elemental subroutine gsw_rho_first_derivatives_wrt_enthalpy (sa, ct, p, &
                                                             rho_sa, rho_h)
! =========================================================================
!
!  Calculates two first-order derivatives of specific volume (v).
!  Note that this method uses the using the computationally-efficient
!  expression for specific volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_SA =  The first derivative of rho with respect to
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  rho_h  =  The first derivative of rho with respect to
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_rho_first_derivatives_wrt_enthalpy (double sa, double ct, double p,
        double *rho_sa, double *rho_h)
{
        double  rec_v2, v_h=0.0, v_sa;

        if ((rho_sa != NULL) && (rho_h != NULL)) {

            gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,&v_h);

        } else if (rho_sa != NULL) {

            gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,NULL);

        } else if (rho_h != NULL) {

            gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,NULL,&v_h);

        }

        rec_v2 = pow(1.0/gsw_specvol(sa,ct,p), 2);

        if (rho_sa != NULL) *rho_sa = -v_sa*rec_v2;

        if (rho_h != NULL) *rho_h = -v_h*rec_v2;
}
/**
==========================================================================
elemental subroutine gsw_specvol_second_derivatives_wrt_enthalpy (sa, ct, &
                                          p, v_sa_sa, v_sa_h, v_h_h, iflag)
! =========================================================================
!
!  Calculates three first-order derivatives of specific volume (v) with
!  respect to enthalpy. Note that this method uses the using the
!  computationally-efficient expression for specific volume
!  (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  v_SA_SA = The second-order derivative of specific volume with respect to
!            Absolute Salinity at constant h & p.       [ J/(kg (g/kg)^2) ]
!  v_SA_h  = The second-order derivative of specific volume with respect to
!            SA and h at constant p.                     [ J/(kg K(g/kg)) ]
!  v_h_h   = The second-order derivative with respect to h at
!            constant SA & p.
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_specvol_second_derivatives_wrt_enthalpy (double sa, double ct, double p,
        double *v_sa_sa, double *v_sa_h, double *v_h_h)
{
        double  h_ct, h_ct_ct, h_sa, h_sa_ct, h_sa_sa, rec_h_ct, v_h_h_part,
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

        if ((v_sa_sa != NULL) || (v_sa_h != NULL)) {

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
method gsw_thermobaric(sa,ct,p)
==========================================================================

!  Calculates the thermobaric coefficient of seawater with respect to
!  Conservative Temperature.  This routine is based on the
!  computationally-efficient expression for specific volume in terms of
!  SA, CT and p (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!
! thermobaric  : thermobaric coefficient with              [1/(K Pa)]
!                    respect to Conservative Temperature (48 term equation)
*/
double TeosSea::gsw_thermobaric(double sa, double ct, double p)
{
        double  v_ct, v_ct_p, v_sa, v_sa_p;

        gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,NULL);

        gsw_specvol_second_derivatives(sa,ct,p,NULL,NULL,NULL,&v_sa_p,&v_ct_p);

        return (gsw_rho(sa,ct,p)*(v_ct_p - (v_ct/v_sa)*v_sa_p));
}
/**
==========================================================================
elemental subroutine gsw_specvol_second_derivatives (sa, ct, p, v_sa_sa, &
                                   v_sa_ct, v_ct_ct, v_sa_p, v_ct_p, iflag)
! =========================================================================
!
!  Calculates five second-order derivatives of specific volume (v).
!  Note that this method uses the computationally-efficient
!  expression for specific volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  v_SA_SA  =  The second derivative of specific volume with respect to
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  v_SA_CT  =  The second derivative of specific volume with respect to
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!  v_CT_CT  =  The second derivative of specific volume with respect to
!              CT at constant SA and p.                      [ J/(kg K^2) ]
!  v_SA_P  =  The second derivative of specific volume with respect to
!              SA and P at constant CT.                  [ J/(kg K(g/kg)) ]
!  v_CT_P  =  The second derivative of specific volume with respect to
!              CT and P at constant SA.                  [ J/(kg K(g/kg)) ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_specvol_second_derivatives (double sa, double ct, double p,
        double *v_sa_sa, double *v_sa_ct, double *v_ct_ct, double *v_sa_p,
        double *v_ct_p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  v_ct_ct_part, v_ct_p_part, v_sa_ct_part, v_sa_p_part,
                v_sa_sa_part, xs, xs2, ys, z;

        xs2 = gtc.gsw_sfac*sa + gtc.offset;
        xs = sqrt(xs2);
        ys = ct*0.025;
        z = p*1e-4;

        if (v_sa_sa != NULL) {

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

        if (v_sa_ct != NULL) {
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

        if (v_ct_ct != NULL) {
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

        if (v_sa_p != NULL) {
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

        if (v_ct_p != NULL) {
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
/**
==========================================================================
method gsw_alpha_on_beta(sa,ct,p)
==========================================================================

!  Calculates alpha divided by beta, where alpha is the thermal expansion
!  coefficient and beta is the saline contraction coefficient of seawater
!  from Absolute Salinity and Conservative Temperature.  This method uses
!  the computationally-efficient expression for specific volume in terms of
!  SA, CT and p (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
!
! alpha_on_beta
!        : thermal expansion coefficient with respect to   [kg g^-1 K^-1]
!          Conservative Temperature divided by the saline
!          contraction coefficient at constant Conservative
!          Temperature
*/
double TeosSea::gsw_alpha_on_beta(double sa, double ct, double p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */
        double  xs, ys, z, v_ct_part, v_sa_part;

        xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys      = ct*0.025;
        z       = p*1e-4;

        v_ct_part = gsvco.a000
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

        return (-(v_ct_part*xs)/(20.0*gtc.gsw_sfac*v_sa_part));
}
/**
==========================================================================
elemental subroutine gsw_enthalpy_second_derivatives (sa, ct, p, h_sa_sa, &
                                                      h_sa_ct, h_ct_ct)
! =========================================================================
!
!  Calculates the following three second-order derivatives of specific
!  enthalpy (h),using the computationally-efficient expression for
!  specific volume in terms of SA, CT and p (Roquet et al., 2014).
!   (1) h_SA_SA, second-order derivative with respect to Absolute Salinity
!       at constant CT & p.
!   (2) h_SA_CT, second-order derivative with respect to SA & CT at
!       constant p.
!   (3) h_CT_CT, second-order derivative with respect to CT at constant SA
!       and p.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_SA_SA  =  The second derivative of specific enthalpy with respect to
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  h_SA_CT  =  The second derivative of specific enthalpy with respect to
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!  h_CT_CT  =  The second derivative of specific enthalpy with respect to
!              CT at constant SA and p.                      [ J/(kg K^2) ]
!--------------------------------------------------------------------------
*/
void TeosSea::gsw_enthalpy_second_derivatives(double sa, double ct, double p,
                double *h_sa_sa, double *h_sa_ct, double *h_ct_ct)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  dynamic_h_ct_ct_part, dynamic_h_sa_ct_part,
                dynamic_h_sa_sa_part, xs, xs2, ys, z;

        xs = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys = ct*0.025;
        z = p*1e-4;

        if (h_sa_sa != NULL) {

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

            *h_sa_sa = 1e8*0.25*gtc.gsw_sfac*gtc.gsw_sfac*dynamic_h_sa_sa_part/
                        pow(xs,3);

        }

        if (h_sa_ct != NULL) {

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

        if (h_ct_ct != NULL) {

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
/**
==========================================================================
method gsw_alpha_wrt_t_exact(sa,t,p)
==========================================================================

! Calculates thermal expansion coefficient of seawater with respect to
! in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : insitu temperature                              [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_alpha_wrt_t_exact : thermal expansion coefficient    [1/K]
!                         wrt (in-situ) temperature
*/
double TeosSea::gsw_alpha_wrt_t_exact(double sa, double t, double p)
{
        int     n0=0, n1=1;

        return (gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}
/**
==========================================================================
method gsw_beta_const_t_exact(sa,t,p)
==========================================================================

! Calculates saline (haline) contraction coefficient of seawater at
! constant in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! beta_const_t_exact : haline contraction coefficient      [kg/g]
*/
double TeosSea::gsw_beta_const_t_exact(double sa, double t, double p)
{
        int     n0=0, n1=1;

        return (-gsw_gibbs(n1,n0,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}
/**
==========================================================================
method gsw_c_from_sp(sp,t,p)
==========================================================================

!  Calculates conductivity, C, from (SP,t,p) using PSS-78 in the range
!  2 < SP < 42.  If the input Practical Salinity is less than 2 then a
!  modified form of the Hill et al. (1986) fomula is used for Practical
!  Salinity.  The modification of the Hill et al. (1986) expression is to
!  ensure that it is exactly consistent with PSS-78 at SP = 2.
!
!  The conductivity ratio returned by this method is consistent with the
!  input value of Practical Salinity, SP, to 2x10^-14 psu over the full
!  range of input parameters (from pure fresh water up to SP = 42 psu).
!  This error of 2x10^-14 psu is machine precision at typical seawater
!  salinities.  This accuracy is achieved by having four different
!  polynomials for the starting value of Rtx (the square root of Rt) in
!  four different ranges of SP, and by using one and a half iterations of
!  a computationally efficient modified Newton-Raphson technique (McDougall
!  and Wotherspoon, 2012) to find the root of the equation.
!
!  Note that strictly speaking PSS-78 (Unesco, 1983) defines Practical
!  Salinity in terms of the conductivity ratio, R, without actually
!  specifying the value of C(35,15,0) (which we currently take to be
!  42.9140 mS/cm).
!
! sp     : Practical Salinity                               [unitless]
! t      : in-situ temperature [ITS-90]                     [deg C]
! p      : sea pressure                                     [dbar]
!
! c      : conductivity                                     [ mS/cm ]
*/
double TeosSea::gsw_c_from_sp(double sp, double t, double p)
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

        double  t68, ft68, x, rtx=0.0, dsp_drtx, sqrty,
                part1, part2, hill_ratio, sp_est,
                rtx_old, rt, aa, bb, cc, dd, ee, ra,r, rt_lc, rtxm,
                sp_hill_raw;

        t68     = t*1.00024e0;
        ft68    = (t68 - 15e0)/(1e0 + gspc.k*(t68 - 15e0));

        x       = sqrt(sp);

    /**-------------------------------------------------------------------------
     ! Finding the starting value of Rtx, the square root of Rt, using four
     ! different polynomials of SP and t68.
     !--------------------------------------------------------------------------
    */

        if (sp >= 9.0) {
            rtx = p0 + x*(p1 + p4*t68 + x*(p3 + p7*t68 + x*(p6
                  + p11*t68 + x*(p10 + p16*t68 + x*p15))))
                  + t68*(p2+ t68*(p5 + x*x*(p12 + x*p17) + p8*x
                  + t68*(p9 + x*(p13 + x*p18)+ t68*(p14 + p19*x + p20*t68))));
        } else if (sp >= 0.25 && sp < 9.0) {
            rtx = q0 + x*(q1 + q4*t68 + x*(q3 + q7*t68 + x*(q6
                  + q11*t68 + x*(q10 + q16*t68 + x*q15))))
                  + t68*(q2+ t68*(q5 + x*x*(q12 + x*q17) + q8*x
                  + t68*(q9 + x*(q13 + x*q18)+ t68*(q14 + q19*x + q20*t68))));
        } else if (sp >= 0.003 && sp < 0.25) {
            rtx =  s0 + x*(s1 + s4*t68 + x*(s3 + s7*t68 + x*(s6
                  + s11*t68 + x*(s10 + s16*t68 + x*s15))))
                  + t68*(s2+ t68*(s5 + x*x*(s12 + x*s17) + s8*x
                  + t68*(s9 + x*(s13 + x*s18)+ t68*(s14 + s19*x + s20*t68))));
        } else if (sp < 0.003) {
            rtx =  u0 + x*(u1 + u4*t68 + x*(u3 + u7*t68 + x*(u6
                  + u11*t68 + x*(u10 + u16*t68 + x*u15))))
                  + t68*(u2+ t68*(u5 + x*x*(u12 + x*u17) + u8*x
                  + t68*(u9 + x*(u13 + x*u18)+ t68*(u14 + u19*x + u20*t68))));
        }

    /**-------------------------------------------------------------------------
     ! Finding the starting value of dSP_dRtx, the derivative of SP with respect
     ! to Rtx.
     !--------------------------------------------------------------------------
    */
        dsp_drtx        =  gspc.a1 + (2e0*gspc.a2 + (3e0*gspc.a3 +
                                (4e0*gspc.a4 + 5e0*gspc.a5*rtx)*rtx)*rtx)*rtx
                          + ft68*(gspc.b1 + (2e0*gspc.b2 + (3e0*gspc.b3 + (4e0*gspc.b4 +
                                5e0*gspc.b5*rtx)*rtx)*rtx)*rtx);

        if (sp < 2.0) {
            x           = 400e0*(rtx*rtx);
            sqrty       = 10.0*rtx;
            part1       = 1e0 + x*(1.5e0 + x);
            part2       = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
            hill_ratio  = gsw_hill_ratio_at_sp2(t);
            dsp_drtx    = dsp_drtx
                          + gspc.a0*800e0*rtx*(1.5e0 + 2e0*x)/(part1*part1)
                          + gspc.b0*ft68*(10e0 + sqrty*(20e0 + 30e0*sqrty))/
                                (part2*part2);
            dsp_drtx    = hill_ratio*dsp_drtx;
        }

    /**-------------------------------------------------------------------------
     ! One iteration through the modified Newton-Raphson method (McDougall and
     ! Wotherspoon, 2012) achieves an error in Practical Salinity of about
     ! 10^-12 for all combinations of the inputs.  One and a half iterations of
     ! the modified Newton-Raphson method achevies a maximum error in terms of
     ! Practical Salinity of better than 2x10^-14 everywhere.
     !
     ! We recommend one and a half iterations of the modified Newton-Raphson
     ! method.
     !
     ! Begin the modified Newton-Raphson method.
     !--------------------------------------------------------------------------
    */
        sp_est  = gspc.a0 + (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4 + gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx
                + ft68*(gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3 + (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);
        if (sp_est <  2.0) {
            x           = 400e0*(rtx*rtx);
            sqrty       = 10e0*rtx;
            part1       = 1e0 + x*(1.5e0 + x);
            part2       = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
            sp_hill_raw = sp_est - gspc.a0/part1 - gspc.b0*ft68/part2;
            hill_ratio  = gsw_hill_ratio_at_sp2(t);
            sp_est      = hill_ratio*sp_hill_raw;
        }

        rtx_old = rtx;
        rtx     = rtx_old - (sp_est - sp)/dsp_drtx;

        rtxm    = 0.5e0*(rtx + rtx_old); /** This mean value of Rtx, Rtxm, is the
                  value of Rtx at which the derivative dSP_dRtx is evaluated. */

        dsp_drtx=  gspc.a1 + (2e0*gspc.a2 + (3e0*gspc.a3 + (4e0*gspc.a4 +
                                5e0*gspc.a5*rtxm)*rtxm)*rtxm)*rtxm
                   + ft68*(gspc.b1 + (2e0*gspc.b2 + (3e0*gspc.b3 + (4e0*gspc.b4 +
                                5e0*gspc.b5*rtxm)*rtxm)*rtxm)*rtxm);
        if (sp_est <  2.0) {
            x   = 400e0*(rtxm*rtxm);
            sqrty       = 10e0*rtxm;
            part1       = 1e0 + x*(1.5e0 + x);
            part2       = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
            dsp_drtx    = dsp_drtx
                          + gspc.a0*800e0*rtxm*(1.5e0 + 2e0*x)/(part1*part1)
                          + gspc.b0*ft68*(10e0 + sqrty*(20e0 + 30e0*sqrty))/
                                (part2*part2);
            hill_ratio  = gsw_hill_ratio_at_sp2(t);
            dsp_drtx    = hill_ratio*dsp_drtx;
        }

    /**------------------------------------------------------------------------
     ! The line below is where Rtx is updated at the end of the one full
     ! iteration of the modified Newton-Raphson technique.
     !-------------------------------------------------------------------------
    */
        rtx     = rtx_old - (sp_est - sp)/dsp_drtx;
    /**------------------------------------------------------------------------
     ! Now we do another half iteration of the modified Newton-Raphson
     ! technique, making a total of one and a half modified N-R iterations.
     !-------------------------------------------------------------------------
    */
        sp_est  = gspc.a0 + (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4 + gspc.a5*rtx)*rtx)*rtx)*rtx)*rtx
                + ft68*(gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3 + (gspc.b4 + gspc.b5*rtx)*rtx)*rtx)*rtx)*rtx);
        if (sp_est <  2.0) {
            x           = 400e0*(rtx*rtx);
            sqrty       = 10e0*rtx;
            part1       = 1e0 + x*(1.5e0 + x);
            part2       = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
            sp_hill_raw = sp_est - gspc.a0/part1 - gspc.b0*ft68/part2;
            hill_ratio  = gsw_hill_ratio_at_sp2(t);
            sp_est      = hill_ratio*sp_hill_raw;
        }
        rtx     = rtx - (sp_est - sp)/dsp_drtx;

    /**
    !--------------------------------------------------------------------------
     ! Now go from Rtx to Rt and then to the conductivity ratio R at pressure p.
     !--------------------------------------------------------------------------
    */
        rt      = rtx*rtx;

        aa      = gspc.d3 + gspc.d4*t68;
        bb      = 1e0 + t68*(gspc.d1 + gspc.d2*t68);
        cc      = p*(gspc.e1 + p*(gspc.e2 + gspc.e3*p));
    /** rt_lc (i.e. rt_lower_case) corresponds to rt as defined in
       the UNESCO 44 (1983) routines. */
        rt_lc   = gspc.c0 + (gspc.c1 + (gspc.c2 + (gspc.c3 + gspc.c4*t68)*t68)*t68)*t68;

        dd      = bb - aa*rt_lc*rt;
        ee      = rt_lc*rt*aa*(bb + cc);
        ra      = sqrt(dd*dd + 4e0*ee) - dd;
        r       = 0.5e0*ra/aa;

    /**
     The dimensionless conductivity ratio, R, is the conductivity input, C,
     ! divided by the present estimate of C(SP=35, t_68=15, p=0) which is
     ! 42.9140 mS/cm (=4.29140 S/m^).
    */
        return (gtc.gsw_c3515*r);
}
/**
==========================================================================
method gsw_ct_from_pt(sa,pt)
==========================================================================

! Calculates Conservative Temperature from potential temperature of seawater
!
! sa      : Absolute Salinity                              [g/kg]
! pt      : potential temperature with                     [deg C]
!           reference pressure of 0 dbar
!
! gsw_ct_from_pt : Conservative Temperature                [deg C]
*/
double TeosSea::gsw_ct_from_pt(double sa, double pt)
{
   /** for GSW_TEOS10_CONSTANTS use gtc */

   double  x2, x, y, pot_enthalpy, rval;

   x2 = gtc.gsw_sfac*sa;
   x = sqrt(x2);
   y = pt*0.025e0;   /*! normalize for F03 and F08 */

   pot_enthalpy    =  61.01362420681071e0 + y*(168776.46138048015e0 +
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

   rval = pot_enthalpy/gtc.gsw_cp0;

   return rval;
}
/**
==========================================================================
subroutine gsw_rho_first_derivatives(sa, ct, p, drho_dsa, drho_dct, drho_dp)
==========================================================================

!  Calculates the three (3) partial derivatives of in situ density with
!  respect to Absolute Salinity, Conservative Temperature and pressure.
!  Note that the pressure derivative is done with respect to pressure in
!  Pa, not dbar.  This method uses the computationally-efficient expression
!  for specific volume in terms of SA, CT and p (Roquet et al., 2014).
!
! sa        : Absolute Salinity                               [g/kg]
! ct        : Conservative Temperature                        [deg C]
! p         : sea pressure                                    [dbar]
! drho_dsa  : partial derivatives of density                  [kg^2/(g m^3)]
!             with respect to Absolute Salinity
! drho_dct  : partial derivatives of density                  [kg/(K m^3)]
!             with respect to Conservative Temperature
! drho_dp   : partial derivatives of density                  [kg/(Pa m^3)]
!             with respect to pressure in Pa
*/
void TeosSea::gsw_rho_first_derivatives(double sa, double ct, double p,
        double *drho_dsa, double *drho_dct, double *drho_dp)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  rho2, v_ct, v_p, v_sa, xs, ys, z, v;

        xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys      = ct*0.025;
        z       = p*1e-4;

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

        rho2 = pow(1.0/v, 2.0);

        if (drho_dsa != NULL) {

            v_sa = gsvco.b000
           + xs*(gsvco.b100 + xs*(gsvco.b200 + xs*(gsvco.b300 + xs*(gsvco.b400 + gsvco.b500*xs))))
           + ys*(gsvco.b010 + xs*(gsvco.b110 + xs*(gsvco.b210 + xs*(gsvco.b310 + gsvco.b410*xs)))
           + ys*(gsvco.b020 + xs*(gsvco.b120 + xs*(gsvco.b220 + gsvco.b320*xs)) + ys*(gsvco.b030
           + xs*(gsvco.b130 + gsvco.b230*xs) + ys*(gsvco.b040 + gsvco.b140*xs + gsvco.b050*ys))))
           + z*(gsvco.b001 + xs*(gsvco.b101 + xs*(gsvco.b201 + xs*(gsvco.b301 + gsvco.b401*xs)))
           + ys*(gsvco.b011 + xs*(gsvco.b111 + xs*(gsvco.b211 + gsvco.b311*xs)) + ys*(gsvco.b021
           + xs*(gsvco.b121 + gsvco.b221*xs) + ys*(gsvco.b031 + gsvco.b131*xs + gsvco.b041*ys)))
           + z*(gsvco.b002 + xs*(gsvco.b102 + xs*(gsvco.b202 + gsvco.b302*xs))+ ys*(gsvco.b012
           + xs*(gsvco.b112 + gsvco.b212*xs) + ys*(gsvco.b022 + gsvco.b122*xs + gsvco.b032*ys))
           + z*(gsvco.b003 + gsvco.b103*xs + gsvco.b013*ys + gsvco.b004*z)));

            *drho_dsa = -rho2*0.5*gtc.gsw_sfac*v_sa/xs;
        }

        if (drho_dct != NULL) {

            v_ct = gsvco.a000
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

            *drho_dct = -rho2*0.025*v_ct;
        }

        if (drho_dp != NULL) {

            v_p = gsvco.c000
          + xs*(gsvco.c100 + xs*(gsvco.c200 + xs*(gsvco.c300 + xs*(gsvco.c400 + gsvco.c500*xs))))
          + ys*(gsvco.c010 + xs*(gsvco.c110 + xs*(gsvco.c210 + xs*(gsvco.c310 + gsvco.c410*xs))) + ys*(gsvco.c020
          + xs*(gsvco.c120 + xs*(gsvco.c220 + gsvco.c320*xs)) + ys*(gsvco.c030 + xs*(gsvco.c130 + gsvco.c230*xs)
          + ys*(gsvco.c040 + gsvco.c140*xs + gsvco.c050*ys)))) + z*(gsvco.c001 + xs*(gsvco.c101 + xs*(gsvco.c201
          + xs*(gsvco.c301 + gsvco.c401*xs))) + ys*(gsvco.c011 + xs*(gsvco.c111 + xs*(gsvco.c211 + gsvco.c311*xs))
          + ys*(gsvco.c021 + xs*(gsvco.c121 + gsvco.c221*xs) + ys*(gsvco.c031 + gsvco.c131*xs + gsvco.c041*ys)))
          + z*(gsvco.c002 + xs*(gsvco.c102 + gsvco.c202*xs) + ys*(gsvco.c012 + gsvco.c112*xs + gsvco.c022*ys)
          + z*(gsvco.c003 + gsvco.c103*xs + gsvco.c013*ys + z*(gsvco.c004 + gsvco.c005*z))));

            *drho_dp = 1e-4*gtc.pa2db*-rho2*v_p;
        }

        return;
}
/**
==========================================================================
method gsw_rho_t_exact(sa,t,p)
==========================================================================

! Calculates in-situ density of seawater from Absolute Salinity and
! in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_rho_t_exact : in-situ density                        [kg/m^3]
*/
double TeosSea::gsw_rho_t_exact(double sa, double t, double p)
{
        int     n0=0, n1=1;

        return (1.0/gsw_gibbs(n0,n0,n1,sa,t,p));
}

/**
==========================================================================
method gsw_pot_rho_t_exact(sa,t,p,p_ref)
==========================================================================

! Calculates the potential density of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
!
! gsw_pot_rho_t_exact : potential density                  [kg/m^3]
*/
double TeosSea::gsw_pot_rho_t_exact(double sa, double t, double p, double p_ref)
{
        double  pt = gsw_pt_from_t(sa,t,p,p_ref);

        return (gsw_rho_t_exact(sa,pt,p_ref));
}

/**
==========================================================================
method gsw_sp_from_sstar(sstar,p,lon,lat)
==========================================================================

! Calculates Practical Salinity, SP, from Preformed Salinity, Sstar.
!
! sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]
*/
double TeosSea::gsw_sp_from_sstar(double sstar, double p, double lon, double lat)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        double  saar, sp_baltic;

    /**
    **! In the Baltic Sea, SA = Sstar.
    */
        sp_baltic       = gsw_sp_from_sa_baltic(sstar,lon,lat);

        if (sp_baltic < cppGSW_ERROR_LIMIT)
            return (sp_baltic);

        saar    = gsw_saar(p,lon,lat);
        if (saar == cppGSW_INVALID_VALUE)
            return (saar);

        return ((sstar/gtc.gsw_ups)/(1.0 - 0.35e0*saar));
}
/**
==========================================================================
method gsw_sp_salinometer(rt,t)
==========================================================================
!  Calculates Practical Salinity SP from a salinometer, primarily using the
!  PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical Salinity
!  is only valid in the range 2 < SP < 42.  If the PSS-78 algorithm
!  produces a Practical Salinity that is less than 2 then the Practical
!  Salinity is recalculated with a modified form of the Hill et al. (1986)
!  formula.  The modification of the Hill et al. (1986) expression is to
!  ensure that it is exactly consistent with PSS-78 at SP = 2.
!
!  A laboratory salinometer has the ratio of conductivities, Rt, as an
!  output, and the present method uses this conductivity ratio and the
!  temperature t of the salinometer bath as the two input variables.
!
!  rt  = C(SP,t_68,0)/C(SP=35,t_68,0)                          [ unitless ]
!  t   = temperature of the bath of the salinometer,
!        measured on the ITS-90 scale (ITS-90)                 [ deg C ]
!
!  gsw_sp_salinometer = Practical Salinity on the PSS-78 scale [ unitless ]
*/
double TeosSea::gsw_sp_salinometer(double rt, double t)
{
  /** for GSW_SP_COEFFICIENTS use gspc */

  double t68, ft68, rtx, sp, hill_ratio,
         x, sqrty, part1, part2, sp_hill_raw;

  if (rt < 0){
    return NAN;
  }

  t68 = t*1.00024;
  ft68 = (t68 - 15)/(1 + gspc.k*(t68 - 15));

  rtx = sqrt(rt);

  sp = gspc.a0 + (gspc.a1 + (gspc.a2 + (gspc.a3 + (gspc.a4 + gspc.a5 * rtx) * rtx) * rtx) * rtx) * rtx
    + ft68 * (gspc.b0 + (gspc.b1 + (gspc.b2+ (gspc.b3 + (gspc.b4 + gspc.b5 * rtx) * rtx) * rtx) * rtx) * rtx);

    /**
     ! The following section of the code is designed for SP < 2 based on the
     ! Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
     ! exactly equal to the PSS-78 algorithm at SP = 2.
    */

   if (sp < 2) {
       hill_ratio  = gsw_hill_ratio_at_sp2(t);
       x           = 400e0*rt;
       sqrty       = 10e0*rtx;
       part1       = 1e0 + x*(1.5e0 + x);
       part2       = 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
       sp_hill_raw = sp - gspc.a0/part1 - gspc.b0*ft68/part2;
       sp          = hill_ratio*sp_hill_raw;
   }

  return sp;

}
/**
==========================================================================
method gsw_sp_from_sr(sr)
==========================================================================

! Calculates Practical Salinity, sp, from Reference Salinity, sr.
!
! sr     : Reference Salinity                              [g/kg]
!
! gsw_sp_from_sr  : Practical Salinity                     [unitless]
*/
double TeosSea::gsw_sp_from_sr(double sr)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        return(sr/gtc.gsw_ups);
}
/**
==========================================================================
method gsw_sp_from_sk(sk)
==========================================================================

! Calculates Practical Salinity, SP, from SK
!
!  SK    : Knudsen Salinity                        [parts per thousand, ppt]
!
! gsw_sp_from_sk  : Practical Salinity                              [unitless]
*/
double TeosSea::gsw_sp_from_sk(double sk)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */

        double  gsw_sp_from_sk_value;


        gsw_sp_from_sk_value = (sk - 0.03e0)*(gtc.gsw_soncl/1.805e0);

        if (gsw_sp_from_sk_value < 0e0)
            gsw_sp_from_sk_value = cppGSW_INVALID_VALUE;

        return (gsw_sp_from_sk_value);
}
/**
==========================================================================
method gsw_sound_speed_t_exact(sa,t,p)
==========================================================================

! Calculates the speed of sound in seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_sound_speed_t_exact : sound speed                    [m/s]
*/
double TeosSea::gsw_sound_speed_t_exact(double sa, double t, double p)
{
        int     n0=0, n1=1, n2=2;
        double  g_tt, g_tp;

        g_tt    = gsw_gibbs(n0,n2,n0,sa,t,p);
        g_tp    = gsw_gibbs(n0,n1,n1,sa,t,p);

        return (gsw_gibbs(n0,n0,n1,sa,t,p) *
                sqrt(g_tt/(g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p))));
}
/**
==========================================================================
method gsw_sound_speed(sa,ct,p)
==========================================================================

!  Calculates the speed of sound in seawater.  This method has inputs of
!  Absolute Salinity and Conservative Temperature.  This method uses the
!  computationally-efficient expression for specific volume in terms of SA,
!  CT and p (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!
! sound_speed  : speed of sound in seawater                [m/s]
*/
double TeosSea::gsw_sound_speed(double sa, double ct, double p)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  v, v_p, xs, ys, z;

        xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys      = ct*0.025;
        z       = p*1e-4;

        v       = gsvco.v000
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

        v_p     = gsvco.c000
        + xs*(gsvco.c100 + xs*(gsvco.c200 + xs*(gsvco.c300 + xs*(gsvco.c400 + gsvco.c500*xs))))
        + ys*(gsvco.c010 + xs*(gsvco.c110 + xs*(gsvco.c210 + xs*(gsvco.c310 + gsvco.c410*xs))) + ys*(gsvco.c020
        + xs*(gsvco.c120 + xs*(gsvco.c220 + gsvco.c320*xs)) + ys*(gsvco.c030 + xs*(gsvco.c130 + gsvco.c230*xs)
        + ys*(gsvco.c040 + gsvco.c140*xs + gsvco.c050*ys)))) + z*(gsvco.c001 + xs*(gsvco.c101 + xs*(gsvco.c201
        + xs*(gsvco.c301 + gsvco.c401*xs))) + ys*(gsvco.c011 + xs*(gsvco.c111 + xs*(gsvco.c211 + gsvco.c311*xs))
        + ys*(gsvco.c021 + xs*(gsvco.c121 + gsvco.c221*xs) + ys*(gsvco.c031 + gsvco.c131*xs + gsvco.c041*ys)))
        + z*(gsvco.c002 + xs*(gsvco.c102 + gsvco.c202*xs) + ys*(gsvco.c012 + gsvco.c112*xs + gsvco.c022*ys)
        + z*(gsvco.c003 + gsvco.c103*xs + gsvco.c013*ys + z*(gsvco.c004 + gsvco.c005*z))));

        return (10000.0*sqrt(-v*v/v_p));
}
/**
==========================================================================
method gsw_sigma0(sa,ct)
==========================================================================

!  Calculates potential density anomaly with reference pressure of 0 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  method has inputs of Absolute Salinity and Conservative Temperature.
!  This method uses the computationally-efficient 48-term expression for
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
!
! gsw_sigma0  : potential density anomaly with reference pressure of 0
!                                                      (48 term equation)
*/
double TeosSea::gsw_sigma0(double sa, double ct)
{
        /** for GSW_TEOS10_CONSTANTS use gtc */
        /** for GSW_SPECVOL_COEFFICIENTS use gsvco */

        double  vp0, xs, ys;

        xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys      = ct*0.025;

        vp0     = gsvco.v000
        + xs*(gsvco.v010 + xs*(gsvco.v020 + xs*(gsvco.v030 + xs*(gsvco.v040 + xs*(gsvco.v050
        + gsvco.v060*xs))))) + ys*(gsvco.v100 + xs*(gsvco.v110 + xs*(gsvco.v120 + xs*(gsvco.v130 + xs*(gsvco.v140
        + gsvco.v150*xs)))) + ys*(gsvco.v200 + xs*(gsvco.v210 + xs*(gsvco.v220 + xs*(gsvco.v230 + gsvco.v240*xs)))
        + ys*(gsvco.v300 + xs*(gsvco.v310 + xs*(gsvco.v320 + gsvco.v330*xs)) + ys*(gsvco.v400 + xs*(gsvco.v410
        + gsvco.v420*xs) + ys*(gsvco.v500 + gsvco.v510*xs + gsvco.v600*ys)))));

        return (1.0/vp0 - 1000.0);
}
/**
==========================================================================
method gsw_sigma1(sa,ct)
==========================================================================

!  Calculates potential density anomaly with reference pressure of 1000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  method has inputs of Absolute Salinity and Conservative Temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
!
! sigma1 : potential density anomaly with reference pressure of 1000
*/
double TeosSea::gsw_sigma1(double sa, double ct)
{
        return (gsw_rho(sa,ct,1000.0) - 1000.0);
}
/**
==========================================================================
method gsw_sigma2(sa,ct)
==========================================================================

!  Calculates potential density anomaly with reference pressure of 2000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  method has inputs of Absolute Salinity and Conservative Temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
!
! sigma2 : potential density anomaly with reference pressure of 2000
*/
double TeosSea::gsw_sigma2(double sa, double ct)
{
        return (gsw_rho(sa,ct,2000.0) - 1000.0);
}
/**
==========================================================================
method gsw_sigma3(sa,ct)
==========================================================================

!  Calculates potential density anomaly with reference pressure of 3000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  method has inputs of Absolute Salinity and Conservative Temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
!
! sigma3 : potential density anomaly with reference pressure of 3000
*/
double TeosSea::gsw_sigma3(double sa, double ct)
{
        return (gsw_rho(sa,ct,3000.0) - 1000.0);
}
/**
==========================================================================
method gsw_sigma4(sa,ct)
==========================================================================

!  Calculates potential density anomaly with reference pressure of 4000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  method has inputs of Absolute Salinity and Conservative Temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
!
! sigma4  : potential density anomaly with reference pressure of 4000
*/
double TeosSea::gsw_sigma4(double sa, double ct)
{
        return (gsw_rho(sa,ct,4000.0) - 1000.0);
}
/**
==========================================================================
method gsw_spiciness0(sa,ct)
==========================================================================
!
!  Calculates spiciness from Absolute Salinity and Conservative
!  Temperature at a pressure of 0 dbar, as described by McDougall and
!  Krzysik (2015).  This routine is based on the computationally-efficient
!  expression for specific volume in terms of SA, CT and p (Roquet et al.,
!  2015).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                        [ deg C ]
!
!  spiciness0  =  spiciness referenced to a pressure of 0 dbar,
!                 i.e. the surface                                 [ kg/m^3 ]
!
*/

double TeosSea::gsw_spiciness0(double sa, double ct)
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
        double  xs, ys, spiciness0;

        xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys      = ct*0.025;

        spiciness0= s01+ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
                +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
                +xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
                +xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
                +xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
                +xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
                +xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))
                ))))));

        return (spiciness0);
}
/**
==========================================================================
method gsw_spiciness1(sa,ct)
==========================================================================
!
!  Calculates spiciness from Absolute Salinity and Conservative
!  Temperature at a pressure of 1000 dbar, as described by McDougall and
!  Krzysik (2015).  This routine is based on the computationally-efficient
!  expression for specific volume in terms of SA, CT and p (Roquet et al.,
!  2015).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  spiciness1  =  spiciness referenced to a pressure of 1000 dbar [ kg/m^3 ]
!
*/
double TeosSea::gsw_spiciness1(double sa, double ct)
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

        double  xs, ys, spiciness1;

        xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys      = ct*0.025;

        spiciness1= s01+ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
                +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
                +xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
                +xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
                +xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
                +xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
                +xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))
                ))))));

        return (spiciness1);
}
/**
==========================================================================
method gsw_spiciness2(sa,ct)
==========================================================================
!
!  Calculates spiciness from Absolute Salinity and Conservative
!  Temperature at a pressure of 2000 dbar, as described by McDougall and
!  Krzysik (2015).  This routine is based on the computationally-efficient
!  expression for specific volume in terms of SA, CT and p (Roquet et al.,
!  2015).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  spiciness2  =  spiciness referenced to a pressure of 2000 dbar [ kg/m^3 ]
!
*/
double TeosSea::gsw_spiciness2(double sa, double ct)
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

        double  xs, ys, spiciness2;

        xs      = sqrt(gtc.gsw_sfac*sa + gtc.offset);
        ys      = ct*0.025;

        spiciness2= s01+ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
                +xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
                +xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
                +xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
                +xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
                +xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
                +xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))
                ))))));

        return (spiciness2);
}

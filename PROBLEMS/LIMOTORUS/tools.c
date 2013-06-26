#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>


ldouble
get_rho(ldouble r,ldouble th)
{
  ldouble rin=150.;
  ldouble rhoin = 2.02e5*pow(1.-pow((th-M_PI/2.)/(M_PI/2.),2.),1.69) * (150./rin);
  ldouble rho = rhoin * (rin/r);

  return rhoCGS2GU(rho);
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//sets rho, uint and velocities according to disk model from Sadowski+2012
int
set_sgradisk(ldouble *pp,ldouble *xx,void *ggg, void* gggBL)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  struct geometry *geomBL
    = (struct geometry *) gggBL;

  // spherical coordinates
  ldouble xx2[4];
  coco_N(xx,xx2,MYCOORDS,BLCOORDS);
  ldouble r=xx2[1];
  ldouble th=xx2[2];

  //rotation here?

  //entropy = fixed at rin
  ldouble xxin[4]={0,get_x(0,0),get_x(NY/2,1),get_x(0,2)};
  coco_N(xxin,xxin,MYCOORDS,BLCOORDS);
  ldouble rin=xxin[1];

  rin=150.;
  ldouble rhoin = 2.02e5*pow(1.-pow((th-M_PI/2.)/(M_PI/2.),2.),1.69) * (150./rin);
  ldouble tempin = pow(10.,9.95+0.24*pow(fabs(th-M_PI/2.),2.93)) * (150./rin);
  //to code units
  rhoin = rhoCGS2GU(rhoin);
  tempin = tempCGS2GU(tempin);
  ldouble uin = calc_PEQ_ufromTrho(tempin,rhoin);
  ldouble K = (GAMMA-1.)*uin / pow(rhoin,GAMMA); //p = K rho^GAMMA

  //empirical fit in cgs
  ldouble rho = rhoin * (rin/r);

  //polytrope G=2
  rho = 1./2./K * 0.7 / r;

  ldouble temp = tempin * (150./r);
  ldouble vphi = pow(10.,9.15-0.24*pow(fabs(th-M_PI/2.),2.04)) * sqrt(150./r);

  //const.fraction of Keplerian:
  vphi=0.7*sqrt(1./r/r/r);

  ldouble chi = 0.1 + 0.31*pow(fabs(th-M_PI/2.),3.89); //pmag/pgas


  //to code units
  //rho = rhoCGS2GU(rho);
  temp = tempCGS2GU(temp);
  vphi = velCGS2GU(vphi);
  chi = chi;

  //temp following constant entropy: - actually should be satisfied already for GAMMA=2
  ldouble u = K * pow(rho,GAMMA) / (GAMMA-1.);
  temp = calc_PEQ_Tfromurho(u,rho); 
 
  pp[0] = rho;
  pp[1] = u;//calc_PEQ_ufromTrho(temp,rho);

  //to add extra magn-related pressure
  //pp[1] *= 1.+chi;

  ldouble ucon[4]={0.,0.,0.,vphi/r};

  conv_vels(ucon,ucon,VEL3,VEL4,geomBL->gg,geomBL->GG);
  trans2_coco(geomBL->xxvec,ucon,ucon,KERRCOORDS,MYCOORDS);
  conv_vels(ucon,ucon,VEL4,VELPRIM,geom->gg,geom->GG);

  pp[2]=ucon[1];
  pp[3]=ucon[2];
  pp[4]=ucon[3];

  return 0;

}


FTYPE    gam = GAMMA;
FTYPE    rin = TORUSRIN;
FTYPE    kappa = TORUSKAPPA;  // AKMARK: entropy constant that appears in EOS
FTYPE    xi = TORUSXI;   // AKMARK: omega is set to this fraction of Keplerian omega
FTYPE    rbreak1 = TORUSRBREAK1;   // AKMARK: locations of breaks in torus angular momentum profile
FTYPE    rbreak2 = TORUSRBREAK2;



void compute_gd( FTYPE r, FTYPE th, FTYPE a, FTYPE *gdtt, FTYPE *gdtp, FTYPE *gdpp ) {
   FTYPE Sigma, tmp;

   Sigma = (pow(a*cos(th),2) + pow(r,2));
   tmp = 2*r*pow(sin(th),2)/Sigma;

   //metric (expressions taken from limotorus4.nb):
   *gdtt = -1 + 2*r/Sigma;
   *gdtp = -a*tmp;
   *gdpp = (pow(r,2) + pow(a,2)*(1+tmp)) * pow(sin(th),2);
}

// Keplerian, equatorial angular momentum density: l = u_phi / u_t
FTYPE lK(FTYPE r, FTYPE a) {
   FTYPE curlyF, curlyG;

   curlyF = 1 - 2*a/pow(r, 1.5) + pow(a/r, 2);
   curlyG = 1 - 2/r + a/pow(r, 1.5);
   return( pow(r,0.5) * curlyF / curlyG );
}

FTYPE l3d(FTYPE lam, FTYPE a, FTYPE lambreak1, FTYPE lambreak2, FTYPE xi) {
   return ( xi * lK( lam<=lambreak1 ? lambreak1 : lam>=lambreak2 ? lambreak2 : lam , a) );
}

FTYPE rtbis(FTYPE (*func)(FTYPE,FTYPE*), FTYPE *parms, FTYPE x1, FTYPE x2, FTYPE xacc)
//Taken from HARM:nrutil.c
{
   int j;
   FTYPE dx,f,fmid,xmid,rtb;
   f=(*func)(x1, parms);
   fmid=(*func)(x2, parms);
   if (f*fmid >= 0.0) {
      printf( "f(%g)=%g f(%g)=%g\n", x1, f, x2, fmid );
      printf("Root must be bracketed for bisection in rtbis\n");
   }
   rtb = (f < 0.0) ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so that f>0 lies at x+dx.
   for (j=1;j<=100;j++) {
      fmid=(*func)(xmid=rtb+(dx *= 0.5),parms); //Bisection loop.
      if (fmid <= 0.0) {
         rtb=xmid;
      }
      if (fabs(dx) < xacc || fmid == 0.0) {
         return rtb;
      }
   }
   printf("Too many bisections in rtbis\n");
   return 0.0; //Never get here.
}

FTYPE lamBL_func(FTYPE lam, FTYPE *parms) {
   FTYPE gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi, l;

   gdtt      = parms[0];
   gdtp      = parms[1];
   gdpp      = parms[2];
   a         = parms[3];
   lambreak1 = parms[4];
   lambreak2 = parms[5];
   xi        = parms[6];

   l = l3d(lam, a, lambreak1, lambreak2, xi);

   return ( lam*lam + l * (l*gdtp + gdpp) / (l*gdtt + gdtp) );
}

// von Zeipel cylinder radius (needs to be calculated iteratively for nonzero spin)
FTYPE lamBL(FTYPE R, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp, FTYPE a, FTYPE lambreak1, FTYPE lambreak2, FTYPE xi) {
   // R = r*sin(th), used as initial guess for lamBL

   FTYPE parms[7];

   //store params in an array before function call
   parms[0] = gdtt;
   parms[1] = gdtp;
   parms[2] = gdpp;
   parms[3] = a;
   parms[4] = lambreak1;
   parms[5] = lambreak2;
   parms[6] = xi;

   //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3)
   //demand accuracy 5x machine prec.
   //in non-rel limit rml = r , use 10x that as the upper limit:
   return( rtbis( &lamBL_func, parms, R, 10*R, 5.*DBL_EPSILON ) );

}

FTYPE omega3d( FTYPE l, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp ) {
   return( -(gdtt*l + gdtp)*pow(gdpp + gdtp*l,-1) );
}

FTYPE compute_Agrav( FTYPE om, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp ){
   return (sqrt(fabs(1./ ( gdtt + 2*om*gdtp + pow(om,2)*gdpp ) )));
}

FTYPE rmidlam( FTYPE x, FTYPE *parms ) {
   FTYPE lamsq, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi;
   FTYPE lam_x, ans;

   lamsq     = parms[0];   // square of target lambda
   gdtt      = parms[1];
   gdtp      = parms[2];
   gdpp      = parms[3];
   a         = parms[4];
   lambreak1 = parms[5];
   lambreak2 = parms[6];
   xi        = parms[7];

   compute_gd(x, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lam_x = lamBL(x, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);   // lambda at current value of x

   ans = lamsq - lam_x*lam_x;

   return(ans);
}

FTYPE limotorus_findrml(FTYPE lam, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp, FTYPE a, FTYPE lambreak1, FTYPE lambreak2, FTYPE xi) {
   FTYPE parms[8];
   FTYPE rml;

   //store params in an array before function call
   parms[0] = lam*lam;
   parms[1] = gdtt;
   parms[2] = gdtp;
   parms[3] = gdpp;
   parms[4] = a;
   parms[5] = lambreak1;
   parms[6] = lambreak2;
   parms[7] = xi;

   //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3)
   //demand accuracy 5x machine prec.
   //in non-rel limit rml = r , use 10x that as the upper limit:
   rml = rtbis( &rmidlam, parms, 6, 10*parms[0], 5.*DBL_EPSILON );

   return( rml );
}

struct lnf_int_params {
   FTYPE a, lambreak1, lambreak2, xi;
};

// Integrand in expression for lnf. This function is called
// by the GSL quadrature routine.
FTYPE lnf_integrand(FTYPE r, void *params) {
   struct lnf_int_params *pars = (struct lnf_int_params *) params;
   FTYPE a = pars->a, lambreak1 = pars->lambreak1, lambreak2 = pars->lambreak2, xi = pars->xi;
   FTYPE r2, r3, a2, lam, lam2, lamroot;
   FTYPE gdtt, gdtp, gdpp;
   FTYPE l, om, dl_dr, dom_dr, integrand;
   FTYPE term1, term2, term3, term4, D, E;
   FTYPE oneplusx, dx_dlam, dlK_dlam, om_numerator, om_denominator;

   // all values below are midplane values
   r2 = pow(r,2);
   r3 = pow(r,3);
   a2 = a*a;
   compute_gd(r, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lam = lamBL(r, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);
   lam2 = lam*lam;

   l = l3d(lam, a, lambreak1, lambreak2, xi);

   term1 = (r-2) * l;
   term2 = 2*a;
   term3 = r3 + a2*(r+2);
   term4 = term2 * l;

   om_numerator = term1 + term2;
   om_denominator = term3 - term4;
   om = om_numerator / om_denominator;

   // derivatives
   if (lam <= lambreak1 || lam >= lambreak2) {
      dl_dr = 0;
   } else {
      oneplusx = 1 - 2/lam + a*pow(lam,-1.5);
      dx_dlam = 2*pow(lam,-2) - 1.5*a*pow(lam,-2.5);
      lamroot = sqrt(lam);
      dlK_dlam = (oneplusx*0.5*pow(lamroot,-1) + (a-lamroot)*dx_dlam) / pow(oneplusx,2);
      D = term3 - 2*term4 - lam2 * (r-2);
      E = l * (3*r2 + a2 - lam2);
      dl_dr = E / ( 2*lam*om_numerator / (xi*dlK_dlam) - D);
   }

   dom_dr = ( om_denominator * (l + (r-2)*dl_dr) - om_numerator * (3*r2+a2+2*a*dl_dr) )
            / pow(om_denominator, 2);

   integrand = -l/(1-om*l) * dom_dr;

   return( integrand );
}

int init_dsandvels_limotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell, FTYPE* omret)
{
   FTYPE hh, eps;
   FTYPE rho, u, ur, uh, up;
   int pl;
   FTYPE R;
   FTYPE lambreak1, lambreak2;
   FTYPE lam, rml, l, om, Agrav, f3d;
   FTYPE lamin, lin, omin, Agravin, f3din;
   FTYPE lnf3d, lnferr;
   FTYPE gdtt, gdtp, gdpp;

   // GSL quadrature stuff
   gsl_integration_workspace *intwork;
   int worksize = 1000;   // workspace size
   gsl_function integrand;
   int intstatus;

   R = r*sin(th);
   if (R < rin) {*rhoout = 0.; *uuout =0.;return(0);}

  

   ///
   /// Computations at break radii
   ///

   // AKMARK: lamBL can calculate lambda at an arbitrary location, but it needs lambreak1,2;
   // how to calculate lambreak1,2 themselves?
   // Solution: note that lambreak1,2 are only needed by "l3d" in order to determine which region of the torus we are in.
   // But lambreak1,2 can be considered to be in region 2, so just need a way to make "l3d" believe that we are in region 2.
   // This can be done by setting lambreak1,2 to very small and very large values respectively.
   compute_gd(rbreak1, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lambreak1 = lamBL(rbreak1, gdtt, gdtp, gdpp, a, 0, 200000, xi);
   compute_gd(rbreak2, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lambreak2 = lamBL(rbreak2, gdtt, gdtp, gdpp, a, lambreak1, 200000, xi);

   ///
   /// Computations at torus inner edge
   ///

   compute_gd(rin, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lamin = lamBL(rin, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);
   lin = l3d(lamin, a, lambreak1, lambreak2, xi);
   omin = omega3d(lin, gdtt, gdtp, gdpp);
   Agravin = compute_Agrav(omin, gdtt, gdtp, gdpp);
   f3din = 1.;   // the way f3d is defined, it equals 1 at r=rin

   ///
   /// Computations at current point: r, th
   ///

   compute_gd(r, th, a, &gdtt, &gdtp, &gdpp);
   lam = lamBL(R, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);
   l = l3d(lam, a, lambreak1, lambreak2, xi);
   om = omega3d(l, gdtt, gdtp, gdpp);
   Agrav = compute_Agrav(om, gdtt, gdtp, gdpp);

   rml = limotorus_findrml( lam, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi );

   // f3d requires numerical evaluation of an integral

   // First, set up things that GSL quadrature routine requires
   struct lnf_int_params pars = {a, lambreak1, lambreak2, xi};
   integrand.function = &lnf_integrand;   // integrand: function
   integrand.params = &pars;   // other parameters to pass to integrand (besides integration variable)
   intwork = gsl_integration_workspace_alloc(worksize);   // integration workspace

   // then perform integration to obtain ln(f3d) ...
   intstatus = gsl_integration_qags(&integrand, rin, rml, 0., 1.e-8, worksize, intwork, &lnf3d, &lnferr);
   gsl_integration_workspace_free( intwork );
   if (intstatus != GSL_SUCCESS) {
      printf("GSL integration failed during limotorus setup at r=%21.15g, th=%21.15g; setting density to zero at this point", r, th);
      lnf3d = GSL_NEGINF;   // cause density to be 0 at the current point
   }

   // ... and finally, f3d
   f3d=exp(lnf3d);

   //hh = w = FA / (FA)_in
   hh = f3d*Agrav / (f3din*Agravin);

#ifdef IMPOSEDRHO
   //rho fixed
   //uint/rho
   if(hh<1.) 
     {
       *uuout = 0.;
       rho= 0.;
     }
   else
     {
       rho = get_rho(r,th);
       eps = (-1 + hh)*pow(gam,-1);
       *uuout = rho * eps;
       //*uuout = rho * 1./gam/(gam-1.)*(pow(hh,gam/(gam-1.))-1.);    
     }

   *rhoout = rho;
   *ell=l;
   *omret=om;
#else
   //entropy constant
   eps = (-1 + hh)*pow(gam,-1);

   if (eps < 0) rho = 0; else rho = pow((-1 + gam)*eps*pow(kappa,-1),pow(-1 + gam,-1));

   *rhoout = rho;
   *ell=l;
   *omret=om;
   *uuout = kappa * pow(rho, gam) / (gam - 1.);
#endif

   return(0);

}



int
donut_analytical_solution(ldouble *pp, ldouble *xxvecBL,ldouble ggBL[][5],ldouble GGBL[][5] )
{
  ldouble r=xxvecBL[1];
  ldouble th=xxvecBL[2];
  
  ldouble rho,uu,ell,om;
  init_dsandvels_limotorus(r, th, BHSPIN, &rho, &uu, &ell,&om);
 

  if(rho<=0.) return -1.; //outside torus
  //printf("%e %e %e\n",rho,uu,uu/pow(rho,gam)); getchar();

  //3-velocity in BL 
  ldouble ucon[4]={0.,0.,0.,om};
  conv_vels(ucon,ucon,VEL3,VELPRIM,ggBL,GGBL);
  
  pp[2]=ucon[1];
  pp[3]=ucon[2];
  pp[4]=ucon[3];

  pp[0]=rho;
  pp[1]=uu;


  return 0;
}

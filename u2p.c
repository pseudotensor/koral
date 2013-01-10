
//KORAL - u2p.c
//conserved to primitives conversion

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//high-level u2p solver
int
u2p(ldouble *uu, ldouble *pp, ldouble gg[][5],ldouble GG[][5],int *corrected)
{
  *corrected=0;
  int verbose=0;
  int hdcorr=0;
  int radcorr=0;

  int u2pret,ret;
  ldouble ppbak[NV];
  for(u2pret=0;u2pret<NV;u2pret++)
    ppbak[u2pret]=pp[u2pret];


  //************************************
  //************************************
  //************************************
  //hydro part
  //************************************
  //************************************
  //************************************

  ldouble u0=pp[1];
  
  //************************************
  //hot hydro - conserving energy
  ret=0;
#ifdef U2P_NUMTEMP
    u2pret=u2p_hot_gsl(uu,pp,gg);  //temporary 3D solver - not perfect - does not work for RADATM!
#else
    u2pret=u2p_hot(uu,pp,gg);  //TODO: to be replaced - catastrophic cancelation!
#endif
  //************************************

  if(u2pret<0) 
    {
      if(verbose>0)
	printf("u2p_hot err  >>> %d <<< %Le %Le\n",u2pret,pp[0],pp[1]);
      
      //************************************
      //entropy solver - conserving entropy
      ret=-1;
      u2pret=u2p_entropy(uu,pp,gg);
      //************************************

      if(verbose>0)
	  printf("u2p_entr     >>> %d <<< %Le > %Le\n",u2pret,u0,pp[1]);

      if(u2pret<0)
	{
	  if(verbose>0)
	    printf("u2p_entr err > %Le %Le\n",pp[0],pp[1]);
	
	  //************************************
	  //cold RHD - assuming u=SMALL
	  ret=-2;
	  u2pret=u2p_cold(uu,pp,gg);
	  //************************************

	  if(u2pret<0)
	    {
	      if(verbose>0)
		printf("u2p_cold err > %Le %Le\n",pp[0],pp[1]);
	
	      //************************************
	      //leaving unchanged primitives - should not happen
	      ret=-3;
	      for(u2pret=0;u2pret<NV;u2pret++)
		pp[u2pret]=ppbak[u2pret];	  
	      //************************************

	    }

	}
      
     }

  //************************************
  //************************************
  //************************************
  //radiation part
  //************************************
  //************************************
  //************************************

#ifdef RADIATION
  int radcor;
  u2p_rad(uu,pp,gg,GG,&radcorr);
#endif
  
  if(radcorr>0 || hdcorr>0) *corrected=1;

  return ret;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver
//'hot grhd' - pure hydro, numerical in 5D
//but currently work only for 1D problems
//TODO: to be replaced by something much better 
int
print_state (int iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3d \n"
	  "x = % .15e  % .15e  % .15e  % .15e  % .15e  "
	  "f(x) = % .15e  % .15e  % .15e  % .15e  % .15e\n",
	  iter,
	  gsl_vector_get (s->x, 0),
	  gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2),
	  gsl_vector_get (s->x, 3),
	  gsl_vector_get (s->x, 4),
	  gsl_vector_get (s->f, 0),
	  gsl_vector_get (s->f, 1),
	  gsl_vector_get (s->f, 2),
	  gsl_vector_get (s->f, 3),
	  gsl_vector_get (s->f, 4));
}

struct u2photpar
{
  ldouble *uuu;
  //  double g[4][5];
};
     
int
f_u2p_hot_gsl(const gsl_vector * x, void *params,
	      gsl_vector * f)
{
  ldouble *uuu = ((struct u2photpar *) params)->uuu;
    
  ldouble rho = (ldouble)gsl_vector_get (x, 0);
  ldouble uu = (ldouble)gsl_vector_get (x, 1);
  ldouble vr = (ldouble)gsl_vector_get (x, 2);
  ldouble vph = 0.*(ldouble)gsl_vector_get (x, 3);
  ldouble vth = 0.*(ldouble)gsl_vector_get (x, 4);
     
  //flat
  ldouble gtt=-1;
  ldouble gtph=1.;
  ldouble gphph=1.;
  ldouble gthth=1.;
  ldouble grr=1.;

  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vph*vph*gphph + vth*vth*gthth);

  if(ut2<0.)
    {
      my_err("ut2.lt.0 in p2u\n"); ut2=0.;
    }

  ldouble ut=sqrtl(ut2);
  ldouble rhout = rho*ut;
  ldouble Sut;

  ldouble Tttt=rhout*(1+ut*(gtt+vph*gtph))+GAMMA*uu*ut2*(gtt+vph*gtph)+uu*(GAMMA-1.);  
  ldouble Ttr=(rho+GAMMA*uu)*ut2*vr*grr;
  ldouble Ttth=(rho+GAMMA*uu)*ut2*vth*gthth;
  ldouble Ttph=(rho+GAMMA*uu)*ut2*(gtph+vph*gphph);

     
  gsl_vector_set (f, 0, (double)(rhout-uuu[0]));
  gsl_vector_set (f, 1, (double)(Tttt-uuu[1]));
  gsl_vector_set (f, 2, (double)(Ttr-uuu[2]));
  //gsl_vector_set (f, 3, (double)(Ttth-uuu[3]));
  //gsl_vector_set (f, 4, (double)(Ttph-uuu[4]));
     
  return GSL_SUCCESS;
}


int
u2p_hot_gsl(ldouble *uuu, ldouble *p, ldouble g[][5])
{

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
     
  int verbose=0;
  int status;
  int iter = 0;
     
  const size_t n = 3;
  struct u2photpar par = {uuu};
  gsl_multiroot_function f = {&f_u2p_hot_gsl, n, &par};
     
  double x_init[5] = {p[0],p[1],p[2],p[3],p[4]};
  gsl_vector *x = gsl_vector_alloc (n);
     
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  gsl_vector_set (x, 2, x_init[2]);
  //gsl_vector_set (x, 3, x_init[3]);
  //gsl_vector_set (x, 4, x_init[4]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);
     
  if(verbose>0) print_state (iter, s);
     
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      
      if(verbose>0) print_state (iter, s);
      
      if (status)   /* check if solver is stuck */
	break;
      
      status =
	gsl_multiroot_test_residual (s->f, 1.e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  if(verbose>0) printf ("status = %s\n", gsl_strerror (status));
  
  ldouble rhout=uuu[0];
  ldouble Tttt=uuu[1];
  ldouble Ttr=uuu[2];
  ldouble Ttth=uuu[3];
  ldouble Ttph=uuu[4];
  ldouble Sut=uuu[5];

  p[0]=gsl_vector_get (s->x,0);
  p[1]=gsl_vector_get (s->x,1);
  p[2]=gsl_vector_get (s->x,2);
  //p[3]=gsl_vector_get (x,3);
  //p[4]=gsl_vector_get (x,4);
  p[3]=p[4]=0.;

  ldouble ut=rhout/p[0];

  p[5]=Sut/ut;

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver
//'hot grhd' - pure hydro, numerical in 1D
//catastrophic cancelation!
int
u2p_hot(ldouble *uuu, ldouble *p, ldouble g[][5])
{
  int verbose=0;

  ldouble gtt=g[0][0];
  ldouble gtph=g[0][3];
  ldouble grr=g[1][1];
  ldouble gthth=g[2][2];
  ldouble gphph=g[3][3];

  ldouble gdet=g[3][4];

  ldouble rhout=uuu[0];
  ldouble Tttt=uuu[1];
  ldouble Ttr=uuu[2];
  ldouble Ttth=uuu[3];
  ldouble Ttph=uuu[4];
  ldouble Sut=uuu[5];
 
  ldouble fff1,fff2,fff3,aaa,bbb,ccc,sqrt1,ut,rho,vph,vr,vth,daaadu,dcccdu,dutdu,dvrdu,dvthdu,dvphdu,ut2,del,S;
  ldouble u,up1,um1,fval,dfval,diffu,utm1;
  int iter=0;ldouble err;
 
  //initial guess
  u=up1=um1=p[1];

  ldouble conv=U2PPREC;
  int itmax=30;
  ldouble fvalmin[2]={0.,-1.};
  ldouble absfval;

  do{
    iter++;
    u=up1;
    
    //this mess is here to make it faster
    fff1=(gtt-gtph*gtph/gphph);

    aaa=GAMMA*u*fff1;
    bbb=rhout*fff1;
    ccc=Ttph*gtph/gphph + rhout + u*(GAMMA-1.) - Tttt;
     
    del=bbb*bbb-4.*aaa*ccc;

    //calculates all terms only if they make sense
    if(del>=0.)
      {	
	sqrt1=sqrtl(del);
	ut=(-bbb-sqrt1)/2./aaa;
	if(ut==0.) 
	  {
	    return -11;
	    my_err("ut zeroed\n");
	    ut=1.;
	  }

	ut2=ut*ut;
	rho=rhout/ut;
	S=Sut/ut;
	fff2=ut2*(u*GAMMA + rho);
    
	vph=-(gtph/gphph) + Ttph/(gphph*fff2);
	vr=Ttr/(grr*fff2);
	vth=Ttth/(gthth*fff2);

	daaadu=GAMMA*fff1;
	dcccdu=GAMMA-1.;

	dutdu=daaadu*(bbb/2./aaa/aaa + ccc/aaa/sqrt1 + sqrt1/2./aaa/aaa) + dcccdu/sqrt1;

	fff3=dutdu*(rhout+2.*GAMMA*u*ut);

	dvrdu=-Ttr/(grr*fff2)/(grr*fff2)*
	  (fff3*grr + GAMMA*ut2*grr);
	dvthdu=-Ttth/(gthth*fff2)/(gthth*fff2)*
	  (fff3*gthth + GAMMA*ut2*gthth);
	dvphdu=-Ttph/(gphph*fff2)/(gphph*fff2)*
	  (fff3*gphph + GAMMA*ut2*gphph);
  
	dfval=-2./ut/ut/ut*dutdu + 2.*grr*vr*dvrdu + 2.*gthth*vth*dvthdu + 2.*gtph*dvphdu + 2.*gphph*vph*dvphdu;
	fval=1./ut/ut+gtt+grr*vr*vr+gthth*vth*vth+2.*gtph*vph+gphph*vph*vph;	
      }

    //handling unphysical input    
    if(del<0. || isnan(dfval) || isnan(fval) || isinf(del))
      {
	if(verbose>0)
	  {
	    printf("unphysical input in u2p at iter: %d\n",iter);
	    printf("u: %Le (%Le, %Le)\n",u,p[1],um1);
	    printf(" %Le %Le %Le %Le\n",del,ut,dfval,fval);
	    getchar();
	  }

	if(iter==1)
	  return -1; //u2p failed due to unphysical u0
	else
	  {
	    return -2; //nan got lateron
	    up1=0.5*(u+um1); //going closer to the previous, successful step
	    continue;
	  }
      }
    
    //input ok, let's carry on

     //absolute error
    if(fval > 0)
      absfval=fval;
    else
      absfval=-fval;

    //putting best value of u to memory
    if(absfval<fvalmin[1] || fvalmin[1]<0.)
      {
	fvalmin[0]=u;
	fvalmin[1]=absfval;
      }

    //Newton
    up1=u-fval/dfval;   

    //difference
    diffu=up1-u;
    
    if(verbose==1)
      {
	printf("       iter %d  > %Le [%Le] %Le >%Le< %Lf\n",iter,u,up1,diffu/u,fval,ut);
	getchar();
      }

    //check
    if(isnan(up1)) 
      {
	printf("  u2p> ======= nan \n");
	if(iter==1)
	  return -3; //up1 nan lateron at iter1
	else
	  {
	    my_err("nan met in u2p\n");
	    return -4; //up1 nan lateron at iter.gt.1

	    up1=0.5*(u+um1); //going closer to the previous, successful step
	    continue;
	  }
      }
    
    err =diffu/u;
    if(err<0.) err*=-1.;
    
    //putting previous step in memory
    um1=u;

  } while(err>conv && iter<itmax);

  //not tried even single value - should not happen - does not change primitives
  if(fvalmin[1]<0.)
    {
      printf("unphysical in u2p at iter: %d\n",iter);
      printf("u: %Le (%Le, %Le)\n",u,p[1],um1);
      printf(" %Le %Le %Le %Le\n",del,ut,dfval,fval);

      my_err("not tried even single value in u2p\n");      
      return -6; //not tried even single value in u2p
    }
 
  //recalculate primitives for u with smallest fval
  if(iter>=itmax)
    {
      return -5; //iteration exceeded
      
      fff1=(gtt-gtph*gtph/gphph);
      aaa=GAMMA*u*fff1;
      bbb=rhout*fff1;
      ccc=Ttph*gtph/gphph + rhout + u*(GAMMA-1.) - Tttt;    
      del=bbb*bbb-4.*aaa*ccc;
      sqrt1=sqrtl(del);
      ut=(-bbb-sqrt1)/2./aaa;
      ut2=ut*ut;
      rho=rhout/ut;
      S=Sut/ut;
      fff2=ut2*(u*GAMMA + rho);
      vph=-(gtph/gphph) + Ttph/(gphph*fff2);
      vr=Ttr/(grr*fff2);
      vth=Ttth/(gthth*fff2);
    }

  if(verbose==1)
    {
      printf("  u/u0   > %Le %Le %Le\n",(u)/p[1],dfval,fval);
      printf("  u,rho,i  > %Le %Le %d\n",u,rho,iter);
      printf("  u2p> =========\n"); getchar();
    }

  if((isnan(u) || isnan(rho)))
    my_err("nan in u2p\n");

  if(u<0. || rho<0.) return -10;

 
  //primitives
  p[0]=rho;
  p[1]=u;
  p[2]=vr;
  p[3]=vth;
  p[4]=vph;
  p[5]=S;

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//auxiliary solver based on the entropy conservation
int
u2p_entropy(ldouble *uuu, ldouble *p, ldouble g[][5])
{
  ldouble gtt=g[0][0];
  ldouble gtph=g[0][3];
  ldouble grr=g[1][1];
  ldouble gthth=g[2][2];
  ldouble gphph=g[3][3];

  ldouble gdet=g[3][4];

  ldouble rhout=uuu[0];
  ldouble Tttt=uuu[1];
  ldouble Ttr=uuu[2];
  ldouble Ttth=uuu[3];
  ldouble Ttph=uuu[4];
  ldouble Sut=uuu[5];

  ldouble rho=p[0]; //initial guess
  ldouble uu=p[1];
  ldouble vr=p[2];
  ldouble vth=p[3];
  ldouble vph=p[4];
  ldouble S=p[5];

  ldouble ut2;
  
  ldouble rhop1=rho;
  ldouble rhom1=rho;  

  int iter=0;
  ldouble err,fff2,fval,dfval,fval1,fval2,diffrho,ut;
  ldouble dudrho,dfdrho,dvphdrho,dvrdrho,dvthdrho;

  ldouble conv=U2PPREC;
  conv=1.e-10;
  int itmax=30;
  ldouble fvalmin[2]={0.,-1.};
  ldouble absfval;

  ldouble ftest[50][3];

  do{
    iter++;
    rho=rhop1;
    ut=rhout/rho;
    ut2=ut*ut;
    S=Sut/ut;
   
    //for some reason powl sometimes breaks
    uu=(ldouble)pow((pow(rho,1./GAMMAM1+1.)*exp(Sut/rhout)),GAMMAM1)/GAMMAM1;
 
    fff2=ut2*(uu*GAMMA + rho);  
    vph=-(gtph/gphph) + Ttph/(gphph*fff2);
    vr=Ttr/(grr*fff2);
    vth=Ttth/(gthth*fff2); 
    fval=1./ut2 + gtt + grr*vr*vr + gthth*vth*vth + gphph*vph*vph + 2*gtph*vph;

    dudrho=(ldouble)GAMMA/GAMMAM1*exp(Sut/rhout)*pow(rho,1./GAMMAM1)*pow(exp(Sut/rhout)*pow(rho,GAMMA/GAMMAM1),GAMMA-2.);
    dfdrho=(rhout*rhout*(rho*(GAMMA*dudrho-1.)-2.*GAMMA*uu))/rho/rho/rho;
    dvrdrho=-Ttr/grr/fff2/fff2*dfdrho;
    dvthdrho=-Ttth/gthth/fff2/fff2*dfdrho;
    dvphdrho=-Ttph/gphph/fff2/fff2*dfdrho;
    dfval=2.*rho/rhout/rhout + 2.*grr*vr*dvrdrho + 2.*gthth*vth*dvthdrho + 2.*gphph*vph*dvphdrho + 2.*gtph*dvphdrho;

    //absolute error
    if(fval > 0)
      absfval=fval;
    else
      absfval=-fval;

    //putting best value of u to memory
    if(absfval<fvalmin[1] || fvalmin[1]<0.)
      {
	fvalmin[0]=ut;
	fvalmin[1]=absfval;
      }

    //Newton
    rhop1=rho-fval/dfval;   

    diffrho=rhop1-rho;
    
    ftest[iter][0]=rho;
    ftest[iter][1]=uu;
    ftest[iter][2]=fval;
    ftest[iter][3]=dfval;

    err =diffrho/rho;
    if(err<0.) err*=-1.;
 
    if(iter>itmax && err>conv) 
      {
	printf("iter exceeded in u2p_entr \n");
	printf(" entr  iter %d> %Le [%Le] %Le >%Le< %Le\n",iter,rho,rhop1,diffrho/rho,fval,err);
	return -1;
      }

    rhom1=rho;

  } while(err>conv);

  rho=rhop1;
  ut=rhout/rho;
  ut2=ut*ut;
  S=Sut/ut;
  uu=(ldouble)pow((pow(rho,1./GAMMAM1+1.)*exp(Sut/rhout)),GAMMAM1)/GAMMAM1;
  fff2=ut2*(uu*GAMMA + rho);  
  vph=-(gtph/gphph) + Ttph/(gphph*fff2);
  vr=Ttr/(grr*fff2);
  vth=Ttth/(gthth*fff2); 

  if(uu<0. || rho<0.)
    {
      printf("u2p_entr didn't work: %Le %Le\n",uu,rho);
      return -1;
    }

  p[0]=rho;
  p[1]=uu;
  p[2]=vr;
  p[3]=vth;
  p[4]=vph;
  p[5]=S;

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//auxiliary solver assuming u=0
int
u2p_cold(ldouble *uuu, ldouble *p, ldouble g[][5])
{
  ldouble gtt=g[0][0];
  ldouble gtph=g[0][3];
  ldouble grr=g[1][1];
  ldouble gthth=g[2][2];
  ldouble gphph=g[3][3];

  ldouble gdet=g[3][4];

  ldouble rhout=uuu[0];
  ldouble Ttt=uuu[1]-rhout;
  ldouble Ttr=uuu[2];
  ldouble Ttth=uuu[3];
  ldouble Ttph=uuu[4];
  ldouble Sut=uuu[5];

  ldouble rho,uu,vr,vth,vph,S,fff,ut;
  
  fff=(Ttt - Ttph*gtph/gphph)/(gtt - gtph*gtph/gphph);
  ut=fff/rhout;
					       
  rho=rhout/ut;
  uu=UFLOOR;
  vr=Ttr/grr/fff;
  vph=Ttph/gphph/fff - gtph/gphph;
  vth=Ttth/gthth/fff;
  S=Sut/ut;

  if(rho<0.)
    {
      rho=RHOFLOOR;
   }

  p[0]=rho;
  p[1]=uu;
  p[2]=vr;
  p[3]=vth;
  p[4]=vph;
  p[5]=S;

  return 0;
}

//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//radiative primitives: (E,\tilde u^i)
//  E - radiative energy density in the rad.rest frame
//  u^i - relative velocity of the rad.rest frame
//**********************************************************************
//**********************************************************************
int
u2p_rad(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5], int *corrected)
{
  //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
  *corrected=0;

  int verbose=0,i;
  ldouble Rij[4][4];

  //conserved - R^t_mu
  ldouble Av[4]={uu[6],uu[7],uu[8],uu[9]};
  //indices up - R^tmu
  indices_12(Av,Av,GG);

  //g_munu R^tmu R^tnu
  ldouble gRR=gg[0][0]*Av[0]*Av[0]+gg[0][1]*Av[0]*Av[1]+gg[0][2]*Av[0]*Av[2]+gg[0][3]*Av[0]*Av[3]+
    gg[1][0]*Av[1]*Av[0]+gg[1][1]*Av[1]*Av[1]+gg[1][2]*Av[1]*Av[2]+gg[1][3]*Av[1]*Av[3]+
    gg[2][0]*Av[2]*Av[0]+gg[2][1]*Av[2]*Av[1]+gg[2][2]*Av[2]*Av[2]+gg[2][3]*Av[2]*Av[3]+
    gg[3][0]*Av[3]*Av[0]+gg[3][1]*Av[3]*Av[1]+gg[3][2]*Av[3]*Av[2]+gg[3][3]*Av[3]*Av[3];
 
  //the quadratic equation for u^t of the radiation rest frame (urf[0])
  //supposed to provide two roots for (u^t)^2 of opposite signs
  ldouble a,b,c,delta,gamma2;
  ldouble urfcon[4],urfcov[4],Erf;
  a=16.*gRR;
  b=8.*(gRR*GG[0][0]+Av[0]*Av[0]);
  c=gRR*GG[0][0]*GG[0][0]-Av[0]*Av[0]*GG[0][0];
  delta=b*b-4.*a*c;
  gamma2=  (-b-sqrtl(delta))/2./a;
  //if unphysical try the other root
  if(gamma2<1.) gamma2=  (-b+sqrtl(delta))/2./a; 

  //cap on u^t
  ldouble gammamax=1000.;
 
  if(gamma2<0 || gamma2>gammamax*gammamax || delta<0.) 
    {
      //top cap
      *corrected=1;
      urfcon[0]=gammamax;
      
      //proper direction for the radiation rest frame, will be normalized later      
      Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);

      ldouble Arad[4];
      for(i=1;i<4;i++)
	{
	  Arad[i]=(Av[i]-1./3.*Erf*GG[0][i])/(4./3.*Erf*gammamax);
	}
      
      //is normalized now
      ldouble Afac;
      c=0.; b=0.;
      for(i=1;i<4;i++)
	{
	  a+=Arad[i]*Arad[i]*gg[i][i];
	  b+=2.*Arad[i]*gg[0][i]*gammamax;
	}
      c=gg[0][0]*gammamax*gammamax;
      delta=b*b-4.*a*c;
      Afac= (-b+sqrtl(delta))/2./a;

      urfcon[0]=gammamax;
      urfcon[1]=Afac*Arad[1];
      urfcon[2]=Afac*Arad[2];
      urfcon[3]=Afac*Arad[3];

      //converting to relative four velocity
      conv_vels(urfcon,urfcon,VEL4,VELR,gg,GG);
    }
  else if(gamma2<1.)
    {
      //low cap
      *corrected=1;

      urfcon[0]=1.;
      urfcon[1]=urfcon[2]=urfcon[3]=0.;

      //radiative energy density in the radiation rest frame
      Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);

      //converting to relative four velocity
      conv_vels(urfcon,urfcon,VEL4,VELR,gg,GG);
    }
  else
    {
      //regular calculation
      urfcon[0]=sqrtl(gamma2);
    
      //radiative energy density in the radiation rest frame
      Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);
      
      //relative velocity
      ldouble alpha=sqrtl(-1./GG[0][0]);
      ldouble gamma=urfcon[0]*alpha;
      for(i=1;i<4;i++)
	{	  
	  urfcon[i]=(3.*Av[i]-Erf*GG[0][i])/(3.*Av[0]-Erf*GG[0][0])/alpha+GG[0][i]/alpha;
	  urfcon[i]*=gamma;
	}
      urfcon[0]=0.;
    }
  
  //new primitives
  pp[6]=Erf;
  pp[7]=urfcon[1];
  pp[8]=urfcon[2];
  pp[9]=urfcon[3];

  return 0;
}

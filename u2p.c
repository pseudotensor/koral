
//KORAL - u2p.c
//conserved to primitives conversion

#include "ko.h"


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates primitives in given cell basing on global array u[]
int
calc_primitives(int ix,int iy,int iz)
{
  int verbose=0;
  int iv,u2pret,u2pretav;
  ldouble uu[NV],uuav[NV],pp[NV],ppav[NV];
  ldouble gg[4][5],GG[4][5], tlo[4][4],tup[4][4];

  pick_g(ix,iy,iz,gg);
  pick_G(ix,iy,iz,GG);

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  //converting to primitives
  int corrected;
  u2pret=u2p(uu,pp,gg,GG,&corrected);

  //update conserved to follow corrections on primitives
  if(corrected!=0)
    {
      if(verbose) {printf("correcting conserved at %d %d %d\n",ix,iy,iz);getchar();}
      p2u(pp,uu,gg,GG);
      for(iv=0;iv<NV;iv++)
	{
	  set_u(u,iv,ix,iy,iz,uu[iv]);
	}
    }

  //sets the flag to mark if hot conversion did not succeed - the entropy will not be updated
  set_cflag(0,ix,iy,iz,u2pret); 
  
  for(iv=0;iv<NV;iv++)    
    set_u(p,iv,ix,iy,iz,pp[iv]);	      

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//high-level u2p solver
int
u2p(ldouble *uu, ldouble *pp, ldouble gg[][5],ldouble GG[][5],int *corrected)
{
  *corrected=0;
  int verbose=1;
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
  //  u2pret=u2p_hot_gsl(uu,pp,gg,GG);  //temporary 1D numerical solver - not perfect
#else
  u2pret=u2p_hot(uu,pp,gg,GG);  //TODO: check cancelation!
#endif
  //************************************

  if(u2pret<0) 
    {
      if(verbose>0)
	printf("u2p_hot err  >>> %d <<< %e %e\n",u2pret,pp[0],pp[1]);
      
      //************************************
      //entropy solver - conserving entropy
      ret=-1;
      u2pret=u2p_entropy(uu,pp,gg,GG);
      //************************************

      if(verbose>0)
	  printf("u2p_entr     >>> %d <<< %e > %e\n",u2pret,u0,pp[1]);

      if(u2pret<0)
	{
	  if(verbose>0)
	    printf("u2p_entr err > %e %e\n",pp[0],pp[1]);
	
	  //************************************
	  //cold RHD - assuming u=SMALL
	  ret=-2;
	  u2pret=u2p_cold(uu,pp,gg,GG);
	  //************************************

	  if(u2pret<0)
	    {
	      if(verbose>0)
		printf("u2p_cold err > %e %e\n",pp[0],pp[1]);
	
	      //************************************
	      //leaving unchanged primitives - should not happen
	      ret=-3;
	      for(u2pret=0;u2pret<NV;u2pret++)
		pp[u2pret]=ppbak[u2pret];	  
	      //************************************
	    }

	}
      
     }

  if(ret<0.)
    hdcorr=1;

  //************************************
  //************************************
  //checking on hd floors
  ret=check_floors_hd(pp,VELPRIM,gg,GG);

  if(ret<0.)
    hdcorr=1;
  //************************************
  //************************************



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
//checks if hydro primitives make sense
//TODO: incorporate structure of state here?
int
check_floors_hd(ldouble *pp, int whichvel,ldouble gg[][5], ldouble GG[][5])
{
  int i,j,k,correct;
  // velocities
  ldouble a,b,c,delta;
  ldouble u1[4],u2[4];
  if(whichvel==VEL4)
    {
      //assumes u^t unknown
      u1[0]=0.;
      u1[1]=pp[2];
      u1[2]=pp[3];
      u1[3]=pp[4];
      a=gg[0][0];
      b=0.;
      c=1.;
      for(i=1;i<4;i++)
	{
	  b+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      c+=u1[i]*u1[j]*gg[i][j];
	    }
	}
      delta=b*b-4.*a*c;
      correct=0;
      if(delta<0.) 
	correct=1;
      else
	{
	  u1[0]=(-b-sqrt(delta))/2./a;
	  if(u1[0]<1.) u1[0]=(-b+sqrt(delta))/2./a;
	  if(u1[0]<1.) 
	    correct=1;
	}
	  
      if(correct==0) return 0; //everything is fine
    }
  else if(whichvel==VEL3)
    {
      //assumes u^t unknown
      u1[0]=0.;
      u1[1]=pp[2];
      u1[2]=pp[3];
      u1[3]=pp[4];
      a=b=0.;
      for(i=1;i<4;i++)
	{
	  a+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      b+=u1[i]*u1[j]*gg[i][j];
	    }
	}

      correct=0;
      if(-1./(gg[0][0]+a+b)<0.)
	correct=1;
     	  
      if(correct==0) return 0; //everything is fine
    }
  else if(whichvel==VELR)
    {
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      ldouble alpha2=-1./GG[0][0];
      
      correct=0;
      if(gamma2/alpha2<0.)
	correct=1;
     	  
      if(correct==0) return 0; //everything is fine

    }
  else
    {
      my_err("VEL not implemented in check_floors_hd()\n");
    }

  //correcting and imposing gammamax keeping the direction given by spatial components
  ldouble Afac;
  ldouble gammamax=100.;

  if(whichvel==VELR)
    {
      //convert spatial u1 to lab-frame
      for(i=1;i<4;i++)
	u1[i]=u1[i]+gammamax*GG[0][i]/GG[0][0];
    }

  //normalizing to u^t=gammamax
  c=0.; b=0.;
  for(i=1;i<4;i++)
    {
      a+=u1[i]*u1[i]*gg[i][i];
      b+=2.*u1[i]*gg[0][i]*gammamax;
    }
  c=gg[0][0]*gammamax*gammamax+1.;
  delta=b*b-4.*a*c;
      
  Afac= (-b+sqrt(delta))/2./a;

  u2[0]=gammamax;
  u2[1]=Afac*u1[1];
  u2[2]=Afac*u1[2];
  u2[3]=Afac*u1[3];

   //converting to whichvel
  conv_vels(u2,u2,VEL4,whichvel,gg,GG);

  /*
  print_4vector(u1);
  print_4vector(u2);
  printf("which: %d det: %e\n",whichvel,dot(u2,u1));
  //  getchar();
  */

  //back to primitives
  pp[2]=u2[1];
  pp[3]=u2[2];
  pp[4]=u2[3];

  //to let the others know to correct conserved
  return -1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//generalized conserved to primitives solver
//'hot grhd' - pure hydro, numerical in 2d
//following Noble+06
ldouble
f_u2p_hot(ldouble W, ldouble* cons)
{
  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];

  return -(Qn+W)*(GAMMA/GAMMAM1)+W*(1.-Qt2/W/W)-D*sqrt(1.-Qt2/W/W);   
}

int
u2p_hot(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5])
{
  int verbose=0;
  int i,j,k;
  ldouble rho,u,p,w,W,gamma,alpha,D;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
  
  if(verbose) {printf("********************\n");print_Nvector(uu,NV);}
  if(verbose) {print_Nvector(pp,NV);}

  //alpha
  alpha=sqrt(-1./GG[0][0]);

  //D
  D=uu[0]*alpha;

  //Q_mu
  Qcov[0]=(uu[1]-uu[0])*alpha;
  Qcov[1]=uu[2]*alpha;
  Qcov[2]=uu[3]*alpha;
  Qcov[3]=uu[4]*alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

  //n_mu = (-alpha, 0, 0, 0)
  ncov[0]=-alpha;
  ncov[1]=ncov[2]=ncov[3]=0.;
  
  //n^mu
  indices_12(ncov,ncon,GG);

  //Q^mu n_mu = -alpha*Q^t
  Qn=Qcon[0] * ncov[0];

  //j^mu_nu=delta^mu_nu +n^mu n_nu
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      jmunu[i][j] = delta(i,j) + ncon[i]*ncov[j];

  //Qtilda^nu = j^nu_mu Q^mu
  for(i=0;i<4;i++)
    {
      Qtcon[i]=0.;
      for(j=0;j<4;j++)
	Qtcon[i]+=jmunu[i][j]*Qcon[j];
    }

  //Qtilda_nu
  indices_21(Qtcon,Qtcov,gg);

  //Qt2=Qtilda^mu Qtilda_mu
  Qt2=dot(Qtcon,Qtcov);

  //initial guess for W = w gamma**2 based on primitives
  rho=pp[0];
  u=pp[1];
  utcon[0]=0.;
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels(utcon,utcon,VELPRIM,VELR,gg,GG);
  ldouble qsq=0.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=utcon[i]*utcon[j]*gg[i][j];
  ldouble gamma2=1.+qsq;
  W=(rho+GAMMA*u)*gamma2;
  if(verbose) printf("initial W:%e\n",W);
 
  //test if does not provide reasonable gamma2
  if(W*W<Qt2)
    {
      W=2.*sqrt(Qt2);
      if(verbose) printf("corrected W:%e\n",W);
    }

  //1d Newton solver
  ldouble CONV=1.e-6;
  ldouble EPS=1.e-8;
  ldouble Wprev=W;
  ldouble f0,f1,dfdW;
  ldouble cons[3]={Qn,Qt2,D};
  if(verbose) printf("in:%e %e %e\n",Qn,Qt2,D);

  int iter=0;
  do
    {
      Wprev=W;
      iter++;
      f0=f_u2p_hot(W,cons);

      f1=f_u2p_hot(W*(1.+EPS),cons);
      dfdW=(f1-f0)/(EPS*W);

      if(verbose) printf("%d %e %e %e %e\n",iter,W,f0,f1,dfdW);

      if(dfdW==0.) {W*=1.1; continue;}
      W-=f0/dfdW;
    }
  while(fabs((W-Wprev)/Wprev)>CONV && iter<50);

  if(iter>=50)
    {
      if(verbose) printf("iter exceeded in u2p_hot\n");
      return -1;
    }
  
  if(isnan(W) || isinf(W)) {if(verbose)printf("nan/inf W: %e\n",W); return -1;}
  if(verbose) {printf("the end: %e\n",W); }

  //W found, let's calculate v2 and the rest
  ldouble v2=Qt2/W/W;
  gamma2=1./(1.-v2);
  gamma=sqrt(gamma2);
  rho=D/gamma;
  u=1./GAMMA*(W/gamma2-rho);
  utcon[0]=0.;
  utcon[1]=gamma/W*Qtcon[1];
  utcon[2]=gamma/W*Qtcon[2];
  utcon[3]=gamma/W*Qtcon[3];

  if(rho<0. || u<0. || gamma2<0. ||isnan(W) || isinf(W)) 
    {
      if(verbose) printf("neg u rho in u2p_hot\n");
      return -1;
    }

  //converting to VELPRIM
  conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  
  //returning new primitives
  pp[0]=rho;
  pp[1]=u;
  pp[2]=utcon[1];
  pp[3]=utcon[2];
  pp[4]=utcon[3];

  //entropy
  ldouble Sut=uu[5];
  ldouble ut=uu[0]/pp[0]; //rhout/rho
  pp[5]=Sut/ut;

  if(verbose) {print_Nvector(pp,NV);}

  return 0;

}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//auxiliary solver based on the entropy conservation
int
u2p_entropy(ldouble *uuu, ldouble *p, ldouble g[][5], ldouble G[][5])
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

  conv_velsinprims(p,VELPRIM,VEL3,g,G);

  if(verbose) print_Nvector(p,NV);

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

    if(rhop1<RHOFLOOR) rhop1=rho/2.;

    if(verbose) printf("%d %e %e %e\n",iter,rhop1,rho,fval);

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
	printf(" entr  iter %d> %e [%e] %e >%e< %e\n",iter,rho,rhop1,diffrho/rho,fval,err);
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

  if(uu<0. || rho<0. || isnan(rho))
    {
      printf("u2p_entr didn't work: %e %e\n",uu,rho); 
      print_Nvector(uuu,NV);
      getchar();
      return -1;
    }

  p[0]=rho;
  p[1]=uu;
  p[2]=vr;
  p[3]=vth;
  p[4]=vph;
  p[5]=S;

  //************************************
  //************************************
  //checking on hd floors
  check_floors_hd(p,VEL3,g,G);
  //************************************
  //************************************

  conv_velsinprims(p,VEL3,VELPRIM,g,G);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//auxiliary solver assuming u=0
int
u2p_cold(ldouble *uuu, ldouble *p, ldouble g[][5], ldouble G[][5])
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

  //************************************
  //************************************
  //checking on hd floors
  check_floors_hd(p,VEL3,g,G);
  //************************************
  //************************************

  conv_velsinprims(p,VEL3,VELPRIM,g,G);

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
  gamma2=  (-b-sqrt(delta))/2./a;
  //if unphysical try the other root
  if(gamma2<0.) gamma2=  (-b+sqrt(delta))/2./a; 

  //cap on u^t
  ldouble gammamax=1000.;

  //gamma in relative velocity definition
  ldouble gammarel2=gamma2/(-GG[0][0]);

   if(gammarel2<0. || gammarel2>gammamax*gammamax || delta<0.) 
    {
      //      printf("top cap\n");
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
      c=gg[0][0]*gammamax*gammamax+1.;
      delta=b*b-4.*a*c;
      Afac= (-b+sqrt(delta))/2./a;

      urfcon[0]=gammamax;
      urfcon[1]=Afac*Arad[1];
      urfcon[2]=Afac*Arad[2];
      urfcon[3]=Afac*Arad[3];

      //converting to relative four velocity
      conv_vels(urfcon,urfcon,VEL4,VELR,gg,GG);
    }
   else if(gammarel2<(-1./GG[0][0]))
    {
      printf("low cap\n");
      //low cap
      *corrected=1;

      //zeros for relative velocity
      urfcon[0]=urfcon[1]=urfcon[2]=urfcon[3]=0.;

      //calculating time component of lab 4-vel
      ldouble gammarel=1.0;
      ldouble urflab[4];
      ldouble alpha = sqrt(-1./GG[0][0]);
      urflab[0]=gammarel/alpha;

      //radiative energy density in the radiation rest frame
      Erf=3.*Av[0]/(4.*urflab[0]*urflab[0]+GG[0][0]);
    }
  else
    {
      //regular calculation
      urfcon[0]=sqrt(gamma2);
    
      //radiative energy density in the radiation rest frame
      Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);

      //relative velocity
      ldouble alpha=sqrt(-1./GG[0][0]);
      ldouble gamma=urfcon[0]*alpha;
      for(i=1;i<4;i++)
	{	  
	  urfcon[i]=(3.*Av[i]-Erf*GG[0][i])/(3.*Av[0]-Erf*GG[0][0])/alpha-GG[0][i]/GG[0][0]/alpha;
	  urfcon[i]*=gamma;
	}
      urfcon[0]=0.;
    }

   conv_vels(urfcon,urfcon,VELR,VELPRIMRAD,gg,GG);
  
   //new primitives
   pp[6]=Erf;
   pp[7]=urfcon[1];
   pp[8]=urfcon[2];
   pp[9]=urfcon[3];

   return 0;
}

//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//**********************************************************************
//**********************************************************************
int
u2p_rad_labfluxes(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5],int *corrected)
{
  int verbose=0;
  ldouble Rij[4][4];

  //R^0mu
  ldouble A[4]={uu[6],uu[7],uu[8],uu[9]};
  //indices up
  indices_12(A,A,GG);

  //covariant formulation
  
  //g_munu R^0mu R^0nu
  ldouble gRR=gg[0][0]*A[0]*A[0]+gg[0][1]*A[0]*A[1]+gg[0][2]*A[0]*A[2]+gg[0][3]*A[0]*A[3]+
    gg[1][0]*A[1]*A[0]+gg[1][1]*A[1]*A[1]+gg[1][2]*A[1]*A[2]+gg[1][3]*A[1]*A[3]+
    gg[2][0]*A[2]*A[0]+gg[2][1]*A[2]*A[1]+gg[2][2]*A[2]*A[2]+gg[2][3]*A[2]*A[3]+
    gg[3][0]*A[3]*A[0]+gg[3][1]*A[3]*A[1]+gg[3][2]*A[3]*A[2]+gg[3][3]*A[3]*A[3];
 
  //the quadratic equation for u^t of the radiation rest frame (urf[0])
  ldouble a,b,c;
  a=16.*gRR;
  b=8.*(gRR*GG[0][0]+A[0]*A[0]);
  c=gRR*GG[0][0]*GG[0][0]-A[0]*A[0]*GG[0][0];
  ldouble delta=b*b-4.*a*c;
  ldouble urf[4],Erf;
  urf[0]=sqrt((-b-sqrt(delta))/2./a);
  if(isnan(urf[0])) 
    {
      my_err("top cap should be imposed\n");
      urf[0]=1.;
    }

  //radiative energy density in the radiation rest frame
  Erf=3.*A[0]/(4.*urf[0]*urf[0]+GG[0][0]);

  //four-velocity of the rest frame
  urf[1]=3./(4.*Erf*urf[0])*(A[1]-1./3.*Erf*GG[0][1]);
  urf[2]=3./(4.*Erf*urf[0])*(A[2]-1./3.*Erf*GG[0][2]);
  urf[3]=3./(4.*Erf*urf[0])*(A[3]-1./3.*Erf*GG[0][3]);

  //converting to three velocity
  //  conv_vels(urf,urf,VEL4,VELR,gg,GG);

  //reading primitives
  pp[6]=Erf;
  pp[7]=urf[1];
  pp[8]=urf[2];
  pp[9]=urf[3];

  return 0;
}


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
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);
 
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  //converting to primitives
  int corrected;
  u2p(uu,pp,gg,GG,tup,tlo,&corrected);
  /*
  if(corrected!=0) 
    { 
      printf("happened at %d %d %d\n",ix,iy,iz);
      getchar();
    }
  */

  //update conserved to follow corrections on primitives
  if(corrected!=0)
    {
      if(verbose) {printf("correcting conserved at %d %d %d\n",ix,iy,iz);}//getchar();}
      p2u(pp,uu,gg,GG);
      for(iv=0;iv<NV;iv++)
	{
	  set_u(u,iv,ix,iy,iz,uu[iv]);
	}
    }

  //sets the flag to mark if hot conversion did not succeed - the entropy will not be updated
   if(corrected!=0)
     set_cflag(ENTROPYFLAG,ix,iy,iz,-1); 
   else 
     set_cflag(ENTROPYFLAG,ix,iy,iz,0); 
  
  for(iv=0;iv<NV;iv++)    
    set_u(p,iv,ix,iy,iz,pp[iv]);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//high-level u2p solver
int
u2p(ldouble *uu, ldouble *pp, ldouble gg[][5],ldouble GG[][5],ldouble tup[][4],ldouble tlo[][4],int *corrected)
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
  u2pret=u2p_hot(uu,pp,gg,GG);  

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
	  //leaving unchanged primitives - should not happen
	  ret=-3;
	  for(u2pret=0;u2pret<NV;u2pret++)
	    pp[u2pret]=ppbak[u2pret];	  
	  //************************************

	  return -3;

	
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
#ifdef EDDINGTON_APR_WRONG
  u2p_rad_onff(uu,pp,gg,GG,tup,tlo,&radcorr);
#else
  u2p_rad(uu,pp,gg,GG,&radcorr);
#endif
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
  ldouble gammamaxhd=GAMMAMAXHD;

  if(whichvel==VELR)
    {
      //convert spatial u1 to lab-frame
      for(i=1;i<4;i++)
	u1[i]=u1[i]+gammamaxhd*GG[0][i]/GG[0][0];
    }

  //normalizing to u^t=gammamaxhd
  c=0.; b=0.;
  for(i=1;i<4;i++)
    {
      a+=u1[i]*u1[i]*gg[i][i];
      b+=2.*u1[i]*gg[0][i]*gammamaxhd;
    }
  c=gg[0][0]*gammamaxhd*gammamaxhd+1.;
  delta=b*b-4.*a*c;
      
  Afac= (-b+sqrt(delta))/2./a;

  u2[0]=gammamaxhd;
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

  //a bit more clear

  ldouble v2 = Qt2/W/W;
  ldouble gamma2 = 1./(1.-v2);
  ldouble gamma = sqrt(gamma2);
  ldouble w = W/gamma2;
  ldouble rho = D/gamma;
  ldouble u = (w - rho) / GAMMA;
  ldouble p = (GAMMA-1)*u;

  return Qn + W - p;
 
  
}

int
u2p_hot(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5])
{
  int verbose=0;
  int i,j,k;
  ldouble rho,u,p,w,W,alpha,D;
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

  //Q_mu n^mu = Q^mu n_mu = -alpha*Q^t
  Qn=Qcon[0] * ncov[0];
  //Qn = dot(Qcov,ncon);

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
  ldouble gamma=sqrt(gamma2);

  //W
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
  ldouble EPS=1.e-6;
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
  ldouble ut2 = Qt2/(W*W - Qt2);
  gamma2=1. + ut2;
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

  if(verbose) {print_Nvector(pp,NV); getchar();}

  return 0;

}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//auxiliary solver based on the entropy conservation
//works for general metric in four velocities
int
u2p_entropy(ldouble *uuu, ldouble *p, ldouble g[][5], ldouble G[][5])
{
  int verbose=0;

  ldouble gtt=g[0][0];
  ldouble gtr=g[0][1];
  ldouble gtth=g[0][2];
  ldouble gtph=g[0][3];

  ldouble grt=g[1][0];
  ldouble grr=g[1][1];
  ldouble grth=g[1][2];
  ldouble grph=g[1][3];

  ldouble gtht=g[2][0];
  ldouble gthr=g[2][1];
  ldouble gthth=g[2][2];
  ldouble gthph=g[2][3];

  ldouble gpht=g[3][0];
  ldouble gphr=g[3][1];
  ldouble gphth=g[3][2];
  ldouble gphph=g[3][3];

  ldouble rhout=uuu[0];
  ldouble Tttt=uuu[1]; //this one unused
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
  ldouble err,W,fval,dfval,fval1,fval2,diffrho,ut;
  ldouble dudrho,dWdrho,dvphdrho,dvrdrho,dvthdrho;

  ldouble conv=1.e-6;
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
 
    W=ut2*(uu*GAMMA + rho);  

    vr=-((gphth*gthr*Ttph - gphr*gthth*Ttph - gphth*gthph*Ttr + gphph*gthth*Ttr + gphr*gthph*Ttth - gphph*gthr*Ttth - gphth*gthr*gtph*W + gphr*gthth*gtph*W + gphth*gthph*gtr*W - gphph*gthth*gtr*W - gphr*gthph*gtth*W + gphph*gthr*gtth*W)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*W));
    vth=-((gphth*grr*Ttph - gphr*grth*Ttph - gphth*grph*Ttr + gphph*grth*Ttr + gphr*grph*Ttth - gphph*grr*Ttth - gphth*grr*gtph*W + gphr*grth*gtph*W + gphth*grph*gtr*W - gphph*grth*gtr*W - gphr*grph*gtth*W + gphph*grr*gtth*W)/((-(gphth*grr*gthph) + gphr*grth*gthph + gphth*grph*gthr - gphph*grth*gthr - gphr*grph*gthth + gphph*grr*gthth)*W));
    vph=(grth*gthr*Ttph - grr*gthth*Ttph - grth*gthph*Ttr + grph*gthth*Ttr + grr*gthph*Ttth - grph*gthr*Ttth - grth*gthr*gtph*W + grr*gthth*gtph*W + grth*gthph*gtr*W - grph*gthth*gtr*W - grr*gthph*gtth*W + grph*gthr*gtth*W)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*W);

    fval=1./ut2 + gtt + grr*vr*vr + gthth*vth*vth + gphph*vph*vph + 2.*gtph*vph + 2.*gtr*vr + 2.*gtth*vth +
      2.*grph*vr*vph + 2.*grth*vr*vth + 2.*gphth*vph*vth;

    dudrho=(ldouble)GAMMA/GAMMAM1*exp(Sut/rhout)*pow(rho,1./GAMMAM1)*pow(exp(Sut/rhout)*pow(rho,GAMMA/GAMMAM1),GAMMA-2.);
    dWdrho=(rhout*rhout*(rho*(GAMMA*dudrho-1.)-2.*GAMMA*uu))/rho/rho/rho;

    dvrdrho=((gphth*gthr*Ttph - gphr*gthth*Ttph - gphth*gthph*Ttr + gphph*gthth*Ttr + gphr*gthph*Ttth - gphph*gthr*Ttth - gphth*gthr*gtph*W + gphr*gthth*gtph*W + gphth*gthph*gtr*W - gphph*gthth*gtr*W - gphr*gthph*gtth*W + gphph*gthr*gtth*W)*dWdrho)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*Power(W,2)) - (-(gphth*gthr*gtph*dWdrho) + gphr*gthth*gtph*dWdrho + gphth*gthph*gtr*dWdrho - gphph*gthth*gtr*dWdrho - gphr*gthph*gtth*dWdrho + gphph*gthr*gtth*dWdrho)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*W);
    dvthdrho=((gphth*grr*Ttph - gphr*grth*Ttph - gphth*grph*Ttr + gphph*grth*Ttr + gphr*grph*Ttth - gphph*grr*Ttth - gphth*grr*gtph*W + gphr*grth*gtph*W + gphth*grph*gtr*W - gphph*grth*gtr*W - gphr*grph*gtth*W + gphph*grr*gtth*W)*dWdrho)/((-(gphth*grr*gthph) + gphr*grth*gthph + gphth*grph*gthr - gphph*grth*gthr - gphr*grph*gthth + gphph*grr*gthth)*Power(W,2)) - (-(gphth*grr*gtph*dWdrho) + gphr*grth*gtph*dWdrho + gphth*grph*gtr*dWdrho - gphph*grth*gtr*dWdrho - gphr*grph*gtth*dWdrho + gphph*grr*gtth*dWdrho)/((-(gphth*grr*gthph) + gphr*grth*gthph + gphth*grph*gthr - gphph*grth*gthr - gphr*grph*gthth + gphph*grr*gthth)*W);
    dvphdrho=-(((grth*gthr*Ttph - grr*gthth*Ttph - grth*gthph*Ttr + grph*gthth*Ttr + grr*gthph*Ttth - grph*gthr*Ttth - grth*gthr*gtph*W + grr*gthth*gtph*W + grth*gthph*gtr*W - grph*gthth*gtr*W - grr*gthph*gtth*W + grph*gthr*gtth*W)*dWdrho)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*Power(W,2))) + (-(grth*gthr*gtph*dWdrho) + grr*gthth*gtph*dWdrho + grth*gthph*gtr*dWdrho - grph*gthth*gtr*dWdrho - grr*gthph*gtth*dWdrho + grph*gthr*gtth*dWdrho)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*W);

    fval=1./ut2 + gtt + grr*vr*vr + gthth*vth*vth + gphph*vph*vph + 2.*gtph*vph + 2.*gtr*vr + 2.*gtth*vth +
      2.*grph*vr*vph + 2.*grth*vr*vth + 2.*gphth*vph*vth;
    dfval=2.*rho/rhout/rhout + 
      2.*grr*vr*dvrdrho + 2.*gthth*vth*dvthdrho + 2.*gphph*vph*dvphdrho + 
      2.*gtph*dvphdrho + 2.*gtr*dvrdrho + 2.*gtth*dvthdrho + 
      2.*grph*dvrdrho*vph + 2.*grph*vr*dvphdrho + 
      2.*grth*dvrdrho*vth + 2.*grth*vr*dvthdrho + 
      2.*gphth*dvphdrho*vth + 2.*gphth*vph*dvthdrho;

    //absolute error
    absfval=fabs(fval);

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

    err =fabs(diffrho/rho);
 
    if(iter>itmax && err>conv) 
      {
	printf("iter exceeded in u2p_entr \n");
	printf(" entr  iter %d> %e [%e] %e >%e< %e\n",iter,rho,rhop1,diffrho/rho,fval,err);
	return -1;
      }

    rhom1=rho;

  } while(err>conv);
  
  if(verbose) {printf("success %d %e %e %e\n",iter,rhop1,rho,fval);getchar();}

  rho=rhop1;
  ut=rhout/rho;
  ut2=ut*ut;
  S=Sut/ut;
  uu=(ldouble)pow((pow(rho,1./GAMMAM1+1.)*exp(Sut/rhout)),GAMMAM1)/GAMMAM1;
  W=ut2*(uu*GAMMA + rho);
  
  vr=-((gphth*gthr*Ttph - gphr*gthth*Ttph - gphth*gthph*Ttr + gphph*gthth*Ttr + gphr*gthph*Ttth - gphph*gthr*Ttth - gphth*gthr*gtph*W + gphr*gthth*gtph*W + gphth*gthph*gtr*W - gphph*gthth*gtr*W - gphr*gthph*gtth*W + gphph*gthr*gtth*W)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*W));
  vth=-((gphth*grr*Ttph - gphr*grth*Ttph - gphth*grph*Ttr + gphph*grth*Ttr + gphr*grph*Ttth - gphph*grr*Ttth - gphth*grr*gtph*W + gphr*grth*gtph*W + gphth*grph*gtr*W - gphph*grth*gtr*W - gphr*grph*gtth*W + gphph*grr*gtth*W)/((-(gphth*grr*gthph) + gphr*grth*gthph + gphth*grph*gthr - gphph*grth*gthr - gphr*grph*gthth + gphph*grr*gthth)*W));
  vph=(grth*gthr*Ttph - grr*gthth*Ttph - grth*gthph*Ttr + grph*gthth*Ttr + grr*gthph*Ttth - grph*gthr*Ttth - grth*gthr*gtph*W + grr*gthth*gtph*W + grth*gthph*gtr*W - grph*gthth*gtr*W - grr*gthph*gtth*W + grph*gthr*gtth*W)/((gphth*grr*gthph - gphr*grth*gthph - gphth*grph*gthr + gphph*grth*gthr + gphr*grph*gthth - gphph*grr*gthth)*W);

  if(uu<0. || rho<0. || isnan(rho))
    {
      printf("u2p_entr didn't work: %e %e\n",uu,rho); 
      //print_Nvector(uuu,NV);
      //getchar();
      return -1;
    }

  p[0]=rho;
  p[1]=uu;
  p[2]=vr;
  p[3]=vth;
  p[4]=vph;
  p[5]=S;
  
  conv_velsinprims(p,VEL3,VELPRIM,g,G);
 
  //************************************
  //************************************
  //checking on hd floors
  check_floors_hd(p,VELPRIM,g,G);
  //************************************
  //************************************

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//auxiliary solver assuming u=0
int
u2p_cold(ldouble *uuu, ldouble *p, ldouble g[][5], ldouble G[][5])
{
  printf("Should not be in u2p_cold() yet - to be generalized.\n");

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
//takes conserved R^t_mu in uu
//**********************************************************************
//**********************************************************************
int
u2p_rad_urf(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5], int *corrected)
{
  //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
  *corrected=0;

  int verbose=0,i,j;
  ldouble Rij[4][4];
  ldouble urfcon[4],urfcov[4],Erf;
  ldouble alpha = sqrt(-1./GG[0][0]);
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
  a=16.*gRR;
  b=8.*(gRR*GG[0][0]+Av[0]*Av[0]);
  c=gRR*GG[0][0]*GG[0][0]-Av[0]*Av[0]*GG[0][0];
  delta=b*b-4.*a*c;
  gamma2=  (-b-sqrt(delta))/2./a;
  //if unphysical try the other root
  if(gamma2<0.) gamma2=  (-b+sqrt(delta))/2./a; 

  /*
  if(isnan(gamma2) || gamma2<0. || 1)
    {
      print_4vector(Av);
      printf("nan gamma2! %e %e\n", (-b-sqrt(delta))/2./a, (-b+sqrt(delta))/2./a); //getchar();
    }
  */

  //cap on u^t
  ldouble gammamax=GAMMAMAXRAD;

  //gamma in relative velocity definition
  ldouble gammarel2=gamma2*alpha*alpha;

  
  /*
  if(delta<0.)
    {
      // can't assume this conditions means large gamma, because if not, then leads to crazy boost of energy.
      Erf=ERADFLOOR;
      //gammarel2=1.0;
      //Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM
      //zeros for relative velocity
      urfcon[0]=urfcon[1]=urfcon[2]=urfcon[3]=0.;

      if(verbose) {printf("topcapbad: gammarel2=%g gamma2=%g delta=%g\n",gammarel2,gamma2,delta);}
    }
  else 
  */

  if(gammarel2>gammamax*gammamax || gamma2<0.) 
    {      
      //top cap
      *corrected=1;
      if(verbose) printf("topcap\n");
			 
      //urfcon[0]=gammamax;
      ldouble gammarel;
      
      gammarel=gammamax;

      gammarel2=gammarel*gammarel;

      //proper direction for the radiation rest frame, will be normalized later      
      //Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);
      Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

      if(Erf<ERADFLOOR)
	{
	  
	  Erf=ERADFLOOR;
	  urfcon[0]=0.;
	  urfcon[1]=0.;
	  urfcon[2]=0.;
	  urfcon[3]=0.;
	  if(verbose) {printf("topcapErfneg: gammarel2=%g gamma2=%g\n",gammarel2,gamma2);}
	}					
      else
	{
	  if(1){ //VELR
	    // lab-frame radiation relative 4-velocity
	    ldouble Aradrel[4];
	    for(i=1;i<4;i++)
	      Aradrel[i] = alpha * (Av[i] + 1./3.*Erf*GG[0][i]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);

	    // skipping as Aradrel[] constructed assuming gammamax
	    // compute \gammarel using this
	    ldouble qsq=0.;
	    for(i=1;i<4;i++)
	      for(j=1;j<4;j++)
		qsq+=Aradrel[i]*Aradrel[j]*gg[i][j];
	    ldouble gammatemp=sqrt(1.+qsq);

	    //MYFUN(gamma_calc_fromuconrel(Aradrel,ptrgeom,&gammatemp,&qsqtemp),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);

	    // now rescale Aradrel[i] so will give desired \gammamax
	    for(i=1;i<4;i++)
	      {
		Aradrel[i] *= (gammamax/gammatemp);
	      }
	    

	    for(i=1;i<4;i++)
	      {
		urfcon[i]=Aradrel[i];
	      }
	    if(verbose) {printf("topcapgamma Erf=%g gammaorg=%g gammatemp=%g gammanow=%g\n",Erf,gamma2/alpha/alpha,gammatemp,gammamax);}
	  }
	  else if(0){ //going through VEL4
	    ldouble Arad[4];
	    for(i=1;i<4;i++)
	      {
		Arad[i]=(Av[i]-1./3.*Erf*GG[0][i])/(4./3.*Erf*gammamax);
	      }
       
	    //is normalized now
	    ldouble Afac;
	    a=0.; c=0.; b=0.;
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
	    conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,gg,GG);

	    if(verbose) {printf("topcapolek: Erf=%g Afac=%g Arad123=%g %g %g\n",Erf,Afac,Arad[1],Arad[2],Arad[3]);}
	  }// end Olek method
	}// end else if Erf>0
    }
  //  else if(gammarel2<(-1./GG[0][0]))
  else if(gammarel2<1. || delta<0.)
    {
      //low cap
     
      //this usually happens when gamma=0.99999999 so enforcing 1. makes sense
       *corrected=0;
       if(verbose) {printf("midcapalt: gamma: %.20e\n",sqrt(gammarel2));}
       if(0)
	{
	  //zeros for relative velocity
	  urfcon[0]=urfcon[1]=urfcon[2]=urfcon[3]=0.;
	  
	  //calculating time component of lab 4-vel
	  ldouble gammarel=1.0;
	  ldouble urflab[4];
	  
	  urflab[0]=gammarel/alpha;
	  
	  //radiative energy density in the radiation rest frame
	  Erf=3.*Av[0]/(4.*urflab[0]*urflab[0]+GG[0][0]);
	}
      else //in VELR
	{
	  gammarel2=1.0;
	  Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM
	  //zeros for relative velocity
	  urfcon[0]=urfcon[1]=urfcon[2]=urfcon[3]=0.;
	}

      if(verbose) {printf("midcapalt: Erf=%g\n",Erf);}
	 
      if(Erf<ERADFLOOR)
	{ 
	  // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
	  if(verbose) {printf("midcapaltnegErf: Erf=%g\n",Erf);}
	  Erf=ERADFLOOR;
	}	
    }
  else
    {
      //regular calculation
      //urfcon[0]=sqrt(gamma2);
      //radiative energy density in the radiation rest frame
      //Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);
      Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

      if(Erf<ERADFLOOR)
	{
	  Erf=ERADFLOOR;
	  urfcon[0]=0.;
	  urfcon[1]=0.;
	  urfcon[2]=0.;
	  urfcon[3]=0.;
	  if(verbose) {printf("nocapbad: gammarel2=%g\n",gammarel2);}
	}					
 
      if(1) //VELR
	{
	  ldouble gammarel=sqrt(gammarel2);
	  for(i=1;i<4;i++)
	    {	  
	      urfcon[i] = alpha * (Av[i] + 1./3.*Erf*GG[0][i]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
	    }
	}
      else
	{
	  //relative velocity
	  ldouble gamma=urfcon[0]*alpha;
	  for(i=1;i<4;i++)
	    {	  
	      urfcon[i]=(3.*Av[i]-Erf*GG[0][i])/(3.*Av[0]-Erf*GG[0][0])/alpha-GG[0][i]/GG[0][0]/alpha;
	      urfcon[i]*=gamma;
	    }
	  urfcon[0]=0.;
	}
    }

   conv_vels(urfcon,urfcon,VELR,VELPRIMRAD,gg,GG);
  
   //new primitives
   pp[6]=Erf;
   pp[7]=urfcon[1];
   pp[8]=urfcon[2];
   pp[9]=urfcon[3];

   return 0;
}

int
u2p_rad(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5], int *corrected)
{
  //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
  *corrected=0;

#ifdef LABRADFLUXES
  //primitives = R^t_mu
  pp[6]=uu[6];
  pp[7]=uu[7];
  pp[8]=uu[8];
  pp[9]=uu[9];
  return 0;
#endif
  
  u2p_rad_urf(uu,pp,gg,GG,corrected);
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//works in ortonormal fluid frame
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
int f_u2prad_num(ldouble *uu,ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble tlo[][4],ldouble *f)
{
  ldouble Rij[4][4];
  ldouble ppp[NV];

  calc_Rij_ff(pp,Rij);  
  trans22_on2cc(Rij,Rij,tlo);  
  boost22_ff2lab(Rij,Rij,pp,gg,GG); 
  indices_2221(Rij,Rij,gg);  

  f[0]=-Rij[0][0]+uu[6];
  f[1]=-Rij[0][1]+uu[7];
  f[2]=-Rij[0][2]+uu[8];
  f[3]=-Rij[0][3]+uu[9];

  return 0;
} 

int
print_state_u2prad_num (int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .3e % .3e % .3e % .3e "
	  "f(x) = % .3e % .3e % .3e % .3e\n",
	  iter,
	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
}

int
u2p_rad_onff(ldouble *uu, ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble tup[][4], ldouble tlo[][4], int *corrected)
{
  ldouble pp0[NV],pporg[NV];
  ldouble J[4][4],iJ[4][4];
  ldouble x[4],f1[4],f2[4],f3[4];
  int i,j,k,iter=0;

  ldouble EPS = 1.e-6;
  ldouble CONV = U2PRADPREC;

  int verbose=0;

  for(i=6;i<NV;i++)
    {
      pporg[i]=pp[i];
    }

  //converting radiative primitives to fluid frame ortonormal
  prad_lab2ff(pp,pp,gg,GG,tup);

  if(verbose!=0)   print_Nvector(uu,NV);
  do
    {
      iter++;
      for(i=6;i<NV;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_u2prad_num(uu,pp,gg,GG,tlo,f1);
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      pp[j+6]=pp[j+6]+EPS*pp[6];
	    
	      f_u2prad_num(uu,pp,gg,GG,tlo,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(EPS*pp[6]);

	      pp[j+6]=pp0[j+6];
	    }
	}

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+6];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}
      if(verbose>0)    print_state_u2prad_num (iter,x,f1); 

      for(i=0;i<4;i++)
	{
	  pp[i+6]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+6]-pp0[i+6]);
	  f3[i]=fabs(f3[i]/pp0[6]);
	}

      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	break;

      if(iter>50)
	{
	  printf("iter exceeded in u2prad_num()\n");
	  
	  for(i=6;i<NV;i++)
	    {
	      pp[i]=pporg[i];
	    }
	  
	  *corrected=1;
	  return -1;

	  break;
	}
     
    }
  while(1);
  
  if(pp[6]<EFLOOR) 
    {
      printf("enegative u2prad()\n");
      pp[6]=EFLOOR;
      *corrected=1;
    }
  
  //converting to lab primitives
  prad_ff2lab(pp,pp,gg,GG,tlo);
  
  if(verbose!=0)   {print_Nvector(pp,NV);}
  if(verbose>0)   {printf("----\n");}

  *corrected=0;
  return 0;

}

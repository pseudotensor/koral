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
  ldouble tlo[4][4],tup[4][4];
  ldouble (*gg)[5],(*GG)[5];

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  //converting to primitives
  int corrected[2], fixups[2];
  u2p(uu,pp,&geom,corrected,fixups);

  //update conserved to follow corrections on primitives
  //should I skip this when going to fixup - if averagin primitives this will have no effect?
  if(corrected[0]!=0 || corrected[1]!=0)
    {
      if(verbose) {printf("correcting conserved at %d %d %d\n",ix,iy,iz);}//getchar();}
      p2u(pp,uu,&geom);
      for(iv=0;iv<NV;iv++)
	{
	  set_u(u,iv,ix,iy,iz,uu[iv]);
	}
    }

  //sets the flags for fixups of unsuccessful cells
  if(fixups[0]>0)
     set_cflag(HDFIXUPFLAG,ix,iy,iz,1); 
  else
    set_cflag(HDFIXUPFLAG,ix,iy,iz,0); 

  if(fixups[1]>0)
     set_cflag(RADFIXUPFLAG,ix,iy,iz,1); 
  else
     set_cflag(RADFIXUPFLAG,ix,iy,iz,0); 

  //sets the flag to mark if hot conversion did not succeed - the entropy will not be updated
   if(corrected[0]!=0)
     {
       //getchar();
       set_cflag(ENTROPYFLAG,ix,iy,iz,-1); 
     }
   else 
     set_cflag(ENTROPYFLAG,ix,iy,iz,0); 
  
  for(iv=0;iv<NV;iv++)    
    set_u(p,iv,ix,iy,iz,pp[iv]);

  //  ldouble u1=get_u(u,5,ix,iy,iz);
  //  printf("%e %e\n",get_u(u,5,ix,iy,iz),get_u(p,5,ix,iy,iz));
  update_entropy(ix,iy,iz,get_cflag(ENTROPYFLAG,ix,iy,iz)); //updates u[5]=Sut(rho,u) when u2p_hot() succeded or calculates p[5]=S from previous Sut if the other case



  // ldouble u2=get_u(u,5,ix,iy,iz);
  // printf("%e %e %d\n",get_u(u,5,ix,iy,iz),get_u(p,5,ix,iy,iz),get_cflag(ENTROPYFLAG,ix,iy,iz));
  //     if(fabs((u1-u2)/u2)>1.e-5 && get_cflag(ENTROPYFLAG,ix,iy,iz)==0)   getchar();

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//high-level u2p solver
int
u2p(ldouble *uu, ldouble *pp,void *ggg,int corrected[2],int fixups[2])
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  corrected[0]=corrected[1]=0;
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
  u2pret=u2p_hot(uu,pp,gg,GG);  

  //************************************

  if(u2pret<0) 
    {
      if(verbose>1)
	printf("u2p_hot err at %d,%d,%d >>> %d <<< %e %e\n",geom->ix,geom->iy,geom->iz,u2pret,pp[0],pp[1]);
      ret=-1;
    }
     
#ifdef ALLOWENTROPYU2P
  if(u2pret<0)
    {
      //************************************
      //entropy solver - conserving entropy
      ret=-1;
      //      print_Nvector(uu,NV);
      //      print_Nvector(ppbak,NV);
      //      print_Nvector(pp,NV);
      u2pret=u2p_entropy(uu,pp,gg,GG);
      //************************************

      if(verbose>1)
	{
	  printf("u2p_entr     >>> %d <<< %e > %e\n",u2pret,u0,pp[1]);
	}
    
      if(u2pret<0)
	{
	  if(verbose>0)
	    {
	      printf("u2p_entr err > %e %e > %e %e > %d %d %d\n",uu[0],uu[1],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
	      getchar();
	    }
	}
    }
#endif

#ifdef ALLOWCOLDU2P

  if(u2pret<0.)
    {
      //***********************************
       //cold RHD - assuming u=SMALL
      ret=-2;
      u2pret=u2p_cold(uu,pp,gg,GG);
      //************************************

      if(u2perr<0)
	  if(verbose>0)
	    {
	      printf("u2p_cold err > %e %e > %e %e > %d %d %d\n",uu[0],uu[1],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
	    }
    }

#endif

  if(u2pret<0)
    {
      //************************************
      //leaving unchanged primitives - should not happen
      ret=-3;
      for(u2pret=0;u2pret<NV;u2pret++)
	pp[u2pret]=ppbak[u2pret];	  
      //************************************
    }
  
  if(ret<0.)
    hdcorr=1;

  //************************************
  //************************************
  if(ret<-1) //do not fix up for entropy solver
    fixups[0]=1;
  else
    fixups[0]=0;
  //************************************
  //************************************


  //************************************
  //************************************
  //checking on hd floors
  
  int floorret;
  floorret=check_floors_hd(pp,VELPRIM,ggg);

  if(floorret<0.)
    {
      hdcorr=1;
      fixups[0]=1;
    }

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
  u2p_rad(uu,pp,geom,&radcorr);
#endif
  
  //************************************
  //************************************
  if(radcorr>0)
    fixups[1]=1;
  else
    fixups[1]=0;
  //************************************
  //************************************

  if(hdcorr>0) corrected[0]=1;
  if(radcorr>0) corrected[1]=1;

  return ret;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//checks if hydro primitives make sense
//TODO: incorporate structure of state here?
int
check_floors_hd(ldouble *pp, int whichvel,void *ggg)
{
  int verbose=1;
  int ret=0;

  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;


  //absolute rho
  if(pp[0]<RHOFLOOR) {pp[0]=RHOFLOOR; ret=-1; if(verbose) printf("hd_floors CASE 1\n");}

  //uint/rho ratios
  ldouble UURHORATIO=1.e-7;
  if(pp[1]<UURHORATIO*pp[0]) {pp[1]=UURHORATIO*pp[0];ret=-1;}//if(verbose) printf("hd_floors CASE 2\n");}
  if(pp[1]>1./UURHORATIO*pp[0]) {pp[1]=1./UURHORATIO*pp[0];ret=-1;}//if(verbose) printf("hd_floors CASE 3\n");}

#ifdef RADIATION
  //EE/rho ratios
  ldouble EERHORATIO=1.e-7;
  ldouble pp2[NV];
  prad_lab2ff(pp, pp2, ggg);
  
  if(pp[6]<EERHORATIO*pp[0]) {pp[6]=EERHORATIO*pp[0];ret=-1;if(verbose) printf("hd_floors CASE R2\n");}
  if(pp[6]>1./EERHORATIO*pp[0]) {pp[6]=1./EERHORATIO*pp[0];ret=-1;if(verbose) printf("hd_floors CASE R3\n");}

  prad_ff2lab(pp2, pp, ggg);
#endif
 

  //velocities
  int i,j,k,correct;
  //velocities
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
    }
  else if(whichvel==VELR)
    {
      u1[0]=0.;
      u1[1]=pp[2];
      u1[2]=pp[3];
      u1[3]=pp[4];

      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      ldouble alpha2=-1./GG[0][0];
      
      correct=0;
      if(gamma2/alpha2<0. || gamma2>GAMMAMAXHD*GAMMAMAXHD)
	{
	  if(verbose) 
	    {
	      printf("hd_floors CASE 4 %e %e %e\n",gamma2/alpha2,gamma2,GAMMAMAXHD*GAMMAMAXHD);
	      print_4vector(u1);//getchar();
	    }
	  correct=1;
	}
    }
  else
    {
      my_err("VEL not implemented in check_floors_hd()\n");
    }

  if(correct==0) return ret;

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
  ret=-1;
  return ret;
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
  
  if(verbose>1) {printf("********************\n");print_Nvector(uu,NV);}
  if(verbose>1) {print_Nvector(pp,NV);}

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

  if(verbose>1) printf("initial W:%e\n",W);
 
  //test if does not provide reasonable gamma2
  if(W*W<Qt2)
    {
      W=2.*sqrt(Qt2);
      if(verbose>0) printf("corrected W:%e\n",W);
    }

  //1d Newton solver
  ldouble CONV=1.e-4;
  ldouble EPS=1.e-6;
  ldouble Wprev=W;
  ldouble f0,f1,dfdW;
  ldouble cons[3]={Qn,Qt2,D};
  if(verbose>1) printf("in:%e %e %e\n",Qn,Qt2,D);

  int iter=0;
  do
    {
      Wprev=W;
      iter++;
      f0=f_u2p_hot(W,cons);

      f1=f_u2p_hot(W*(1.+EPS),cons);
      dfdW=(f1-f0)/(EPS*W);

      if(verbose>1) printf("%d %e %e %e %e\n",iter,W,f0,f1,dfdW);

      if(dfdW==0.) {W*=1.1; continue;}

      ldouble Wnew=W-f0/dfdW;

      //test if does produce nan and damp solution if so
      if(Wnew*Wnew<Qt2)
	{
	  int idump=0;
	  ldouble dumpfac=1.;
	  do
	    {
	      idump++;
	      dumpfac/=2.;
	      Wnew=W-dumpfac*f0/dfdW;
	    }
	  while(Wnew*Wnew<Qt2 && idump<100);
	  
	  if(idump>=100) return -1;

	  if(verbose>0) printf("damped successfuly\n");
	}

      W=Wnew; 
    }
  while(fabs((W-Wprev)/Wprev)>CONV && iter<50);

  if(iter>=50)
    {
      if(verbose>0) printf("iter exceeded in u2p_hot\n");
      return -1;
    }
  
  if(isnan(W) || isinf(W)) {if(verbose>0)printf("nan/inf W: %e\n",W); getchar();return -1;}
  if(verbose>1) {printf("the end: %e\n",W); }

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
      if(verbose>0) printf("neg u rho in u2p_hot %e %e %e %e\n",rho,u,gamma2,W);//getchar();
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
  //  ldouble Sut=uu[5];
  //  ldouble ut=uu[0]/pp[0]; //rhout/rho
  //  pp[5]=Sut/ut;

  pp[5]=calc_Sfromu(rho,u);

  if(verbose>1) {print_Nvector(pp,NV); getchar();}

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
 
  ldouble rho=p[0]; //initial guess
  ldouble rho0=p[0]; //initial guess
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
	if(verbose || 1) printf("iter exceeded in u2p_entr \n");
	if(verbose || 1) printf(" entr  iter %d> %e [%e] %e >%e< %e\n",iter,rho,rhop1,diffrho/rho,fval,err);
	if(verbose || 1) print_Nvector(p,NV);
	if(verbose || 1) print_Nvector(uuu,NV);

	//	getchar();
	return -1;
      }

    rhom1=rho;

  } while(err>conv);
  
  if(verbose) {printf("success %d %e %e %e\n",iter,rhop1,rho,fval);
    
  }

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
      if(verbose || 1) printf("iter didn't work in u2p_entr \n");
      if(verbose || 1) printf(" entr  iter %d> %e [%e] %e >%e< %e\n",iter,rho,rhop1,diffrho/rho,fval,err);
      if(verbose || 1) print_Nvector(p,NV);
      if(verbose || 1) print_Nvector(uuu,NV);
      //      getchar();
      return -1;
    }

  p[0]=rho;
  p[1]=uu;
  p[2]=vr;
  p[3]=vth;
  p[4]=vph;
  p[5]=S;
  
  if(conv_velsinprims(p,VEL3,VELPRIM,g,G)!=0) 
    {
      print_Nvector(p,NV);
      printf("conv vels in _entropy failed %e ut\n",ut);
      //      getchar();
      return -1;
    }
 
  //************************************
  //************************************
  //checking on hd floors
  //check_floors_hd(p,VELPRIM,g,G);
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
  int verbose=1;

  my_err("cold");

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
  ldouble Tttt=uuu[1]; 
  ldouble Ttr=uuu[2];
  ldouble Ttth=uuu[3];
  ldouble Ttph=uuu[4];
  ldouble Sut=uuu[5];
  ldouble Ttt=Tttt+rhout;

  //  conv_velsinprims(p,VELPRIM,VEL3,g,G);


  ldouble rho,uu,vr,vth,vph,S,fff,ut;
  
  fff=(Ttt - Ttph*gtph/gphph)/(gtt - gtph*gtph/gphph);
  ut=fff/rhout;
					       
  rho=rhout/ut;
  uu=UUFLOOR;
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
  //check_floors_hd(p,VEL3,g,G);
  //************************************
  //***********************************

  /*
 if(verbose) print_Nvector(uuu,NV);
 if(verbose) print_Nvector(p,NV);
 if(verbose) getchar();
  */
  conv_velsinprims(p,VEL3,VELPRIM,g,G);

  return 0;
}


// get's gamma assuming fixed E rather than using original R^t_t that we assume is flawed near floor regions.  We want to preserve R^t_i (i.e conserved momentum)
static int get_m1closure_gammarel2_cold(int verbose,
					void *ggg, 
					FTYPE *Avcon, 
					FTYPE *Avcov, 
					FTYPE *gammarel2return, 
					FTYPE *deltareturn, 
					FTYPE *numeratorreturn, 
					FTYPE *divisorreturn, 
					FTYPE *Erfreturn, 
					FTYPE *urfconrel)
{
   struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  FTYPE gamma2,gammarel2,delta;
  FTYPE Erf;
  FTYPE alpha=geom->alpha;
  int jj;

  static FTYPE gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=GG[0][0];
  gn12=GG[0][1];
  gn13=GG[0][2];
  gn14=GG[0][3];
  gn22=GG[1][1];
  gn23=GG[1][2];
  gn24=GG[1][3];
  gn33=GG[2][2];
  gn34=GG[2][3];
  gn44=GG[3][3]; 

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];

  // choose gamma
  if(gammarel2return==NULL){
    FTYPE gammamax=GAMMAMAXRAD;
    gammarel2=gammamax*gammamax;
  }
  else gammarel2=*gammarel2return; // feed in desired gammarel2

  FTYPE utsq=gammarel2/(alpha*alpha);

  FTYPE Avcovorig[NDIM],Avconorig[NDIM];
  DLOOPA(jj) Avcovorig[jj]=Avcov[jj];
  DLOOPA(jj) Avconorig[jj]=Avcon[jj];

  // get new Avcov[0]=R^t_t

  // NOTEMARK: Note that Sqrt() is only ever negative when gammarel2<0, so never has to be concern.
  Avcov[0]=(0.25*(-4.*(gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*utsq*(gn11 + utsq) + 
                        Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                              1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                       gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2))))/(gn11*utsq*(gn11 + utsq));

  Erf=(0.75*Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                           1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                    gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2)))/(utsq*(gn11 + utsq)*(gn11 + 4.*utsq));

  if(verbose) printf("NOR SOL: Avcov0new=%g Avcov0old=%g Erf=%g :: %g %g %g\n",Avcov[0],Avcovorig[0],Erf,Rtx,Rty,Rtz);

  //modify Avcon
  indices_12(Avcov,Avcon,GG);
  if(verbose) DLOOPA(jj) printf("jj=%d Avconorig=%g Avconnew=%g\n",jj,Avconorig[jj],Avcon[jj]);

 
  delta=0; // not yet

  *gammarel2return=gammarel2;
  *deltareturn=delta;

  // get good relative velocity
  FTYPE gammarel=sqrt(gammarel2);

  // get relative 4-velocity
  if(Erf>0.0) SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*GG[0][jj]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
  else SLOOPA(jj) urfconrel[jj] = 0.0;

  *Erfreturn=Erf; // pass back new Erf to pointer

  return(0);
}




// get's gamma^2 for lab-frame gamma  using Rd and gcon
static int get_m1closure_gammarel2(int verbose,void *ggg, ldouble *Avcon, ldouble *Avcov, ldouble *gammarel2return,ldouble *deltareturn, ldouble *numeratorreturn, ldouble *divisorreturn)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble gamma2,gammarel2,delta,numerator,divisor;
  ldouble gamma2a,gamma2b;

  // mathematica solution that avoids catastrophic cancellation when Rtt very small (otherwise above gives gamma2=1/2 oddly when gamma2=1) -- otherwise same as above
  // well, then had problems for R~1E-14 for some reason when near BH.  Couldn't quickly figure out, so use no replacement of gv11.
  // see u2p_inversion.nb
  static ldouble gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=GG[0][0];
  gn12=GG[0][1];
  gn13=GG[0][2];
  gn14=GG[0][3];
  gn22=GG[1][1];
  gn23=GG[1][2];
  gn24=GG[1][3];
  gn33=GG[2][2];
  gn34=GG[2][3];
  gn44=GG[3][3];

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];

  gamma2a=(-0.25*(2.*Power(gn11,2)*Power(Rdtt,2) + (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*
        (gn12*Rdtx + gn13*Rdty + gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
            gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
               8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))) + 
       gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 
          2.*gn34*Rdty*Rdtz + gn44*Power(Rdtz,2) + Rdtt*
           (4.*gn13*Rdty + 4.*gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
               gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                  8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))))))/
   (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
    2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));


  if( gamma2a<GAMMASMALLLIMIT || !isfinite(gamma2a) ){
    gamma2b=(0.25*(-2.*Power(gn11,2)*Power(Rdtt,2) - 1.*gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 
                                                              Rdty*(4.*gn13*Rdtt + 2.*gn23*Rdtx + gn33*Rdty) + 2.*(2.*gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2)) + 
                   gn11*Rdtt*Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                  gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                        8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))) + 
                   (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*(-1.*gn12*Rdtx - 1.*gn13*Rdty - 1.*gn14*Rdtz + 
                                                        Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                                             gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                                                   8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))))))/
      (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
       2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));
    gamma2=gamma2b;
  }
  else{
    // choose
    gamma2=gamma2a;
  }

  ////////////////////////
  //
  //cap on u^t
  //
  ///////////////////////
  ldouble alpha=geom->alpha;


  // get relative 4-velocity, that is always >=1 even in GR
  gammarel2 = gamma2*alpha*alpha;

  // check for machine error away from 1.0 that happens sometimes
  if(gammarel2>GAMMASMALLLIMIT && gammarel2<1.0){
    if(verbose) printf("Hit machine error of gammarel2=%27.20g fixed to be 1.0\n",gammarel2);
    gammarel2=1.0;
  }

  *gammarel2return=gammarel2; 
  *deltareturn=delta=0;
  *numeratorreturn=numerator=0;
  *divisorreturn=divisor=0;
  return(0);
}


// get Erf
static int get_m1closure_Erf(void *ggg, ldouble *Avcon, ldouble gammarel2, ldouble *Erfreturn)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble alpha=geom->alpha;

  ////////////
  //
  // get initial attempt for Erf
  // If delta<0, then gammarel2=nan and Erf<RADLIMIT check below will fail as good.
  //
  ////////////
  *Erfreturn = 3.*Avcon[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

  //  printf("ee %d %d %d %e\n",geom->ix,geom->iy,geom->iz,*Erfreturn);

  return(0);
}

// get contravariant relative 4-velocity in lab frame
static int get_m1closure_urfconrel(int verbose, 
				   void *ggg, 
				   ldouble *pp, 
				   ldouble *Avcon, 
				   ldouble *Avcov, 
				   ldouble gammarel2, 
				   ldouble delta, 
				   ldouble numerator, 
				   ldouble divisor, 
				   ldouble *Erfreturn,
				   ldouble *urfconrel, 
				   int *corflag)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble alpha=geom->alpha;

  ldouble Erf=*Erfreturn; // get initial Erf
  ldouble gammamax=GAMMAMAXRAD;
  int ii,jj,kk;

  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  // NOTE: gammarel2 just below 1.0 already fixed to be =1.0
  int nonfailure=gammarel2>=1.0 && Erf>ERADLIMIT && gammarel2<=gammamax*gammamax/GAMMASMALLLIMIT/GAMMASMALLLIMIT;
  // falilure1 : gammarel2 normal, but already Erf<ERADLIMIT (note for M1 that gammarel2>=1/4 for any reasonable chance for correct non-zero Erf
  int failure1=Avcon[0]<0.0 || (gammarel2>0.0 && gammarel2<=0.25 && delta>=0.0 && divisor!=0.0) || numerator==0.0 || gammarel2>=1.0 && delta>=0.0 && divisor!=0.0 && Erf<ERADLIMIT;
  // gamma probably around 1
  int failure2=gammarel2<1.0 && gammarel2>0.0 && delta>=0.0;
  // i.e. all else, so not really used below.
  int failure3=gammarel2>gammamax*gammamax && Erf>=ERADLIMIT || gammarel2<0.0 || delta<0.  || divisor==0.0 && numerator==0.0 || divisor==0.0 && numerator!=0.0;

  // any failure
  int failure=!nonfailure || !isfinite(gammarel2) || !isfinite(Erf);

  if(failure && (failure1==0 && failure2==0 && failure3==0)){
    printf("Undetected failure, now considered\n");
  }


  if(nonfailure){
    // get good relative velocity
    ldouble gammarel=sqrt(gammarel2);
 
    for(ii=1;ii<4;ii++)
      {	  
	urfconrel[ii] = alpha * (Avcon[ii] + 1./3.*Erf*GG[0][ii]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
      }

    *Erfreturn=Erf; // pass back new Erf to pointer
    *corflag=0;
    return(0);
  }
  else{
    if(verbose) 
      {
	printf("failure: %d %d %d %d %d\n",nonfailure,failure1,failure2,failure3,failure);
	printf("%e %e %e\n",Erf,gammarel2,Avcov[0]);
	print_4vector(Avcon);
	print_4vector(Avcov);
	getchar();
      }
    // get \gammarel=1 case
    ldouble gammarel2slow=pow(1.0+10.0*NUMEPSILON,2.0);
    ldouble Avconslow[4],Avcovslow[4],Erfslow,urfconrelslow[4];
    for(jj=0;jj<4;jj++)
      {
	Avconslow[jj]=Avcon[jj];
	Avcovslow[jj]=Avcov[jj];
      }
    Erfslow=Erf;
    get_m1closure_gammarel2_cold(verbose,ggg,Avconslow,Avcovslow,&gammarel2slow,&delta,&numerator,&divisor,&Erfslow,urfconrelslow);

    // get \gammarel=gammamax case
    ldouble gammarel2fast=gammamax*gammamax;
    ldouble Avconfast[4],Avcovfast[4],Erffast,urfconrelfast[4];
    for(jj=0;jj<4;jj++)
      {
	Avconfast[jj]=Avcon[jj];
	Avcovfast[jj]=Avcov[jj];
      }
    Erffast=Erf;
    get_m1closure_gammarel2_cold(verbose,ggg,Avconfast,Avcovfast,&gammarel2fast,&delta,&numerator,&divisor,&Erffast,urfconrelfast);

    int usingfast=1;
    // choose by which Avcov[0] is closest to original
    if( fabs(Avcovslow[0]-Avcov[0])>fabs(Avcovfast[0]-Avcov[0])){
      usingfast=1;
      for(jj=0;jj<4;jj++)
	{
	  Avcon[jj]=Avconfast[jj];
	  Avcov[jj]=Avcovfast[jj];
	  urfconrel[jj]=urfconrelfast[jj];
	}
      gammarel2=gammarel2fast;
      Erf=Erffast;
    }
    else{
      usingfast=0;
      for(jj=0;jj<4;jj++)
	{
	  Avcon[jj]=Avconslow[jj];
	  Avcov[jj]=Avcovslow[jj];
	  urfconrel[jj]=urfconrelslow[jj];
	}
      gammarel2=gammarel2slow;
      Erf=Erfslow;
    }    

    // report
    *corflag=1;
    if(verbose) printf("CASEGEN: gammarel>gammamax (cold, usingfast=%d): gammarel2=%g Erf=%g : i=%d j=%d k=%d\n",usingfast,gammarel2,Erf,geom->ix,geom->iy,geom->iz);
  }


  *Erfreturn=Erf; // pass back new Erf to pointer

  if(!isfinite(Erf) || !isfinite(gammarel2) || !isfinite(urfconrel[0])|| !isfinite(urfconrel[1])|| !isfinite(urfconrel[2])|| !isfinite(urfconrel[3]) ){
      printf("JONNAN: ijk=%d %d %d :  %g %g : %g %g %g : %d %d %d %d : %g %g %g %g\n",geom->ix,geom->iy,geom->iz,Erf,gammarel2,urfconrel[1],urfconrel[2],urfconrel[3],failure1,failure2,failure3,failure,Avcon[0],Avcon[1],Avcon[2],Avcon[3]);
  }

  return(0);
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
u2p_rad_urf(ldouble *uu, ldouble *pp,void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  int irf,verbose=0;

  //multi-rad-fluids
  for(irf=0;irf<NRF;irf++)
    {
      //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
      *corrected=0;

      ldouble Rij[4][4];
      ldouble urfcon[4],urfcov[4],Erf;
      //conserved - R^t_mu
      ldouble Avcov[4]={uu[EE(irf)],uu[FX(irf)],uu[FY(irf)],uu[FZ(irf)]};
      ldouble Avcon[4];
      //indices up - R^tmu
      indices_12(Avcov,Avcon,GG);

      ldouble gammarel2,delta,numerator,divisor;

      // get \gamma^2 for relative 4-velocity
      get_m1closure_gammarel2(verbose,ggg,Avcon,Avcov,&gammarel2,&delta,&numerator,&divisor);

      // get E in radiation frame
      get_m1closure_Erf(ggg,Avcon,gammarel2,&Erf);

      // get relative 4-velocity
      get_m1closure_urfconrel(verbose,ggg,pp,Avcon,Avcov,gammarel2,delta,numerator,divisor,&Erf,urfcon,corrected);

      //new primitives
      pp[EE(irf)]=Erf;
      pp[FX(irf)]=urfcon[1];
      pp[FY(irf)]=urfcon[2];
      pp[FZ(irf)]=urfcon[3];
    }

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
u2p_rad_urf_old(ldouble *uu, ldouble *pp,void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  int irf;

  //multi-rad-fluids
  for(irf=0;irf<NRF;irf++)
    {

      //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
      *corrected=0;

      int verbose=0,debug=0;
      int i,j;
      ldouble Rij[4][4];
      ldouble urfcon[4],urfcov[4],Erf;
      ldouble alpha = sqrt(-1./GG[0][0]);
      //conserved - R^t_mu
      ldouble Av[4]={uu[EE(irf)],uu[FX(irf)],uu[FY(irf)],uu[FZ(irf)]};
      //indices up - R^tmu
      indices_12(Av,Av,GG);

#if(0)
      //ortonormal classic formulation
      ldouble Avlen2=Av[1]*Av[1]+Av[2]*Av[2]+Av[3]*Av[3];
      if(Avlen2>Av[0]*Av[0])
	{
	  ldouble AE=sqrt(Avlen2/Av[0]/Av[0]);
	  Av[1]/=1.00001*AE;
	  Av[2]/=1.00001*AE;
	  Av[3]/=1.00001*AE;
	  //printf("F/E: %e %e %e f:%e\n",Av[1]/Av[0],Av[2]/Av[0],Av[3]/Av[0],sqrt(Av[1]*Av[1]+Av[2]*Av[2]+Av[3]*Av[3])/Av[0]);getchar();
	}

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

      //cap on u^t
      ldouble gammamax=GAMMAMAXRAD;

      //gamma in relative velocity definition
      ldouble gammarel2=gamma2*alpha*alpha;
      if((gammarel2-1.0)>-1E-13 && gammarel2<1.0) gammarel2=1.0;
      
      debug=0;
      Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

      if(Erf<EEFLOOR)
	{
	  Erf=EEFLOOR;      printf("imposing EEFLOOR 1\n");

	  urfcon[0]=0.;
	  urfcon[1]=0.;
	  urfcon[2]=0.;
	  urfcon[3]=0.;
	  if(verbose) {printf("nocapbad: gammarel2=%g\n",gammarel2);}
	}				
      else
	{ 
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


#else


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

      //cap on u^t
      ldouble gammamax=GAMMAMAXRAD;

      //gamma in relative velocity definition
      ldouble gammarel2=gamma2*alpha*alpha;

      //to avoid the machine precision error
      if((gammarel2-1.0)>-1E-13 && gammarel2<1.0) gammarel2=1.0;

      if((gammarel2>1.0*gammamax*gammamax || gammarel2<.9 || delta<0.) && verbose) 
	{
	  debug=1;

	  printf("ixy: %d %d irf: %d\n",geom->ix,geom->iy,irf);
	  print_Nvector(uu,NV);

	  print_4vector(Av);
	  printf("F/E: %e %e %e f:%e\n",Av[1]*sqrt(gg[1][1])/Av[0],Av[2]*sqrt(gg[2][2])/Av[0],sqrt(gg[3][3])*Av[3]/Av[0],sqrt(Av[1]*Av[1]+Av[2]*Av[2]+Av[3]*Av[3])/Av[0]);
	  printf("gamma2: %.20e %e\n", (-b-sqrt(delta))/2./a, (-b+sqrt(delta))/2./a); 
	  printf("delta: %e gRR: %e Erf: %e\n",delta,gRR,3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0));
      
	  //	  if(geom->ix==0&&geom->iy==0&&irf==1)getchar();//getchar();
	  
	}

      if(gammarel2>1.0*gammamax*gammamax || gammarel2<0. || delta<0.) 
	{      
      
	  //top cap
	  *corrected=1;
	  if(verbose) printf("topcap %e %e %e \n",gammarel2,gammamax*gammamax,delta);
			 
	  //urfcon[0]=gammamax;
	  ldouble gammarel;
      
	  gammarel=gammamax;

	  gammarel2=gammarel*gammarel;

	  //proper direction for the radiation rest frame, will be normalized later      
	  //Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);
	  Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

	  if(Erf<EEFLOOR)
	    {
	  
	      Erf=EEFLOOR;      printf("imposing EEFLOOR 2\n");

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
		if(verbose) {printf("topcapgamma Erf=%g gammaorg=%.7e gammatemp=%g gammanow=%g\n",Erf,gamma2/alpha/alpha,gammatemp,gammamax);}
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
      else if(gammarel2<1.)
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
	 
	  if(Erf<EEFLOOR)
	    { 
	      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
	      if(verbose) {printf("midcapaltnegErf: Erf=%g\n",Erf);}
	      Erf=EEFLOOR;      printf("imposing EEFLOOR 3\n");

	    }	
	}
      else
	{
	  //regular calculation
	  //urfcon[0]=sqrt(gamma2);
	  //radiative energy density in the radiation rest frame
	  //Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);
	  Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

	  if(Erf<EEFLOOR)
	    {
	      Erf=EEFLOOR;      printf("imposing EEFLOOR 4\n");

	      urfcon[0]=0.;
	      urfcon[1]=0.;
	      urfcon[2]=0.;
	      urfcon[3]=0.;
	      if(verbose) {printf("nocapbad: gammarel2=%g\n",gammarel2);}
	    }	
	  else
	    {
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
	}
    
#endif

      conv_vels(urfcon,urfcon,VELR,VELPRIMRAD,gg,GG);
  
      /*
   if(debug==1 && verbose)
     {
#ifndef LABRADFLUXES
       ldouble Rij[4][4];
       calc_Rij(pp,ggg,Rij);
       trans22_cc2on(Rij,Rij,tup);
       ldouble Avv[4]={Rij[0][0],Rij[0][1], Rij[0][2],Rij[0][3]};
       print_4vector(Avv);
       printf("F/E: %e %e %e f:%e\n",Avv[1]/Avv[0],Avv[2]/Avv[0],Avv[3]/Avv[0],sqrt(Avv[1]*Avv[1]+Avv[2]*Avv[2]+Avv[3]*Avv[3])/Avv[0]);
#endif

       getchar();
       }
      */
      
      //new primitives
      pp[EE(irf)]=Erf;
      pp[FX(irf)]=urfcon[1];
      pp[FY(irf)]=urfcon[2];
      pp[FZ(irf)]=urfcon[3];
   
      if(debug==1 && verbose)
	{
	  ldouble Rij[4][4];
	  calc_Rij(pp,ggg,Rij);
	  trans22_cc2on(Rij,Rij,tup);
	  ldouble Avv[4]={Rij[0][0],Rij[0][1], Rij[0][2],Rij[0][3]};
	  print_4vector(Avv);
	  printf("F/E: %e %e %e f:%e\n",Avv[1]/Avv[0],Avv[2]/Avv[0],Avv[3]/Avv[0],sqrt(Avv[1]*Avv[1]+Avv[2]*Avv[2]+Avv[3]*Avv[3])/Avv[0]);
 
	  getchar();
	}
    }

  return 0;
}

int
u2p_rad(ldouble *uu, ldouble *pp, void *ggg, int *corrected)
{
  //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
  *corrected=0;

#ifdef EDDINGTON_APR
  int irf,ii;
  ldouble pp10[10],uu10[10];
  for(ii=0;ii<NVHD;ii++)
    {
      pp10[ii]=pp[ii];
      uu10[ii]=uu[ii];
    }
  
  for(irf=0;irf<NRF;irf++)
    {
      pp10[EE(0)]=pp[EE(irf)];
      pp10[FX(0)]=pp[FX(irf)];
      pp10[FY(0)]=pp[FY(irf)];
      pp10[FZ(0)]=pp[FZ(irf)];
      uu10[EE(0)]=uu[EE(irf)];
      uu10[FX(0)]=uu[FX(irf)];
      uu10[FY(0)]=uu[FY(irf)];
      uu10[FZ(0)]=uu[FZ(irf)];

      if(uu10[EE(0)]<EEFLOsOR)
	u2p_rad_onff(uu10,pp10, ggg, corrected);
      else
	{
	  p2u_rad(pp10,uu10,ggg);
	}

      pp[EE(irf)]=pp10[EE(0)];
      pp[FX(irf)]=pp10[FX(0)];
      pp[FY(irf)]=pp10[FY(0)];
      pp[FZ(irf)]=pp10[FZ(0)];
    }
  return 0;
#endif

#ifdef LABRADFLUXES
  //primitives = R^t_mu
  pp[6]=uu[6];
  pp[7]=uu[7];
  pp[8]=uu[8];
  pp[9]=uu[9];
  return 0;
#endif
  
  u2p_rad_urf(uu,pp,ggg,corrected);
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//works in ortonormal fluid frame
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
int f_u2prad_num(ldouble *uu,ldouble *pp, void* ggg,ldouble *f)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  ldouble Rij[4][4];

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
u2p_rad_onff(ldouble *uu, ldouble *pp, void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble pp0[NV],pporg[NV];
  ldouble J[4][4],iJ[4][4];
  ldouble x[4],f1[4],f2[4],f3[4];
  int i,j,k,iter=0;

  ldouble EPS = 1.e-6;
  ldouble CONV = 1.e-6;

  int verbose=0;

  for(i=6;i<NV;i++)
    {
      pporg[i]=pp[i];
    }

  //converting radiative primitives to fluid frame ortonormal
  //ad_lab2ff(pp,pp,geom);

  if(verbose!=0)   print_Nvector(uu,NV);
  do
    {
      iter++;
      for(i=6;i<NV;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_u2prad_num(uu,pp,geom,f1);
 
      //calculating approximate Jacobian
      for(j=0;j<4;j++)
	{
	  pp[j+6]=pp[j+6]+EPS*pp[6];
	    
	  if(verbose>0)    print_Nvector(pp,NV);
 	  f_u2prad_num(uu,pp,geom,f2);
	  if(verbose>0)    print_state_u2prad_num (iter,x,f2); 
     
	  for(i=0;i<4;i++)
	    {
	      J[i][j]=(f2[i] - f1[i])/(EPS*pp[6]);
	    }

	  pp[j+6]=pp0[j+6];
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

      if(iter>20)
	{
	  printf("iter exceeded in u2prad_num() %d %d %d\n",geom->ix,geom->iy,geom->iz);getchar();

	  pp[6]=pporg[6];
	  pp[7]=pp[8]=pp[9]=0.;
	  
	  *corrected=1;
	  return -1;

	  break;
	}
     
    }
  while(1);
  
  if(pp[6]<EEFLOOR) 
    {
      printf("enegative u2prad() %d %d %d\n",geom->ix,geom->iy,geom->iz); getchar();
      pp[6]=EEFLOOR;
      pp[7]=pp[8]=pp[9]=0.;
      *corrected=1;
      return -1;
    }
  
  //converting to lab primitives
  //no - in EDD_APR I use fluid frame fluxes as primitives
  //prad_ff2lab(pp,pp,geom);
  
  if(verbose!=0)   {print_Nvector(pp,NV);}
  if(verbose>0)   {printf("----\n");}//getchar();}

  *corrected=0;
  return 0;

}

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
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;
  gdet=geom.gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int corrected[2]={0,0}, fixups[2]={0,0};


  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  //converting to primitives
  u2p(uu,pp,&geom,corrected,fixups);

  //imposing floors

  //************************************
  //************************************
  //checking on hd floors  
  int floorret;
  floorret=check_floors_hd(pp,VELPRIM,&geom);

  if(floorret<0.)
    {
      corrected[0]=1;
      fixups[0]=1;
    }
  //************************************
  //************************************
  //checking on rad floors
    
  floorret=check_floors_rad(pp,VELPRIMRAD,&geom);
  if(floorret<0.)
    {
      corrected[1]=1;
      fixups[1]=1;
    }

  //************************************
  //************************************
  
  //************************************
  //update conserved to follow corrections on primitives
  //or to be sane
  //if(corrected[0]!=0 || corrected[1]!=0)

  p2u(pp,uu,&geom);
  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz,uu[iv]);
      set_u(p,iv,ix,iy,iz,pp[iv]);
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

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates primitives in given cell basing on global array u[]
//does not update global arrays, ignores returns 
//returns primitives from ix,iy,iz
int
calc_primitives_local(int ix,int iy,int iz,ldouble *pp)
{
  int verbose=0;
  int iv,u2pret,u2pretav;
  ldouble uu[NV],uuav[NV],ppav[NV];
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

  int verbose=1;
  int hdcorr=0;
  int radcorr=0;
  corrected[0]=corrected[1]=0;
  fixups[0]=fixups[1]=0;

  int u2pret,ret;
  ldouble ppbak[NV];
  for(u2pret=0;u2pret<NV;u2pret++)
    ppbak[u2pret]=pp[u2pret];


  //************************************
  //************************************
  //************************************
  //magneto-hydro part
  //************************************
  //************************************
  //************************************

  ldouble u0=pp[1];  

  //************************************
  //************************************
  //************************************
  //hot hydro - conserving energy
  ret=0;
#ifdef ENFORCEENTROPY
  u2pret=-1;//u2p_entropy(uu,pp,ggg);
#else
  //u2pret=u2p_hot(uu,pp,ggg);  
  u2pret=u2p_solver(uu,pp,ggg,U2P_HOT);  
  //************************************
  if(u2pret<0) 
    {
      if(verbose>1)
	printf("u2p_hot err at %d,%d,%d >>> %d <<< %e %e\n",geom->ix,geom->iy,geom->iz,u2pret,pp[0],pp[1]);

      if(u2pret==-105) 
	//negative rho but everything else right (u2pret==-105)
	{
	  /*
	  if(1)
	    {
	      printf("neg rho: %e (%d)\n",pp[0],u2pret);
	      print_Nvector(uu,NV);
	      print_Nvector(pp,NV);
	      getchar();
	    }
	  */
	  pp[0]=RHOFLOOR; 
	  ret=-1; //to ask for conserved update
	  u2pret=0;
	}
    }
#endif

  
  if(ALLOWENTROPYU2P)
    if(u2pret<0)
      {
	ret=-1;

	//u2p_entropy cannot handle negative rhos - correcting
	if(uu[0]<GAMMAMAXHD*RHOFLOOR) 
	  {
	    /*
	    printf("at %d %d %d neg uu[0] - imposing RHOFLOOR and other floors\n",geom->ix,geom->iy,geom->iz);
	    u2pret=u2p_hot(uu,pp,geom);
	    printf("u2p_hot out at %d,%d,%d >>> %d <<< %e %e\n",geom->ix,geom->iy,geom->iz,u2pret,pp[0],pp[1]);
	    getchar();
	    */

	    //using old state to estimate the correction
	    pp[0]=RHOFLOOR;
	    check_floors_hd(pp,VELPRIM,&geom);
	    pp[5]=calc_Sfromu(pp[0],pp[1]);
	    p2u(pp,uu,&geom);
	   
	    u2pret=0;
	  }
	else //regular entropy
	  {

	    //************************************
	    //************************************
	    //************************************
	    //entropy solver - conserving entropy
	    //u2pret=u2p_entropy(uu,pp,ggg);
	    u2pret=u2p_solver(uu,pp,ggg,U2P_ENTROPY);  

	    //************************************

	    if(verbose>1)
	      {
		printf("u2p_entr     >>> %d <<< %e > %e\n",u2pret,u0,pp[1]);
	      }
    
	    if(u2pret<0)
	      {
		if(verbose>0 && u2pret!=-103 && u2pret!=-107)
		  {
		    printf("u2p_entr err No. %d > %e %e %e > %e %e > %d %d %d\n",u2pret,uu[0],uu[1],uu[5],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
		  }

		if(u2pret==-107)
		  //solver converged but p2u(u2p()).neq.1
		  //requesting fixup
		  {
		    ret=-2;
		  }

		if(u2pret==-103) 
		  //solver went rho->D meaning entropy too large 
		  //imposing URHOLIMIT 
		  {		
		    //u2pret=u2p_hotmax(uu,pp,ggg);
		    u2pret=u2p_solver(uu,pp,ggg,U2P_HOTMAX);
		    if(u2pret<0)
		      {
			if(verbose>0)
			  {
			    printf("u2p_hotmax err No. %d > %e %e %e > %e %e > %d %d %d\n",u2pret,uu[0],uu[1],uu[5],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
			  }
			
			//should not happen but if happens use the old state to impose URHOLIMIT
			pp[1]=UURHORATIOMAX*pp[0];
			check_floors_hd(pp,VELPRIM,&geom);
			pp[5]=calc_Sfromu(pp[0],pp[1]);
			p2u(pp,uu,&geom);	

			//no need for another entropy solver - p2u does its job
			u2pret=0;
		      }		    
		  }
	      }	
	  }
      }

  if(ALLOWCOLDU2P)
    if(u2pret<0.)
      {
	//***********************************
	//cold RHD - assuming u=SMALL
	ret=-2;
	//u2pret=u2p_cold(uu,pp,ggg);
	u2pret=u2p_solver(uu,pp,ggg,U2P_COLD);
	//************************************

	if(u2pret<0)
	  if(verbose>0)
	    {
	      printf("u2p_cold err > %e %e > %e %e > %d %d %d\n",uu[0],uu[1],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
	    }
      }

  if(u2pret<0)
    {
      //************************************
      //leaving unchanged primitives - should not happen
      ret=-3;
      for(u2pret=0;u2pret<NV;u2pret++)
	pp[u2pret]=ppbak[u2pret];	  
      //************************************
    }
  
  
  //************************************
  //************************************
  if(ret<0) //to update conserved
    hdcorr=1;
  //************************************
  //************************************


  //************************************
  //************************************
  if(ret<-1) //request fixup when entropy failed
    fixups[0]=1;
  else
    fixups[0]=0;
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
  int radcor,radret;
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
  //************************************

  if(hdcorr>0) corrected[0]=1;
  if(radcorr>0) corrected[1]=1;

  return ret;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//checks if hydro primitives make sense
int
check_floors_hd(ldouble *pp, int whichvel,void *ggg)
{
  int verbose=0;
  int ret=0;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

#ifdef TRACER
  if(pp[TRA]<0.) {pp[0]=0.; ret=-1; if(verbose) printf("hd_floors CASE TRA 1\n");}
#endif

  //absolute rho
  if(pp[0]<RHOFLOOR) {pp[0]=RHOFLOOR; ret=-1; if(verbose) printf("hd_floors CASE 1\n");}

  //uint/rho ratios  
  if(pp[1]<UURHORATIOMIN*pp[0]) 
    {
      if(verbose) {printf("hd_floors CASE 2 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],pp[1]);}//getchar();}
      pp[1]=UURHORATIOMIN*pp[0]; //increasing uint
      ret=-1;

    }

  if(pp[1]>UURHORATIOMAX*pp[0]) 
    {
      pp[1]=UURHORATIOMAX*pp[0]; //decreasing uint
      //pp[0]=pp[1]/UURHORATIOMAX; //increasing rho
      ret=-1;      
      if(verbose) printf("hd_floors CASE 3 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],pp[1]);
    }

  //Rho to rho max
  //TODO: calculate rho max
  //ldouble rhomax=100.;
  //if(pp[0]<RHORHOMAXRATIOMIN*rhomax) {pp[0]=RHORHOMAXRATIOMIN*rhomax; ret=-1; if(verbose) printf("hd_floors CASE 5\n");}

  //updates entropy after floor corrections
  pp[5]=calc_Sfromu(pp[0],pp[1]);

  return ret;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//routines for various types of inversions
//**********************************************************************
//**********************************************************************
//**********************************************************************

/********************************************
Harm u2p_hot
********************************************/

static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq)
{
  FTYPE W=Wp+D;

  return( (GAMMA - 1.) * (1. - vsq) /  GAMMA ) ;

}

// 1 / (d(u+p)/dp)
FTYPE compute_idwmrho0dp(FTYPE wmrho0)
{
  return(GAMMAM1/GAMMA);
}


// 1 / (drho0/dp) holding wmrho0 fixed
FTYPE compute_idrho0dp(FTYPE wmrho0)
{
  return(0.0);
}


int
f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df)
{

  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  
  ldouble W=Wp+D;

  FTYPE W3,X3,Ssq,Wsq,X,X2; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X2 = X*X;
  X3 = X2*X;
  //  return -(Qn+W)*(GAMMA/GAMMAM1)+W*(1.-Qt2/W/W)-D*sqrt(1.-Qt2/W/W);   

  //a bit more clear

  ldouble v2 = Qt2/W/W;
  ldouble vsq=v2;
  ldouble gamma2 = 1./(1.-v2);
  ldouble gammasq=gamma2;
  ldouble gamma = sqrt(gamma2);
  ldouble w = W/gamma2;
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0;
  ldouble u = (w - rho0) / GAMMA;
  ldouble p = (GAMMA-1)*u;

  //*f= Qn + W - p;

  *f = Qn + W - p + 0.5*Bsq*(1.+vsq) - QdotBsq/2./Wsq;

  // dp/dW = dp/dW + dP/dv^2 dv^2/dW
    
  ldouble dvsq=(-2.0/X3 * ( Qtsq  +  QdotBsq * (3.0*W*X + Bsq*Bsq)/W3));
  ldouble dp1 = dpdWp_calc_vsq(Wp, D, vsq ); // vsq can be unphysical

  ldouble idwmrho0dp=compute_idwmrho0dp(wmrho0);
  ldouble dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp);

  ldouble drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma
  ldouble idrho0dp = compute_idrho0dp(wmrho0);

  ldouble dp2 =   drho0dvsq *idrho0dp  +   dwmrho0dvsq *idwmrho0dp;

  ldouble dpdW = dp1  + dp2*dvsq; // dp/dW = dp/dWp

  *df=1.-dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;

  return 0;  
}

/********************************************
Harm u2p_entropy
********************************************/

// p(rho0, w-rho0 = u+p)
FTYPE pressure_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0)
{
  return(IGAMMAR*wmrho0) ;
}


// local aux function
FTYPE compute_inside_entropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE pressure,indexn,insideentropy;

  pressure=pressure_wmrho0_idealgas(rho0,wmrho0);
  indexn=1.0/GAMMAM1;

  // Don't limit rho0 and pressure since this is used for iterative scheme that requires to know if beyond valid domain or not.  Na winll be terminated during inversion.
  //  if(rho0<SMALL) rho0=SMALL;
  //  if(pressure<SMALL) pressure=SMALL;
  
  insideentropy=pow(pressure,indexn)/pow(rho0,indexn+1.0);

  return(insideentropy);
}

// specific entropy as function of rho0 and internal energy (u)
// Ss(rho0,\chi=u+p)
// specific entropy = \ln( p^n/\rho^{n+1} )
FTYPE compute_specificentropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE insideentropy,specificentropy;

  insideentropy=compute_inside_entropy_wmrho0_idealgas(rho0, wmrho0);
  
  specificentropy=log(insideentropy);

  return(specificentropy);

}

// used for utoprim_jon when doing entropy evolution
// dSspecific/d\chi
FTYPE compute_dspecificSdwmrho0_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdchi;

  dSdchi = 1.0/(GAMMAM1*wmrho0);

  // Again, GAMMA->1 means dSdchi->\infty unless \chi->0 or rho0->0

  return(dSdchi);

}

// dSspecific/drho0
FTYPE compute_dspecificSdrho_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dSdrho;
  
  dSdrho=GAMMA/((1.0-GAMMA)*rho0);

  return(dSdrho);
}

/* evaluate dv^2/dW */
// does NOT depend on EOS
// Note that this does NOT use Qdotn or Qdotnp (from energy equation) so works for entropy evolution too
static FTYPE dvsq_dW(FTYPE W, FTYPE *wglobal,FTYPE Bsq,FTYPE QdotB,FTYPE QdotBsq,FTYPE Qtsq,FTYPE Qdotn,FTYPE Qdotnp,FTYPE D,FTYPE Sc, int whicheos, FTYPE *EOSextra)
{
  FTYPE W3,X3,Ssq,Wsq,X;
 
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X3 = X*X*X;

  // if(fabs(Bsq)==0.0) Ssq=0.0;
  // else Ssq = QdotBsq / Bsq;

  //return( -2.*( Ssq * ( 1./W3 - 1./X3 )  +  Qtsq / X3 ) ); 
  // return( -2.*( W3*Qtsq + QdotBsq * ( 3*W*X + Bsq*Bsq ) ) / ( W3 * X3 )   );

  return( -2.0/X3 * ( Qtsq  +  QdotBsq * (3.0*W*X + Bsq*Bsq)/W3  )  ); // RIGHT (avoids catastrophic cancellation with W^3 term in numerator)

  // return( -2.*( Qtsq/X3  +  QdotBsq/Bsq * (1.0/W3 - 1.0/X3)  )  ); // RIGHT (said was WRONG!)


}


//**********************************************************************
//**********************************************************************
//**********************************************************************
int
f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df)
{
  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Sc=cons[5];
 
  ldouble W=Wp+D;

  ldouble v2 = Qt2/W/W;
  ldouble gamma2 = 1./(1.-v2);
  ldouble gammasq = gamma2;
  ldouble gamma = sqrt(gamma2);
  ldouble w = W/gamma2;
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0;
  ldouble u = wmrho0 / GAMMA;
  ldouble p = (GAMMA-1)*u;

  ldouble Ssofchi=compute_specificentropy_wmrho0_idealgas(rho0,wmrho0);

  *f= -Sc + D*Ssofchi;

  FTYPE dSsdW,dSsdvsq,dSsdWp,dScprimedWp,dSsdrho,dSsdchi;
  FTYPE dvsq,dwmrho0dW,drho0dW;
  FTYPE dwmrho0dvsq,drho0dvsq;

  dSsdrho=compute_dspecificSdrho_wmrho0_idealgas(rho0,wmrho0);
  dSsdchi=compute_dspecificSdwmrho0_wmrho0_idealgas(rho0,wmrho0);

  dwmrho0dW = 1.0/gammasq; // holding utsq fixed
  drho0dW = 0.0; // because \rho=D/\gamma and holding utsq fixed
  dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp); // holding Wp fixed
  drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma and holding Wp fixed

  FTYPE W3,X3,Ssq,Wsq,X; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X3 = X*X*X;

  dvsq=(-2.0/X3 * ( Qtsq  +  QdotBsq * (3.0*W*X + Bsq*Bsq)/W3));

  dSsdW =   drho0dW   *dSsdrho +   dwmrho0dW   *dSsdchi; // dSs/dW' holding utsq fixed
  dSsdvsq = drho0dvsq *dSsdrho +   dwmrho0dvsq *dSsdchi;
  dSsdWp = dSsdW  + dSsdvsq*dvsq; // dSs/dW = dSs/dWp [total derivative]

  dScprimedWp = D*dSsdWp;

  *df = dScprimedWp;
  
  return 0;
 
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
int
f_u2p_cold(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df)
{
  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble Sc=cons[3];

  ldouble W=Wp+D;

  ldouble v2 = Qt2/W/W;
  ldouble gamma2 = 1./(1.-v2);
  ldouble gammasq = gamma2;
  ldouble gamma = sqrt(gamma2);
  ldouble w = W/gamma2;
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0;
  ldouble u = wmrho0 / GAMMA;
  ldouble p = (GAMMA-1)*u;

  *f = u - UURHORATIOU2PCOLD*rho0;

  ldouble dWdWp=1.;
  ldouble dv2dWp=Qt2*(-2.)/W/W/W*dWdWp;
  ldouble dgammadWp=.5*1./sqrt(1./(1.-v2))*(-1.)/(1.-v2)/(1.-v2)*(-1.)*dv2dWp;
  ldouble drho0dWp=-D/gamma2*dgammadWp;
  ldouble dwdWp = dWdWp/gamma2 - 2. *W/gamma/gamma/gamma*dgammadWp;
  ldouble dudWp = 1./GAMMA*(dwdWp - drho0dWp);

  *df = dudWp - UURHORATIOU2PCOLD*drho0dWp;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//uses D,T^t_i to get solution with u=UURHORATIOMAX*rho
int
f_u2p_hotmax(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df)
{
  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];

  ldouble W=Wp+D;

  ldouble v2 = Qt2/W/W;
  ldouble gamma2 = 1./(1.-v2);
  ldouble gammasq = gamma2;
  ldouble gamma = sqrt(gamma2);
  ldouble w = W/gamma2;
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0;
  ldouble u = wmrho0 / GAMMA;
  ldouble p = (GAMMA-1)*u;

  *f = u - UURHORATIOMAX*rho0;

  ldouble dWdWp=1.;
  ldouble dv2dWp=Qt2*(-2.)/W/W/W*dWdWp;
  ldouble dgammadWp=.5*1./sqrt(1./(1.-v2))*(-1.)/(1.-v2)/(1.-v2)*(-1.)*dv2dWp;
  ldouble drho0dWp=-D/gamma2*dgammadWp;
  ldouble dwdWp = dWdWp/gamma2 - 2. *W/gamma/gamma/gamma*dgammadWp;
  ldouble dudWp = 1./GAMMA*(dwdWp - drho0dWp);

  *df = dudWp - UURHORATIOMAX*drho0dWp;
  
  return 0;
 
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//Newton-Rapshon solver 
//Etype == 0 -> hot inversion (uses D,Ttt,Tti)
//Etype == 1 -> entropy inversion (uses D,S,Tti)
//Etype == 2 -> hotmax inversion (uses D,Tti,u over rho max ratio)
//Etype == 3 -> cold inversion (uses D,Tti,u over rho min ratio)

int
u2p_solver(ldouble *uu, ldouble *pp, void *ggg,int Etype)
{
  int verbose=0;
  int i,j,k;
  ldouble rho,u,p,w,W,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
  ldouble QdotB,QdotBsq,Bcon[4],Bcov[4],Bsq;

  /****************************/
  //prepare geometry
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  /****************************/

  /****************************/
  //equations choice
  int (*f_u2p)(ldouble,ldouble*,ldouble*,ldouble*);
 if(Etype==U2P_HOT) 
   f_u2p=&f_u2p_hot;
 if(Etype==U2P_ENTROPY) 
   f_u2p=&f_u2p_entropy;
 if(Etype==U2P_HOTMAX) 
   f_u2p=&f_u2p_hotmax;
 if(Etype==U2P_COLD) 
   f_u2p=&f_u2p_cold;
  /****************************/
 
  
  if(verbose>1) {printf("********************\n");print_Nvector(uu,NV);}
  if(verbose>1) {print_Nvector(pp,NV);}

  /****************************/
  //conserved quantities etc
  
  //alpha
  alpha=sqrt(-1./GG[0][0]);

  //D
  D=uu[0]/gdetu*alpha; //uu[0]=gdetu rho ut

  //conserved entropy "S u^t"
  Sc=uu[5]/gdetu*alpha; 

  //Q_mu
  Qcov[0]=(uu[1]/gdetu-uu[0]/gdetu)*alpha;
  Qcov[1]=uu[2]/gdetu*alpha;
  Qcov[2]=uu[3]/gdetu*alpha;
  Qcov[3]=uu[4]/gdetu*alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

#ifdef MAGNFIELD
  //B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1]/gdet*alpha;
  Bcon[2]=uu[B2]/gdet*alpha;
  Bcon[3]=uu[B3]/gdet*alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
#else
  Bsq=QdotB=QdotBsq=0.;
#endif  

  //n_mu = (-alpha, 0, 0, 0)
  ncov[0]=-alpha;
  ncov[1]=ncov[2]=ncov[3]=0.;
  
  //n^mu
  indices_12(ncov,ncon,GG);

  //Q_mu n^mu = Q^mu n_mu = -alpha*Q^t
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
  FTYPE Qtsq = Qt2;

  /****************************/
  
  //initial guess for W = w gamma**2 based on primitives
  rho=pp[0];
  u=pp[1];
  utcon[0]=0.;
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  //conv_vels(utcon,ucon,VELPRIM,VEL4,gg,GG);
  //indices_21(ucon,ucov,gg);
  conv_vels(utcon,utcon,VELPRIM,VELR,gg,GG);

  /*
  ldouble bcon[4],bcov[4],bsq;
#ifdef MAGNFIELD
  bcon_calc(pp,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
#else
  bcon[0]=bcon[1]=bcon[2]=bcon[3]=0.;
  bsq=0.;
#endif
  */

  ldouble qsq=0.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=utcon[i]*utcon[j]*gg[i][j];
  ldouble gamma2=1.+qsq;
  ldouble gamma=sqrt(gamma2);

  //W
  W=(rho+GAMMA*u)*gamma2;

  if(verbose>1) printf("initial W:%e\n",W);

  /****************************/
  
  //test if does not provide reasonable gamma2
  // Make sure that W is large enough so that v^2 < 1 : 
  int i_increase = 0;
  ldouble f0,f1,dfdW;
  ldouble CONV=1.e-8;
  ldouble EPS=1.e-4;
  ldouble Wprev=W;
  ldouble cons[6]={Qn,Qt2,D,QdotBsq,Bsq,Sc};
 
  do
    {
      f0=dfdW=0.;
      if(Etype!=U2P_HOT) //entropy-like solvers require this additional check
	(*f_u2p)(W-D,cons,&f0,&dfdW);
      
      if( ((( W*W*W * ( W + 2.*Bsq ) 
	    - QdotBsq*(2.*W + Bsq) ) <= W*W*(Qtsq-Bsq*Bsq))
	  || isinf(f0) || isnan(f0)
	  || isinf(dfdW) || isnan(dfdW))	  
	  && (i_increase < 10)) //if not enough will complain later returnin negative number
	{
	  W *= 10.;
	  i_increase++;
	  continue;
	}
      else
	break;    
    }
  while(1);


  //1d Newton solver
  if(verbose>1) printf("in:%e %e %e %e %e\n",Qn,Qt2,D,QdotBsq,Bsq);

  int iter=0,fu2pret;
  do
    {
      Wprev=W;
      iter++;
     
      fu2pret=(*f_u2p)(W-D,cons,&f0,&dfdW);
     
      if(verbose>1) printf("%d %e %e %e %e\n",iter,W,f0,f1,dfdW);

      if(dfdW==0.) {W*=1.1; continue;}

      ldouble Wnew=W-f0/dfdW;
      int idump=0;
      ldouble dumpfac=1.;

      //test if goes out of bounds and damp solution if so
      //if(Wnew*Wnew<Qt2 || isnan(Wnew))
      
      do
	{
	  ldouble f0tmp,dfdWtmp;
	  f0tmp=dfdWtmp=0.;
	  if(Etype!=U2P_HOT) //entropy-like solvers require this additional check
	    (*f_u2p)(Wnew-D,cons,&f0tmp,&dfdWtmp);

	  if( ((( Wnew*Wnew*Wnew * ( Wnew + 2.*Bsq ) 
		  - QdotBsq*(2.*Wnew + Bsq) ) <= Wnew*Wnew*(Qtsq-Bsq*Bsq))
	       || isinf(f0tmp) || isnan(f0tmp)
	       || isinf(dfdWtmp) || isnan(dfdWtmp))
	      && (idump<100))
	    {
	      idump++;
	      dumpfac/=2.;
	      Wnew=W-dumpfac*f0/dfdW;
	      continue;
	    }
	  else
	    break;
	}
      while(1);
	  
      if(idump>=100) {if(verbose>0) printf("damped unsuccessfuly\n");return -101;}

      if(verbose>0) printf("damped successfuly\n");
	
      W=Wnew; 

      if(fabs(W)>BIG) 
	{
	  if(verbose>1) printf("W has gone out of bounds at %d,%d,%d\n",geom->ix,geom->iy,geom->iz); 
	  return -103;
	}
    }
  while(fabs((W-Wprev)/Wprev)>CONV && iter<50);

  if(iter>=50)
    {
      if(verbose>0) printf("iter exceeded in u2p_hot\n"); //getchar();
      return -102;
    }
  
  if(isnan(W) || isinf(W)) {if(verbose) printf("nan/inf W in u2p_solver with Etype: %d\n",Etype); return -103;}
  if(verbose>1) {printf("the end: %e\n",W); }

  //W found, let's calculate v2 and the rest
  //ldouble v2=Qt2/W/W;

  ldouble Wsq,Xsq,v2;
	
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);  
  v2 = ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);

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

  if(u<0. || gamma2<0. ||isnan(W) || isinf(W)) 
    {
      if(verbose>0) printf("neg u in u2p_hot %e %e %e %e\n",rho,u,gamma2,W);//getchar();
      return -104;
    }

  if(rho<0.) 
    {
      if(verbose>0) printf("neg rho in u2p_hot %e %e %e %e\n",rho,u,gamma2,W);//getchar();
      return -105;
    }

  //converting to VELPRIM
  conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  
  //returning new primitives
  pp[0]=rho;
  pp[1]=u;
  pp[2]=utcon[1];
  pp[3]=utcon[2];
  pp[4]=utcon[3];

  //entropy based on UU[1]
  pp[5]=calc_Sfromu(rho,u);

#ifdef TRACER
  ldouble Dtr=uu[TRA]/gdetu*alpha; //uu[0]=gdetu rho ut
  pp[TRA]=Dtr/gamma;
#endif


  // test the inversion
  ldouble uu2[NV];
  int iv;
  int lostprecision=0;
  p2u(pp,uu2,ggg);
  for(iv=0;iv<NVMHD;iv++)
    {
      if(Etype==U2P_HOT) if(iv==5) continue;
      if(Etype==U2P_ENTROPY) if(iv==1) continue;
      if(Etype==U2P_COLD || Etype==U2P_HOTMAX) if(iv==1 || iv==5) continue;
      if(((iv==0 || iv==1) && fabs(uu2[iv]-uu[iv])/fabs(uu[iv]+uu2[iv])>1.e-3))
	//|| ((iv>1) && fabs(uu2[iv]-uu[iv])/fabs(uu[iv])>1.e-6 && fabs(uu[iv])>1.e-6))
	lostprecision=1;
    }     

  if(lostprecision)
    {
      if(verbose>0) printf("u2p_hot lost precision:\n");
      //print_Nvector(uu,NV);
      //print_Nvector(uu2,NV);  
      return -106;
      //getchar();
    }

  if(verbose>1) {print_Nvector(pp,NV); getchar();}

  return 0; //ok

}

//**********************************************************************
//**********************************************************************
//**********************************************************************

//no - longer used
//compactified into u2p_solver

//**********************************************************************
//**********************************************************************
//**********************************************************************
//uses D,Ttt,Tti
int
u2p_hot(ldouble *uu, ldouble *pp, void *ggg)
{
   struct geometry *geom
   = (struct geometry *) ggg;

   ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int verbose=0;
  int i,j,k;
  ldouble rho,u,p,w,W,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
  ldouble QdotB,QdotBsq,Bcon[4],Bcov[4],Bsq;
  
  if(verbose>1) {printf("********************\n");print_Nvector(uu,NV);}
  if(verbose>1) {print_Nvector(pp,NV);}

  //alpha
  alpha=sqrt(-1./GG[0][0]);

  //D
  D=uu[0]/gdetu*alpha; //uu[0]=gdetu rho ut

  //conserved entropy "S u^t"
  Sc=uu[5]/gdetu*alpha; 

  //Q_mu
  Qcov[0]=(uu[1]/gdetu-uu[0]/gdetu)*alpha;
  Qcov[1]=uu[2]/gdetu*alpha;
  Qcov[2]=uu[3]/gdetu*alpha;
  Qcov[3]=uu[4]/gdetu*alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

#ifdef MAGNFIELD
  //B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1]/gdet*alpha;
  Bcon[2]=uu[B2]/gdet*alpha;
  Bcon[3]=uu[B3]/gdet*alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
#else
  Bsq=QdotB=QdotBsq=0.;
#endif
  

  //n_mu = (-alpha, 0, 0, 0)
  ncov[0]=-alpha;
  ncov[1]=ncov[2]=ncov[3]=0.;
  
  //n^mu
  indices_12(ncov,ncon,GG);

  //Q_mu n^mu = Q^mu n_mu = -alpha*Q^t
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
  ldouble gamma=sqrt(gamma2);

  //W
  W=(rho+GAMMA*u)*gamma2;

  if(verbose>1) printf("initial W:%e\n",W);
 
  //test if does not provide reasonable gamma2
  //TODO!
  if(W*W<Qt2)
    {
      W=1.001*sqrt(Qt2);
      if(verbose>1) printf("corrected W:%e\n",W);
    }

  //1d Newton solver
  ldouble CONV=1.e-8;
  ldouble EPS=1.e-4;
  ldouble Wprev=W;
  ldouble f0,f1,dfdW;
  ldouble cons[6]={Qn,Qt2,D,QdotBsq,Bsq,Sc};
  if(verbose>1) printf("in:%e %e %e %e %e\n",Qn,Qt2,D,QdotBsq,Bsq);

  int iter=0;
  do
    {
      Wprev=W;
      iter++;
      f_u2p_hot(W,cons,&f0,&dfdW);

      //f_u2p_hot(W*(1.+EPS),cons,&f1,&dfdW);
      //dfdW=(f1-f0)/(EPS*W);

      if(verbose>1) printf("%d %e %e %e %e\n",iter,W,f0,f1,dfdW);

      if(dfdW==0.) {W*=1.1; continue;}

      ldouble Wnew=W-f0/dfdW;

      //test if does produce nan and damp solution if so
      if(Wnew*Wnew<Qt2 || isnan(Wnew))
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
	  
	  if(idump>=100) {if(verbose>0) printf("damped unsuccessfuly\n");return -101;}

	  if(verbose>0) printf("damped successfuly\n");
	}

      W=Wnew; 
    }
  while(fabs((W-Wprev)/Wprev)>CONV && iter<50);

  if(iter>=50)
    {
      if(verbose>0) printf("iter exceeded in u2p_hot\n"); //getchar();
      return -102;
    }
  
  if(isnan(W) || isinf(W)) {printf("nan/inf W in u2p_hot: %e\n",W); exit(0);return -103;}
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

  if(u<0. || gamma2<0. ||isnan(W) || isinf(W)) 
    {
      if(verbose>0) printf("neg u in u2p_hot %e %e %e %e\n",rho,u,gamma2,W);//getchar();
      return -104;
    }

  if(rho<0.) 
    {
      if(verbose>0) printf("neg rho in u2p_hot %e %e %e %e\n",rho,u,gamma2,W);//getchar();
      return -105;
    }

  //converting to VELPRIM
  conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  
  //returning new primitives
  pp[0]=rho;
  pp[1]=u;
  pp[2]=utcon[1];
  pp[3]=utcon[2];
  pp[4]=utcon[3];

  //entropy based on UU[1]
  pp[5]=calc_Sfromu(rho,u);

#ifdef TRACER
  ldouble Dtr=uu[TRA]/gdetu*alpha; //uu[0]=gdetu rho ut
  pp[TRA]=Dtr/gamma;
#endif


  // test the inversion
  ldouble uu2[NV];
  int iv;
  int lostprecision=0;
  p2u(pp,uu2,ggg);
  for(iv=0;iv<NVMHD;iv++)
    {
      if(iv==5) continue;
       if(((iv==0 || iv==1) && fabs(uu2[iv]-uu[iv])/fabs(uu[iv]+uu2[iv])>1.e-3))
	//|| ((iv>1) && fabs(uu2[iv]-uu[iv])/fabs(uu[iv])>1.e-6 && fabs(uu[iv])>1.e-6))
	lostprecision=1;
    }     

  if(lostprecision)
    {
      if(verbose>0) printf("u2p_hot lost precision:\n");
      //print_Nvector(uu,NV);
      //print_Nvector(uu2,NV);  
      return -106;
      //getchar();
    }

  if(verbose>1) {print_Nvector(pp,NV); getchar();}

  return 0; //ok

}

int
u2p_entropy(ldouble *uu, ldouble *pp, void *ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int verbose=1;
  
  //if(uu[0]<1.e-49) verbose=2;
  int superverbose=0;

  if(superverbose)
    {
      printf("start of u2p_entropy()\n");
      printf("in:\n");
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);

    }

  int i,j,k;
  ldouble rho,u,p,w,W,Wp,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
  ldouble QdotB,QdotBsq,Bcon[4],Bcov[4],Bsq;
 
  if(verbose>1 && !superverbose) {printf("********************\n");print_Nvector(uu,NV);}
  if(verbose>1 && !superverbose) {print_Nvector(pp,NV);}

 
  //alpha
  alpha=sqrt(-1./GG[0][0]);
  
  //conserved entropy "S u^t"
  Sc=uu[5]/gdetu*alpha; 

  //D
  D=uu[0]/gdetu*alpha;

  //Q_mu
  Qcov[0]=(uu[1]/gdetu-uu[0]/gdetu)*alpha;
  Qcov[1]=uu[2]/gdetu*alpha;
  Qcov[2]=uu[3]/gdetu*alpha;
  Qcov[3]=uu[4]/gdetu*alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

#ifdef MAGNFIELD
  //B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1]/gdet*alpha;
  Bcon[2]=uu[B2]/gdet*alpha;
  Bcon[3]=uu[B3]/gdet*alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
#else
  Bsq=QdotB=QdotBsq=0.;
#endif

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

  //will solve for Wp = W - D
  
  //initial guess
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

  //Wp
  Wp=W-D; 

  //w
  w = W/gamma2;

  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0; //u+p

  //1d Newton solver
  ldouble CONV=1.e-8;
  ldouble EPS=1.e-6;
  ldouble Wpprev=Wp,Wpprev2=Wp;
  ldouble f0,f1,dfdWp,Wpnew,v2,ut2;
  ldouble cons[6]={Qn,Qt2,D,QdotBsq,Bsq,Sc};
  //if(verbose>1) printf("in:%e %e %e\n",Qn,Qt2,D);
  int idump;
  int iter=0;

  //testing if initial guess works
  f_u2p_entropy(Wp,cons,&f0,&dfdWp);
  
  //initial guess wrong - leading to negative pressure
  //arbitrary look up
  if(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp))
    {
      W=1.00001*sqrt(D*D+Qt2);
      Wp=W-D;
      v2 = Qt2/W/W;
      gamma2 = 1./(1.-v2);
      gamma = sqrt(gamma2);
      w = W/gamma2;
      rho0 = D/gamma;
      wmrho0 = w - rho0;

      f_u2p_entropy(Wp,cons,&f0,&dfdWp);
      if(verbose>1) printf(">>> %e %e\n",f0,wmrho0);      
    }
    
  //still not enough, iterating
  idump=0;  
  if(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp))
    {
      W=sqrt(D*D+Qt2);
      do
	{
	  idump++;
	  W*=1.1;
	  Wp=W-D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  f_u2p_entropy(Wp,cons,&f0,&dfdWp);
	    
	  if(verbose>1) 
	    {
	      printf("alala %e %e %e %e %e %e\n",W,Wp,wmrho0,f0,compute_specificentropy_wmrho0_idealgas(rho0,wmrho0),cons[3]);
	    }

	  if(idump>50)
	    {
	      print_Nvector(uu,NV);
	      printf("%e %e %e %e %e %e %e %e\n",W,Wp,D,rho0,wmrho0,f0,compute_specificentropy_wmrho0_idealgas(rho0,wmrho0),cons[3]);
	      printf("idump exceeded before cons: %e %e %e %e\n",cons[0],cons[1],cons[2],cons[3]);
	      return -101;//getchar();
	    }
	}
      while(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp));
    }

  if(verbose>1) printf("initial:%e %e %e %e\n",Wp,rho0,wmrho0,compute_inside_entropy_wmrho0_idealgas(rho0,wmrho0));

  ldouble dampfac;
  ldouble f0temp,dfdWptemp;
  do
    {
      Wpprev=Wp;
      iter++;
      dampfac=1.;

      f_u2p_entropy(Wp,cons,&f0,&dfdWp);

      do
	{
	  Wp=Wpprev-dampfac*f0/dfdWp;
	  
	  W=Wp+D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  u = wmrho0 / GAMMA;

	  f_u2p_entropy(Wp,cons,&f0temp,&dfdWptemp);

	  if(verbose>1) printf("substep: %f %e %e %e %e\n",dampfac,Wp,f0temp,dfdWptemp,u);
	  dampfac/=2.;

	  if(dampfac<1.e-7) break;
	}
      while(u<0. || isnan(f0temp) || isinf(f0temp) || isnan(dfdWptemp) || isinf(dfdWptemp));

      if(dampfac<1.e-7) {printf("damping unsuccesful at entropy\n");return -102; getchar(); }
      else  if(verbose>1) {printf("damping succesful\n"); }

      if(verbose>1) 
	{
	  W=Wp+D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  u = wmrho0 / GAMMA;

	  printf("%d %e %e %e -> wmrho %e\n",iter,Wp,f0,dfdWp,u);
	}

      if(fabs(Wp)>BIG) 
	{
	  if(verbose>1) printf("Wp has gone out of bounds at %d,%d,%d\n",geom->ix,geom->iy,geom->iz); 

	  //means that entropy to big
	  //imposing uint/rho = URHOMAXRATIO amd calling u2p_hotmax()

	  return -103;

	  /*
	  print_Nvector(uu,NV);	  
	  ldouble ddd;
	  for(ddd=1.00001;ddd<30000.;ddd=1.+(ddd-1.)*2.)
	    {
	      
	      W=ddd*sqrt(D*D+Qt2);
	      Wp=W-D;
	      v2 = Qt2/W/W;
	      gamma2 = 1./(1.-v2);
	      gamma = sqrt(gamma2);
	      w = W/gamma2;
	      rho0 = D/gamma;
	      wmrho0 = w - rho0;
	      ldouble Ssofchi=compute_specificentropy_wmrho0_idealgas(rho0,wmrho0);

	      f_u2p_entropy(Wp,cons,&f0,&dfdWp);
	      printf("%f %e %e %e %e %e %e -> %e %e -> %e %e\n",ddd,W,Wp,gamma,wmrho0,rho0,alpha,f0,dfdWp,D*Ssofchi,Sc);
	    }
	  getchar();
	  */

	}

      if(isnan(Wp) || isinf(Wp)) {printf("nan/inf Wp: %e %e %e %e\n",Wp,f0,dfdWp,W);  return -104;getchar();}
    }
  //  while(( fabs((Wp-Wpprev)/Wpprev)>CONV || u<0. || rho<0.) && iter<50);
  while(( fabs((Wp-Wpprev)/Wpprev)>CONV) && iter<50);
  
  if(iter>=50)
    { 
      if(fabs(1.-gamma)<1.e-10) //going out of bounds
	return -103; 
      printf("iter exceeded in u2p_entropy at %d %d %d\n",geom->ix,geom->iy,geom->iz);
      return -105;
    }
  

  if(verbose>1) {printf("the end: %e\n",Wp); }//getchar();}

  //Wp found, let's calculate W, v2 and the rest
  W=Wp+D;
  v2=Qt2/W/W;
  gamma2=1./(1.-v2);
  ut2 = Qt2/(W*W - Qt2);
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
      printf("neg u rho in u2p_entropy %e %e %e %e\n",rho,u,gamma2,W);
      //getchar();
      return -106;
    }

  //converting to VELPRIM
  conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  
  //returning new primitives
  pp[0]=rho;
  pp[1]=u;
  pp[2]=utcon[1];
  pp[3]=utcon[2];
  pp[4]=utcon[3];

  //entropy corresponding to UU[5]
  pp[5]=calc_Sfromu(rho,u);

#ifdef TRACER
  ldouble Dtr=uu[TRA]/gdetu*alpha; //uu[0]=gdetu rho ut
  pp[TRA]=Dtr/gamma;
#endif

  //precision test
  ldouble uu2[NV];
  int iv;
  int lostprecision=0;
  p2u(pp,uu2,ggg);
  for(iv=0;iv<NVMHD;iv++)
    {
      if(iv==1) continue;
      if(iv==2 || iv==3 || iv==4) continue; //do not check momenta which could be zero
      if(fabs(uu2[iv]-uu[iv])/fabs(uu[iv]+uu2[iv])>1.e-3 && fabs(uu[iv])>SMALL) lostprecision=1;
    }
     
  if(lostprecision)
    {
      //printf("u2p_entropy lost precision at %d %d %d:\n",geom->ix,geom->iy,geom->iz);
      return -107;
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);
      print_Nvector(uu2,NV); 
      getchar();
    }

  if(superverbose)
    {
      printf("superverbose in u2p_entropy:\n");
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);
      print_Nvector(uu2,NV);  
    }

  if(verbose>1 && !superverbose) {print_Nvector(pp,NV); }//getchar();}

  return 0; //ok!

}

int
u2p_cold(ldouble *uu, ldouble *pp, void *ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int verbose=2;
  int superverbose=1;

  if(superverbose)
    {
      printf("start of u2p_cold()\n");
      printf("in:\n");
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);

    }

  int i,j,k;
  ldouble rho,u,p,w,W,Wp,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
  ldouble QdotB,QdotBsq,Bcon[4],Bcov[4],Bsq;
 
  if(verbose>1 && !superverbose) {printf("********************\n");print_Nvector(uu,NV);}
  if(verbose>1 && !superverbose) {print_Nvector(pp,NV);}

 
  //alpha
  alpha=sqrt(-1./GG[0][0]);
  
  //conserved entopy "S u^t"
  Sc=uu[5]/gdetu*alpha; //alpha?

  //D
  D=uu[0]/gdetu*alpha;

  //Q_mu
  Qcov[0]=(uu[1]/gdetu-uu[0]/gdetu)*alpha;
  Qcov[1]=uu[2]/gdetu*alpha;
  Qcov[2]=uu[3]/gdetu*alpha;
  Qcov[3]=uu[4]/gdetu*alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

#ifdef MAGNFIELD
  //B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1]/gdet*alpha;
  Bcon[2]=uu[B2]/gdet*alpha;
  Bcon[3]=uu[B3]/gdet*alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
#else
  Bsq=QdotB=QdotBsq=0.;
#endif

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

  //will solve for Wp = W - D
  
  //initial guess
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

  //Wp
  Wp=W-D; 

  //w
  w = W/gamma2;

  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0; //u+p

  //1d Newton solver
  ldouble CONV=1.e-8;
  ldouble EPS=1.e-6;
  ldouble Wpprev=Wp,Wpprev2=Wp;
  ldouble f0,f1,dfdWp,Wpnew,v2,ut2;
  ldouble cons[4]={Qn,Qt2,D,Sc};
  //if(verbose>1) printf("in:%e %e %e\n",Qn,Qt2,D);
  int idump;
  int iter=0;

  //testing if initial guess works
  f_u2p_cold(Wp,cons,&f0,&dfdWp);
  if(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp))
    {
      W=1.01*sqrt(D*D+Qt2);
      Wp=W-D;
      v2 = Qt2/W/W;
      gamma2 = 1./(1.-v2);
      gamma = sqrt(gamma2);
      w = W/gamma2;
      rho0 = D/gamma;
      wmrho0 = w - rho0;

      f_u2p_cold(Wp,cons,&f0,&dfdWp);
      if(verbose>1) printf(">>> %e %e\n",f0,wmrho0);      
    }
  
  idump=0;
  if(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp))
    {
      W=sqrt(D*D+Qt2);
      do
	{
	  idump++;
	  W*=1.001;
	  Wp=W-D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  f_u2p_cold(Wp,cons,&f0,&dfdWp);
	    
	  if(verbose>1 || 1) 
	    {
	      printf("alala %e %e %e %e %e %e\n",W,Wp,wmrho0,f0,compute_specificentropy_wmrho0_idealgas(rho0,wmrho0),cons[3]);
	    }

	  if(idump>50)
	    {
	      print_Nvector(uu,NV);
	      printf("%e %e %e %e %e %e\n",W,Wp,wmrho0,f0,compute_specificentropy_wmrho0_idealgas(rho0,wmrho0),cons[3]);
	      printf("idump exceeded before cons: %e %e %e %e\n",cons[0],cons[1],cons[2],cons[3]);
	      return -1;//getchar();
	    }
	}
      while(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp));
    }

  if(verbose>1) printf("initial:%e %e %e %e\n",Wp,rho0,wmrho0,compute_inside_entropy_wmrho0_idealgas(rho0,wmrho0));

  ldouble dampfac;
  ldouble f0temp,dfdWptemp;
  do
    {
      Wpprev=Wp;
      iter++;
      dampfac=1.;

      f_u2p_cold(Wp,cons,&f0,&dfdWp);
      //f_u2p_cold(Wp*(1.+EPS),cons,&f1,&dfdWp);
      //dfdWp=(f1-f0)/(EPS*Wp);

      do
	{
	  Wp=Wpprev-dampfac*f0/dfdWp;
	  
	  W=Wp+D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  u = wmrho0 / GAMMA;

	  f_u2p_cold(Wp,cons,&f0temp,&dfdWptemp);

	  if(verbose>1) printf("substep: %f %e %e %e %e\n",dampfac,Wp,f0temp,dfdWptemp,u);
	  dampfac/=2.;

	  if(dampfac<1.e-7) break;
	}
      while(u<0. || isnan(f0temp) || isinf(f0temp) || isnan(dfdWptemp) || isinf(dfdWptemp));

      if(dampfac<1.e-7) {printf("damping unsuccesful at entropy\n");return -1; getchar(); }
      else  if(verbose>1) {printf("damping succesful\n"); }

      if(verbose>1) 
	{
	  W=Wp+D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  u = wmrho0 / GAMMA;

	  printf("%d %e %e %e -> wmrho %e\n",iter,Wp,f0,dfdWp,u);
	}

      //      if(fabs(dfdWp)<SMALL ) {if(verbose) printf("derivative zero. asssuming found solution\n"); return -1;getchar(); break;}
      if(fabs(Wp)>BIG) {if(verbose) printf("Wp has gone out of bounds\n"); return -1;getchar(); break;}

      if(isnan(Wp) || isinf(Wp)) {printf("nan/inf Wp: %e %e %e %e\n",Wp,f0,dfdWp,W);  return -1;getchar();}
    }
  //  while(( fabs((Wp-Wpprev)/Wpprev)>CONV || u<0. || rho<0.) && iter<50);
  while(( fabs((Wp-Wpprev)/Wpprev)>CONV) && iter<50);

  if(iter>=50)
    {
      if(verbose>0 || 1) printf("iter exceeded in u2p_cold at %d %d %d\n",geom->ix,geom->iy,geom->iz);
      //getchar();
      return -1;
    }
  

  if(verbose>1) {printf("the end: %e\n",Wp); }//getchar();}

  //Wp found, let's calculate W, v2 and the rest
  W=Wp+D;
  v2=Qt2/W/W;
  gamma2=1./(1.-v2);
  ut2 = Qt2/(W*W - Qt2);
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
      if(verbose>0 || 1) printf("neg u rho in u2p_cold %e %e %e %e\n",rho,u,gamma2,W);
      //getchar();
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
  pp[5]=calc_Sfromu(rho,u);
#ifdef TRACER
  ldouble Dtr=uu[TRA]/gdetu*alpha; //uu[0]=gdetu rho ut
  pp[TRA]=Dtr/gamma;
#endif

  ldouble uu2[NV];
  int iv;
  int lostprecision=0;
  p2u(pp,uu2,ggg);
  for(iv=0;iv<NVMHD;iv++)
    {
      if(iv==1) continue;
      //      if(fabs(uu2[iv]-uu[iv])/fabs(uu[iv]+uu2[iv])>1.e-3 && fabs(uu[iv])>SMALL) lostprecision;
    }
     

  if(superverbose)
    {
      printf("superverbose or u2p_cold lost precision:\n");
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);
      print_Nvector(uu2,NV);  
      getchar();
    }

  if(verbose>1 && !superverbose) {print_Nvector(pp,NV); }//getchar();}

  return 0;

}

int
u2p_hotmax(ldouble *uu, ldouble *pp, void *ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int verbose=0;
  int superverbose=0;

  if(superverbose)
    {
      printf("start of u2p_hotmax()\n");
      printf("in:\n");
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);

    }

  int i,j,k;
  ldouble rho,u,p,w,W,Wp,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
  ldouble QdotB,QdotBsq,Bcon[4],Bcov[4],Bsq;
   
  if(verbose>1 && !superverbose) {printf("********************\n");print_Nvector(uu,NV);}
  if(verbose>1 && !superverbose) {print_Nvector(pp,NV);}

 
  //alpha
  alpha=sqrt(-1./GG[0][0]);
  
  //conserved entopy "S u^t"
  Sc=uu[5]/gdetu*alpha; //alpha?

  //D
  D=uu[0]/gdetu*alpha;

  //Q_mu
  Qcov[0]=(uu[1]/gdetu-uu[0]/gdetu)*alpha;
  Qcov[1]=uu[2]/gdetu*alpha;
  Qcov[2]=uu[3]/gdetu*alpha;
  Qcov[3]=uu[4]/gdetu*alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

 #ifdef MAGNFIELD
  //B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1]/gdet*alpha;
  Bcon[2]=uu[B2]/gdet*alpha;
  Bcon[3]=uu[B3]/gdet*alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
#else
  Bsq=QdotB=QdotBsq=0.;
#endif

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

  //will solve for Wp = W - D
  
  //initial guess
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

  //Wp
  Wp=W-D; 

  //w
  w = W/gamma2;

  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0; //u+p

  //1d Newton solver
  ldouble CONV=1.e-8;
  ldouble EPS=1.e-6;
  ldouble Wpprev=Wp,Wpprev2=Wp;
  ldouble f0,f1,dfdWp,Wpnew,v2,ut2;
  ldouble cons[4]={Qn,Qt2,D,Sc};
  //if(verbose>1) printf("in:%e %e %e\n",Qn,Qt2,D);
  int idump;
  int iter=0;

  //testing if initial guess works
  f_u2p_hotmax(Wp,cons,&f0,&dfdWp);
  if(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp))
    {
      W=1.01*sqrt(D*D+Qt2);
      Wp=W-D;
      v2 = Qt2/W/W;
      gamma2 = 1./(1.-v2);
      gamma = sqrt(gamma2);
      w = W/gamma2;
      rho0 = D/gamma;
      wmrho0 = w - rho0;

      f_u2p_hotmax(Wp,cons,&f0,&dfdWp);
      if(verbose>1) printf(">>> %e %e\n",f0,wmrho0);      
    }
  
  idump=0;
  if(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp))
    {
      W=sqrt(D*D+Qt2);
      do
	{
	  idump++;
	  W*=1.001;
	  Wp=W-D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  f_u2p_hotmax(Wp,cons,&f0,&dfdWp);
	    
	  if(verbose>1 || 1) 
	    {
	      printf("alala %e %e %e %e %e %e\n",W,Wp,wmrho0,f0,compute_specificentropy_wmrho0_idealgas(rho0,wmrho0),cons[3]);
	    }

	  if(idump>50)
	    {
	      print_Nvector(uu,NV);
	      printf("%e %e %e %e %e %e\n",W,Wp,wmrho0,f0,compute_specificentropy_wmrho0_idealgas(rho0,wmrho0),cons[3]);
	      printf("idump exceeded before cons: %e %e %e %e\n",cons[0],cons[1],cons[2],cons[3]);
	      return -1;//getchar();
	    }
	}
      while(isnan(f0)|| isnan(dfdWp) || isinf(f0) || isinf(dfdWp));
    }

  if(verbose>1) printf("initial:%e %e %e %e\n",Wp,rho0,wmrho0,compute_inside_entropy_wmrho0_idealgas(rho0,wmrho0));

  ldouble dampfac;
  ldouble f0temp,dfdWptemp;
  do
    {
      Wpprev=Wp;
      iter++;
      dampfac=1.;

      f_u2p_hotmax(Wp,cons,&f0,&dfdWp);
      f_u2p_hotmax(Wp*(1.+EPS),cons,&f1,&dfdWp);
      //printf("%e %e\n",dfdWp,(f1-f0)/(EPS*Wp));getchar();
      //dfdWp=(f1-f0)/(EPS*Wp);


      do
	{
	  Wp=Wpprev-dampfac*f0/dfdWp;
	  
	  W=Wp+D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  u = wmrho0 / GAMMA;

	  f_u2p_hotmax(Wp,cons,&f0temp,&dfdWptemp);

	  if(verbose>1) printf("substep: %f %e %e %e %e\n",dampfac,Wp,f0temp,dfdWptemp,u);
	  dampfac/=2.;

	  if(dampfac<1.e-7) break;
	}
      while(u<0. || isnan(f0temp) || isinf(f0temp) || isnan(dfdWptemp) || isinf(dfdWptemp));

      if(dampfac<1.e-7) {printf("damping unsuccesful at hotmax\n");return -1; getchar(); }
      else  if(verbose>1) {printf("damping succesful\n"); }

      if(verbose>1) 
	{
	  W=Wp+D;
	  v2 = Qt2/W/W;
	  gamma2 = 1./(1.-v2);
	  gamma = sqrt(gamma2);
	  w = W/gamma2;
	  rho0 = D/gamma;
	  wmrho0 = w - rho0;
	  u = wmrho0 / GAMMA;

	  printf("%d %e %e %e -> u %e\n",iter,Wp,f0,dfdWp,u);
	}

      //      if(fabs(dfdWp)<SMALL ) {if(verbose) printf("derivative zero. asssuming found solution\n"); return -1;getchar(); break;}
      if(fabs(Wp)>BIG) {printf("Wp in hotmax has gone out of bounds\n"); return -1;}

      if(isnan(Wp) || isinf(Wp)) {printf("nan/inf Wp in hotmax: %e %e %e %e\n",Wp,f0,dfdWp,W);  return -1;}
    }
  //  while(( fabs((Wp-Wpprev)/Wpprev)>CONV || u<0. || rho<0.) && iter<50);
  while(( fabs((Wp-Wpprev)/Wpprev)>CONV) && iter<50);

  if(iter>=50)
    {
      printf("iter exceeded in u2p_hotmax at %d %d %d\n",geom->ix,geom->iy,geom->iz);
      //getchar();
      return -1;
    }
  

  if(verbose>1) {printf("the end: %e\n",Wp); }//getchar();}

  //Wp found, let's calculate W, v2 and the rest
  W=Wp+D;
  v2=Qt2/W/W;
  gamma2=1./(1.-v2);
  ut2 = Qt2/(W*W - Qt2);
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
      if(verbose>0 || 1) printf("neg u rho in u2p_hotmax %e %e %e %e\n",rho,u,gamma2,W);
      //getchar();
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
  pp[5]=calc_Sfromu(rho,u);

#ifdef TRACER
  ldouble Dtr=uu[TRA]/gdetu*alpha; //uu[0]=gdetu rho ut
  pp[TRA]=Dtr/gamma;
#endif

  ldouble uu2[NV];
  int iv;
  int lostprecision=0;
  p2u(pp,uu2,ggg);
  for(iv=0;iv<NVMHD;iv++)
    {
      if(iv==1) continue;
      //      if(fabs(uu2[iv]-uu[iv])/fabs(uu[iv]+uu2[iv])>1.e-3 && fabs(uu[iv])>SMALL) lostprecision;
    }
     

  if(superverbose)
    {
      printf("superverbose or u2p_hotmax lost precision:\n");
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);
      print_Nvector(uu2,NV);  
      getchar();
    }

  if(verbose>1 && !superverbose) {print_Nvector(pp,NV); }//getchar();}

  return 0;

}

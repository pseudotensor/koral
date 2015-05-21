//KORAL - u2p.c
//conserved to primitives conversion

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates primitives in given cell basing on global array u[]
//type: 0 - regular
//1 - extra check for decreasing entropy
int
calc_primitives(int ix,int iy,int iz,int type,int setflags)
{
  int verbose=0;
  int iv,u2pret,u2pretav;
  ldouble uu[NV],uuav[NV],pp[NV],ppav[NV];
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

  int corrected[3]={0,0,0}, fixups[2]={0,0};

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }


  //aux
  if(setflags)
    {
      set_cflag(ENTROPYFLAG,ix,iy,iz,0); 
      set_cflag(ENTROPYFLAG2,ix,iy,iz,0); 
    }
 
  //converting to primitives

  if(is_cell_corrected_polaraxis(ix,iy,iz))
    {
      u2p_solver_Bonly(uu,pp,&geom); //invert only the magnetic field, the rest will be overwritten
    }
  else
    //    u2p(uu,pp,&geom,corrected,fixups,type); //regular inversion
    u2p(uu,pp,&geom,corrected,fixups,0); //regular inversion

  if(corrected[0]==1 && setflags) //hd correction - entropy solver
    set_cflag(ENTROPYFLAG,ix,iy,iz,1); 
  if(corrected[2]==1 && setflags) //borrowing energy from radiation didn't work
    set_cflag(ENTROPYFLAG2,ix,iy,iz,1); 

  //imposing floors

  //************************************
  //************************************
  //checking on hd floors  
  int floorret=0;

  if(is_cell_active(ix,iy,iz) && !is_cell_corrected_polaraxis(ix,iy,iz))
    floorret=check_floors_mhd(pp,VELPRIM,&geom);

  if(floorret<0.)
    {
      corrected[0]=1;
    }
  //************************************
  //************************************
  //checking on rad floors
    
  #ifdef RADIATION
  floorret=0;
  if(is_cell_active(ix,iy,iz) &&  !is_cell_corrected_polaraxis(ix,iy,iz))
    floorret=check_floors_rad(pp,VELPRIMRAD,&geom);
  if(floorret<0.)
    {
      corrected[1]=1;
    }
  #endif

  //************************************
  //************************************
  
  //************************************
  //update conserved to follow corrections on primitives

if(!is_cell_corrected_polaraxis(ix,iy,iz))
  p2u(pp,uu,&geom);

  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz,uu[iv]);
      set_u(p,iv,ix,iy,iz,pp[iv]);
    }

  //sets the flags for fixups of unsuccessful cells
  if(setflags)
    {
      if(fixups[0]>0)
	{
	  set_cflag(HDFIXUPFLAG,ix,iy,iz,1); 
	  global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++;      
	}
      else
	set_cflag(HDFIXUPFLAG,ix,iy,iz,0); 

  
      if(fixups[1]>0)
	{
	  set_cflag(RADFIXUPFLAG,ix,iy,iz,-1); 
	  global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]++;      
	}
      else
	set_cflag(RADFIXUPFLAG,ix,iy,iz,0); 
    }
  

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
  int corrected[3], fixups[2];
  u2p(uu,pp,&geom,corrected,fixups,0);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//high-level u2p solver
int
u2p(ldouble *uu0, ldouble *pp,void *ggg,int corrected[3],int fixups[2],int type)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble uu[NV];
  int iv;
  PLOOP(iv) uu[iv]=uu0[iv];

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int verbose=0;
  int hdcorr=0;
  int radcor=0;
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
  
  //test
  ldouble ppold[NV];
  PLOOP(iv)
    ppold[iv]=pp[iv];

  //negative uu[0] = rho u^t
  if(uu[0]/gdetu<0.)
    {
       int gix,giy,giz;
       mpi_local2globalidx(geom->ix,geom->iy,geom->iz,&gix,&giy,&giz);
       printf("%4d > %4d %4d %4d > NEGUU  > neg uu[0] - requesting fixup\n",PROCID,gix,giy,giz);
      pp[0]=RHOFLOOR; //used when not fixin up
      ret=-2; //to request fixup
      u2pret=0; //to skip solvers
    }
  else
    {
      int method=U2P_HOT;


#ifdef ENFORCEENTROPY
      u2pret=-1; //skip hot
#else

    

      u2pret=u2p_solver(uu,pp,ggg,method,0); 
  
      //check if u2p_hot went mad by making entropy decrease
      //this check performed only when type==1
      //which is used after the advection operator
      if(u2pret==0 && type==1)
	{
	  //by comparing the Lagrangian uu[ENTR] value and the one from u2p_hot
	  ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};
	  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
	  ldouble s1=exp(uu[ENTR]/gdetu/ucon[0]/pp[RHO]);
	  ldouble s2=exp(pp[ENTR]/pp[RHO]);
	  
	  //TODO: double check on this - why in ORSZAG and RADTUBE 0.999 often exceeded?
	  if(s2/s1 < 0.1)
	    {  
	      //go to entropy
	      if(verbose || 1) printf("enforcing entr at %4d %4d %4d\n",geom->ix,geom->iy,geom->iz);
	      u2pret=-1;
	    }
	}

#endif
    }
 
  if(ALLOWENTROPYU2P)
    if(u2pret<0)
      {
	ret=-1;
	
	if(verbose>2)
	  {
	    printf("u2p_entr     >>> %d %d <<< %d >>> %e > %e\n",geom->ix + TOI, geom->iy + TOJ,u2pret,u0,pp[1]);
	  }

	//************************************
	//************************************
	//************************************
	//entropy solver - conserving entropy
	u2pret=u2p_solver(uu,pp,ggg,U2P_ENTROPY,0);  

	

	if(u2pret<0)
	  {
	    ret=-2;
	    
	    if(verbose>1)
	      {
		printf("u2p_entr err No. %4d > %e %e %e > %e %e > %4d %4d %4d\n",u2pret,uu[0],uu[1],uu[5],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
	      }

	    //mark: u2perror
	    //test, to print it if error and stop
	    //u2pret=u2p_solver(uu,pp,ggg,U2P_HOT,2);  
	    //getch();
	    //u2pret=u2p_solver(uu,pp,ggg,U2P_ENTROPY,2);  
	    //getchar();

	    /*
	    if(u2pret==-103 && 0)  //TODO: work out hotmax
	      //solver went rho->D meaning entropy too large 
	      //imposing URHOLIMIT 
	      {		
		//u2pret=u2p_hotmax(uu,pp,ggg);
		u2pret=u2p_solver(uu,pp,ggg,U2P_HOTMAX,0);
		if(u2pret<0)
		  {
		    if(verbose>0)
		      {
			printf("u2p_hotmax err No. %4d > %e %e %e > %e %e > %4d %4d %4d\n",u2pret,uu[0],uu[1],uu[5],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
		      }
			
		    //should not happen but if happens use the old state to impose URHOLIMIT
		    pp[1]=UURHORATIOMAX*pp[0];
		    check_floors_mhd(pp,VELPRIM,&geom);
		    pp[5]=calc_Sfromu(pp[0],pp[1]);
		    p2u(pp,uu,&geom);	

		    //no need for another entropy solver - p2u does its job
		    u2pret=0;
		  }		    
	      }
	    */
	  }	
      }

  /*
  if(ALLOWCOLDU2P)
    if(u2pret<0.)
      {
	//cold RHD - assuming u=SMALL
	ret=-2;
	u2pret=u2p_solver(uu,pp,ggg,U2P_COLD,0);

	if(u2pret<0)
	  if(verbose>0)
	    {
	      printf("u2p_cold err > %e %e > %e %e > %4d %4d %4d\n",uu[0],uu[1],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
	    }
      }
  */

  if(u2pret<0)
    {
      //************************************
      //leaving unchanged primitives - should not happen
      if(verbose>1 || 1)
	{
	  printf("%4d > %4d %4d %4d > MHDU2PFAIL > u2p prim. unchanged > %d \n",PROCID,geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,u2pret);
	}
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

corrected[2]=0;
  
#ifdef RADIATION
  #ifdef BALANCEENTROPYWITHRADIATION
 //trying to balance gain of energy because of entropy inversion
  //by borrowing from the radiation field

 if(ret==-1) //entropy u2p was used
   {            
     ldouble uunew[NV],ppnew[NV];
     PLOOP(iv) {uunew[iv]=uu[iv]; ppnew[iv]=pp[iv];}
     p2u_mhd(pp,uunew,geom); 
     ldouble dugas = uunew[UU]-uu[UU]; //this much energy was introduced
     if(fabs(dugas)<0.1*fabs(uunew[EE0]))//correction relatively small - is this general enough?
       {
	 uunew[EE0]-=dugas; //balancing with radiation      
	 u2p_rad(uunew,ppnew,geom,&radcor);
       }
     else
       radcor=1;
   
	
     if(radcor==0) //there was enough energy to borrow from and uunew inverts with hot
       {
	 PLOOP(iv) 
	   uu[iv]=uunew[iv];
	 //printf("entropy correction did work at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
       }
     else
       {
	 corrected[2]=1; //entropy correction didn't work
	 //printf("entropy correction didn't work at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
       }
   }
  
#endif

 u2p_rad(uu,pp,geom,&radcor);

#ifdef BALANCERADCORRWITHGAS
 if(radcor!=0) //some type of radiative correction applied
   {
     ldouble uunew[NV],ppnew[NV];
     PLOOP(iv) 
     {uunew[iv]=uu[iv];ppnew[iv]=pp[iv];}
     p2u(pp,uunew,geom);
     ldouble durad = uunew[EE0]-uu[EE0];
     uunew[UU]-=durad;
     u2pret=u2p_solver(uunew,ppnew,geom,U2P_HOT,0); 
      
      if(u2pret==0) //there was enough energy to borrow from
	{
	  PLOOP(iv) 
	  {
	    uu[iv]=uunew[iv];
	    pp[iv]=ppnew[iv];
	  }
	  printf("radcorr correction did work at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
	}
      else
	{
	  corrected[2]=1; //entropy correction didn't work
	  printf("radcorr correction didn't work at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
	}
    }

  #endif


#endif
  
  //************************************
  //************************************

  //commented out - rad fixups only for critical failure in implicit
#if (DORADFIXUPS==1)
  if(radcor>0)     
    fixups[1]=1;
  else
    fixups[1]=0;
#endif

  //************************************
  //************************************
  //************************************

  if(hdcorr>0) corrected[0]=1;
  if(radcor>0) corrected[1]=1;

  return ret;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//checks if hydro primitives make sense
int
check_floors_mhd(ldouble *pp, int whichvel,void *ggg)
{
  //return 0;

  int verbose=0;
  int ret=0;
  int iv;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble uu[NV];
  p2u_mhd(pp,uu,ggg);
 
#ifdef TRACER
  if(pp[TRA]<0.) {pp[TRA]=0.; ret=-1; if(verbose) printf("hd_floors CASE TRA 1\n");}
#endif

  //**********************************************************************
  //rho too small
  if(pp[0]<RHOFLOOR) 
    {
      if(verbose) printf("hd_floors CASE 1 at %d %d %d | %d %d %d (%e) | tijk: %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,geom->ix,geom->iy,geom->iz,pp[0],TI,TJ,TK);
      pp[0]=RHOFLOOR; 
      ret=-1; 
    }

#ifdef VXFLOOR
  if(fabs(pp[VX])<1.e-10) 
    { 
      struct geometry geomBL;
      fill_geometry_arb(geom->ix,geom->iy,geom->iz,&geomBL,BLCOORDS);
      ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};
      conv_vels(ucon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);
      trans2_coco(geom->xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
      if(fabs(ucon[1])<VXFLOOR)
	{
	  ucon[1]=my_sign(ucon[1])*VXFLOOR;
	}
      trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
      conv_vels(ucon,ucon,VEL4,VELPRIM,geom->gg,geom->GG);
     
      pp[VX]=ucon[1];
       
      ret=-1;
    }
#endif

  //**********************************************************************
  //rho too small, BH-disk like
#ifdef RHOFLOOR_BH
  ldouble xxBL[4];
  coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
  ldouble rr = xxBL[1] / rhorizonBL;
  ldouble rhofloor = RHOFLOOR_BH_NORM / sqrt(rr*rr*rr);
  if(pp[0]<rhofloor) 
    {
      if(verbose) printf("hd_floors BH CASE 1 at %d %d (%e)\n",geom->ix+TOI,geom->iy+TOJ,pp[0]);
      pp[0]=rhofloor;
      ret=-1; 
    }
#endif

  //**********************************************************************
  //too cold
  if(pp[1]<UURHORATIOMIN*pp[0]) 
    {
      if(verbose) {printf("hd_floors CASE 2 at (%d,%d,%d | %d,%d,%d): %e %e | tijk: %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,geom->ix,geom->iy,geom->iz,pp[0],pp[1],TI,TJ,TK);}//getchar();}
      pp[1]=UURHORATIOMIN*pp[0]; //increasing uint
      ret=-1;
    }

  //**********************************************************************
  //to hot
  if(pp[1]>UURHORATIOMAX*pp[0]) 
    {
      pp[1]=UURHORATIOMAX*pp[0]; //decreasing uint
      //pp[0]=pp[1]/UURHORATIOMAX; //increasing rho
      ret=-1;      
      if(verbose) printf("hd_floors CASE 3 at (%d,%d,%d): %e %e\n",geom->ix+TOI,geom->iy+TOJ,geom->iz,pp[0],pp[1]);
    }

  
  //**********************************************************************
  //too magnetized
#ifdef MAGNFIELD
  ldouble ucond[4],ucovd[4];
  ldouble bcond[4],bcovd[4],magpre;  
  ldouble etacon[4],etarel[4];
  for(iv=1;iv<4;iv++)
    ucond[iv]=pp[1+iv];
  conv_vels(ucond,ucond,VELPRIM,VEL4,gg,GG);
  indices_21(ucond,ucovd,gg);
  calc_bcon_4vel(pp,ucond,ucovd,bcond);
  indices_21(bcond,bcovd,gg); 
  magpre = dot(bcond,bcovd)/2.;
  calc_normalobs_4vel(GG,etacon);
  conv_vels_ut(etacon,etarel,VEL4,VELPRIM,gg,GG);

  if(magpre>B2RHORATIOMAX*pp[RHO]) 
    {
      if(verbose) printf("mag_floors CASE 2 at (%d,%d,%d): %e %e\n",geom->ix+TOI,geom->iy+TOJ,geom->iz,pp[RHO],magpre);
      ldouble f=magpre/(B2RHORATIOMAX*pp[RHO]);

#if (B2RHOFLOORFRAME==ZAMOFRAME) //new mass in ZAMO

      ldouble dpp[NV],duu[NV];
      ldouble drho=pp[RHO]*(f-1.);
    
      for(iv=0;iv<NVMHD;iv++)
	dpp[iv]=0.0;

      dpp[RHO]=drho;
      dpp[UU] = pp[UU]*drho/pp[RHO];
      dpp[VX] = etarel[1];
      dpp[VY] = etarel[2];
      dpp[VZ] = etarel[3];
      dpp[ENTR] = 0.;
      dpp[B1] = dpp[B2] = dpp[B3] = 0.;

      p2u_mhd(dpp,duu,geom);
 
      for(iv=0;iv<NVMHD;iv++)
	uu[iv]+=duu[iv];

      int rettemp=0;
      rettemp=u2p_solver(uu,pp,geom,U2P_HOT,0); 
      if(rettemp<0)
	rettemp=u2p_solver(uu,pp,geom,U2P_ENTROPY,0); 
      
      if(rettemp<0) 
	{
#ifdef BHDISK_PROBLEMTYPE
	  if(geom->ix+TOI>5) //report only outside horizon
#endif
	    printf("u2p failed after imposing bsq over rho floors at %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK);
	  
	}
    
#elif(B2RHOFLOORFRAME==FFFRAME) //new mass in fluid frame

      pp[RHO]*=f;
      pp[UU]*=f;

#endif
      
      ret=-1;      
    }
  
  /*
  //independent check on ugas
  if(magpre>B2UURATIOMAX*pp[UU]) 
    {
      if(verbose) printf("mag_floors CASE 3 at (%d,%d,%d): %e %e\n",geom->ix+TOI,geom->iy+TOJ,geom->iz,pp[UU],magpre);
      pp[UU]*=magpre/(B2UURATIOMAX*pp[UU]);
      ret=-1;      
    }
  */

#endif

  //**********************************************************************
  //too fast
  
  if(VELPRIM==VELR) 
    {
      ldouble qsq=0.;
      int i,j;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=pp[UU+i]*pp[UU+j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      if(gamma2>GAMMAMAXHD*GAMMAMAXHD)
	{
	  ldouble qsqmax=GAMMAMAXHD*GAMMAMAXHD-1.;
	  ldouble A=sqrt(qsqmax/qsq);
	  for(j=1;j<4;j++)
	    pp[UU+j]*=A;
	  ret=-1;
	  if(verbose)
	    {
	      printf("hd_floors CASE 4 at (%d,%d,%d): %e",geom->ix+TOI,geom->iy+TOJ,geom->iz,sqrt(gamma2));
	      qsq=0.;
	      for(i=1;i<4;i++)
		for(j=1;j<4;j++)
		  qsq+=pp[UU+i]*pp[UU+j]*gg[i][j];
	      gamma2=1.+qsq;
	      printf(" -> %e\n",sqrt(gamma2));
	    }
	}
    }
  //TODO: implement checks for other VELPRIM


  //updates entropy after floor corrections
  if(ret<0)
    pp[5]=calc_Sfromu(pp[RHO],pp[UU]);

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
f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df,ldouble *err)
{

  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Qdotnp=cons[6];
  
  ldouble W=Wp+D;

  FTYPE W3,X3,Ssq,Wsq,X,X2,Xsq; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X2 = X*X;
  Xsq = X2;
  X3 = X2*X;
  //  return -(Qn+W)*(GAMMA/GAMMAM1)+W*(1.-Qt2/W/W)-D*sqrt(1.-Qt2/W/W);   

  //a bit more clear

  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
  ldouble gamma2 = 1./(1.-v2);

  //test
  /*
  ldouble S = sqrt(QdotBsq);
  ldouble SoW = S/W;
  ldouble utsqtop = Qtsq + (Bsq+2.0*W)*SoW*SoW;
  ldouble utsqbottom = (Bsq*Bsq-Qtsq) + 2.0*Bsq*W + W*W - SoW*SoW*(Bsq+2.0*W);
  ldouble utsq = utsqtop/utsqbottom;
  gamma2=1.0+utsq;
  */

  ldouble gamma = sqrt(gamma2);
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);
  ldouble u = wmrho0 / GAMMA;
  ldouble p = (GAMMA-1)*u;

  //original:
#if (U2P_EQS==U2P_EQS_NOBLE)
  *f = Qn + W - p + 0.5*Bsq*(1.+v2) - QdotBsq/2./Wsq;
 *err = fabs(*f) / (fabs(Qn) + fabs(W) + fabs(p) + fabs(0.5*Bsq*(1.+v2)) + fabs(QdotBsq/2./Wsq));
#endif

  //JONS:
#if (U2P_EQS==U2P_EQS_JON)
 *f = Qdotnp + Wp - p + 0.5*Bsq + (Bsq*Qtsq - QdotBsq)/X2;
 *err = fabs(*f) / (fabs(Qdotnp) + fabs(Wp) + fabs(p) + fabs(0.5*Bsq) + fabs((Bsq*Qtsq - QdotBsq)/X2));
#endif

  // dp/dWp = dp/dW + dP/dv^2 dv^2/dW
    
  ldouble dvsq=(-2.0/X3 * ( Qtsq  +  QdotBsq * (3.0*W*X + Bsq*Bsq)/W3));
  ldouble dp1 = dpdWp_calc_vsq(Wp, D, v2 ); // vsq can be unphysical

  ldouble idwmrho0dp=compute_idwmrho0dp(wmrho0);
  ldouble dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp);

  ldouble drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma
  ldouble idrho0dp = compute_idrho0dp(wmrho0);

  ldouble dp2 =   drho0dvsq *idrho0dp  +   dwmrho0dvsq *idwmrho0dp;

  ldouble dpdW = dp1  + dp2*dvsq; // dp/dW = dp/dWp

  //original:
  #if (U2P_EQS==U2P_EQS_NOBLE)
  *df=1.-dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;
  #endif

  //JONs:
  #if (U2P_EQS==U2P_EQS_JON)
  *df=1. -dpdW + (Bsq*Qtsq - QdotBsq)/X3*(-2.0);
  #endif

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
f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df, ldouble *err)
{
  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Sc=cons[5];
 
  ldouble W=Wp+D;

  FTYPE W3,X3,Ssq,Wsq,X,X2,Xsq; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X2 = X*X;
  Xsq = X2;
  X3 = X2*X;

  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
  ldouble gamma2 = 1./(1.-v2);
  ldouble gamma = sqrt(gamma2);
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);
  ldouble u = wmrho0 / GAMMA;
  ldouble p = (GAMMA-1)*u;

  ldouble Ssofchi=compute_specificentropy_wmrho0_idealgas(rho0,wmrho0);

  *f= -Sc + D*Ssofchi;

  *err = fabs(*f) / (fabs(Sc) + fabs(D*Ssofchi));

  //printf("%e %e %e %e\n",*err,*f, Sc , D*Ssofchi);

  FTYPE dSsdW,dSsdvsq,dSsdWp,dScprimedWp,dSsdrho,dSsdchi;
  FTYPE dvsq,dwmrho0dW,drho0dW;
  FTYPE dwmrho0dvsq,drho0dvsq;

  dSsdrho=compute_dspecificSdrho_wmrho0_idealgas(rho0,wmrho0);
  dSsdchi=compute_dspecificSdwmrho0_wmrho0_idealgas(rho0,wmrho0);

  dwmrho0dW = 1.0/gamma; // holding utsq fixed
  drho0dW = 0.0; // because \rho=D/\gamma and holding utsq fixed
  dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp); // holding Wp fixed
  drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma and holding Wp fixed

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
f_u2p_cold(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df, ldouble *err)
{
  //my_err("Think f_u2p_cold over in terms of mhd\n");

  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Sc=cons[5];
 
  ldouble W=Wp+D;

  FTYPE W3,X3,Ssq,Wsq,X,X2,Xsq; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X2 = X*X;
  Xsq = X2;
  X3 = X2*X;
 
  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
   ldouble gamma2 = 1./(1.-v2);
  ldouble gammasq = gamma2;
  ldouble gamma = sqrt(gamma2);
  ldouble w = W/gamma2;
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0;
  ldouble u = wmrho0 / GAMMA;
  ldouble p = (GAMMA-1)*u;

  *f = u - UURHORATIOU2PCOLD*rho0;

  *err = fabs(*f) / (fabs(u) + fabs(UURHORATIOU2PCOLD*rho0));

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
f_u2p_hotmax(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df, ldouble *err)
{
  my_err("Think f_u2p_hotmax over in terms of mhd\n");

  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Sc=cons[5];

  ldouble W=Wp+D;

  FTYPE W3,X3,Ssq,Wsq,X,X2,Xsq; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X2 = X*X;
  Xsq = X2;
  X3 = X2*X; 

  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
  ldouble gamma2 = 1./(1.-v2);
  ldouble gammasq = gamma2;
  ldouble gamma = sqrt(gamma2);
  ldouble w = W/gamma2;
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = w - rho0;
  ldouble u = wmrho0 / GAMMA;
  ldouble p = (GAMMA-1)*u;

  *f = u - UURHORATIOMAX*rho0;

  *err = fabs(*f) / (fabs(u) + fabs(UURHORATIOMAX*rho0));

  ldouble dWdWp=1.;
  ldouble dv2dWp=Qt2*(-2.)/W/W/W*dWdWp;
  ldouble dgammadWp=.5*1./sqrt(1./(1.-v2))*(-1.)/(1.-v2)/(1.-v2)*(-1.)*dv2dWp;
  ldouble drho0dWp=-D/gamma2*dgammadWp;
  ldouble dwdWp = dWdWp/gamma2 - 2. *W/gamma/gamma/gamma*dgammadWp;
  ldouble dudWp = 1./GAMMA*(dwdWp - drho0dWp);

  *df = dudWp - UURHORATIOMAX*drho0dWp;
  
  return 0;
 
}

double
fWplim (double Wp, void *params)
{
  ldouble *cons
    = (double *) params;
  
  ldouble Qn=cons[0];
  ldouble Qtsq=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Qdotnp=cons[6];
  
  return Wp*(1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
       (Power(D + Wp,2)*Power(Bsq + D + Wp,2))) - 
   (D*(Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp))))/
    (Power(D + Wp,2)*Power(Bsq + D + Wp,2)*
      (1 + Sqrt(1/
          (1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
	   (Power(D + Wp,2)*Power(Bsq + D + Wp,2))))));
}

double
fWplim_deriv (double Wp, void *params)
{
  ldouble *cons
    = (double *) params;
  
  ldouble Qn=cons[0];
  ldouble Qtsq=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Qdotnp=cons[6];

  return 1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
    (Power(D + Wp,2)*Power(Bsq + D + Wp,2)) + 
   (2*Wp*(Power(Bsq,2)*QdotBsq + 3*Bsq*QdotBsq*(D + Wp) + 
        Power(D + Wp,2)*(3*QdotBsq + Qtsq*(D + Wp))))/
    (Power(D + Wp,3)*Power(Bsq + D + Wp,3)) - 
   (D*(Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))*
      (Power(Bsq,2)*QdotBsq + 3*Bsq*QdotBsq*(D + Wp) + 
        Power(D + Wp,2)*(3*QdotBsq + Qtsq*(D + Wp)))*
      Power(1/(1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
           (Power(D + Wp,2)*Power(Bsq + D + Wp,2))),1.5))/
    (Power(D + Wp,5)*Power(Bsq + D + Wp,5)*
      Power(1 + Sqrt(1/
          (1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
             (Power(D + Wp,2)*Power(Bsq + D + Wp,2)))),2)) - 
   (2*D*(QdotBsq + Qtsq*(D + Wp)))/
    (Power(D + Wp,2)*Power(Bsq + D + Wp,2)*
      (1 + Sqrt(1/
          (1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
             (Power(D + Wp,2)*Power(Bsq + D + Wp,2)))))) + 
   (2*D*(Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp))))/
    (Power(D + Wp,2)*Power(Bsq + D + Wp,3)*
      (1 + Sqrt(1/
          (1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
             (Power(D + Wp,2)*Power(Bsq + D + Wp,2)))))) + 
   (2*D*(Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp))))/
    (Power(D + Wp,3)*Power(Bsq + D + Wp,2)*
      (1 + Sqrt(1/
          (1 - (Bsq*QdotBsq + (D + Wp)*(2*QdotBsq + Qtsq*(D + Wp)))/
	   (Power(D + Wp,2)*Power(Bsq + D + Wp,2))))));
}

void
fWplim_fdf (double Wp, void *params, 
               double *y, double *dy)
{
  *y = fWplim (Wp,params);
  *dy =  fWplim_deriv (Wp,params);
}


int
find_Wplim(ldouble *Wp,ldouble *cons)
{
  int verbose=0;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x = *Wp, x0;
  gsl_function_fdf FDF;

  FDF.f = &fWplim;
  FDF.df = &fWplim_deriv;
  FDF.fdf = &fWplim_fdf;
  FDF.params = cons;

  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);

  if(verbose) printf ("using %s method\n", 
          gsl_root_fdfsolver_name (s));

  if(verbose) printf ("%-5s %10s %10s\n",
          "iter", "root", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-3);

      if(isnan(x))
	return -1;

      //forbid  Wp<0
      while(x<=0)
	{
	  x=0.5*(x0+x);
	}
	

      if (verbose && status == GSL_SUCCESS)
        printf ("Converged:\n");

      if(verbose) printf ("%5d %10.7e %10.7e\n",
              iter, x, x - x0);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fdfsolver_free (s);
  *Wp = x;

  return 0;
}

//*** wrapper ***/
int
u2p_solver(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  #ifdef NONRELMHD
  return u2p_solver_nonrel(uu,pp,ggg,Etype,verbose);
  #endif

  int (*solver)(ldouble*,ldouble*,void*,int,int);
#if (U2P_SOLVER==U2P_SOLVER_WP)
  solver = & u2p_solver_Wp;
#endif
#if (U2P_SOLVER==U2P_SOLVER_W)
  solver = & u2p_solver_W;
#endif

  return (*solver)(uu,pp,ggg,Etype,verbose);
}
 

//**********************************************************************
//**********************************************************************
//**********************************************************************
//non-relativistic, analytical, solver

int
u2p_solver_nonrel(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
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
  //density
  ldouble rho=uu[RHO]/gdetu;
  pp[RHO]=rho;

  //velocities u_i
  ldouble ucov[4],ucon[4],vcov[4];
  ucov[0]=-1.;
  ucov[1]=uu[VX]/rho/gdetu;
  ucov[2]=uu[VY]/rho/gdetu;
  ucov[3]=uu[VZ]/rho/gdetu;
  fill_utinucov(ucov,gg,GG);

  indices_12(ucov,ucon,GG);

  ldouble v2=dot3nr(ucon,ucov);

  pp[VX]=ucon[1];
  pp[VY]=ucon[2];
  pp[VZ]=ucon[3];

  ldouble bcon[4],bcov[4],bsq=0.;

#ifdef MAGNFIELD
 
  bcon[0]=0.;
  bcon[1]=uu[B1]/gdetu;
  bcon[2]=uu[B2]/gdetu;
  bcon[3]=uu[B3]/gdetu;

  pp[B1]=bcon[1];
  pp[B2]=bcon[2];
  pp[B3]=bcon[3];

  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
#endif

  
  ldouble uint;
  if(Etype==U2P_HOT)
    {
      uint = -uu[UU]/gdetu-bsq/2. - rho*v2/2.;
      if(uint<NONRELMHDENTROPYCUT*rho) return -1;
      pp[UU]=uint;
    }
  else if(Etype==U2P_ENTROPY)
    {
      ldouble S=uu[ENTR]/gdetu;
      uint= calc_ufromS(S,rho);
      pp[UU]=uint;
    }

  //entropy 
  pp[ENTR]=calc_Sfromu(rho,uint);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//Newton-Rapshon solver 
//upgraded - uses Wp
//Etype == 0 -> hot inversion (uses D,Ttt,Tti)
//Etype == 1 -> entropy inversion (uses D,S,Tti)
//Etype == 2 -> hotmax inversion (uses D,Tti,u over rho max ratio)
//Etype == 3 -> cold inversion (uses D,Tti,u over rho min ratio)

int
u2p_solver_Wp(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  int i,j,k;
  ldouble rho,uint,w,W,Wp,Wpprev,alpha,D,Sc,alphasq,betasqoalphasq;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],Qconp[4],Qcovp[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn,Qdotnp;
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
  int (*f_u2p)(ldouble,ldouble*,ldouble*,ldouble*,ldouble*);
 if(Etype==U2P_HOT) 
   f_u2p=&f_u2p_hot;
 if(Etype==U2P_ENTROPY) 
   f_u2p=&f_u2p_entropy;
 if(Etype==U2P_HOTMAX) 
   f_u2p=&f_u2p_hotmax;
 if(Etype==U2P_COLD) 
   f_u2p=&f_u2p_cold;
  /****************************/
 
  
  if(verbose>1) {printf("********************\n");print_conserved(uu);print_primitives(pp);}

  /****************************/
  //conserved quantities etc
  
  //alpha
  alpha=geom->alpha;
  alphasq=alpha*alpha;

  //D
  D=uu[0]/gdetu*alpha; //uu[0]=gdetu rho ut

  //conserved entropy "S u^t"
  Sc=uu[5]/gdetu*alpha; 

  //Qp_mu=alpha T^t_mu 
  Qcovp[0]=uu[1]/gdetu*alpha;
  Qcovp[1]=uu[2]/gdetu*alpha;
  Qcovp[2]=uu[3]/gdetu*alpha;
  Qcovp[3]=uu[4]/gdetu*alpha;

  //Qp^mu
  indices_12(Qcovp,Qconp,GG);

  //Q_mu=alpha (T^t_mu - rho u^t delta(t,mu)) - avoid this one
  Qcov[0]=(uu[1]/gdetu-uu[0]/gdetu)*alpha;
  Qcov[1]=uu[2]/gdetu*alpha;
  Qcov[2]=uu[3]/gdetu*alpha;
  Qcov[3]=uu[4]/gdetu*alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

#ifdef MAGNFIELD
  //curly B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1]/gdetu*alpha;
  Bcon[2]=uu[B2]/gdetu*alpha;
  Bcon[3]=uu[B3]/gdetu*alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
#else
  Bsq=QdotB=QdotBsq=0.;
  Bcon[0]=Bcon[1]=Bcon[2]=Bcon[3]=0.;
#endif  

  //n_mu = (-alpha, 0, 0, 0)
  ncov[0]=-alpha;
  ncov[1]=ncov[2]=ncov[3]=0.;
  
  //n^mu
  indices_12(ncov,ncon,GG);

  //Q_mu n^mu = Q^mu n_mu = -alpha*Q^t
  Qn=Qcon[0] * ncov[0];

  //\beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
  betasqoalphasq=gg[0][1]*GG[0][1] + gg[0][2]*GG[0][2] + gg[0][3]*GG[0][3]; 
  
  //Qdotnp=-E'=-E+D
  ldouble Dfactor = (-geom->gttpert + alphasq*betasqoalphasq)/(alphasq+alpha);
  Qdotnp = Qconp[0]*ncov[0] + D*(Dfactor) ; // -Qdotn-W = -Qdotnp-Wp

 

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
  
  //initial guess for Wp = w gamma**2 based on primitives
  rho=pp[0];
  uint=pp[1];
  utcon[0]=0.;
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  //conv_vels(utcon,ucon,VELPRIM,VEL4,gg,GG);
  //indices_21(ucon,ucov,gg);
  conv_vels(utcon,utcon,VELPRIM,VELR,gg,GG);

  ldouble qsq=0.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=utcon[i]*utcon[j]*gg[i][j];
  ldouble gamma2=1.+qsq;
  ldouble gamma=sqrt(gamma2);

  //Wp
  Wp=(GAMMA*uint)*gamma2;

  ldouble Wpinit, Winit;
  if(verbose>1) printf("initial Wp:%e\n",Wp);
  Wpinit=Wp;
  Winit=Wpinit+D;

  /****************************/
  
  //test if does not provide reasonable gamma2
  // Make sure that W is large enough so that v^2 < 1 and w-rho > 0 : 
  int i_increase = 0;
  ldouble f0,f1,dfdW,err;
  ldouble CONV=U2PCONV; 
  ldouble cons[7]={Qn,Qt2,D,QdotBsq,Bsq,Sc,Qdotnp};
  int iter=0,fu2pret;
  
  /*
  if(find_Wplim(&Wp,cons)!=0) //does not modify if error
    {
      if(verbose) printf("err in find_Wplim\n");
    }
  */

  do
    {
      W=Wp+D;
      f0=dfdW=0.;

      FTYPE Wsq,Xsq,X; 
      X = Bsq + W;
      Xsq = X*X;
      Wsq = W*W;

      ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
      ldouble gamma2 = 1./(1.-v2);
      ldouble gamma = sqrt(gamma2);
      ldouble rho0 = D/gamma;
      ldouble wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);


      //if(Etype!=U2P_HOT) 
      (*f_u2p)(Wp,cons,&f0,&dfdW,&err);
      
      if((gamma2<0. || Wp<0. || wmrho0<0.|| !isfinite(f0) || !isfinite(dfdW)) && (i_increase < 50))
	//if((gamma2<0. || Wp<0. || wmrho0<0.) && (i_increase < 50))
	{
	  if(verbose>0) printf("init Wp : %e - %e %e %e %e\n",Wp,v2,wmrho0,f0,dfdW);
	  Wp *= 2.;
	  i_increase++;
	  continue;
	}
      else
	break;    
    }
  while(1);

  

  if(i_increase>=50)
    {
      if(verbose>0) 
	{printf("failed to find initial W for Etype: %d\n",Etype);
	  printf("at %d %d\n",geom->ix+TOI,geom->iy+TOJ);}
      return -150;

      print_NVvector(uu);
      print_NVvector(pp);
      getchar();
    }



  //1d Newton solver
 
  do
    {
      Wpprev=Wp;
      iter++;
     
      fu2pret=(*f_u2p)(Wp,cons,&f0,&dfdW,&err);

      //numerical derivative
      //      double EPS=1.e-8;
      //fu2pret=(*f_u2p)((1.+EPS)*W-D,cons,&f1,&dfdW,&err);
      //dfdW=(f1-f0)/(EPS*W);

      if(verbose>1) printf("%d %e %e %e %e\n",iter,Wp,f0,dfdW,err);
 
      //convergence test
      //if(err<CONV)
      //break;
      
      if(dfdW==0.) {Wp*=1.1; continue;}

      ldouble Wpnew=Wp-f0/dfdW;
      
      //not to jump zero over
      //test
      Wpnew = my_max(Wpnew, Wp/100.);

      ldouble Wnew=Wpnew+D;
      int idump=0;
      ldouble dumpfac=1.;

      //test if goes out of bounds and damp step if so
      int itmaxdamp=50;
      do
	{
	  ldouble f0tmp,dfdWtmp,errtmp;
	  Wnew=Wpnew+D;
	  f0tmp=dfdWtmp=0.;

	  FTYPE Wsq,Xsq,X; 
	  X = Bsq + Wnew;
	  Xsq = X*X;
	  Wsq = Wnew*Wnew;

	  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*Wnew)) / (Wsq*Xsq);
	  ldouble gamma2 = 1./(1.-v2);
	  ldouble gamma = sqrt(gamma2);
	  ldouble rho0 = D/gamma;
	  ldouble wmrho0 = Wpnew/gamma2 - D*v2/(1.+gamma);

	  //if(Etype!=U2P_HOT) 
	  //(*f_u2p)(Wpnew,cons,&f0tmp,&dfdWtmp,&errtmp);

	  if(verbose>1) printf("sub (%d) :%d %e %e %e %e %e %e\n",idump,iter,Wpnew,f0tmp,dfdWtmp,v2,gamma2,wmrho0);

	  //if((gamma2<0. || Wpnew<0. || wmrho0<0. || !isfinite(f0tmp) || !isfinite(dfdWtmp)) && (idump<itmaxdamp))
	  if((gamma2<0. || Wpnew<0. || wmrho0<=0.) && (idump<itmaxdamp))
	    {
	      idump++;
	      dumpfac/=2.;
	      Wpnew=Wp-dumpfac*f0/dfdW;
	      continue;
	    }
	  else
	    break;
	}
      while(1);



      if(idump>=itmaxdamp) 
	{
	  if(verbose>0) printf("damped unsuccessfuly\n");
	  return -101;
	}
	
      Wp=Wpnew; 

      if(fabs(W)>BIG) 
	{
	  if(verbose>1) printf("W has gone out of bounds at %d,%d,%d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz); 
	  return -103;
	}


      //convergence test:
      if(err<CONV || (Etype==U2P_HOT && fabs((Wp-Wpprev)/Wpprev)<CONV && err<(sqrt(CONV)))
	 || (Etype==U2P_ENTROPY && fabs((Wp-Wpprev)/Wpprev)<CONV && err<0.99)) break;
      //if(err<CONV || fabs((Wp-Wpprev)/Winit)<CONV) break;
    }
  while(iter<50);
   
  if(iter>=50)
    {
      if(verbose>0) printf("iter exceeded in u2p_solver with Etype: %d\n",Etype); //getchar();
      return -102;
    }


  if(!isfinite(Wp) || !isfinite(Wp)) {if(verbose) printf("nan/inf W in u2p_solver with Etype: %d\n",Etype); return -103;}
 
  if(verbose>1) 
    {
      fu2pret=(*f_u2p)(Wp,cons,&f0,&dfdW,&err);
      printf("end: %d %e %e %e %e\n",iter,Wp,f0,dfdW,err);
    }

  //W found, let's calculate v2 and the rest
  //ldouble v2=Qt2/W/W;

  W=Wp+D;
  ldouble Wsq,Xsq,v2,wmrho0;
	
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);  
  v2 = ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);

  gamma2=1./(1.-v2);
  gamma=sqrt(gamma2);
  rho=D/gamma;
  // w-\rho_0 = (u+p) = W'/\gamma^2 - D v^2/(1+\gamma)
  wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);
  uint=1./GAMMA*wmrho0;
  utcon[0]=0.;
  utcon[1]=gamma/(W+Bsq)*(Qtcon[1]+QdotB*Bcon[1]/W);
  utcon[2]=gamma/(W+Bsq)*(Qtcon[2]+QdotB*Bcon[2]/W);
  utcon[3]=gamma/(W+Bsq)*(Qtcon[3]+QdotB*Bcon[3]/W);

  if(!isfinite(utcon[1]))
    {
      //print_4vector(utcon);
      return -120;
    }


  if(uint<0. || gamma2<0. ||isnan(Wp) || !isfinite(Wp)) 
    {
      if(verbose>0) printf("neg u in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);//getchar();
      return -104;
    }

  //converting to VELPRIM
  conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  
  //returning new primitives
  pp[RHO]=rho;
  pp[UU]=uint;
  pp[VX]=utcon[1];
  pp[VY]=utcon[2];
  pp[VZ]=utcon[3];


  if(rho<0.) 
    {
      if(verbose>0) printf("neg rho in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);//getchar();
      return -105;
    }


  //entropy based on Etype
  pp[ENTR]=calc_Sfromu(rho,uint);

#ifdef MAGNFIELD
  //magnetic conserved=primitives
  pp[B1]=uu[B1]/gdetu;
  pp[B2]=uu[B2]/gdetu;
  pp[B3]=uu[B3]/gdetu;
#endif

#ifdef TRACER
  ldouble Dtr=uu[TRA]/gdetu*alpha; //uu[0]=gdetu rho ut
  pp[TRA]=Dtr/gamma;
#endif

  if(verbose) print_primitives(pp);

  //verify uunew agains uuoriginal
  //seems unnecessary - if this point reached, it already did its best and its best to stay like that
  if(Etype==U2P_HOT)
    {
      ldouble uunew[NV];
      p2u(pp,uunew,geom);
      
      ldouble errinv,maxerrinv=-1.;
      int iv;
      
      //do we recover rho properly
      iv=RHO;
      errinv = fabs((uunew[iv]-uu[iv])/uu[iv]);
      if(errinv > maxerrinv) maxerrinv=errinv;
      //internal energy
      if(Etype==U2P_HOT)
	{
	  iv=UU;
	  errinv = fabs((uunew[iv]-uu[iv])/uu[iv]);
	  if(errinv > maxerrinv) maxerrinv=errinv;
	}
      if(Etype==U2P_ENTROPY)
	{
	  iv=ENTR;
	  errinv = fabs((uunew[iv]-uu[iv])/uu[iv]);
	  if(errinv > maxerrinv) maxerrinv=errinv;
	}
      
      double inverr=1.e-2;
      
      if(Etype==U2P_ENTROPY) inverr=0.999;
      
      if(Etype==U2P_HOT) inverr=0.1;
      
      
      if(maxerrinv>inverr)// && verbose>0) 
    {
      
      if(Etype==U2P_ENTROPY && 0) { 
	printf("verify u2p (%d) failed: %e || ",Etype,maxerrinv);
	printf("%e %e | %e %e | %e %e\n",uunew[RHO],uu[RHO],uunew[ENTR],uu[ENTR],uunew[VX],uu[VX]);
	print_conserved(uu);
	print_conserved(uunew);
	print_primitives(pp);
	getch();
      }
      if(Etype==U2P_HOT && 0) { 
	printf("verify u2p (%d) failed: %e || ",Etype,maxerrinv);
	printf("%e %e | %e %e | %e %e\n",uunew[RHO],uu[RHO],uunew[UU],uu[UU],uunew[VX],uu[VX]);
	print_conserved(uu);
	print_conserved(uunew);
	print_primitives(pp);
	getch();
      }
      return -200;
    }
    }

  if(verbose>0) printf("u2p_solver returns 0\n");
  return 0; //ok

}

int count_entropy(int *n, int *n2)
{
  int nentr=0,nentrloc=0,ii,ix,iy,iz;
  int nentr2=0,nentrloc2=0;

  //counting the number of entropy inversions

  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  nentrloc+=get_cflag(ENTROPYFLAG,ix,iy,iz); 
	  nentrloc2+=get_cflag(ENTROPYFLAG2,ix,iy,iz); 
	}
#ifdef MPI
  MPI_Allreduce(&nentrloc, &nentr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  
  MPI_Allreduce(&nentrloc2, &nentr2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  
#else
  nentr=nentrloc;
  nentr2=nentrloc2;
#endif
    
  *n = nentr;
  *n2 = nentr2;
  return 0;
}

//backups entropy count to spit it into a silo file
int copy_entropycount()
{
  int ii,ix,iy,iz;
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2];
      ldouble val=get_cflag(ENTROPYFLAG,ix,iy,iz);
      set_cflag(ENTROPYFLAG3,ix,iy,iz,val);
    }

  return 0;
}

//tests
int
test_inversion()
{
  ldouble pp[NV],pp2[NV],uu[NV],ucon[4]={0.,0.,0.,0.};
  struct geometry geom,geomBL;
  int iv;

  fill_geometry(NX-2,0,0,&geom);
  fill_geometry_arb(NX-2,0,0,&geomBL,BLCOORDS);
  ucon[1]=1.e-12;
  conv_vels(ucon,ucon,VEL4,VEL4,geomBL.gg,geomBL.GG);
  trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
  conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
  pp[RHO]=1.;
  pp[UU]=1.;
  pp[VX]=ucon[1];
  pp[VY]=pp[VZ]=0.;

  print_primitives(pp);
  p2u(pp,uu,&geom);

  pp[VX]*=100000.*M_PI;
  pp[RHO]*=0.001245325124;
  pp[UU]*=23.124124214421124;
  
  PLOOP(iv) pp2[iv]=pp[iv];
  

  
  print_conserved(uu);
  printf("gdet = %e\n",geom.gdet);
  u2p_solver(uu,pp,&geom,U2P_HOT,0); 
  print_primitives(pp);

  return 0;

}


//tests
int
test_inversion_nonrel()
{
  ldouble pp[NV],pp2[NV],uu[NV],ucon[4]={0.,0.,0.,0.};
  struct geometry geom,geomBL;
  int iv;

  fill_geometry(NX-2,NY/2,0,&geom);

  print_metric(geom.gg);

  ucon[1]=0.001;
  ucon[2]=0.001;
  ucon[3]=0.00001;
  conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

  pp[RHO]=100.;
  pp[UU]=1.e-5;
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);
  pp[VX]=ucon[1];
  pp[VY]=ucon[2];
  pp[VZ]=ucon[3];

#ifdef MAGNFIELD
  pp[B1]=pp[B2]=pp[B3]=0.;
  pp[B1]=1.e-5;
#endif

#ifdef RADIATION
  pp[EE]=pp[UU];

  pp[FX]=pp[FY]=pp[FZ]=0.;
#endif


  print_primitives(pp);
  p2u(pp,uu,&geom);

  print_conserved(uu);
  printf("gdet = %e\n",geom.gdet);

  u2p_solver(uu,pp,&geom,U2P_HOT,0); 
  print_primitives(pp);

  return 0;

}



//**********************************************************************
//**********************************************************************
//**********************************************************************
//old Newton-Rapshon solver 
//iterates W, not Wp
//Etype == 0 -> hot inversion (uses D,Ttt,Tti)
//Etype == 1 -> entropy inversion (uses D,S,Tti)
//Etype == 2 -> hotmax inversion (uses D,Tti,u over rho max ratio)
//Etype == 3 -> cold inversion (uses D,Tti,u over rho min ratio)

int
u2p_solver_W(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  int i,j,k;
  ldouble rho,uint,w,W,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],Qconp[4],Qcovp[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
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
  int (*f_u2p)(ldouble,ldouble*,ldouble*,ldouble*,ldouble*);
 if(Etype==U2P_HOT) 
   f_u2p=&f_u2p_hot;
 if(Etype==U2P_ENTROPY) 
   f_u2p=&f_u2p_entropy;
 if(Etype==U2P_HOTMAX) 
   f_u2p=&f_u2p_hotmax;
 if(Etype==U2P_COLD) 
   f_u2p=&f_u2p_cold;
  /****************************/
 
  
  if(verbose>1) {printf("********************\n");print_conserved(uu);print_primitives(pp);}

  /****************************/
  //conserved quantities etc
  
  //alpha
  alpha=sqrt(-1./GG[0][0]);

  //D
  D=uu[0]/gdetu*alpha; //uu[0]=gdetu rho ut

  //conserved entropy "S u^t"
  Sc=uu[5]/gdetu*alpha; 

  //Q_mu=alpha T^t_mu
  Qcov[0]=(uu[1]/gdetu-uu[0]/gdetu)*alpha;
  Qcov[1]=uu[2]/gdetu*alpha;
  Qcov[2]=uu[3]/gdetu*alpha;
  Qcov[3]=uu[4]/gdetu*alpha;

  //Qp_mu=alpha T^t_mu 
  Qcovp[0]=uu[1]/gdetu*alpha;
  Qcovp[1]=uu[2]/gdetu*alpha;
  Qcovp[2]=uu[3]/gdetu*alpha;
  Qcovp[3]=uu[4]/gdetu*alpha;


  //Qp^mu
  indices_12(Qcovp,Qconp,GG);

  //Q^mu
  indices_12(Qcov,Qcon,GG);

#ifdef MAGNFIELD
  //curly B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1]/gdetu*alpha;
  Bcon[2]=uu[B2]/gdetu*alpha;
  Bcon[3]=uu[B3]/gdetu*alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
#else
  Bsq=QdotB=QdotBsq=0.;
  Bcon[0]=Bcon[1]=Bcon[2]=Bcon[3]=0.;
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

  //\beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
  ldouble betasqoalphasq=gg[0][1]*GG[0][1] + gg[0][2]*GG[0][2] + gg[0][3]*GG[0][3]; 
  ldouble alphasq=alpha*alpha;
  //Qdotnp=-E'=-E+D
  ldouble Dfactor = (-geom->gttpert + alphasq*betasqoalphasq)/(alphasq+alpha);
  ldouble Qdotnp = Qconp[0]*ncov[0] + D*(Dfactor) ; // -Qdotn-W = -Qdotnp-Wp


  /****************************/
  
  //initial guess for W = w gamma**2 based on primitives
  rho=pp[0];
  uint=pp[1];
  utcon[0]=0.;
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  //conv_vels(utcon,ucon,VELPRIM,VEL4,gg,GG);
  //indices_21(ucon,ucov,gg);
  conv_vels(utcon,utcon,VELPRIM,VELR,gg,GG);

  ldouble qsq=0.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=utcon[i]*utcon[j]*gg[i][j];
  ldouble gamma2=1.+qsq;
  ldouble gamma=sqrt(gamma2);

  //W
  W=(rho+GAMMA*uint)*gamma2;

  if(verbose>1) printf("initial W:%e\n",W);

  /****************************/
  
  //test if does not provide reasonable gamma2
  // Make sure that W is large enough so that v^2 < 1 : 
  int i_increase = 0;
  ldouble f0,f1,dfdW,err;
  ldouble CONV=U2PCONV; 
  ldouble EPS=1.e-4;
  ldouble Wprev=W;
  ldouble cons[7]={Qn,Qt2,D,QdotBsq,Bsq,Sc,Qdotnp};
  
  do
    {
      f0=dfdW=0.;

      //now invoked for all solvers:
      //if(Etype!=U2P_HOT) //entropy-like solvers require this additional check
      (*f_u2p)(W-D,cons,&f0,&dfdW,&err);

      if( ((( W*W*W * ( W + 2.*Bsq ) 
	    - QdotBsq*(2.*W + Bsq) ) <= W*W*(Qtsq-Bsq*Bsq))
	  || !isfinite(f0) || !isfinite(f0)
	  || !isfinite(dfdW) || !isfinite(dfdW))	  
	  && (i_increase < 50)) 
	{
	  if(verbose>0) printf("init W : %e -> %e (%e %e)\n",W,100.*W,f0,dfdW);
	  W *= 10.;
	  i_increase++;
	  continue;
	}
      else
	break;    
    }
  while(1);

  if(i_increase>=50)
    {
      return -150;
      printf("failed to find initial W for Etype: %d\n",Etype);
      printf("at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
      print_NVvector(uu);
      print_NVvector(pp);
      getchar();
    }

  //1d Newton solver
  int iter=0,fu2pret;
  do
    {
      Wprev=W;
      iter++;
     
      fu2pret=(*f_u2p)(W-D,cons,&f0,&dfdW,&err);

      //numerical derivative
      //fu2pret=(*f_u2p)((1.+EPS)*W-D,cons,&f1,&dfdW,&err);
      //dfdW=(f1-f0)/(EPS*W);

      if(verbose>1) printf("%d %e %e %e %e\n",iter,W,f0,dfdW,err);
 
      //convergence test
      if(err<CONV)
	break;
      
      if(dfdW==0.) {W*=1.1; continue;}

      ldouble Wnew=W-f0/dfdW;
      int idump=0;
      ldouble dumpfac=1.;

      //test if goes out of bounds and damp solution if so
      do
	{
	  ldouble f0tmp,dfdWtmp,errtmp;
	  f0tmp=dfdWtmp=0.;
	  //now for all solvers
	  //if(Etype!=U2P_HOT) //entropy-like solvers require this additional check
	  (*f_u2p)(Wnew-D,cons,&f0tmp,&dfdWtmp,&errtmp);
	  if(verbose>1) printf("sub (%d) :%d %e %e %e %e\n",idump,iter,Wnew,f0tmp,dfdWtmp,errtmp);
	  if( ((( Wnew*Wnew*Wnew * ( Wnew + 2.*Bsq ) 
		  - QdotBsq*(2.*Wnew + Bsq) ) <= Wnew*Wnew*(Qtsq-Bsq*Bsq))
	       || !isfinite(f0tmp) || !isfinite(f0tmp)
	       || !isfinite(dfdWtmp) || !isfinite(dfdWtmp))
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
	  
      if(idump>=100) 
	{
	  if(verbose>0) printf("damped unsuccessfuly\n");
	  return -101;
	}
	
      W=Wnew; 

      if(fabs(W)>BIG) 
	{
	  if(verbose>1) printf("W has gone out of bounds at %d,%d,%d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz); 
	  return -103;
	}


      if(fabs((W-Wprev)/Wprev)<CONV && err<1.e-1) break;
      //if(fabs((W-Wprev)/Wprev)<CONV && err<sqrt(CONV)) break;
    }
  while(iter<50);

 
  if(iter>=50)
    {
      if(verbose>0) printf("iter exceeded in u2p_solver with Etype: %d\n",Etype); //getchar();
      return -102;
    }


  if(!isfinite(W) || !isfinite(W)) {if(verbose) printf("nan/inf W in u2p_solver with Etype: %d\n",Etype); return -103;}
 
  if(verbose>1) 
    {
      fu2pret=(*f_u2p)(W-D,cons,&f0,&dfdW,&err);
      printf("end: %d %e %e %e %e\n",iter,W,f0,dfdW,err);
    }

  //W found, let's calculate v2 and the rest
  //ldouble v2=Qt2/W/W;

  ldouble Wsq,Xsq,v2;
	
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);  
  v2 = ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);

  gamma2=1./(1.-v2);
  gamma=sqrt(gamma2);
  rho=D/gamma;
  uint=1./GAMMA*(W/gamma2-rho);
  utcon[0]=0.;
  utcon[1]=gamma/(W+Bsq)*(Qtcon[1]+QdotB*Bcon[1]/W);
  utcon[2]=gamma/(W+Bsq)*(Qtcon[2]+QdotB*Bcon[2]/W);
  utcon[3]=gamma/(W+Bsq)*(Qtcon[3]+QdotB*Bcon[3]/W);

  if(!isfinite(utcon[1]))
    {
      //print_4vector(utcon);
      return -120;
    }


  if(uint<0. || gamma2<0. ||isnan(W) || !isfinite(W)) 
    {
      if(verbose>0) printf("neg u in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);//getchar();
      return -104;
    }

  //converting to VELPRIM
  conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  
  //returning new primitives
  pp[RHO]=rho;
  pp[UU]=uint;
  pp[VX]=utcon[1];
  pp[VY]=utcon[2];
  pp[VZ]=utcon[3];


  if(rho<0.) 
    {
      if(verbose>0) printf("neg rho in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);//getchar();
      return -105;
    }


  //entropy based on Etype
  pp[ENTR]=calc_Sfromu(rho,uint);

#ifdef MAGNFIELD
  //magnetic conserved=primitives
  pp[B1]=uu[B1]/gdetu;
  pp[B2]=uu[B2]/gdetu;
  pp[B3]=uu[B3]/gdetu;
#endif

#ifdef TRACER
  ldouble Dtr=uu[TRA]/gdetu*alpha; //uu[0]=gdetu rho ut
  pp[TRA]=Dtr/gamma;
#endif


  if(verbose) print_primitives(pp);

  // test the inversion
  /*
  ldouble uu2[NV];
  int iv;
  int lostprecision=0;
  p2u(pp,uu2,ggg);

  //if(verbose) print_NVvector(uu2);
  //if(verbose) print_NVvector(uu);
  
 
  for(iv=0;iv<NVMHD;iv++)
    {
      if(Etype==U2P_HOT) if(iv==5) continue;
      if(Etype==U2P_ENTROPY) if(iv==1) continue;
      if(Etype==U2P_COLD || Etype==U2P_HOTMAX) if(iv==1 || iv==5) continue;
      if(((iv==0 || iv==1 || iv==5) && fabs(uu2[iv]-uu[iv])/fabs(uu[iv]+uu2[iv])>1.e-1))
	lostprecision=1;
    }     
  
  if(lostprecision)
    {
      if(verbose>0 || 1)
	{
	  print_Nvector(uu,NV);
	  print_Nvector(uu2,NV);  
	  printf("u2p_solver lost precision:\n");      
	  //getchar();
	}
      
      //test
      return -106;
    }
  */

  if(verbose>0) printf("u2p_solver returns 0\n");
  return 0; //ok

}

//recovers only magnetic field primitives - used when correcting polar axis
int
u2p_solver_Bonly(ldouble *uu, ldouble *pp, void *ggg)
{
  //prepare geometry
  struct geometry *geom
    = (struct geometry *) ggg;

 ldouble gdet,gdetu;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  
#ifdef MAGNFIELD
  //magnetic conserved=primitives
  pp[B1]=uu[B1]/gdetu;
  pp[B2]=uu[B2]/gdetu;
  pp[B3]=uu[B3]/gdetu;
#endif

  return 0; //always succeeds
}

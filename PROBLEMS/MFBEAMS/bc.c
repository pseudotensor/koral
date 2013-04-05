//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

/**********************/
  
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

gdet_bc=get_g(g,3,4,ix,iy,iz);  
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
ldouble xx=get_x(ix,0);
ldouble rho,E,Fx,Fy,Fz,uint;

if(ix<0 ) 
  {
    ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;
    
    //flat gas profiles
    Tgas=T_AMB;
    rho=RHO_AMB;
    uint=calc_PEQ_ufromTrho(Tgas,rho);
    Fz=Fy=Fx=0.;
    pp[0]=rho;
    pp[1]=uint;
    pp[2]=0.;
    pp[3]=0.;
    pp[4]=0.;
    pp[5]=calc_Sfromu(rho,uint);

#ifdef RADIATION
    E=EEBEAM1;
    Fx=FRATIO1*E;
    pp[6]=E;
    pp[7]=Fx;
    pp[8]=Fy;
    pp[9]=Fz; 

    int irf;
#ifdef MULTIRADFLUID

    for(irf=1;irf<NRF;irf++)
      {
	pp[FX(irf)]=pp[FY(irf)]=pp[FZ(irf)]=0.;
	pp[EE(irf)]=SMALL;
      }

#endif
#endif

    prad_ff2lab(pp,pp,&geom);
    p2u(pp,uu,gg,GG);	 

#ifdef MULTIRADFLUID

    redistribute_radfluids(pp,uu,&geom);
    u2p_rad(uu,pp,&geom,&irf);
#endif						
    
    return 0;
  }

if(ix>=NX) 
  {
    ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;
    
    //flat gas profiles
    Tgas=T_AMB;
    rho=RHO_AMB;
    uint=calc_PEQ_ufromTrho(Tgas,rho);
    Fz=Fy=Fx=0.;
    pp[0]=rho;
    pp[1]=uint;
    pp[2]=0.;
    pp[3]=0.;
    pp[4]=0.;
    pp[5]=calc_Sfromu(rho,uint);

#ifdef RADIATION
    E=EEBEAM2;
    Fx=-FRATIO2*E;


    pp[6]=E;
    pp[7]=Fx;
    pp[8]=Fy;
    pp[9]=Fz; 

    int irf;
#ifdef MULTIRADFLUID

    for(irf=1;irf<NRF;irf++)
      {
	pp[FX(irf)]=pp[FY(irf)]=pp[FZ(irf)]=0.;
	pp[EE(irf)]=SMALL;
      }

#endif
#endif

    prad_ff2lab(pp,pp,&geom);
    p2u(pp,uu,gg,GG);	 

#ifdef MULTIRADFLUID

    redistribute_radfluids(pp,uu,&geom);
#endif						
    
    return 0;
  }


    //periodic
    while(iiz<0)    iiz+=NZ;
    while(iiz>=NZ)    iiz-=NZ; 
    //periodic
    while(iiy<0)    iiy+=NY;
    while(iiy>=NY)    iiy-=NY; 


    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);      
      }
    p2u(pp,uu,gg,GG);

  
    return 0;
  

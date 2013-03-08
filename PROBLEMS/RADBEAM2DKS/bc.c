//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

  gdet_bc=get_g(g,3,4,ix,iy,iz); 
  ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4],tup[4][4],tlo[4][4];
  pick_g(ix,iy,iz,gg);
  pick_T(emuup,ix,iy,iz,eup);
  pick_T(emulo,ix,iy,iz,elo);
  pick_G(ix,iy,iz,GG);
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);
ldouble xxvec[4];
get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvec,MYCOORDS,BLCOORDS);
//BL
ldouble xx=xxvec[1];

  //radius
  if(iz<0 && xx>BEAML && xx<BEAMR && IFBEAM) //hot boundary
    {
      ldouble Fx,Fy,Fz,rho,E,uint,vx;
      iix=ix;
      iiy=iy;
      iiz=0;
      rho=RHOAMB;
      E=calc_LTE_EfromT(TLEFT);
      uint=get_u(p,1,iix,iiy,iiz);
      rho=get_u(p,0,iix,iiy,iiz);
      vx=get_u(p,2,iix,iiy,iiz);
      
      Fz=NLEFT*E;
      Fy=Fx=0.;     

      pp[0]=rho;
      pp[1]=uint;
      pp[2]=vx;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz;

      //working in BL
      ldouble ggBL[4][5],GGBL[4][5];
      calc_g_arb(xxvec,ggBL,KERRCOORDS);
      calc_G_arb(xxvec,GGBL,KERRCOORDS);
      ldouble eupBL[4][4],eloBL[4][4];
      ldouble tupBL[4][4],tloBL[4][4];
      calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
      calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);
      prad_zamo2ff(pp,pp,ggBL,GGBL,eupBL);
      prad_ff2lab(pp,pp,ggBL,GGBL,tloBL);

      //to transform radiative primitives from BL to MYCOORDS
      trans_prad_coco(pp, pp, KERRCOORDS, MYCOORDS,xxvec,ggBL,GGBL,gg,GG);
     
      p2u(pp,uu,gg,GG);
      
      return 0.;
    }
  else if(iz<0 )
    {
      /*
      ldouble Fx,Fy,Fz,rho,E,uint,vx,T;
      rho=RHOAMB;
      T=TAMB;
      E=calc_LTE_EfromT(T);
      uint=calc_PEQ_ufromTrho(T,rho);

      Fx=Fy=Fz=0.*E;
      pp[0]=rho;
      pp[1]=uint;
      pp[2]=0;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz;
      p2u(pp,uu,gg,eup,elo);
      return 0;
      */
      
      
      iiz=0;
      iiy=iy;
      iix=ix;
      for(iv=0;iv<NV;iv++)
	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      //flux azimuthal only
      pp[7]=0.;
      pp[9]=0.;

      p2u(pp,uu,gg,GG);


      return 0;
      
    }
  
  else if(iz>=NZ ) //copy
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  

      iix=ix;
      iiy=iy;
      iiz=NZ-1;
      for(iv=0;iv<NV;iv++)
	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      p2u(pp,uu,gg,GG);
      return 0;
    }
  else if(ix<0 && 1) //copy
    {
      ldouble gdet=get_g(g,3,4,ix,iy,iz);  

      iix=0;
      iiy=iy;
      iiz=iz;
      for(iv=0;iv<NV;iv++)
	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      
      ldouble Fx,Fy,Fz,rho,E,uint,vx;
      //       rho=RHOAMB;
       //      E=calc_LTE_EfromT(TAMB);
      //pp[6]=E;


      p2u(pp,uu,gg,GG);


      return 0;
    }
#ifdef FLATBACKGROUND
   else if(ix>NX-1 || ix<0) //copy
    {

      ldouble gdet=get_g(g,3,4,ix,iy,iz);  

      iix=NX-1;
      iiy=iy;
      iiz=iz;
      for(iv=0;iv<NV;iv++)
 	{ 
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      ldouble Fx,Fy,Fz,rho,E,uint,vx;
      //       rho=RHOAMB;
       //      E=calc_LTE_EfromT(TAMB);
      //pp[6]=E;
      p2u(pp,uu,gg,GG);
      return 0;
    }
#endif
   else if(ix>NX-1 || 1 ) //fixed outern boundary
    {
      //TODO
    }
 


  iix=ix;
  iiz=iz;
  iiy=iy;
  //copy
  while(iix<0)    iix=0.;//iix+=NX;
  while(iix>=NX)    iix=NX-1;//iix-=NX; 
  //periodic
  while(iiy<0)    iiy+=NY;
  while(iiy>=NY)    iiy-=NY; 
 
  for(iv=0;iv<NV;iv++)
    {
      
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }
p2u(pp,uu,gg,GG);
return 0;
  
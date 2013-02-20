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
  ldouble xx=get_x(ix,0);

  //radius
  if(iz<0 && xx>BEAML && xx<BEAMR && IFBEAM) //hot boundary
    {
      ldouble Fx,Fy,Fz,rho,E,uint,vx;
      iix=ix;
      iiy=iy;
      iiz=0;
     

      E= 1000.*RADBEAMFLAT_ERAD;

      uint=get_u(p,1,iix,iiy,iiz);
      rho=get_u(p,0,iix,iiy,iiz);
      
      Fz=RADBEAMFLAT_FRATIO*E;
      Fy=Fx=0.;     

      pp[0]=rho;
      pp[1]=uint;
      pp[2]=0.;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz;

     
      //      prad_zamo2ff(pp,pp,gg,GG,eup);
      //      prad_ff2lab(pp,pp,gg,GG,tlo);
     
      p2u(pp,uu,gg,GG);

      
      return 0.;
    }
  else if(iz<0 )
    {
       ldouble Fx,Fy,Fz,rho,E,uint,vx;
      
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

      E= RADBEAMFLAT_ERAD;

      uint=get_u(p,1,iix,iiy,iiz);
      rho=get_u(p,0,iix,iiy,iiz);
      
      Fz=0.;
      Fy=Fx=0.;     

      pp[0]=rho;
      pp[1]=uint;
      pp[2]=0.;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz;

     
      //      prad_zamo2ff(pp,pp,gg,GG,eup);
      //      prad_ff2lab(pp,pp,gg,GG,tlo);
 

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
  else if(ix<0) //copy
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
   else if(ix>NX-1) //copy
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
  

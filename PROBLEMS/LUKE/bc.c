//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

gdet_bc=get_g(g,3,4,ix,iy,iz);  
gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],ggsrc[4][5],eup[4][4],elo[4][4],GG[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
ldouble xx=get_x(ix,0);
ldouble yy=get_x(iy,1);
ldouble zz=get_x(iz,2);

/**********************/


#ifdef SPHERICAL
  if(ix>=NX) //imposing velocity inflow
    {
      if(zz>Pi/2.) //outflow
	{
	  iiy=iy;
	  iiz=iz;
	  iix=NX-1;
	  for(iv=0;iv<NV;iv++)
	    {
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	    }
	}
      else
	{
	  pp[0]=1.;
	  pp[1]=UINTFRAC;
	}
      
      pp[2]=VELINX*cos(zz)*sin(yy)/sqrt(gg[1][1]);
      pp[3]=VELINX*cos(zz)*cos(yy)/sqrt(gg[2][2]);
      pp[4]=-VELINX*sin(zz)/sqrt(gg[3][3]);

      pp[5]=calc_Sfromu(pp[0],pp[1]);	      

      //converting from 3vel to VELPRIM
      conv_velsinprims(pp,VEL3,VELPRIM,gg,GG);
	      
      p2u(pp,uu,gg,GG);

      return 0.;
    }
  else if(ix<0) //bullet like wall
    {
      //reflection
      iix=-ix-1;
      iiz=iz;
      iiy=iy;
      for(iv=0;iv<NV;iv++)
	{
	  if(iv==2)
	    pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	  else
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      
      //converting from 3vel to VELPRIM
      conv_velsinprims(pp,VEL3,VELPRIM,gg,GG);

      p2u(pp,uu,gg,GG);
      return 0;
    }

  //theta 
  if(iy<0.) //spin axis, transmissive
    {      
      iiy=-iy-1;
      iiz=NZ-1-iz;
      iix=ix;
      for(iv=0;iv<NV;iv++)
	{
	   if(iv==4)
	      pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	   else
	     pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      p2u(pp,uu,gg,GG);
      return 0;
     }
  if(iy>=NY) //equatorial plane, reflection
    {
      iiy=NY-(iy-NY)-1;
      iiz=iz;
      iix=ix;
    	  
      for(iv=0;iv<NV;iv++)
	  {
	    if(iv==3)
	      pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	    else
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }
      p2u(pp,uu,gg,GG); 
      return 0; 
    }

  //reflections in phi
  if(iz<0.) //velocity axis
    {      
      iiz=-iz-1;
      iiy=iy;
      iix=ix;
      for(iv=0;iv<NV;iv++)
	{
	  if(iv==4)
	    pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	  else
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      p2u(pp,uu,gg,GG);
      return 0;
     }
  if(iz>=NZ) 
    {
      iiz=NZ-(iz-NZ)-1;
      iiy=iy;
      iix=ix;
    	  
      for(iv=0;iv<NV;iv++)
	  {
	    if(iv==4)
	      pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	    else
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }
      p2u(pp,uu,gg,GG); 
      return 0; 
    }
#endif

#ifdef CYLINDRICAL
  if(ix>=NX) //imposing velocity inflow
    {
      pp[0]=1.;
      pp[1]=UINTFRAC;
      pp[2]=pp[3]=pp[4]=0.;

      pp[2]=VELINX*cos(zz);
      pp[4]=-VELINX*sin(zz)/sqrt(gg[3][3]);


      pp[5]=calc_Sfromu(pp[0],pp[1]);	      

      //converting from 3vel to VELPRIM
      conv_velsinprims(pp,VEL3,VELPRIM,gg,GG);
	      
      p2u(pp,uu,gg,GG);

      return 0.;
    }
  else if(ix<0) //bullet like wall
    {
      //reflection
      iix=-ix-1;
      iiz=iz;
      iiy=iy;
      for(iv=0;iv<NV;iv++)
	{
	  if(iv==2)
	    pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	  else
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	}
      
      //converting from 3vel to VELPRIM
      conv_velsinprims(pp,VEL3,VELPRIM,gg,GG);

      p2u(pp,uu,gg,GG);
      return 0;
    }

iix=ix;
iiy=iy;
iiz=iz;

  //periodic in z (NY=1)
  while(iiy<0)    iiy+=NY;
  while(iiy>=NY)    iiy-=NY;   
   
  //periodic in phi
  while(iiz<0)    iiz+=NZ;
  while(iiz>=NZ)    iiz-=NZ; 

#endif

 for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,iix,iiy,iiz);
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }
  
  return 0;
  

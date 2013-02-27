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

/**********************/

  if(ix>=NX) //imposing velocity inflow
    {
      pp[0]=1.;
      pp[1]=UINTFRAC;
      pp[2]=pp[3]=pp[4]=0.;
      pp[2]=VELINX;
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

  //reflections in theta 
  if(iy<0.) //spin axis
    {      
      iiy=-iy-1;
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
  if(iy>=NY) //equatorial plane
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
   
  //symmetric in phi:
  while(iz<0)    iz+=NZ;
  while(iz>=NZ)    iz-=NZ; 

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,iix,iiy,iiz);
      pp[iv]=get_u(p,iv,iix,iiy,iiz);      
    }
  
  return 0;
  

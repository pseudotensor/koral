//returns problem specific BC
//ix, iy, iz hold the ghost cell indices, e.g., (-1,0,0)
//BCtype gives the type of the boundary

int iix,iiy,iiz,iv;

/***********************************************/
//structure of geometry
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
/***********************************************/

//atmosphere
if(BCtype==ZBCLO) 
  {
   /*
   pp[RHO]=RHOAMB; 
   pp[UU]=UUAMB; 
   pp[VZ]=0.;
   pp[VY]=0.;
   pp[VX]=0.;
   pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

   p2u(pp,uu,&geom);
    */
   return 0;  
  }


//reflection:
if(BCtype==XBCLO) 
  {
    pp[RHO]=RHOAMB; 
    pp[UU]=UUAMB; 
    pp[VZ]=0.;
    pp[VY]=0.;
    pp[VX]=WINDVX;
    pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

    p2u(pp,uu,&geom);
    return 0;  
  }

//periodic:
if(BCtype==YBCLO) 
  {
    iiy=iy+NY;
    iix=ix;
    iiz=iz;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    p2u(pp,uu,&geom);
    return 0;  
  }

//outflow:
if(BCtype==XBCHI) 
  {
    //iix=NX - (iz-NX) -1;
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    //pp[VX]=-get_u(p,VX,iix,iiy,iiz);

    p2u(pp,uu,&geom);
    return 0;  
  }

//periodic:
if(BCtype==YBCHI) 
  {
    iiy=iy-NY;
    iix=ix;
    iiz=iz;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    p2u(pp,uu,&geom);
    return 0;  
  }

//outflow:
if(BCtype==ZBCHI) 
  {
    iiz=NZ-1;
    iiy=iy;
    iix=ix;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    p2u(pp,uu,&geom);
    return 0;  
  }

//and that is all
 
return 0;


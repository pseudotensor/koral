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
   pp[RHO]=RHOAMB; 
   pp[UU]=UUAMB; 
   pp[VZ]=0.;
   pp[VY]=0.;
   pp[VX]=0.;
   pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

   p2u(pp,uu,&geom);
   return 0;  
  }


//reflection:
if(BCtype==XBCLO) 
  {
    iix=-ix-1;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    pp[VX]=-get_u(p,VX,iix,iiy,iiz);

    p2u(pp,uu,&geom);
    return 0;  
  }

//reflection:
if(BCtype==YBCLO) 
  {
    iiy=-iy-1;
    iix=ix;
    iiz=iz;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    pp[VY]=-get_u(p,VY,iix,iiy,iiz);

    p2u(pp,uu,&geom);
    return 0;  
  }

//outflow:
if(BCtype==XBCHI) 
  {
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    p2u(pp,uu,&geom);
    return 0;  
  }

//outflow:
if(BCtype==YBCHI) 
  {
    iiy=NY-1;
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


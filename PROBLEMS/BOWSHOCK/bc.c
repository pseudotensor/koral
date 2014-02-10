//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
/***********************************************/

//wind:
if(BCtype==ZBCLO) 
  {
   pp[RHO]=RHOAMB; 
   pp[UU]=UUAMB; 
   pp[VZ]=VELAMB;
   pp[VY]=0.;
   pp[VX]=0.;
   pp[5]=calc_Sfromu(pp[0],pp[1]);

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

//outflows:
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


for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//and that is all
 
return 0;


//returns problem specific BC
//ix, iy, iz hold the ghost cell indices, e.g., (-1,0,0)
//BCtype gives the type of the boundary

int iix,iiy,iiz,iv;

/***********************************************/
//structure of geometry
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
/***********************************************/


if(BCtype==XBCHI) 
  {
    //outflow
    iix=NX-1;
    iiy=iy;
    iiz=iz;

    
    for(iv=0;iv<NV;iv++) //NV - number of variables
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

 //check for inflow
    if(pp[VX]<0.)
      {
	pp[VX]=0.;
      }
    p2u(pp,uu,&geom); //convert to conserved
   return 0;  
  }

if(BCtype==XBCLO) 
  {
    //outflow
    iix=0;
    iiy=iy;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++) //NV - number of variables
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

 //check for inflow
    if(pp[VX]>0.)
      {
	pp[VX]=0.;
      }

    p2u(pp,uu,&geom); //convert to conserved
   return 0;  
  }

//outflow:
if(BCtype==YBCLO) 
  {
    iix=ix;
    iiy=0;
    iiz=iz;
    
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

 //check for inflow
    if(pp[VY]>0.)
      {
	pp[VY]=0.;
      }
   

    p2u(pp,uu,&geom);


    return 0;  
  }


//outflow:                                                                                                   
if(BCtype==YBCHI)
  {
    iix=ix;
    iiy=NY-1;
    iiz=iz;

    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }    

    //check for inflow
    if(pp[VY]<0.)
      {
	pp[VY]=0.;
      }

    p2u(pp,uu,&geom);
    return 0;
  }
 
//periodic in z, used under OMP
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }


//and that is all

return 0;

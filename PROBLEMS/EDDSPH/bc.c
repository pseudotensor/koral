//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
ldouble rsph=geomBL.xx;
/***********************************************/

/***********************************************/
/***********************************************/
if(ix>=NX) //analytical solution within the torus and atmosphere outside
  {
    /*
    //hydro atmosphere
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,2);
    //rad atmosphere
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,1);
#endif
    */

    //outflow BC
    for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,NX-1,iy,iz);
       }

    //zero radial velocity
    if(pp[VX]<0.) pp[VX]=0.;
    //pp[RHO]=0.;
    
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom);
 
    return 0;
  }
 
/***********************************************/
/***********************************************/
if(ix<0) //star surface
   {
     for(iv=0;iv<NVHD;iv++)
       {
	 pp[iv]=get_u(p,iv,0,iy,iz);
       }

     ldouble PSTAR=calc_PEQ_ufromTrho(TEMPSTAR,RHOSTAR)*GAMMAM1;
     ldouble KKK = PSTAR / pow(RHOSTAR,GAMMA);
     ldouble C0 = (GAMMA/GAMMAM1*PSTAR/RHOSTAR+1.)*sqrt(1.-2./RSTAR);
     ldouble pre = pow(KKK,-1./GAMMAM1)*pow(GAMMAM1/GAMMA*(C0*sqrt(rsph/(rsph-2.))-1.),GAMMA/GAMMAM1);
     ldouble rho = pow(pre/KKK,1./GAMMA);
     ldouble temp=calc_PEQ_Tfromurho(pre/GAMMAM1,rho);

     pp[RHO]=rho;
     pp[UU]=calc_PEQ_ufromTrho(temp,pp[RHO]);//density
    
     pp[VX]=pp[VY]=pp[VZ]=0.;

     //flux
     ldouble Fr = FLUXBETA / KAPPA_ES_COEFF / RSTAR / (RSTAR-2.);

     //ortho-normal
     pp[7]=Fr; ////R^tr
     pp[6]=2.*Fr; //R^tt
     pp[8]=pp[9]=0.; //R^tinne

     //     ldouble Rijlab[4][4];
     //     calc_Rij_ff(pp,Rijlab);  
     //     trans22_on2cc(Rijlab,Rijlab,geom.tlo);  

     //converts upper row (R^tmu) from orthormal to coord. basis
     prad_on2lab(pp,pp,&geom);
 
     //testing if interpolated primitives make sense
     check_floors_hd(pp,VELPRIM,&geom);
     //end of floor section

     
     //     printf("%d | rho %e | pre %e | R(tr) %e | R(tt) %e | Rtr %e | Rtt %e | kes %e\n",ix,pp[RHO],GAMMAM1*pp[UU],Fr,2.*Fr,Rijlab[0][1],Rijlab[1][1],KAPPA_ES_COEFF); getchar();

     p2u(pp,uu,&geom);
     return 0;
   }

/***********************************************/
/***********************************************/
//periodic in phi:

if(iy<0 || iy>NY-1)
  {
    iiz=iz;
    iiy=iy;
    iix=ix;
    if(iy<0) iiy=iy+NY;
    if(iy>NY-1) iiy=iy-NY;

    for(iv=0;iv<NV;iv++)
      {
	uu[iv]=get_u(u,iv,iix,iiy,iiz);
	pp[iv]=get_u(p,iv,iix,iiy,iiz);      
      }
  }

/***********************************************/
/***********************************************/
//periodic in phi:
if(iz<0 || iz>NZ-1)
  {
    iiz=iz;
    iiy=iy;
    iix=ix;
    if(iz<0) iiz=iz+NZ;
    if(iz>NZ-1) iiz=iz-NZ;

    for(iv=0;iv<NV;iv++)
      {
	uu[iv]=get_u(u,iv,iix,iiy,iiz);
	pp[iv]=get_u(p,iv,iix,iiy,iiz);      
      }
  }

//and that is all
 
return 0;


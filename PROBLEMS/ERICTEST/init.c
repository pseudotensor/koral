
ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/
//hydro atmosphere
ldouble d2 = (geom.xx*geom.xx + geom.yy*geom.yy + geom.zz*geom.zz); 
pp[RHO]=RHOAMB + (RHOCLOUD - RHOAMB)* exp(-d2/2./SIGMACLOUD/SIGMACLOUD); 
pp[UU]=UUAMB + (UUCLOUD - UUAMB)* exp(-d2/2./SIGMACLOUD/SIGMACLOUD); 
pp[VZ]=(RHOAMB*VELAMB + (RHOCLOUD*0. - RHOAMB*VELAMB)* exp(-d2/2./SIGMACLOUD/SIGMACLOUD))/pp[RHO]; 
pp[VY]=0.;
pp[VX]=0.;

/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
//to conserved
p2u(pp,uu,&geom);

/*
print_primitives(pp);
printf("%e\n",(RHOAMB*VELAMB + (RHOCLOUD*0. - RHOAMB*VELAMB)* exp(-d2/2./SIGMACLOUD/SIGMACLOUD))/pp[RHO]);
print_conserved(uu);
getch();
*/
/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);

ldouble uu[NV];
ldouble pp[NV];

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomcart;
fill_geometry_arb(ix,iy,iz,&geomcart,MINKCOORDS);

/***********************************************/

ldouble r=geom.xx;
ldouble th=geom.yy;
ldouble ph=geom.zz;


ldouble x=geomcart.xx;
ldouble y=geomcart.yy;
ldouble z=geomcart.zz;


pp[RHO]=RHO_AMB + RHO_BLOB*exp(-((y-5.)*(y-5.) + (z-20.)*(z-20.) + (x*x))/10.);
pp[UU]=U_AMB*RHO_AMB/pp[RHO];

pp[VX]=0.;
pp[VY]=0.;
pp[VZ]=0.;

//velocity in cartesian
ldouble ucart[3]={0.,-.1,0.};

/*
ldouble carttetrad[3][3];
carttetrad[0][0]=sin(th)*cos(ph);
carttetrad[0][1]=sin(th)*sin(ph);
carttetrad[0][2]=cos(th);
	 
carttetrad[1][0]=cos(th)*cos(ph);
carttetrad[1][1]=cos(th)*sin(ph);
carttetrad[1][2]=-sin(th);

carttetrad[2][0]=-sin(ph);
carttetrad[2][1]=cos(ph);
carttetrad[2][2]=0.;


//velocity in spherical
usph[1]=dot(ucart,carttetrad[0]);
usph[2]=dot(ucart,carttetrad[1])/r;
usph[3]=dot(ucart,carttetrad[2])/r/sin(th);
*/
ldouble usph[4];
ldouble vx,vy,vz;
vx=vz=0.;
vy=-.1;

ldouble vr,vth,vph;
ldouble cosph,sinth,costh,sinph;
sinth=sin(th);
costh=cos(th);
sinph=sin(ph);
cosph=cos(ph);


vr =-((-(cosph*sinth*vx) - sinph*sinth*vy - Power(cosph,2)*costh*vz - 
          costh*Power(sinph,2)*vz)/
      ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
vth = -((-(cosph*costh*vx) - costh*sinph*vy + Power(cosph,2)*sinth*vz + 
          Power(sinph,2)*sinth*vz)/
        ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
vph = -((sinph*vx - cosph*vy)/(Power(cosph,2) + Power(sinph,2)));

vth /= r;
vph /= r*sinth;

usph[1]=vr;
usph[2]=vth;
usph[3]=vph;

if (pp[RHO]/RHO_AMB > 2.)
  {   
    pp[VX]=vr;
    pp[VY]=vth;
    pp[VZ]=vph;
    /*
    ldouble vx,vy,vz;
    usph[2]*=r;
    usph[3]*=r*sin(th);
		  
    vx = sin(th)*cos(ph)*usph[1] 
      + cos(th)*cos(ph)*usph[2]
      - sin(ph)*usph[3];

    vy = sin(th)*sin(ph)*usph[1] 
      + cos(th)*sin(ph)*usph[2]
      + cos(ph)*usph[3];


    vz = cos(th)*usph[1] 
      - sin(th)*usph[2];

		  printf("%f %f %f > %e %e %e > %e %e %e\n",r,th,ph,pp[VX],pp[VY],pp[VZ],vx,vy,vz);
    */
  }

/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
p2u(pp,uu,&geom);
/***********************************************/

//print_primitives(pp); getch();

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

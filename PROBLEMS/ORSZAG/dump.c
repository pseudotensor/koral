
//calculate div. B
ldouble divB;
if(ix>=0 && ix<NX && iy>=0 && iy<NY && iz==0)
  {
    divB = (get_u(p,B1,ix,iy,iz) + get_u(p,B1,ix,iy-1,iz) - get_u(p,B1,ix-1,iy,iz) - get_u(p,B1,ix-1,iy-1,iz))/(2.*(get_x(ix+1,0)-get_x(ix,0)))
      + (get_u(p,B2,ix,iy,iz) + get_u(p,B2,ix-1,iy,iz) - get_u(p,B2,ix,iy-1,iz) - get_u(p,B2,ix-1,iy-1,iz))/(2.*(get_x(iy+1,1)-get_x(iy,1)));
  }
v2=divB;

v1=calc_divB(ix,iy,iz);

//printf("%e %e\n",v1,v2);getch();      

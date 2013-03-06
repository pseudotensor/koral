//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

ldouble minx,miny,maxx,maxy;

#ifdef myMKS1COORDS
minx= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxx= .66*(exp(get_xb(NX,0))+MKS1R0);
miny= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxy= .66*(exp(get_xb(NX,0))+MKS1R0);
#else
minx= -.02*get_xb(NX,0);
maxx= 1.02*get_xb(NX,0);
miny= -.02*get_xb(NX,0);
maxy= 1.02*get_xb(NX,0);
#endif

  fprintf(fgnu,
	  "set table \"table.gp\"\n"
	  "set contour base\n"
	  "unset surface\n"
	  "set log z\n"
	  "set cntrparam level discrete .01,.05,.1,.3,.5,.7,.9,1.1,1.3,1.5,1.7 \n"
	  //	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):24 w l\n"
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 3\n"
	  "set style line 11 lt 1 lw 2 lc 6\n"
	  "set style line 2 lt 1 lw 2 lc 2\n"
	  "set style line 3 lt 1 lw 2 lc 2\n"
	  "set style line 21 lt 3 lw 1 lc -1\n"
	  //	  "set style arrow 1 head nofilled size screen 0.002,35 ls 3\n"
	  "set style arrow 1 ls 1\n"
	  
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 23,28,3\n"
	  "set palette model RGB rgbformulae 7,8,9\n"
	  "set palette model RGB rgbformulae 6,3,21\n"
	  "set palette model RGB rgbformulae 35,3,9\n"

	  "unset log cb\n"
	  "set term gif large size 1000,700\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set lmargin at screen 0.08\n"
	  "set rmargin at screen 0.43\n"
	  "set bmargin at screen .45\n"
	  "set tmargin at screen .92\n"
	  "set size .5,1\n"
	  "set autoscale\n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"
	  "unset log cb\n"
	  //	  "set cbrange [0.005:.35]\n"
	  "set ylabel \"z\"\n"
	  "set cblabel \"\"\n"
	  //	  "set log cb\n"
	  "set title \"rho / velocity\"\n"
	  
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($14) ti \"\" w l ls 1\n"

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset log cb\n"
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):(($16*sin($2)+$17*cos($2))/(%f)):((-$17*sin($2)+$16*cos($2))/(%f)) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  //"plot \"table.gp\" w l ls 1\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"

	  "set table \"table.gp\"\n"
	  "set contour base\n"
	  "set log z\n"
	  "set cntrparam level discrete .01,.05,.1,.3,.5,.7,.9,1.1,1.3,1.5,1.7 \n"
	  //	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):14 w l\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"

	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  //"plot \"table.gp\" w l ls 11\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"	  

	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .45\n"
	  "set tmargin at screen .92\n"
	  "set cblabel \"\"\n"
	  "set ylabel \"y\"\n"
	  "set title \"radiative E / flux\"\n"
	  "unset log cb\n"
	  //"set log cb\n"
	  
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($20) ti \"\" w l ls 1\n"

          "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset log cb\n"
	  //	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):(($21*sin($2)+$22*cos($2))/(%f)):((-$22*sin($2)+$21*cos($2))/(%f)) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):(($21*sin($2)+$22*cos($2))/($20/%f)):((-$22*sin($2)+$21*cos($2))/($20/%f)) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	  "unset pm3d\n"	
	  "set table \"table.gp\"\n"
	  "set contour base\n"
	  "set log z\n"
	  "set cntrparam level discrete .01,.05,.1,.3,.5,.7,.9,1.1,1.3,1.5,1.7 \n"
	  //	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):24 w l\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"

	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  //"plot \"table.gp\" w l ls 1\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"	  
	  
	  "unset log cb\n"
	  "unset colorbox\n"
	  "set lmargin at screen 0.08\n"
	  "set rmargin at screen 0.43\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .40\n"
	  "set title \"\"\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"density\"\n"
	  "set autoscale\n"
	  "unset log y\n"
	  "plot \"%s\" u 1:14 ti \"\" w l ls 11  \n"

	  "set lmargin at screen 0.55\n"
	  "set rmargin at screen 0.9\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .40\n"
	  "set title \"\"\n"
	  "set yrange [0.09:.5]\n"
	  "set ylabel \"radiative en.density\"\n"
	  "set autoscale\n"
	  "unset log y\n"
	  "set label 1 \"t=%f\" at screen .45, .015\n"
	  "plot \"%s\" u 1:($20) ti \"\" w l ls 2 \n"
	  ,fname2,
	  minx,
	  maxx,
	  miny,
	  maxy,
	  fname,
	  fname,1./(get_xb(NX,0)/11)*.5,1./(get_xb(NX,0)/11)*.5,NX/11+1,NY/11+1,
	  fname,
	  fname,3.,3.,NX/11+1,NY/11+1,
	  fname,t,fname);  
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/

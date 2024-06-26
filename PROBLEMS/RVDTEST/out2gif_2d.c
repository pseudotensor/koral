//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

ldouble minx,miny,maxx,maxy;

#ifdef myMKS1COORDS
minx= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxx= 1.06*(exp(get_xb(NX,0))+MKS1R0);
miny= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxy= 1.06*(exp(get_xb(NX,0))+MKS1R0);
#else
minx= -.02*get_xb(NX,0);
maxx= 1.02*get_xb(NX,0);
miny= -.02*get_xb(NX,0);
maxy= 1.02*get_xb(NX,0);
#endif
/*
#ifndef RADIATION
  fprintf(fgnu,
          

	  "set term gif large size 1100,500\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"

	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 3\n"
	  "set style line 11 lt 1 lw 2 lc 6\n"
	  "set style line 2 lt 1 lw 2 lc 2\n"
	  "set style line 3 lt 1 lw 2 lc 2\n"
	  "set style line 21 lt 3 lw 1 lc -1\n"
	 
	  "set style arrow 1 ls 2\n"
	  
	  "set palette model RGB rgbformulae 7,5,15\n"
	 
	  "set autoscale\n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"

	  "set lmargin at screen 0.04\n"
	  "set rmargin at screen 0.41\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .95\n"
	  "set log cb\n"
	  "set ylabel \"z\"\n"
	  "set xlabel \"x\" offset 0,1\n"
	  "set cblabel \"\"\n"
	  "set title \"rho\" offset 0,-1\n"
	  "set cbrange [1e-3:1000]\n"
	  
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($14) ti \"\" w l ls 1\n"

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "set xlabel \"\"\n"
	  "unset tics\n"
	  "unset title\n"
	  "unset border\n"
	  "unset log cb\n"
	  //	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  //	  "(($16*sin($2)+$17*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5*%f):"
	  //	  "((-$17*sin($2)+$16*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	  "set border\n"
	  "set tics\n"
	  "set lmargin at screen 0.53\n"
	  "set rmargin at screen 0.90\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .95\n"
	  "unset log cb\n"
	  "set xlabel \"x\" offset 0,1\n"
	  "set ylabel \"z\"\n"
	  "set cblabel \"\"\n"
	  "set title \"uu\" offset 0,-1\n"
	  "set cbrange [0:1]\n"
	  
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($15) ti \"\" w l ls 1\n"

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "set xlabel \"\"\n"
	  "unset tics\n"
	  "unset title\n"
	  "unset border\n"
	  "unset log cb\n"
	  //	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  //	  "(($16*sin($2)+$17*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5*%f):"
	  //	  "((-$17*sin($2)+$16*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  "(($16*sin($2)+$17*cos($2))*%f):"
	  "((-$17*sin($2)+$16*cos($2))*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

 	  ,fname2,
	  minx,
	  maxx,
	  miny,
	  maxy,
	  fname,
	  //	  fname,1.,1.,NX/21+1,NY/21+1,
	  fname,
	  //	  fname,1.,1.,NX/21+1,NY/21+1
	  fname,3.,3.,NX/21+1,NY/21+1
	  );  
#else
*/
  fprintf(fgnu,
"set table \"table.gp\"\n"
	  "set contour base\n"
	  "unset surface\n"
	  "set log z\n"
	  "set cntrparam levels discrete 1 \n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):25 w l\n"
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

	  "set term gif large size 1100,500\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"

	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 3\n"
	  "set style line 11 lt 1 lw 2 lc 6\n"
	  "set style line 2 lt 1 lw 2 lc 2\n"
	  "set style line 3 lt 1 lw 2 lc 2\n"
	  "set style line 21 lt 3 lw 1 lc -1\n"
	 
	  "set style arrow 1 ls 2\n"
	  
	  "set palette model RGB rgbformulae 7,5,15\n"
	  /*
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 23,28,3\n"
	  "set palette model RGB rgbformulae 7,8,9\n"
	  //"set palette model RGB rgbformulae 35,3,9\n"
	  //	  "set palette model RGB rgbformulae 6,3,21\n"
	  */
	  "set autoscale\n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"

	  "set lmargin at screen 0.04\n"
	  "set rmargin at screen 0.41\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .95\n"
	  "set log cb\n"
	  "set ylabel \"z\"\n"
	  "set xlabel \"x\" offset 0,1\n"
	  "set cblabel \"\"\n"
	  "set title \"rho\" offset 0,-1\n"
	  //	  "set cbrange [1e-3:1000]\n"
	  
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($14) ti \"\" w l ls 1\n"

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "set xlabel \"\"\n"
	  "unset tics\n"
	  "unset title\n"
	  "unset border\n"
	  "unset log cb\n"

	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  "(($16*sin($2)+$17*cos($2))*%f):"
	  "((-$17*sin($2)+$16*cos($2))*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l ls 1\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"

	  "set border\n"
	  "set tics\n"
	  "set lmargin at screen 0.53\n"
	  "set rmargin at screen 0.90\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .95\n"
	  "set log cb\n"
	  "set xlabel \"x\" offset 0,1\n"
	  "set ylabel \"z\"\n"
	  "set cblabel \"\"\n"
	  //	  "set cbrange [1.e-5:100]\n"
	  //"set autoscale cb\n"

#ifdef RADIATION	  
	  "set title \"E\" offset 0,-1\n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($20) ti \"\" w l ls 1\n"
#else
	  "set title \"u_int\" offset 0,-1\n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($15) ti \"\" w l ls 1\n"
#endif

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "set xlabel \"\"\n"
	  "unset tics\n"
	  "unset title\n"
	  "unset border\n"
	  "unset log cb\n"
#ifdef RADIATION	  
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  "(($21*sin($2)+$22*cos($2))/(($21*sin($2)+$22*cos($2))**2+(-$22*sin($2)+$21*cos($2))**2)**.5*%f):"
	  "((-$22*sin($2)+$21*cos($2))/(($21*sin($2)+$22*cos($2))**2+(-$22*sin($2)+$21*cos($2))**2)**.5*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
#else
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  "(($16*sin($2)+$17*cos($2))*%f):"
	  "((-$17*sin($2)+$16*cos($2))*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
#endif
	  //	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  //	  "(($21*sin($2)+$22*cos($2))*%f):"
	  //	  "((-$22*sin($2)+$21*cos($2))*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l ls 1\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"
	  ,fname
 	  ,fname2,
	  minx,
	  maxx,
	  miny,
	  maxy,
	  fname,
	  fname,3.,3.,NX/21+1,NY/21+1,
	  fname,
#ifdef RADIATION
	  fname,1.,1.,NX/21+1,NY/21+1
#else
	  fname,3.,3.,NX/21+1,NY/21+1
#endif
	  );  
//#endif	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/

//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

double minx,maxx,miny,maxy;

#ifdef CYLINDRICAL
minx= -1.01*get_xb(NX+1,0);
maxx= 1.01*get_xb(NX+1,0);
miny=minx;
maxy=maxx;
#endif
#ifdef SPHERICAL 
minx= -1.1*get_xb(NX+1,0);
maxx= 1.1*get_xb(NX+1,0);
miny=-0.05*maxx;
  maxy=maxx;
#endif


  fprintf(fgnu,
	  "set view map\n"
	  "set pm3d\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 2\n"

	  "set style arrow 1 ls 1\n"
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 35,3,9\n"



	  "unset surface\n"
#ifdef CYLINDRICAL
	  "set term gif large size 800,700\n"
#endif
#ifdef SPHERICAL
	  "set term gif large size 1200,550\n"
#endif
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set lmargin at screen 0.1\n"
	  "set rmargin at screen 0.8\n"
	  "set bmargin at screen .09\n"
	  "set tmargin at screen .92\n"
	  "set size .5,1\n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"density / velocity \" offset 0,0\n"
#ifdef CYLINDRICAL
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):14 ti \"\" w l \n"
#endif
#ifdef SPHERICAL
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):14 ti \"\" w l \n"
#endif
	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
#ifdef CYLINDRICAL
	  "plot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($16*cos($3)-$18*$1*sin($3))*(%e)):(($18*$1*cos($3)+$16*sin($3))*%e) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
,fname2,minx,maxx,miny,maxy,fname, fname, 20.,20., (int)(NX/10),(int)(NZ/20)
#endif
#ifdef SPHERICAL
	  "plot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($16*cos($3)-$18*$1*sin($3))*(%e)):(($18*$1*cos($3)+$16*sin($3))*%e) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
,fname2,minx,maxx,miny,maxy,fname, fname, 15.,15., (int)(NX/10),(int)(NZ/20)
#endif

);
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/

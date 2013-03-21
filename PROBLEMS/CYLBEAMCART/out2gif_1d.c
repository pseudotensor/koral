//KORAL - bc.c
//gives the gnuplot script used for plotting gif files on the go - 1 dimension
//**********************
//included in:
//int convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
//from misc.c
//**********************

/****************************************/
/****************************************/
/****************************************/

ldouble minx,maxx,miny,maxy;
#ifdef myMCYL1COORDS
minx= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxx= 1.02*(exp(get_xb(NX,0))+MKS1R0);
miny= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxy= 1.02*(exp(get_xb(NX,0))+MKS1R0);
#else
minx= -.02*get_xb(NX,0);
maxx= 1.02*get_xb(NX,0);
miny= -.02*get_xb(NX,0);
maxy= 1.02*get_xb(NX,0);
#endif


fprintf(fgnu,
	"set style line 1 lw 2 lc 1 lt 1\n"
	"set style line 2 lw 2 lc 2 lt 1\n"
	"set style line 3 lw 2 lc 3 lt 1\n"
	"set style line 4 lw 2 lc 4 lt 1\n"
	"set style line 10 lw 2 lc 3 lt 3\n"
	"set style line 11 lw 2 lc 4 lt 3\n"
	"set style line 12 lw 2 lc 5 lt 3\n"
	"set style line 13 lw 2 lc 0 lt 3\n"
	"set style line 14 lw 2 lc 7 lt 3\n"
	"set term gif large size 1200,600\n"
	"set output \"%s\"\n"
	"set size 1,1\n"
	"set origin 0,0\n"
	"set multiplot\n"
	"set autoscale\n" 
	"set label \"t=%.2e (%.2e s)\" at screen .48, .98\n"

	"set lmargin at screen 0.07\n"
	"set rmargin at screen 0.33\n"
	"set bmargin at screen .07\n"
	"set tmargin at screen .42\n"
	"set autoscale\n"
	"set xrange [%f:%f]\n"
	//	  "unset log y\n"
	"set format x \"%%.1f\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:($27) w lp ls 1 pt 7 ps .5 ti \"rad fixup\"\n"

	"set lmargin at screen 0.40\n"
	"set rmargin at screen 0.66\n"
	"set bmargin at screen .07\n"
	"set tmargin at screen .42\n"
	//	  "unset log y\n"
	"set format x \"%%.1f\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	//	  "plot \"%s\" u 1:27 w lp ls 2 pt 7 ps .5  ti \"tau_abs\", \"%s\" u 1:26 w lp ls 3 pt 7 ps .5  ti \"tau_tot\"\n"
	"plot \"%s\" u 1:20 w lp ls 2 pt 7 ps .5  ti \"E_rad\"\n"

	"set lmargin at screen 0.73\n"
	"set rmargin at screen 0.99\n"
	"set bmargin at screen .07\n"
	"set tmargin at screen .42\n"
	//	  "unset log y\n"
	"set format x \"%%.1f\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:($21+1.e-80) w lp ls 2 pt 7 ps .5  ti \"Fx\"\n"

	"set lmargin at screen 0.07\n"
	"set rmargin at screen 0.33\n"
	"set bmargin at screen .5\n"
	"set tmargin at screen .9\n"
	//	  "unset log y\n"
	"set format x \"\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	//	  "set log y\n"
	"plot \"%s\" u 1:15 w lp ls 3 pt 7 ps .5 ti \"u_int\", \"%s\" u 1:20 w lp ls 7 pt 7 ps .5 ti \"E_rad\"\n"
	//	  "unset log y\n"

	"set lmargin at screen 0.40\n"
	"set rmargin at screen 0.66\n"
	"set bmargin at screen .5\n"
	"set tmargin at screen .9\n"
	"set format x \"\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:26 w p ls 4 pt 7 ti \"f\"\n"

	"set lmargin at screen 0.73\n"
	"set rmargin at screen 0.99\n"
	"set bmargin at screen .5\n"
	"set tmargin at screen .9\n"
	"set format x \"\"\n"
	"set format y \"%%.1e\"\n" 
	"set xlabel \"\"\n"
	"set ylabel \"\"\n"
	"plot \"%s\" u 1:($23) w lp ls 4 pt 7 ti \"Fph\""
	,fname2,t,t/CCC,minx,maxx,fname,fname,fname,fname,fname,fname,fname);

/****************************************/
/****************************************/
/****************************************/


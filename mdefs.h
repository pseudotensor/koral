/*************************************************************************

                        Mathematica source file

        Copyright 1986 through 1999 by Wolfram Research Inc.


*************************************************************************/

/* C language definitions for use with Mathematica output */


#define Power(x, y)	(powl((ldouble)(x), (ldouble)(y)))
#define Sqrt(x)		(sqrtl((ldouble)(x)))
#define Sqrtl(x)        (sqrtl((ldouble)(x)))

#define Abs(x)		(fabsl((ldouble)(x)))

#define Exp(x)		(expl((ldouble)(x)))
#define Log(x)		(logl((ldouble)(x)))

#define Sin(x)		(sinl((ldouble)(x)))
#define Cos(x)		(cosl((ldouble)(x)))
#define Tan(x)		(tanl((ldouble)(x)))

#define ArcSin(x)       (asinl((ldouble)(x)))
#define ArcCos(x)       (acosl((ldouble)(x)))
#define ArcTan(x)       (atanl((ldouble)(x)))

#define Sinh(x)          (sinhl((ldouble)(x)))
#define Cosh(x)          (coshl((ldouble)(x)))
#define Tanh(x)          (tanhl((ldouble)(x)))

#define Cot(x)          (1./tanl((ldouble)(x)))
#define Csc(x)          (1./sinl((ldouble)(x)))




/** Could add definitions for Random(), SeedRandom(), etc. **/



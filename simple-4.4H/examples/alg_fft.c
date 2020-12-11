#include <stdio.h>
#include <math.h>

#include "alg_2dfft.h"

#define REAL double

/*
** FFT and FHT routines
**  Copyright 1988, 1993; Ron Mayer
**  
**  fht(fz,n);
**      Does a hartley transform of "n" points in the array "fz".
**  fft(n,real,imag)
**      Does a fourier transform of "n" points of the "real" and
**      "imag" arrays.
**  ifft(n,real,imag)
**      Does an inverse fourier transform of "n" points of the "real"
**      and "imag" arrays.
**  realfft(n,real)
**      Does a real-valued fourier transform of "n" points of the
**      "real" and "imag" arrays.  The real part of the transform ends
**      up in the first half of the array and the imaginary part of the
**      transform ends up in the second half of the array.
**  realifft(n,real)
**      The inverse of the realfft() routine above.
**      
**      
** NOTE: This routine uses at least 2 patented algorithms, and may be
**       under the restrictions of a bunch of different organizations.
**       Although I wrote it completely myself; it is kind of a derivative
**       of a routine I once authored and released under the GPL, so it
**       may fall under the free software foundation's restrictions;
**       it was worked on as a Stanford Univ project, so they claim
**       some rights to it; it was further optimized at work here, so
**       I think this company claims parts of it.  The patents are
**       held by R. Bracewell (the FHT algorithm) and O. Buneman (the
**       trig generator), both at Stanford Univ.
**       If it were up to me, I'd say go do whatever you want with it;
**       but it would be polite to give credit to the following people
**       if you use this anywhere:
**           Euler     - probable inventor of the fourier transform.
**           Gauss     - probable inventor of the FFT.
**           Hartley   - probable inventor of the hartley transform.
**           Buneman   - for a really cool trig generator
**           Mayer(me) - for authoring this particular version and
**                       including all the optimizations in one package.
**       Thanks,
**       Ron Mayer; mayer@acuson.com
**
*/


#define FHT_SWAP(a,b,t) {(t)=(a);(a)=(b);(b)=(t);}
#define TRIG_VARS                                                \
      int t_lam=0;
#define TRIG_INIT(k,c,s)                                         \
     {                                                           \
      int i;                                                     \
      for (i=2 ; i<=k ; i++)                                     \
	  {coswrk[i]=costab[i];sinwrk[i]=sintab[i];}             \
      t_lam = 0;                                                 \
      c = 1;                                                     \
      s = 0;                                                     \
     }
#define TRIG_NEXT(k,c,s)                                         \
     {                                                           \
	 int i,j;                                                \
	 (t_lam)++;                                              \
	 for (i=0 ; !((1<<i)&t_lam) ; i++);                      \
	 i = k-i;                                                \
	 s = sinwrk[i];                                          \
	 c = coswrk[i];                                          \
	 if (i>1)                                                \
	    {                                                    \
	     for (j=k-i+2 ; (1<<j)&t_lam ; j++);                 \
	     j         = k - j;                                  \
	     sinwrk[i] = halsec[i] * (sinwrk[i-1] + sinwrk[j]);  \
	     coswrk[i] = halsec[i] * (coswrk[i-1] + coswrk[j]);  \
	    }                                                    \
     }
#define TRIG_RESET(k,c,s)

static REAL halsec[20]=
    {
     0,
     0,
     .54119610014619698439972320536638942006107206337801,
     .50979557910415916894193980398784391368261849190893,
     .50241928618815570551167011928012092247859337193963,
     .50060299823519630134550410676638239611758632599591,
     .50015063602065098821477101271097658495974913010340,
     .50003765191554772296778139077905492847503165398345,
     .50000941253588775676512870469186533538523133757983,
     .50000235310628608051401267171204408939326297376426,
     .50000058827484117879868526730916804925780637276181,
     .50000014706860214875463798283871198206179118093251,
     .50000003676714377807315864400643020315103490883972,
     .50000000919178552207366560348853455333939112569380,
     .50000000229794635411562887767906868558991922348920,
     .50000000057448658687873302235147272458812263401372
    };
static REAL costab[20]=
    {
     .00000000000000000000000000000000000000000000000000,
     .70710678118654752440084436210484903928483593768847,
     .92387953251128675612818318939678828682241662586364,
     .98078528040323044912618223613423903697393373089333,
     .99518472667219688624483695310947992157547486872985,
     .99879545620517239271477160475910069444320361470461,
     .99969881869620422011576564966617219685006108125772,
     .99992470183914454092164649119638322435060646880221,
     .99998117528260114265699043772856771617391725094433,
     .99999529380957617151158012570011989955298763362218,
     .99999882345170190992902571017152601904826792288976,
     .99999970586288221916022821773876567711626389934930,
     .99999992646571785114473148070738785694820115568892,
     .99999998161642929380834691540290971450507605124278,
     .99999999540410731289097193313960614895889430318945,
     .99999999885102682756267330779455410840053741619428
    };
static REAL sintab[20]=
    {
     1.0000000000000000000000000000000000000000000000000,
     .70710678118654752440084436210484903928483593768846,
     .38268343236508977172845998403039886676134456248561,
     .19509032201612826784828486847702224092769161775195,
     .09801714032956060199419556388864184586113667316749,
     .04906767432741801425495497694268265831474536302574,
     .02454122852291228803173452945928292506546611923944,
     .01227153828571992607940826195100321214037231959176,
     .00613588464915447535964023459037258091705788631738,
     .00306795676296597627014536549091984251894461021344,
     .00153398018628476561230369715026407907995486457522,
     .00076699031874270452693856835794857664314091945205,
     .00038349518757139558907246168118138126339502603495,
     .00019174759731070330743990956198900093346887403385,
     .00009587379909597734587051721097647635118706561284,
     .00004793689960306688454900399049465887274686668768
    };
static REAL coswrk[20]=
    {
     .00000000000000000000000000000000000000000000000000,
     .70710678118654752440084436210484903928483593768847,
     .92387953251128675612818318939678828682241662586364,
     .98078528040323044912618223613423903697393373089333,
     .99518472667219688624483695310947992157547486872985,
     .99879545620517239271477160475910069444320361470461,
     .99969881869620422011576564966617219685006108125772,
     .99992470183914454092164649119638322435060646880221,
     .99998117528260114265699043772856771617391725094433,
     .99999529380957617151158012570011989955298763362218,
     .99999882345170190992902571017152601904826792288976,
     .99999970586288221916022821773876567711626389934930,
     .99999992646571785114473148070738785694820115568892,
     .99999998161642929380834691540290971450507605124278,
     .99999999540410731289097193313960614895889430318945,
     .99999999885102682756267330779455410840053741619428
    };
static REAL sinwrk[20]=
    {
     1.0000000000000000000000000000000000000000000000000,
     .70710678118654752440084436210484903928483593768846,
     .38268343236508977172845998403039886676134456248561,
     .19509032201612826784828486847702224092769161775195,
     .09801714032956060199419556388864184586113667316749,
     .04906767432741801425495497694268265831474536302574,
     .02454122852291228803173452945928292506546611923944,
     .01227153828571992607940826195100321214037231959176,
     .00613588464915447535964023459037258091705788631738,
     .00306795676296597627014536549091984251894461021344,
     .00153398018628476561230369715026407907995486457522,
     .00076699031874270452693856835794857664314091945205,
     .00038349518757139558907246168118138126339502603495,
     .00019174759731070330743990956198900093346887403385,
     .00009587379909597734587051721097647635118706561284,
     .00004793689960306688454900399049465887274686668768
    };

char fht_version[] = "Brcwl-Hrtly-Ron-dbld";

#define SQRT2_2   0.70710678118654752440084436210484
#define SQRT2   2*0.70710678118654752440084436210484

void
fht(fz,n)
int n;
REAL *fz;
{
 REAL a,b;
 REAL c1,s1,s2,c2,s3,c3,s4,c4;
 REAL f0,g0,f1,g1,f2,g2,f3,g3;
 int i,k,k1,k2,k3,k4,kx;
 REAL *fi,*fn,*gi;
 TRIG_VARS;

 for (k1=1,k2=0;k1<n;k1++)
  {
     REAL a;
     for (k=n>>1; (!((k2^=k)&k)); k>>=1);
     if (k1>k2)
	{
	     a=fz[k1];fz[k1]=fz[k2];fz[k2]=a;
	}
    }
 for ( k=0 ; (1<<k)<n ; k++ );
 k  &= 1;
 if (k==0)
    {
	 for (fi=fz,fn=fz+n;fi<fn;fi+=4)
	    {
	     REAL f0,f1,f2,f3;
	     f1     = fi[0 ]-fi[1 ];
	     f0     = fi[0 ]+fi[1 ];
	     f3     = fi[2 ]-fi[3 ];
	     f2     = fi[2 ]+fi[3 ];
	     fi[2 ] = (f0-f2);  
	     fi[0 ] = (f0+f2);
	     fi[3 ] = (f1-f3);  
	     fi[1 ] = (f1+f3);
	    }
    }
 else
    {
	 for (fi=fz,fn=fz+n,gi=fi+1;fi<fn;fi+=8,gi+=8)
	    {
	     REAL s1,c1,s2,c2,s3,c3,s4,c4,g0,f0,f1,g1,f2,g2,f3,g3;
	     c1     = fi[0 ] - gi[0 ];
	     s1     = fi[0 ] + gi[0 ];
	     c2     = fi[2 ] - gi[2 ];
	     s2     = fi[2 ] + gi[2 ];
	     c3     = fi[4 ] - gi[4 ];
	     s3     = fi[4 ] + gi[4 ];
	     c4     = fi[6 ] - gi[6 ];
	     s4     = fi[6 ] + gi[6 ];
	     f1     = (s1 - s2);        
	     f0     = (s1 + s2);
	     g1     = (c1 - c2);        
	     g0     = (c1 + c2);
	     f3     = (s3 - s4);        
	     f2     = (s3 + s4);
	     g3     = SQRT2*c4;         
	     g2     = SQRT2*c3;
	     fi[4 ] = f0 - f2;
	     fi[0 ] = f0 + f2;
	     fi[6 ] = f1 - f3;
	     fi[2 ] = f1 + f3;
	     gi[4 ] = g0 - g2;
	     gi[0 ] = g0 + g2;
	     gi[6 ] = g1 - g3;
	     gi[2 ] = g1 + g3;
	    }
    }
 if (n<16) return;

 do
    {
     REAL s1,c1;
     k  += 2;
     k1  = 1  << k;
     k2  = k1 << 1;
     k4  = k2 << 1;
     k3  = k2 + k1;
     kx  = k1 >> 1;
	 fi  = fz;
	 gi  = fi + kx;
	 fn  = fz + n;
	 do
	    {
	     REAL g0,f0,f1,g1,f2,g2,f3,g3;
	     f1      = fi[0 ] - fi[k1];
	     f0      = fi[0 ] + fi[k1];
	     f3      = fi[k2] - fi[k3];
	     f2      = fi[k2] + fi[k3];
	     fi[k2]  = f0         - f2;
	     fi[0 ]  = f0         + f2;
	     fi[k3]  = f1         - f3;
	     fi[k1]  = f1         + f3;
	     g1      = gi[0 ] - gi[k1];
	     g0      = gi[0 ] + gi[k1];
	     g3      = SQRT2  * gi[k3];
	     g2      = SQRT2  * gi[k2];
	     gi[k2]  = g0         - g2;
	     gi[0 ]  = g0         + g2;
	     gi[k3]  = g1         - g3;
	     gi[k1]  = g1         + g3;
	     gi     += k4;
	     fi     += k4;
	    } while (fi<fn);
     TRIG_INIT(k,c1,s1);
     for (i=1;i<kx;i++)
	{
	 REAL c2,s2;
	 TRIG_NEXT(k,c1,s1);
	 c2 = c1*c1 - s1*s1;
	 s2 = 2*(c1*s1);
	     fn = fz + n;
	     fi = fz +i;
	     gi = fz +k1-i;
	     do
		{
		 REAL a,b,g0,f0,f1,g1,f2,g2,f3,g3;
		 b       = s2*fi[k1] - c2*gi[k1];
		 a       = c2*fi[k1] + s2*gi[k1];
		 f1      = fi[0 ]    - a;
		 f0      = fi[0 ]    + a;
		 g1      = gi[0 ]    - b;
		 g0      = gi[0 ]    + b;
		 b       = s2*fi[k3] - c2*gi[k3];
		 a       = c2*fi[k3] + s2*gi[k3];
		 f3      = fi[k2]    - a;
		 f2      = fi[k2]    + a;
		 g3      = gi[k2]    - b;
		 g2      = gi[k2]    + b;
		 b       = s1*f2     - c1*g3;
		 a       = c1*f2     + s1*g3;
		 fi[k2]  = f0        - a;
		 fi[0 ]  = f0        + a;
		 gi[k3]  = g1        - b;
		 gi[k1]  = g1        + b;
		 b       = c1*g2     - s1*f3;
		 a       = s1*g2     + c1*f3;
		 gi[k2]  = g0        - a;
		 gi[0 ]  = g0        + a;
		 fi[k3]  = f1        - b;
		 fi[k1]  = f1        + b;
		 gi     += k4;
		 fi     += k4;
		} while (fi<fn);
	}
     TRIG_RESET(k,c1,s1);
    } while (k4<n);
}

void
fht_c(n, fz, plane)
int n;
complex_t *fz;
int plane;
{
 REAL a,b;
 REAL c1,s1,s2,c2,s3,c3,s4,c4;
 REAL f0,g0,f1,g1,f2,g2,f3,g3;
 int i,k,k1,k2,k3,k4,kx;
 complex_t *fi,*fn,*gi;
 TRIG_VARS;

 for (k1=1,k2=0;k1<n;k1++) {
   REAL a;
   for (k=n>>1; (!((k2^=k)&k)); k>>=1);
   if (k1>k2) {
     if (plane==0) {
       a=fz[k1].r;fz[k1].r=fz[k2].r;fz[k2].r=a;
     }
     else {
       a=fz[k1].i;fz[k1].i=fz[k2].i;fz[k2].i=a;
     }
   }
 }
 for ( k=0 ; (1<<k)<n ; k++ );
 k  &= 1;
 if (k==0) {
   for (fi=fz,fn=fz+n;fi<fn;fi+=4) {
     REAL f0,f1,f2,f3;
     if (plane==0) {
       f1     = fi[0 ].r-fi[1 ].r;
       f0     = fi[0 ].r+fi[1 ].r;
       f3     = fi[2 ].r-fi[3 ].r;
       f2     = fi[2 ].r+fi[3 ].r;
       fi[2 ].r = (f0-f2);  
       fi[0 ].r = (f0+f2);
       fi[3 ].r = (f1-f3);  
       fi[1 ].r = (f1+f3);
     }
     else {
       f1     = fi[0 ].i-fi[1 ].i;
       f0     = fi[0 ].i+fi[1 ].i;
       f3     = fi[2 ].i-fi[3 ].i;
       f2     = fi[2 ].i+fi[3 ].i;
       fi[2 ].i = (f0-f2);  
       fi[0 ].i = (f0+f2);
       fi[3 ].i = (f1-f3);  
       fi[1 ].i = (f1+f3);
     }
   }
 }
 else {
   for (fi=fz,fn=fz+n,gi=fi+1;fi<fn;fi+=8,gi+=8) {
     REAL s1,c1,s2,c2,s3,c3,s4,c4,g0,f0,f1,g1,f2,g2,f3,g3;
     if (plane==0) {
       c1     = fi[0 ].r - gi[0 ].r;
       s1     = fi[0 ].r + gi[0 ].r;
       c2     = fi[2 ].r - gi[2 ].r;
       s2     = fi[2 ].r + gi[2 ].r;
       c3     = fi[4 ].r - gi[4 ].r;
       s3     = fi[4 ].r + gi[4 ].r;
       c4     = fi[6 ].r - gi[6 ].r;
       s4     = fi[6 ].r + gi[6 ].r;
       f1     = (s1 - s2);        
       f0     = (s1 + s2);
       g1     = (c1 - c2);        
       g0     = (c1 + c2);
       f3     = (s3 - s4);        
       f2     = (s3 + s4);
       g3     = SQRT2*c4;         
       g2     = SQRT2*c3;
       fi[4 ].r = f0 - f2;
       fi[0 ].r = f0 + f2;
       fi[6 ].r = f1 - f3;
       fi[2 ].r = f1 + f3;
       gi[4 ].r = g0 - g2;
       gi[0 ].r = g0 + g2;
       gi[6 ].r = g1 - g3;
       gi[2 ].r = g1 + g3;
     }
     else {
       c1     = fi[0 ].i - gi[0 ].i;
       s1     = fi[0 ].i + gi[0 ].i;
       c2     = fi[2 ].i - gi[2 ].i;
       s2     = fi[2 ].i + gi[2 ].i;
       c3     = fi[4 ].i - gi[4 ].i;
       s3     = fi[4 ].i + gi[4 ].i;
       c4     = fi[6 ].i - gi[6 ].i;
       s4     = fi[6 ].i + gi[6 ].i;
       f1     = (s1 - s2);        
       f0     = (s1 + s2);
       g1     = (c1 - c2);        
       g0     = (c1 + c2);
       f3     = (s3 - s4);        
       f2     = (s3 + s4);
       g3     = SQRT2*c4;         
       g2     = SQRT2*c3;
       fi[4 ].i = f0 - f2;
       fi[0 ].i = f0 + f2;
       fi[6 ].i = f1 - f3;
       fi[2 ].i = f1 + f3;
       gi[4 ].i = g0 - g2;
       gi[0 ].i = g0 + g2;
       gi[6 ].i = g1 - g3;
       gi[2 ].i = g1 + g3;
     }
   }
 }
 if (n<16) return;

 do {
   REAL s1,c1;
   k  += 2;
   k1  = 1  << k;
   k2  = k1 << 1;
   k4  = k2 << 1;
   k3  = k2 + k1;
   kx  = k1 >> 1;
   fi  = fz;
   gi  = fi + kx;
   fn  = fz + n;
   do {
     REAL g0,f0,f1,g1,f2,g2,f3,g3;
     if (plane==0) {
       f1      = fi[0 ].r - fi[k1].r;
       f0      = fi[0 ].r + fi[k1].r;
       f3      = fi[k2].r - fi[k3].r;
       f2      = fi[k2].r + fi[k3].r;
       fi[k2].r  = f0         - f2;
       fi[0 ].r  = f0         + f2;
       fi[k3].r  = f1         - f3;
       fi[k1].r  = f1         + f3;
       g1      = gi[0 ].r - gi[k1].r;
       g0      = gi[0 ].r + gi[k1].r;
       g3      = SQRT2  * gi[k3].r;
       g2      = SQRT2  * gi[k2].r;
       gi[k2].r  = g0         - g2;
       gi[0 ].r  = g0         + g2;
       gi[k3].r  = g1         - g3;
       gi[k1].r  = g1         + g3;
       gi     += k4;
       fi     += k4;
     }
     else {
       f1      = fi[0 ].i - fi[k1].i;
       f0      = fi[0 ].i + fi[k1].i;
       f3      = fi[k2].i - fi[k3].i;
       f2      = fi[k2].i + fi[k3].i;
       fi[k2].i  = f0         - f2;
       fi[0 ].i  = f0         + f2;
       fi[k3].i  = f1         - f3;
       fi[k1].i  = f1         + f3;
       g1      = gi[0 ].i - gi[k1].i;
       g0      = gi[0 ].i + gi[k1].i;
       g3      = SQRT2  * gi[k3].i;
       g2      = SQRT2  * gi[k2].i;
       gi[k2].i  = g0         - g2;
       gi[0 ].i  = g0         + g2;
       gi[k3].i  = g1         - g3;
       gi[k1].i  = g1         + g3;
       gi     += k4;
       fi     += k4;
     }
   } while (fi<fn);
   TRIG_INIT(k,c1,s1);
   for (i=1;i<kx;i++) {
     REAL c2,s2;
     TRIG_NEXT(k,c1,s1);
     c2 = c1*c1 - s1*s1;
     s2 = 2*(c1*s1);
     fn = fz + n;
     fi = fz +i;
     gi = fz +k1-i;
     do {
       REAL a,b,g0,f0,f1,g1,f2,g2,f3,g3;
       if (plane ==0) {
	 b       = s2*fi[k1].r - c2*gi[k1].r;
	 a       = c2*fi[k1].r + s2*gi[k1].r;
	 f1      = fi[0 ].r    - a;
	 f0      = fi[0 ].r    + a;
	 g1      = gi[0 ].r    - b;
	 g0      = gi[0 ].r    + b;
	 b       = s2*fi[k3].r - c2*gi[k3].r;
	 a       = c2*fi[k3].r + s2*gi[k3].r;
	 f3      = fi[k2].r    - a;
	 f2      = fi[k2].r    + a;
	 g3      = gi[k2].r    - b;
	 g2      = gi[k2].r    + b;
	 b       = s1*f2     - c1*g3;
	 a       = c1*f2     + s1*g3;
	 fi[k2].r  = f0        - a;
	 fi[0 ].r  = f0        + a;
	 gi[k3].r  = g1        - b;
	 gi[k1].r  = g1        + b;
	 b       = c1*g2     - s1*f3;
	 a       = s1*g2     + c1*f3;
	 gi[k2].r  = g0        - a;
	 gi[0 ].r  = g0        + a;
	 fi[k3].r  = f1        - b;
	 fi[k1].r  = f1        + b;
	 gi     += k4;
	 fi     += k4;
       }
       else {
	 b       = s2*fi[k1].i - c2*gi[k1].i;
	 a       = c2*fi[k1].i + s2*gi[k1].i;
	 f1      = fi[0 ].i    - a;
	 f0      = fi[0 ].i    + a;
	 g1      = gi[0 ].i    - b;
	 g0      = gi[0 ].i    + b;
	 b       = s2*fi[k3].i - c2*gi[k3].i;
	 a       = c2*fi[k3].i + s2*gi[k3].i;
	 f3      = fi[k2].i    - a;
	 f2      = fi[k2].i    + a;
	 g3      = gi[k2].i    - b;
	 g2      = gi[k2].i    + b;
	 b       = s1*f2     - c1*g3;
	 a       = c1*f2     + s1*g3;
	 fi[k2].i  = f0        - a;
	 fi[0 ].i  = f0        + a;
	 gi[k3].i  = g1        - b;
	 gi[k1].i  = g1        + b;
	 b       = c1*g2     - s1*f3;
	 a       = s1*g2     + c1*f3;
	 gi[k2].i  = g0        - a;
	 gi[0 ].i  = g0        + a;
	 fi[k3].i  = f1        - b;
	 fi[k1].i  = f1        + b;
	 gi     += k4;
	 fi     += k4;
       }
     } while (fi<fn);
   }
   TRIG_RESET(k,c1,s1);
 } while (k4<n);
}

void
ifft(n,real,imag)
int n;
double *real,*imag;
{
 double a,b,c,d;
 double q,r,s,t;
 int i,j,k;
 fht(real,n);
 fht(imag,n);
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i]; b = real[j];  q=a+b; r=a-b;
  c = imag[i]; d = imag[j];  s=c+d; t=c-d;
  imag[i] = (s+r)*0.5;  imag[j] = (s-r)*0.5;
  real[i] = (q-t)*0.5;  real[j] = (q+t)*0.5;
 }
}

void
realfft(n,real)
int n;
double *real;
{
 double a,b,c,d;
 int i,j,k;
 fht(real,n);
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i];
  b = real[j];
  real[j] = (a-b)*0.5;
  real[i] = (a+b)*0.5;
 }
}

void
fft(n,real,imag)
int n;
double *real,*imag;
{
 double a,b,c,d;
 double q,r,s,t;
 int i,j,k;
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i]; b = real[j];  q=a+b; r=a-b;
  c = imag[i]; d = imag[j];  s=c+d; t=c-d;
  real[i] = (q+t)*.5; real[j] = (q-t)*.5;
  imag[i] = (s-r)*.5; imag[j] = (s+r)*.5;
 }
 fht(real,n);
 fht(imag,n);
}

void
fft_complex(n,data)
int n;
complex_t *data;
{
 double a,b,c,d;
 double q,r,s,t;
 int i,j,k;
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = data[i].r; b = data[j].r;  q=a+b; r=a-b;
  c = data[i].i; d = data[j].i;  s=c+d; t=c-d;
  data[i].r = (q+t)*.5; data[j].r = (q-t)*.5;
  data[i].i = (s-r)*.5; data[j].i = (s+r)*.5;
 }
 fht_c(n, data, 0);
 fht_c(n, data, 1);
}

void
realifft(n,real)
int n;
double *real;
{
 double a,b,c,d;
 int i,j,k;
 for (i=1,j=n-1,k=n/2;i<k;i++,j--) {
  a = real[i];
  b = real[j];
  real[j] = (a-b);
  real[i] = (a+b);
 }
 fht(real,n);
}



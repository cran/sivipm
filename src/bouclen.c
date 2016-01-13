/* ---------------------------------------------------------------
  sivipm R package
  Copyright INRA 2016
  INRA, UR1404, Research Unit MaIAGE
  F78352 Jouy-en-Josas, France.
 
  URL: http://cran.r-project.org/web/packages/sivipm
 
  This file is part of sivipm R package.
  sivipm is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  See the GNU General Public License at:
  http://www.gnu.org/licenses/
 
----------- --------------------------------------------------- */

#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> // to allow user interrupt
#include "string.h" // for memcpy

/* ---------------------------------------------------------------
 Programmation in C of the regression PLS without missing values
 NOTE: in C, the values are stored by column
----------- --------------------------------------------------- */


/* multpp
FUNCTION
 Multiplication of a  matrix t(mat1) by the first column of mat2
with excluding row number 'indice' of mat1.
Elements of mat2 are divided by 'somme'.
 mat1 should have one more row as mat2.
INPUT
mat1: nrow X ncol1
mat2: (nrow-1) X ignored
OUTPUT
res:  ncol1 X 1
RETURN
sum(res**2)
*/
  double multpp(double *mat1, double *mat2, double somme,
	      int nrow, int ncol1,
	      int indice,	double *res) {
  int i,j,l=0 ;
  double som2=0.0;

  for (j=0; j< ncol1; j++) {
    res[j]=0.0;
    l=0;
    for (i=0; i< nrow; i++) {
      if (i != indice) {
	res[j] += (mat1[nrow * j+i] * (mat2[l++]/somme));
      }
    }
    som2 += (res[j]*res[j]);
  } /* fin j */
  return(som2);
} /* fin multpp */
/* --------------------------------------------------- */
/* multpps
FUNCTION
 Multiplication of a  matrix t(mat1) by the first column of mat2
Elements of mat2 are divided by 'somme'.
 mat1 should have the same number of rows as mat2.
Same as 'multpp', but with all the rows of mat1.
INPUT
mat1: nrow X ncol1
mat2: nrow X ignored
OUTPUT
res:  ncol1 X 1
RETURN
sum(res**2)
*/
  double multpps(double *mat1, double *mat2, double somme,
	      int nrow, int ncol1,
	      double *res) {
  int i,j,l=0 ;
  double som2=0.0;

  for (j=0; j< ncol1; j++) {
    res[j]=0.0;
    l=0;
    for (i=0; i< nrow; i++) {
	res[j] += (mat1[nrow * j+i] * (mat2[l++]/somme));
    }
    som2 += (res[j]*res[j]);
  } /* fin j */
  return(som2);
} /* fin multpps */
/* --------------------------------------------------- */
/* multppsc
FUNCTION
 Multiplication of a  matrix t(mat1) by the first column of mat2
Elements of mat2 are divided by 'somme'.
 mat1 should have the same number of rows as mat2.
Same as 'multpps', but return nothing
INPUT
mat1: nrow X ncol1
mat2: nrow X ignored
OUTPUT
res:  ncol1 X 1
*/
  void multppsc(double *mat1, double *mat2, double somme,
	      int nrow, int ncol1,
	      double *res) {
  int i,j,l=0 ;

  for (j=0; j< ncol1; j++) {
    res[j]=0.0;
    l=0;
    for (i=0; i< nrow; i++) {
	res[j] += (mat1[nrow * j+i] * (mat2[l++]/somme));
    }
  } /* fin j */

} /* fin multppsc */
/* --------------------------------------------------- */
/* colSums
FUNCTION
Add the elements of each column of a matrix
INPUT
mat: nrow X ncol
OUTPUT
res: ncol
*/

void colSums(double *mat,  int nrow, int ncol, 	double *res) {
  int i,j;

  for (j=0; j < ncol; j++) {
    res[j] = 0.0;
    for (i=0; i< nrow; i++) {
      res[j] += mat[j*nrow + i];
    }
  }
} /* fin colSums */


/* --------------------------------------------------- */
/*  multmat
FUNCTION
Multiplication of a  matrix (mat1) by the first column of t(mat2)
with excluding the row number 'indice' of mat1.
The number of columns of mat1 should be equal to the number 
of rows of mat2.
INPUT
mat1: nrow X ncol
mat2: nrow X ignored
OUTPUT
res:  ncol1 X 1
RETURN
somme(res*res)
*/

  double multmat(double *mat1, double *mat2,
	      int nrow, int ncol,
	      int indice,	double *res) {
  int i,j,l=0 ;
  double som2=0.0;

  for (i=0;  i<  nrow; i++) {
     if (i != indice) {
       res[l] =0.0;
       for (j=0; j< ncol; j++) {
	 /*	 Rprintf(" %g, %g ;", mat1[nrow * j+i] , mat2[j]); */
	 res[l]  += (mat1[nrow * j+i] * mat2[j]);
       }/* fin j */
       /*       Rprintf(" R[%d]=%g\n", l, res[l]); */
       som2 += (res[l]*res[l]);
       l++;
     } /* fin (i != indice) */
  } /* fin i */
  return(som2);
} /* fin multmat */




/* --------------------------------------------------- */
/* multmat2
Same as multmat, where mat2 is divided by 'somme'.
*/

double multmat2(double *mat1, double *mat2, double somme, 
	      int nrow, int ncol,
	      int indice,	double *res) {
  int i,j,l=0 ;
  double som2=0.0;
 
  for (i=0;  i<  nrow; i++) {
     if (i != indice) {
       res[l] =0.0;
       for (j=0; j< ncol; j++) {
	 /*	 Rprintf(" %g, %g ;", mat1[nrow * j+i] , mat2[j]); */
	 res[l]  += (mat1[nrow * j+i] * (mat2[j]/somme));
       }/* fin j */
       /*       Rprintf(" R[%d]=%g\n", l, res[l]); */
       som2 += ( res[l]* res[l]);
       l++;
     } /* fin (i != indice) */
  } /* fin i */
  return(som2);

} /* fin multmat2 */

/* --------------------------------------------------- */
/*
bouclen
FUNCTION
Code the loop over n in the function R reglps2nomissing
INPUT
Xold: matrix nlig X ncolX
YYold: matrix nlig X ncolYY
nlig: number of rows of Xold ans YYold
ncolX: number of columns of Xold
ncolYY: number of columns of YYold
WORKING
whsi, whsiold, wsidif:  vectors of length  ncolX 
thsi, uhsi:  vectors of length nlig-1
chsi, YYhatsi :  vectors of length ncolYY
OUTPUT
pression: matrix nlig X  ncolYY
*/

 void bouclen( double *Xold, double * YYold, 
	       int *nlig, int *ncolX,
	       int *ncolYY,
	       double *whsi, double *whsiold, double *wsidif,
	       double *thsi, double *uhsi,
	       double *chsi, double *YYhatsi,
	       double *pression) {
   int i,j,l, indice;
   int iloop; 
   double somuhsi2, somwhsi2, somthsi2, somchsi2, somwsidif2, som, a;


  /*
  The allocations are done before the call
  whsi = (double *) S_alloc( *ncolX,  sizeof(double));
  whsiold = (double *) S_alloc( *ncolX,  sizeof(double));
  wsidif = (double *) S_alloc( *ncolX,  sizeof(double));
  thsi = (double *) S_alloc( (*nlig-1),  sizeof(double));
  uhsi = (double *) S_alloc( (*nlig-1),  sizeof(double));
  chsi  = (double *) S_alloc( *ncolYY,  sizeof(double));
  YYhatsi =  (double *) S_alloc( *ncolYY,  sizeof(double));
  pression = (double *) S_alloc( (*nlig) *(*ncolYY),  sizeof(double));
  */


  for (indice=0;  indice <*nlig; indice++) {
     for (i=0; i<*ncolX; i++) {
       whsiold[i] =1.0;
     }
     l=0;
     somuhsi2 = 0.0;
     for (i=0; i< *nlig; i++) {
       if (i != indice) {
	 somuhsi2 += (YYold[i]*YYold[i]);
	 uhsi[l++] = YYold[i];
       }
     }

     iloop=0;
    while(1) {
      iloop++; /* pour s'assurer que le while se termine */
      /* whsi = t(Xold[-indice,])*(uhsi/somuhsi2) */
      somwhsi2 = multpp(Xold, uhsi, somuhsi2, *nlig, *ncolX, indice, whsi);


      /* somwhsi2 = sqrt(somme(whsi**2)) */
      somwhsi2 = R_pow(somwhsi2, (double)0.5);

      for (i=0; i<*ncolX; i++) {
	whsi[i] /= somwhsi2;
      }

      somthsi2=multmat(Xold, whsi, *nlig, *ncolX, indice, thsi);

      /* somthsi2=somme(thsi**2) */
      /* chsi = t(YYold[-indice,])*(thsi/somthsi2) */
      /* somchsi2=somme(chsi**2) */
      somchsi2=multpp(YYold, thsi, somthsi2, *nlig, *ncolYY, indice, chsi);

      /* uhsi = YYold[-indice,] * (chsi/somchsi2) */
      somuhsi2 = multmat2(YYold, chsi, somchsi2, *nlig, *ncolYY, indice, uhsi);


      /* wsidif = whsi - whsiold; whsiold=whsi */
      somwsidif2 = 0.0;
      for (i=0; i< *ncolX; i++) {
	wsidif[i] =  whsi[i] - whsiold[i];
	somwsidif2 += (wsidif[i]*wsidif[i]);
	whsiold[i] = whsi[i];
      }

      if ( (somwsidif2 < 1.0e-12)|| (iloop >3000))	{
	break;
      }

      som=0.0;
      for (i=0; i< *ncolX; i++) {
	 som += (Xold[ (*nlig) *i + indice]* whsi[i]);
      }

      for (i=0; i < *ncolYY; i++) {
	 YYhatsi[i] = chsi[i] * som;
	 a= YYold[ (*nlig)*i+indice] - YYhatsi[i];
	 pression[ (*nlig)*i+indice] = a*a;
      }
      R_CheckUserInterrupt(); // check User interrupt


    } /* fin while */

  } /* fin indice */

 
} /* fin bouclen */

/* --------------------------------------------------- */
/* boucler
FUNCTION
Code the inner repeat-loop of  regpls2nomissing
INPUT
Xold: matrix nlig X ncolX
YYold: matrix nlig X ncolYY
unew: vector of length nlig
wold:  vector of length  ncolX 
OUTPUT
tnew: vector of length nlig
RETURN
somtnew2: sum(tnew**2)
*/

double boucler (double *Xold, double * YYold,
	      double *unew, double *wold,
	      int *nlig, int *ncolX,
	      int *ncolYY,
		double *wnew, double  *wdif, double  *cnew,
	      double *tnew) {

  int i, indice=-1;
  double somwnew2, somcnew2, somunew2, somtnew2, somwdif2;

  int iloop=0;

  /*
  The allocations are done before the call
  wnew= (double *) S_alloc( *ncolX,  sizeof(double));
  wdif= (double *) S_alloc( *ncolX,  sizeof(double));
  tnew= (double *) S_alloc( *nlig,  sizeof(double)); 
  cnew = (double *) S_alloc( *ncolYY,  sizeof(double));
*/
  somunew2= 0.0;
  for (i=0; i<*nlig; i++) {
    somunew2 += (unew[i] * unew[i]);
  }

  while (1) {
    iloop++; /* pour s'assurer que le while se termine */
  /* w.new <- t(X.old) %*% u.new/sum(u.new^2) */

      somwnew2 = multpps(Xold, unew, somunew2, *nlig, *ncolX,   wnew);
  /* w.new <- w.new/sqrt(sum(w.new^2)) */
  somwnew2 = R_pow(somwnew2, (double)0.5);

  for (i=0; i<*ncolX; i++) {
	wnew[i] /= somwnew2;
      }
  /* t.new <- X.old %*% w.new */
  somtnew2 = multmat(Xold, wnew,  *nlig, *ncolX, indice, tnew);

  /*  c.new <- t(YY.old) %*% t.new/sum(t.new^2) */
  somcnew2 = multpps(YYold, tnew, somtnew2, *nlig, *ncolYY,  cnew);

  /* u.new <- YY.old %*% c.new/sum(c.new^2) */
  somunew2 =  multmat2(YYold, cnew, somcnew2, *nlig, *ncolYY, indice, unew);

  /* w.dif <- w.new - w.old ; w.old <- w.new */
  somwdif2=0.0;
  for (i=0; i< *ncolX; i++) {
    wdif[i] = wnew[i] - wold[i];
    somwdif2 += ( wdif[i]* wdif[i]);
    wold[i] = wnew[i];
  }


  if ( (somwdif2 < 1.0e-12)|| (iloop >3000))	{ 
    break;
  }
  R_CheckUserInterrupt(); // check User interrupt
  } /* fin while */
  return(somtnew2);
} /* fin boucler */
/* --------------------------------------------------- */
/* multsom
FUNCTION
Compute:
RSS[, h + 1 ] <- colSums((YY.old - t.new %*% t(c.new))^2)
YYold: matrix nlig X ncolYY
tnew : vector of length nlig
cnew :  vector of length ncolYY
nlig: number of rows of YYold
ncolYY: number of columns of YYold
nc: number of components
h: current component= current row in RSS
OUTPUT
(allocated before the call)
RSS: vector of length ncolYY. In fact, pointer to the column h+1 
of the matrix RSS (nc+1 X ncolYY)
*/

    void multsom( double *YYold, double *tnew, double *cnew,
		  int nlig, int ncolYY, 
		  double *RSS) {

     int i,j; 
    double b;


    for (j=0; j< ncolYY; j++) {
      /*     RSS[h+1, j] =0.0; */
      RSS[ j ] = 0.0;
      for (i=0; i < nlig; i++) {
	/* b = YYold[i,j] - ( tnew[i] * cnew[j]); */
	b= YYold[ (nlig) * j +i]  - ( tnew[i] * cnew[j]);
	/* RSS[h+1, j] += (b*b); */
	 RSS[  j ] += (b*b);
      } /* fin i */
    } /* fin j */
    } /* fin multsom */

/* --------------------------------------------------- */
/* matmoins
FUNCTION
Compute: mat - tn %*% t(pn)
INPUT
mat: matrix nlig X ncol
tn: vector nlig
pn: vector ncol
nlig
ncol
RETURN
The result is in mat
*/

void matmoins( double *mat, double *tn, double *pn,
	       int nlig, int ncol) {
  int i,j ;
  for (i=0; i < nlig; i++) {
    for (j=0; j < ncol; j++) {
      mat[nlig*j +i] -= tn[i]*pn[j];
    }
  }
} /* fin matmoins */


/* --------------------------------------------------- */
/* boucle
FUNCTION
The main program: NIPALS loop
INPUT 
Xold: matrix nlig X ncolX
YYold: matrix nlig X ncolYY
nlig: number of rows of Xold ans YYold
ncolX: number of columns of Xold
ncolYY: number of columns of YYold
nc: number of components
WORKING
(allocated before the call)
wold, wdif:  vectors of length  ncolX 
whsi, whsiold, wsidif:  vectors of length  ncolX 
thsi, uhsi:  vectors of length nlig-1
chsi, YYhatsi, cnew :  vectors of length ncolYY
pression: matrix nlig X  ncolYY
OUTPUT
(allocated before the call)
RSS: matrix   nc+1 X ncolYY 
( RSS[, 1] : initialized before the call)
P: matrix ncolX X nc
C, PRESS, Q2:  matrix ncolYY  X nc 
U, T: matrix nlig X nc
W: matrix ncolX X nc

Notations R versus C
n = nlig
p= ncolX
q= ncolYY

RETURN
h : should be equal to nc
------------------------------------------------------- */
int boucle (double *Xold, double * YYold,
	    int *nlig, int *ncolX,
	    int *ncolYY, int *nc,
	    double *wold, double *wnew, double *wdif,
	    double *whsi, double *whsiold, double *wsidif,
            double *thsi, double *uhsi,
	    double *chsi,  double *YYhatsi, 
	    double *pression,
	    double *RSS, double *C, double *PRESS, 
	    double *Q2, double *P, double *W, double *T, 
	    double *U) {

  int i,j, l, h;
  int decalUT, decalPW, decalC;
  double b, somtnew2 = 0.0;


  for (h=0; h < *nc; h++) {
    decalUT = (*nlig)*h; /* gap in the matrix T and U */
    decalPW = (*ncolX)*h; /* gap in the matrix P and W */
    decalC = (*ncolYY)*h; /* gap in the matrix C */

    /* u.new <- YY.old[, 1] */
  memcpy(U+decalUT, YYold, (size_t)(*nlig*sizeof(double)));

    /* w.old <- rep(1, p) */
    for (i=0; i< *ncolX; i++)
      wold[i] = 1.0;
    /* tnew <- rep(0, n) */
    for (i=0; i< *nlig; i++)
      T[decalUT+i] = 0.0;
      /* tnew[i] = 0.0; */

    somtnew2 = boucler(Xold, YYold, U+decalUT,  wold, nlig, ncolX, ncolYY,
		       W+decalPW, wdif, C+decalC,T+decalUT);

    /* p.new <- t(X.old) %*% t.new/sum(t.new^2) */
    multppsc(Xold, T+decalUT, somtnew2, *nlig, *ncolX, P+decalPW);
    /* c.new <- t(YY.old) %*% t.new/sum(t.new^2) */
    multppsc(YYold, T+decalUT, somtnew2, *nlig, *ncolYY, C+decalC);
    /*  RSS[, h + 1 ] <- colSums((YY.old - t.new %*% t(c.new))^2) */


    multsom( YYold, T+decalUT, C+decalC, *nlig, *ncolYY, 
	     (RSS+ (*ncolYY) *(h+1)));


    /* compute pression */   
    bouclen(Xold, YYold, nlig, ncolX, ncolYY,
	    whsi,  whsiold,wsidif,
	       thsi, uhsi,
	       chsi, YYhatsi,
	    pression);


    /*  PRESS[h, ] <- colSums(pression) */
    colSums(pression,  *nlig, *ncolYY, PRESS+(*ncolYY)*h );
    /* Q2[h, ] <- 1 - PRESS[h, ]/RSS[h, ] */
    l= (*ncolYY) * h;
    for (i=0; i< *ncolYY; i++) {
      Q2[l] = 1.0 - PRESS[l]/RSS[l];
      l++;
    }

    /* X.old <- X.old - t.new %*% t(p.new) */
    matmoins( Xold, T+decalUT, P+decalPW, *nlig, *ncolX);
    matmoins( YYold, T+decalUT, C+decalC, *nlig, *ncolYY);

    R_CheckUserInterrupt(); // check User interrupt
 
  } /* fin for h */
  return(h);
} /* fin boucle */


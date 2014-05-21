/*
   FDHS-DSSM.c -- this file is part of sparseHessianFD, a contributed package
   for the R statistical programming platform.

   Please see the LICENSE file for details on the package license


     THIS SOFTWARE IS PROVIDED WITH PERMISSION OF THE ASSOCIATION
     FOR COMPUTING MACHINERY, WHICH HOLDS THE COPYRIGHT.
     PLEASE REVIEW THE TERMS OF THE ACM SOFTWARE LICENSE AGREEMENT,
     WHICH IS IN THE LICENSE FILE OF THE TRUSTOPTIM PACKAGE.

     THE CODE HAS BEEN MODIFIED TO SUPPORT DOUBLE PRECISION TYPES,
     AND HAD BEEN CONVERTED FROM THE ORGINAL FORTRAN USING F2C.

     THANKS TO THE AUTHORS OF THE ORIGINAL ALGORITHM:
     THOMAS F. COLEMAN, BURTON S. GARBOW AND JORGE J. MORE 

     ALGORITHM 636 COLLECTED ALGORITHMS FROM ACM.
     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.11, NO. 4,
     DEC., 1985, P. 378.

*/


#include <stdlib.h>
#define Extern extern

#define GMAX(a,b) ((a) > (b) ? (a):(b))
#define GMIN(a,b) ((a) < (b) ? (a):(b))

/* Table of constant values */

static int c__1 = 1;
static int c_n1 = -1;


/* Subroutine */ 
#ifdef __cplusplus 
extern "C"
#endif
 int dssm_(int *n, int *npairs, int *indrow, 
	int *indcol, int *method, int *listp, int *ngrp, 
	int *maxgrp, int *mingrp, int *info, int *ipntr, 
	int *jpntr, int *iwa, int *liwa)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Local variables */
    static int i__, j, k, jp, ir;
    extern /* Subroutine */ int ido_(int *, int *, int *, int 
	    *, int *, int *, int *, int *, int *, int 
	    *, int *, int *, int *), seq_(int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *), slo_(int *, int *, int *, int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *);
    static int nnz;
    extern /* Subroutine */ int degr_(int *, int *, int *, 
	    int *, int *, int *, int *), idog_(int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *), slog_(int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *, int *, int *), sdpt_(int *, 
	    int *, int *, int *, int *, int *, int *, 
	    int *, int *), setr_(int *, int *, int *, 
	    int *, int *, int *, int *);
    static int maxid, maxvd, maxclq;
    extern /* Subroutine */ int srtdat_(int *, int *, int *, 
	    int *, int *, int *);
    static int numgrp;
    extern /* Subroutine */ int numsrt_(int *, int *, int *, 
	    int *, int *, int *, int *);

/*     ********** */

/*     SUBROUTINE DSSM */

/*     GIVEN THE SPARSITY PATTERN OF A SYMMETRIC MATRIX A OF ORDER N, */
/*     THIS SUBROUTINE DETERMINES A SYMMETRIC PERMUTATION OF A AND A */
/*     PARTITION OF THE COLUMNS OF A CONSISTENT WITH THE DETERMINATION */
/*     OF A BY A LOWER TRIANGULAR SUBSTITUTION METHOD. */

/*     THE SPARSITY PATTERN OF THE MATRIX A IS SPECIFIED BY THE */
/*     ARRAYS INDROW AND INDCOL. ON INPUT THE INDICES FOR THE */
/*     NON-ZERO ELEMENTS IN THE LOWER TRIANGULAR PART OF A ARE */

/*           (INDROW(K),INDCOL(K)), K = 1,2,...,NPAIRS. */

/*     THE (INDROW(K),INDCOL(K)) PAIRS MAY BE SPECIFIED IN ANY ORDER. */
/*     DUPLICATE INPUT PAIRS ARE PERMITTED, BUT THE SUBROUTINE */
/*     ELIMINATES THEM. THE SUBROUTINE REQUIRES THAT ALL THE DIAGONAL */
/*     ELEMENTS BE PART OF THE SPARSITY PATTERN AND REPLACES ANY PAIR */
/*     (INDROW(K),INDCOL(K)) WHERE INDROW(K) IS LESS THAN INDCOL(K) */
/*     BY THE PAIR (INDCOL(K),INDROW(K)). */

/*     THE DIRECT METHOD (METHOD = 1) FIRST DETERMINES A PARTITION */
/*     OF THE COLUMNS OF A SUCH THAT TWO COLUMNS IN A GROUP HAVE A */
/*     NON-ZERO ELEMENT IN ROW K ONLY IF COLUMN K IS IN AN EARLIER */
/*     GROUP. USING THIS PARTITION, THE SUBROUTINE THEN COMPUTES A */
/*     SYMMETRIC PERMUTATION OF A CONSISTENT WITH THE DETERMINATION */
/*     OF A BY A LOWER TRIANGULAR SUBSTITUTION METHOD. */

/*     THE INDIRECT METHOD FIRST COMPUTES A SYMMETRIC PERMUTATION OF A */
/*     WHICH MINIMIZES THE MAXIMUM NUMBER OF NON-ZERO ELEMENTS IN ANY */
/*     ROW OF L, WHERE L IS THE LOWER TRIANGULAR PART OF THE PERMUTED */
/*     MATRIX. THE SUBROUTINE THEN PARTITIONS THE COLUMNS OF L INTO */
/*     GROUPS SUCH THAT COLUMNS OF L IN A GROUP DO NOT HAVE A NON-ZERO */
/*     IN THE SAME ROW POSITION. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE DSSM(N,NPAIRS,INDROW,INDCOL,METHOD,LISTP,NGRP, */
/*                       MAXGRP,MINGRP,INFO,IPNTR,JPNTR,IWA,LIWA) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE ORDER OF A. */

/*       NPAIRS IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF (INDROW,INDCOL) PAIRS USED TO DESCRIBE THE SPARSITY */
/*         PATTERN OF A. */

/*       INDROW IS AN INT ARRAY OF LENGTH NPAIRS. ON INPUT INDROW */
/*         MUST CONTAIN THE ROW INDICES OF THE NON-ZERO ELEMENTS IN */
/*         THE LOWER TRIANGULAR PART OF A. ON OUTPUT INDROW IS */
/*         PERMUTED SO THAT THE CORRESPONDING COLUMN INDICES ARE IN */
/*         NON-DECREASING ORDER. THE COLUMN INDICES CAN BE RECOVERED */
/*         FROM THE ARRAY JPNTR. */

/*       INDCOL IS AN INT ARRAY OF LENGTH NPAIRS. ON INPUT INDCOL */
/*         MUST CONTAIN THE COLUMN INDICES OF THE NON-ZERO ELEMENTS */
/*         IN THE LOWER TRIANGULAR PART OF A. ON OUTPUT INDCOL IS */
/*         PERMUTED SO THAT THE CORRESPONDING ROW INDICES ARE IN */
/*         NON-DECREASING ORDER. THE ROW INDICES CAN BE RECOVERED */
/*         FROM THE ARRAY IPNTR. */

/*       METHOD IS AN INT INPUT VARIABLE. IF METHOD = 1, THE */
/*         DIRECT METHOD IS USED TO DETERMINE THE PARTITION AND */
/*         SYMMETRIC PERMUTATION. OTHERWISE, THE INDIRECT METHOD IS */
/*         USED TO DETERMINE THE SYMMETRIC PERMUTATION AND PARTITION. */

/*       LISTP IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE SYMMETRIC PERMUTATION OF THE MATRIX A. ELEMENT (I,J) */
/*         OF A IS THE (LISTP(I),LISTP(J)) ELEMENT OF THE PERMUTED */
/*         MATRIX. */

/*       NGRP IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE PARTITION OF THE COLUMNS OF A. COLUMN J BELONGS TO */
/*         GROUP NGRP(J). */

/*       MAXGRP IS AN INT OUTPUT VARIABLE WHICH SPECIFIES THE */
/*         NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF A. */

/*       MINGRP IS AN INT OUTPUT VARIABLE WHICH SPECIFIES A LOWER */
/*         BOUND FOR THE NUMBER OF GROUPS IN ANY PARTITION OF THE */
/*         COLUMNS OF A CONSISTENT WITH THE DETERMINATION OF A BY A */
/*         LOWER TRIANGULAR SUBSTITUTION METHOD. */

/*       INFO IS AN INT OUTPUT VARIABLE SET AS FOLLOWS. FOR */
/*         NORMAL TERMINATION INFO = 1. IF N OR NPAIRS IS NOT */
/*         POSITIVE OR LIWA IS LESS THAN 6*N, THEN INFO = 0. IF THE */
/*         K-TH ELEMENT OF INDROW OR THE K-TH ELEMENT OF INDCOL IS */
/*         NOT AN INT BETWEEN 1 AND N, OR IF THE K-TH DIAGONAL */
/*         ELEMENT IS NOT IN THE SPARSITY PATTERN, THEN INFO = -K. */

/*       IPNTR IS AN INT OUTPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS IN THE LOWER TRIANGULAR PART OF THE MATRIX A. */

/*       JPNTR IS AN INT OUTPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS IN THE LOWER TRIANGULAR PART OF THE MATRIX A. */

/*       IWA IS AN INT WORK ARRAY OF LENGTH LIWA. */

/*       LIWA IS A POSITIVE INT INPUT VARIABLE NOT LESS THAN 6*N. */

/*     SUBPROGRAMS CALLED */

/*       MINPACK-SUPPLIED ... DEGR,IDO,IDOG,NUMSRT,SDPT,SEQ,SETR, */
/*                            SLO,SLOG,SRTDAT */

/*       FORTRAN-SUPPLIED ... MAX,MIN */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. DECEMBER 1984. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     CHECK THE INPUT DATA. */

    /* Parameter adjustments */
    --jpntr;
    --ipntr;
    --ngrp;
    --listp;
    --indcol;
    --indrow;
    --iwa;

    /* Function Body */
    *info = 0;
    if (*n < 1 || *npairs < 1 || *liwa < *n * 6) {
	return 0;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	iwa[k] = 0;
/* L10: */
    }
    i__1 = *npairs;
    for (k = 1; k <= i__1; ++k) {
	*info = -k;
	if (indrow[k] < 1 || indrow[k] > *n || indcol[k] < 1 || indcol[k] > *
		n) {
	    return 0;
	}
	if (indrow[k] == indcol[k]) {
	    iwa[indrow[k]] = 1;
	}
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	*info = -k;
	if (iwa[k] != 1) {
	    return 0;
	}
/* L30: */
    }
    *info = 1;

/*     GENERATE THE SPARSITY PATTERN FOR THE LOWER */
/*     TRIANGULAR PART OF A. */

    i__1 = *npairs;
    for (k = 1; k <= i__1; ++k) {
	i__ = indrow[k];
	j = indcol[k];
	indrow[k] = GMAX(i__,j);
	indcol[k] = GMIN(i__,j);
/* L40: */
    }

/*     SORT THE DATA STRUCTURE BY COLUMNS. */

    srtdat_(n, npairs, &indrow[1], &indcol[1], &jpntr[1], &iwa[1]);

/*     COMPRESS THE DATA AND DETERMINE THE NUMBER OF NON-ZERO */
/*     ELEMENTS IN THE LOWER TRIANGULAR PART OF A. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwa[i__] = 0;
/* L50: */
    }
    nnz = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = nnz;
	i__2 = jpntr[j + 1] - 1;
	for (jp = jpntr[j]; jp <= i__2; ++jp) {
	    ir = indrow[jp];
	    if (iwa[ir] != j) {
		++nnz;
		indrow[nnz] = ir;
		iwa[ir] = j;
	    }
/* L60: */
	}
	jpntr[j] = k + 1;
/* L70: */
    }
    jpntr[*n + 1] = nnz + 1;

/*     EXTEND THE DATA STRUCTURE TO ROWS. */

    setr_(n, n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[1]);

/*     DETERMINE THE SMALLEST-LAST ORDERING OF THE VERTICES OF THE */
/*     ADJACENCY GRAPH OF A, AND FROM IT DETERMINE A LOWER BOUND */
/*     FOR THE NUMBER OF GROUPS. */

    slog_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[1], &maxclq, &
	    maxvd, &iwa[*n + 1], &iwa[(*n << 1) + 1], &iwa[*n * 3 + 1]);
    *mingrp = maxvd + 1;

/*     USE THE SELECTED METHOD. */

    if (*method == 1) {

/*        DIRECT METHOD. DETERMINE A PARTITION OF THE COLUMNS */
/*        OF A BY THE POWELL-TOINT METHOD. */

	sdpt_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &ngrp[1], 
		maxgrp, &iwa[*n + 1], &iwa[(*n << 1) + 1]);

/*        DEFINE A SYMMETRIC PERMUTATION OF A ACCORDING TO THE */
/*        ORDERING OF THE COLUMN GROUP NUMBERS IN THE PARTITION. */

	numsrt_(n, maxgrp, &ngrp[1], &c__1, &iwa[1], &iwa[(*n << 1) + 1], &
		iwa[*n + 1]);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    listp[iwa[i__]] = i__;
/* L80: */
	}
    } else {

/*        INDIRECT METHOD. DETERMINE THE INCIDENCE DEGREE ORDERING */
/*        OF THE VERTICES OF THE ADJACENCY GRAPH OF A AND, TOGETHER */
/*        WITH THE SMALLEST-LAST ORDERING, DEFINE A SYMMETRIC */
/*        PERMUTATION OF A. */

	idog_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &listp[1], &
		maxclq, &maxid, &iwa[*n + 1], &iwa[(*n << 1) + 1], &iwa[*n * 
		3 + 1]);
	if (maxid > maxvd) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		listp[i__] = iwa[i__];
/* L90: */
	    }
	}

/*        GENERATE THE SPARSITY PATTERN FOR THE LOWER */
/*        TRIANGULAR PART L OF THE PERMUTED MATRIX. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = jpntr[j + 1] - 1;
	    for (jp = jpntr[j]; jp <= i__2; ++jp) {
		i__ = indrow[jp];
/* Computing MAX */
		i__3 = listp[i__], i__4 = listp[j];
		indrow[jp] = GMAX(i__3,i__4);
/* Computing MIN */
		i__3 = listp[i__], i__4 = listp[j];
		indcol[jp] = GMIN(i__3,i__4);
/* L100: */
	    }
/* L110: */
	}

/*        SORT THE DATA STRUCTURE BY COLUMNS. */

	srtdat_(n, &nnz, &indrow[1], &indcol[1], &jpntr[1], &iwa[1]);

/*        EXTEND THE DATA STRUCTURE TO ROWS. */

	setr_(n, n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[1]);

/*        DETERMINE THE DEGREE SEQUENCE FOR THE INTERSECTION */
/*        GRAPH OF THE COLUMNS OF L. */

	degr_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[*n * 5 + 
		1], &iwa[*n + 1]);

/*        COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF L */
/*        WITH THE SMALLEST-LAST (SL) ORDERING. */

	slo_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[*n * 5 + 1]
		, &iwa[(*n << 2) + 1], &maxclq, &iwa[1], &iwa[*n + 1], &iwa[(*
		n << 1) + 1], &iwa[*n * 3 + 1]);
	seq_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[(*n << 2) 
		+ 1], &iwa[1], maxgrp, &iwa[*n + 1]);
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    ngrp[j] = iwa[listp[j]];
/* L120: */
	}

/*        EXIT IF THE SMALLEST-LAST ORDERING IS OPTIMAL. */

	if (*maxgrp == maxclq) {
	    goto L140;
	}

/*        COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF L */
/*        WITH THE INCIDENCE DEGREE (ID) ORDERING. */

	ido_(n, n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[*n * 5 
		+ 1], &iwa[(*n << 2) + 1], &maxclq, &iwa[1], &iwa[*n + 1], &
		iwa[(*n << 1) + 1], &iwa[*n * 3 + 1]);
	seq_(n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[(*n << 2) 
		+ 1], &iwa[1], &numgrp, &iwa[*n + 1]);

/*        RETAIN THE BETTER OF THE TWO ORDERINGS. */

	if (numgrp < *maxgrp) {
	    *maxgrp = numgrp;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		ngrp[j] = iwa[listp[j]];
/* L130: */
	    }
	}
L140:

/*        GENERATE THE SPARSITY PATTERN FOR THE LOWER */
/*        TRIANGULAR PART OF THE ORIGINAL MATRIX. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    iwa[listp[j]] = j;
/* L150: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = jpntr[j + 1] - 1;
	    for (jp = jpntr[j]; jp <= i__2; ++jp) {
		i__ = indrow[jp];
/* Computing MAX */
		i__3 = iwa[i__], i__4 = iwa[j];
		indrow[jp] = GMAX(i__3,i__4);
/* Computing MIN */
		i__3 = iwa[i__], i__4 = iwa[j];
		indcol[jp] = GMIN(i__3,i__4);
/* L160: */
	    }
/* L170: */
	}

/*        SORT THE DATA STRUCTURE BY COLUMNS. */

	srtdat_(n, &nnz, &indrow[1], &indcol[1], &jpntr[1], &iwa[1]);

/*        EXTEND THE DATA STRUCTURE TO ROWS. */

	setr_(n, n, &indrow[1], &jpntr[1], &indcol[1], &ipntr[1], &iwa[1]);
    }
    return 0;

/*     LAST CARD OF SUBROUTINE DSSM. */

} /* dssm_ */

/* Subroutine */ int idog_(int *n, int *nghbrp, int *npntrp, 
	int *nghbrs, int *npntrs, int *listp, int *maxclq, 
	int *maxid, int *iwa1, int *iwa2, int *iwa3)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int i__, j, k, ncomp, maxdeg, maxinc, numdeg, numinc, numord, 
	    maxlst;
    extern /* Subroutine */ int numsrt_(int *, int *, int *, 
	    int *, int *, int *, int *);

/*     ********** */

/*     SUBROUTINE IDOG */

/*     GIVEN A LOOPLESS GRAPH G = (V,E), THIS SUBROUTINE DETERMINES */
/*     THE INCIDENCE DEGREE ORDERING OF THE VERTICES OF G. */

/*     THE INCIDENCE DEGREE ORDERING IS DETERMINED RECURSIVELY BY */
/*     LETTING LIST(K), K = 1,...,N BE A VERTEX WITH MAXIMAL */
/*     INCIDENCE TO THE SUBGRAPH SPANNED BY THE ORDERED VERTICES. */
/*     AMONG ALL THE VERTICES OF MAXIMAL INCIDENCE, A VERTEX OF */
/*     MAXIMAL DEGREE IS CHOSEN. THIS SUBROUTINE DETERMINES THE */
/*     INVERSE OF THE INCIDENCE DEGREE ORDERING, THAT IS, AN ARRAY */
/*     LISTP SUCH THAT LISTP(LIST(K)) = K FOR K = 1,2,...,N. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE IDOG(N,NGHBRP,NPNTRP,NGHBRS,NPNTRS,LISTP, */
/*                       MAXCLQ,MAXID,IWA1,IWA2,IWA3) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF VERTICES OF G. */

/*       NGHBRP IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         PREDECESSOR ADJACENCY LISTS FOR THE GRAPH G. */

/*       NPNTRP IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE PREDECESSOR ADJACENCY */
/*         LISTS IN NGHBRP. THE VERTICES PRECEDING AND ADJACENT */
/*         TO VERTEX J ARE */

/*               NGHBRP(K), K = NPNTRP(J),...,NPNTRP(J+1)-1. */

/*         NOTE THAT NPNTRP(N+1)-1 IS THEN THE NUMBER OF VERTICES */
/*         PLUS EDGES OF THE GRAPH G. */

/*       NGHBRS IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         SUCCESSOR ADJACENCY LISTS FOR THE GRAPH G. */

/*       NPNTRS IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE SUCCESSOR ADJACENCY */
/*         LISTS IN NGHBRS. THE VERTICES SUCCEEDING AND ADJACENT */
/*         TO VERTEX J ARE */

/*               NGHBRS(K), K = NPNTRS(J),...,NPNTRS(J+1)-1. */

/*         NOTE THAT NPNTRS(N+1)-1 IS THEN THE NUMBER OF VERTICES */
/*         PLUS EDGES OF THE GRAPH G. */

/*       LISTP IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE INVERSE OF THE INCIDENCE DEGREE ORDERING OF THE */
/*         VERTICES. VERTEX J IS IN POSITION LISTP(J) OF THIS ORDERING. */

/*       MAXCLQ IS AN INT OUTPUT VARIABLE SET TO THE SIZE */
/*         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING. */

/*       MAXID IS AN INT OUTPUT VARIABLE SET TO THE MAXIMUM */
/*         INCIDENCE DEGREE FOUND DURING THE ORDERING. */

/*       IWA1,IWA2, AND IWA3 ARE INT WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       MINPACK-SUPPLIED ... NUMSRT */

/*       FORTRAN-SUPPLIED ... MAX */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. DECEMBER 1984. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa3;
    --iwa2;
    --listp;
    --npntrs;
    --npntrp;
    --nghbrp;
    --nghbrs;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	listp[j] = npntrp[j + 1] - npntrp[j] - 1 + (npntrs[j + 1] - npntrs[j] 
		- 1);
/* L10: */
    }
    maxlst = (npntrp[*n + 1] + npntrs[*n + 1]) / *n;

/*     SORT THE DEGREE SEQUENCE. */

    i__1 = *n - 1;
    numsrt_(n, &i__1, &listp[1], &c__1, iwa1, &iwa2[1], &iwa3[1]);

/*     CREATE A DOUBLY-LINKED LIST TO ACCESS THE INCIDENCES OF THE */
/*     VERTICES. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS. */

/*     EACH UN-ORDERED VERTEX I IS IN A LIST (THE INCIDENCE LIST) */
/*     OF VERTICES WITH THE SAME INCIDENCE. */

/*     IWA1(NUMINC) IS THE FIRST VERTEX IN THE NUMINC LIST */
/*     UNLESS IWA1(NUMINC) = 0. IN THIS CASE THERE ARE */
/*     NO VERTICES IN THE NUMINC LIST. */

/*     IWA2(I) IS THE VERTEX BEFORE I IN THE INCIDENCE LIST */
/*     UNLESS IWA2(I) = 0. IN THIS CASE I IS THE FIRST */
/*     VERTEX IN THIS INCIDENCE LIST. */

/*     IWA3(I) IS THE VERTEX AFTER I IN THE INCIDENCE LIST */
/*     UNLESS IWA3(I) = 0. IN THIS CASE I IS THE LAST */
/*     VERTEX IN THIS INCIDENCE LIST. */

/*     IF I IS AN UN-ORDERED VERTEX, THEN -LISTP(I) IS THE */
/*     INCIDENCE OF I TO THE GRAPH INDUCED BY THE ORDERED */
/*     VERTICES. IF J IS AN ORDERED VERTEX, THEN LISTP(J) */
/*     IS THE INCIDENCE DEGREE ORDER OF VERTEX J. */

    maxinc = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__ = iwa1[j - 1];
	iwa1[j - 1] = 0;
	iwa2[i__] = 0;
	iwa3[i__] = iwa1[0];
	if (iwa1[0] > 0) {
	    iwa2[iwa1[0]] = i__;
	}
	iwa1[0] = i__;
	listp[j] = 0;
/* L20: */
    }
    *maxclq = 0;
    *maxid = 0;
    numord = 1;

/*     BEGINNING OF ITERATION LOOP. */

L30:

/*        CHOOSE A VERTEX J OF MAXIMAL DEGREE AMONG THE */
/*        VERTICES OF MAXIMAL INCIDENCE MAXINC. */

L40:
    k = iwa1[maxinc];
    if (k > 0) {
	goto L50;
    }
    --maxinc;
    goto L40;
L50:
    maxdeg = -1;
    i__1 = maxlst;
    for (i__ = 1; i__ <= i__1; ++i__) {
	numdeg = npntrp[k + 1] - npntrp[k] - 1 + (npntrs[k + 1] - npntrs[k] - 
		1);
	if (numdeg > maxdeg) {
	    maxdeg = numdeg;
	    j = k;
	}
	k = iwa3[k];
	if (k <= 0) {
	    goto L70;
	}
/* L60: */
    }
L70:
    listp[j] = numord;
    *maxid = GMAX(*maxid,maxinc);

/*        UPDATE THE SIZE OF THE LARGEST CLIQUE */
/*        FOUND DURING THE ORDERING. */

    if (maxinc == 0) {
	ncomp = 0;
    }
    ++ncomp;
    if (maxinc + 1 == ncomp) {
	*maxclq = GMAX(*maxclq,ncomp);
    }

/*        TERMINATION TEST. */

    ++numord;
    if (numord > *n) {
	goto L100;
    }

/*        DELETE VERTEX J FROM THE MAXINC LIST. */

    if (iwa2[j] == 0) {
	iwa1[maxinc] = iwa3[j];
    } else {
	iwa3[iwa2[j]] = iwa3[j];
    }
    if (iwa3[j] > 0) {
	iwa2[iwa3[j]] = iwa2[j];
    }

/*        DETERMINE ALL THE NEIGHBORS OF VERTEX J WHICH PRECEDE J */
/*        IN THE SUBGRAPH SPANNED BY THE UN-ORDERED VERTICES. */

    i__1 = npntrp[j + 1] - 1;
    for (k = npntrp[j]; k <= i__1; ++k) {
	i__ = nghbrp[k];

/*           UPDATE THE POINTERS TO THE CURRENT INCIDENCE LISTS. */

	numinc = -listp[i__];
	if (numinc >= 0) {
	    --listp[i__];
/* Computing MAX */
	    i__2 = maxinc, i__3 = -listp[i__];
	    maxinc = GMAX(i__2,i__3);

/*              DELETE VERTEX I FROM THE NUMINC LIST. */

	    if (iwa2[i__] == 0) {
		iwa1[numinc] = iwa3[i__];
	    } else {
		iwa3[iwa2[i__]] = iwa3[i__];
	    }
	    if (iwa3[i__] > 0) {
		iwa2[iwa3[i__]] = iwa2[i__];
	    }

/*              ADD VERTEX I TO THE NUMINC+1 LIST. */

	    iwa2[i__] = 0;
	    iwa3[i__] = iwa1[numinc + 1];
	    if (iwa1[numinc + 1] > 0) {
		iwa2[iwa1[numinc + 1]] = i__;
	    }
	    iwa1[numinc + 1] = i__;
	}
/* L80: */
    }

/*        DETERMINE ALL THE NEIGHBORS OF VERTEX J WHICH SUCCEED J */
/*        IN THE SUBGRAPH SPANNED BY THE UN-ORDERED VERTICES. */

    i__1 = npntrs[j + 1] - 1;
    for (k = npntrs[j]; k <= i__1; ++k) {
	i__ = nghbrs[k];

/*           UPDATE THE POINTERS TO THE CURRENT INCIDENCE LISTS. */

	numinc = -listp[i__];
	if (numinc >= 0) {
	    --listp[i__];
/* Computing MAX */
	    i__2 = maxinc, i__3 = -listp[i__];
	    maxinc = GMAX(i__2,i__3);

/*              DELETE VERTEX I FROM THE NUMINC LIST. */

	    if (iwa2[i__] == 0) {
		iwa1[numinc] = iwa3[i__];
	    } else {
		iwa3[iwa2[i__]] = iwa3[i__];
	    }
	    if (iwa3[i__] > 0) {
		iwa2[iwa3[i__]] = iwa2[i__];
	    }

/*              ADD VERTEX I TO THE NUMINC+1 LIST. */

	    iwa2[i__] = 0;
	    iwa3[i__] = iwa1[numinc + 1];
	    if (iwa1[numinc + 1] > 0) {
		iwa2[iwa1[numinc + 1]] = i__;
	    }
	    iwa1[numinc + 1] = i__;
	}
/* L90: */
    }

/*        END OF ITERATION LOOP. */

    goto L30;
L100:
    return 0;

/*     LAST CARD OF SUBROUTINE IDOG. */

} /* idog_ */

/* Subroutine */ int sdpt_(int *n, int *nghbrp, int *npntrp, 
	int *nghbrs, int *npntrs, int *ngrp, int *maxgrp, 
	int *iwa1, int *iwa2)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int j, k, l, jp, kp, numv, maxdeg, numdeg;

/*     ********** */

/*     SUBROUTINE SDPT */

/*     GIVEN A LOOPLESS GRAPH G = (V,E), THIS SUBROUTINE DETERMINES */
/*     A SYMMETRIC COLORING OF G BY THE POWELL-TOINT DIRECT METHOD. */

/*     THE POWELL-TOINT METHOD ASSIGNS THE K-TH COLOR BY EXAMINING */
/*     THE UN-COLORED VERTICES U(K) IN ORDER OF NON-INCREASING DEGREE */
/*     AND ASSIGNING COLOR K TO VERTEX V IF THERE ARE NO PATHS OF */
/*     LENGTH 1 OR 2 (IN THE GRAPH INDUCED BY U(K)) BETWEEN V AND */
/*     SOME K-COLORED VERTEX. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SDPT(N,NGHBRP,NPNTRP,NGHBRS,NPNTRS,NGRP,MAXGRP, */
/*                       IWA1,IWA2) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF VERTICES OF G. */

/*       NGHBRP IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         PREDECESSOR ADJACENCY LISTS FOR THE GRAPH G. */

/*       NPNTRP IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE PREDECESSOR ADJACENCY */
/*         LISTS IN NGHBRP. THE VERTICES PRECEDING AND ADJACENT */
/*         TO VERTEX J ARE */

/*               NGHBRP(K), K = NPNTRP(J),...,NPNTRP(J+1)-1. */

/*         NOTE THAT NPNTRP(N+1)-1 IS THEN THE NUMBER OF VERTICES */
/*         PLUS EDGES OF THE GRAPH G. */

/*       NGHBRS IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         SUCCESSOR ADJACENCY LISTS FOR THE GRAPH G. */

/*       NPNTRS IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE SUCCESSOR ADJACENCY */
/*         LISTS IN NGHBRS. THE VERTICES SUCCEEDING AND ADJACENT */
/*         TO VERTEX J ARE */

/*               NGHBRS(K), K = NPNTRS(J),...,NPNTRS(J+1)-1. */

/*         NOTE THAT NPNTRS(N+1)-1 IS THEN THE NUMBER OF VERTICES */
/*         PLUS EDGES OF THE GRAPH G. */

/*       NGRP IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE SYMMETRIC COLORING OF G. VERTEX J IS COLORED WITH */
/*         COLOR NGRP(J). */

/*       MAXGRP IS AN INT OUTPUT VARIABLE WHICH SPECIFIES THE */
/*         NUMBER OF COLORS IN THE SYMMETRIC COLORING OF G. */

/*       IWA1 AND IWA2 ARE INT WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... MAX */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. DECEMBER 1984. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. NUMV IS THE CURRENT NUMBER OF UN-COLORED */
/*     VERTICES, MAXDEG IS THE MAXIMUM INDUCED DEGREE OF THESE */
/*     VERTICES, AND MAXGRP IS THE CURRENT GROUP NUMBER (COLOR). */

    /* Parameter adjustments */
    --iwa2;
    --ngrp;
    --npntrs;
    --npntrp;
    --nghbrp;
    --nghbrs;

    /* Function Body */
    numv = *n;
    maxdeg = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ngrp[j] = npntrp[j] - npntrp[j + 1] + 1 + (npntrs[j] - npntrs[j + 1] 
		+ 1);
/* Computing MAX */
	i__2 = maxdeg, i__3 = -ngrp[j];
	maxdeg = GMAX(i__2,i__3);
	iwa2[j] = -j;
/* L10: */
    }
    *maxgrp = 0;

/*     BEGINNING OF ITERATION LOOP. */

L20:

/*        SORT THE LIST OF UN-COLORED VERTICES SO THAT THEIR */
/*        INDUCED DEGREES ARE IN NON-DECREASING ORDER. */

    i__1 = maxdeg;
    for (numdeg = 0; numdeg <= i__1; ++numdeg) {
	iwa1[numdeg] = 0;
/* L30: */
    }
    i__1 = numv;
    for (l = 1; l <= i__1; ++l) {
	numdeg = -ngrp[-iwa2[l]];
	++iwa1[numdeg];
/* L40: */
    }
    k = 1;
    for (numdeg = maxdeg; numdeg >= 0; --numdeg) {
	l = iwa1[numdeg];
	iwa1[numdeg] = k;
	k += l;
/* L50: */
    }
    k = 1;
L60:
    j = iwa2[k];
    if (j > 0) {
	k = iwa1[-ngrp[j]];
    } else {
	numdeg = -ngrp[-j];
	l = iwa1[numdeg];
	iwa2[k] = iwa2[l];
	iwa2[l] = -j;
	++iwa1[numdeg];
    }
    if (k <= numv) {
	goto L60;
    }
    ++(*maxgrp);

/*        DETERMINE THE VERTICES IN GROUP MAXGRP. */

    i__1 = numv;
    for (l = 1; l <= i__1; ++l) {
	j = iwa2[l];

/*           EXAMINE EACH VERTEX K PRECEDING VERTEX J AND ALL */
/*           THE NEIGHBORS OF VERTEX K TO DETERMINE IF VERTEX */
/*           J CAN BE CONSIDERED FOR GROUP MAXGRP. */

	i__2 = npntrp[j + 1] - 1;
	for (jp = npntrp[j]; jp <= i__2; ++jp) {
	    k = nghbrp[jp];
	    if (ngrp[k] == *maxgrp) {
		goto L150;
	    }
	    if (ngrp[k] <= 0) {
		i__3 = npntrp[k + 1] - 1;
		for (kp = npntrp[k]; kp <= i__3; ++kp) {
		    if (ngrp[nghbrp[kp]] == *maxgrp) {
			goto L150;
		    }
/* L70: */
		}
		i__3 = npntrs[k + 1] - 1;
		for (kp = npntrs[k]; kp <= i__3; ++kp) {
		    if (ngrp[nghbrs[kp]] == *maxgrp) {
			goto L150;
		    }
/* L80: */
		}
	    }
/* L90: */
	}

/*           EXAMINE EACH VERTEX K SUCCEEDING VERTEX J AND ALL */
/*           THE NEIGHBORS OF VERTEX K TO DETERMINE IF VERTEX */
/*           J CAN BE ADDED TO GROUP MAXGRP. */

	i__2 = npntrs[j + 1] - 1;
	for (jp = npntrs[j]; jp <= i__2; ++jp) {
	    k = nghbrs[jp];
	    if (ngrp[k] == *maxgrp) {
		goto L150;
	    }
	    if (ngrp[k] <= 0) {
		i__3 = npntrp[k + 1] - 1;
		for (kp = npntrp[k]; kp <= i__3; ++kp) {
		    if (ngrp[nghbrp[kp]] == *maxgrp) {
			goto L150;
		    }
/* L100: */
		}
		i__3 = npntrs[k + 1] - 1;
		for (kp = npntrs[k]; kp <= i__3; ++kp) {
		    if (ngrp[nghbrs[kp]] == *maxgrp) {
			goto L150;
		    }
/* L110: */
		}
	    }
/* L120: */
	}

/*           ADD VERTEX J TO GROUP MAXGRP AND REMOVE VERTEX J */
/*           FROM THE LIST OF UN-COLORED VERTICES. */

	ngrp[j] = *maxgrp;
	iwa2[l] = 0;

/*           UPDATE THE DEGREES OF THE NEIGHBORS OF VERTEX J. */

	i__2 = npntrp[j + 1] - 1;
	for (jp = npntrp[j]; jp <= i__2; ++jp) {
	    k = nghbrp[jp];
	    if (ngrp[k] < 0) {
		++ngrp[k];
	    }
/* L130: */
	}
	i__2 = npntrs[j + 1] - 1;
	for (jp = npntrs[j]; jp <= i__2; ++jp) {
	    k = nghbrs[jp];
	    if (ngrp[k] < 0) {
		++ngrp[k];
	    }
/* L140: */
	}
L150:
/* L160: */
	;
    }

/*        COMPRESS THE UPDATED LIST OF UN-COLORED VERTICES. */
/*        RESET NUMV AND RECOMPUTE MAXDEG. */

    k = 0;
    maxdeg = 0;
    i__1 = numv;
    for (l = 1; l <= i__1; ++l) {
	if (iwa2[l] != 0) {
	    ++k;
	    iwa2[k] = -iwa2[l];
/* Computing MAX */
	    i__2 = maxdeg, i__3 = -ngrp[iwa2[l]];
	    maxdeg = GMAX(i__2,i__3);
	}
/* L170: */
    }
    numv = k;

/*        END OF ITERATION LOOP. */

    if (numv > 0) {
	goto L20;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SDPT. */

} /* sdpt_ */

/* Subroutine */ int slog_(int *n, int *nghbrp, int *npntrp, 
	int *nghbrs, int *npntrs, int *listp, int *maxclq, 
	int *maxvd, int *iwa1, int *iwa2, int *iwa3)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int i__, j, k, mindeg, numdeg, numord;

/*     ********** */

/*     SUBROUTINE SLOG */

/*     GIVEN A LOOPLESS GRAPH G = (V,E), THIS SUBROUTINE DETERMINES */
/*     THE SMALLEST-LAST ORDERING OF THE VERTICES OF G. */

/*     THE SMALLEST-LAST ORDERING IS DETERMINED RECURSIVELY BY */
/*     LETTING LIST(K), K = N,...,1 BE A VERTEX WITH LEAST DEGREE */
/*     IN THE SUBGRAPH SPANNED BY THE UN-ORDERED VERTICES. */
/*     THIS SUBROUTINE DETERMINES THE INVERSE OF THE SMALLEST-LAST */
/*     ORDERING, THAT IS, AN ARRAY LISTP SUCH THAT LISTP(LIST(K)) = K */
/*     FOR K = 1,2,...,N. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SLOG(N,NGHBRP,NPNTRP,NGHBRS,NPNTRS,LISTP, */
/*                       MAXCLQ,MAXVD,IWA1,IWA2,IWA3) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF VERTICES OF G. */

/*       NGHBRP IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         PREDECESSOR ADJACENCY LISTS FOR THE GRAPH G. */

/*       NPNTRP IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE PREDECESSOR ADJACENCY */
/*         LISTS IN NGHBRP. THE VERTICES PRECEDING AND ADJACENT */
/*         TO VERTEX J ARE */

/*               NGHBRP(K), K = NPNTRP(J),...,NPNTRP(J+1)-1. */

/*         NOTE THAT NPNTRP(N+1)-1 IS THEN THE NUMBER OF VERTICES */
/*         PLUS EDGES OF THE GRAPH G. */

/*       NGHBRS IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         SUCCESSOR ADJACENCY LISTS FOR THE GRAPH G. */

/*       NPNTRS IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE SUCCESSOR ADJACENCY */
/*         LISTS IN NGHBRS. THE VERTICES SUCCEEDING AND ADJACENT */
/*         TO VERTEX J ARE */

/*               NGHBRS(K), K = NPNTRS(J),...,NPNTRS(J+1)-1. */

/*         NOTE THAT NPNTRS(N+1)-1 IS THEN THE NUMBER OF VERTICES */
/*         PLUS EDGES OF THE GRAPH G. */

/*       LISTP IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE INVERSE OF THE SMALLEST-LAST ORDERING OF THE VERTICES. */
/*         VERTEX J IS IN POSITION LISTP(J) OF THIS ORDERING. */

/*       MAXCLQ IS AN INT OUTPUT VARIABLE SET TO THE SIZE */
/*         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING. */

/*       MAXVD IS AN INT OUTPUT VARIABLE SET TO THE MAXIMUM */
/*         VERTEX DEGREE FOUND DURING THE ORDERING. */

/*       IWA1,IWA2, AND IWA3 ARE INT WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... MAX,MIN */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. DECEMBER 1984. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa3;
    --iwa2;
    --listp;
    --npntrs;
    --npntrp;
    --nghbrp;
    --nghbrs;

    /* Function Body */
    mindeg = *n;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	iwa1[j - 1] = 0;
	listp[j] = npntrp[j] - npntrp[j + 1] + 1 + (npntrs[j] - npntrs[j + 1] 
		+ 1);
/* Computing MIN */
	i__2 = mindeg, i__3 = -listp[j];
	mindeg = GMIN(i__2,i__3);
/* L10: */
    }

/*     CREATE A DOUBLY-LINKED LIST TO ACCESS THE DEGREES OF THE */
/*     VERTICES. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS. */

/*     EACH UN-ORDERED VERTEX I IS IN A LIST (THE DEGREE LIST) */
/*     OF VERTICES WITH THE SAME DEGREE. */

/*     IWA1(NUMDEG) IS THE FIRST VERTEX IN THE NUMDEG LIST */
/*     UNLESS IWA1(NUMDEG) = 0. IN THIS CASE THERE ARE */
/*     NO VERTICES IN THE NUMDEG LIST. */

/*     IWA2(I) IS THE VERTEX BEFORE I IN THE DEGREE LIST */
/*     UNLESS IWA2(I) = 0. IN THIS CASE I IS THE FIRST */
/*     VERTEX IN THIS DEGREE LIST. */

/*     IWA3(I) IS THE VERTEX AFTER I IN THE DEGREE LIST */
/*     UNLESS IWA3(I) = 0. IN THIS CASE I IS THE LAST */
/*     VERTEX IN THIS DEGREE LIST. */

/*     IF I IS AN UN-ORDERED VERTEX, THEN -LISTP(I) IS THE */
/*     DEGREE OF I IN THE GRAPH INDUCED BY THE UN-ORDERED */
/*     VERTICES. IF J IS AN ORDERED VERTEX, THEN LISTP(J) */
/*     IS THE SMALLEST-LAST ORDER OF VERTEX J. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	numdeg = -listp[j];
	iwa2[j] = 0;
	iwa3[j] = iwa1[numdeg];
	if (iwa1[numdeg] > 0) {
	    iwa2[iwa1[numdeg]] = j;
	}
	iwa1[numdeg] = j;
/* L20: */
    }
    *maxclq = 0;
    *maxvd = 0;
    numord = *n;

/*     BEGINNING OF ITERATION LOOP. */

L30:

/*        CHOOSE A VERTEX J OF MINIMAL DEGREE MINDEG. */

L40:
    j = iwa1[mindeg];
    if (j > 0) {
	goto L50;
    }
    ++mindeg;
    goto L40;
L50:
    listp[j] = numord;
    *maxvd = GMAX(*maxvd,mindeg);

/*        MARK THE SIZE OF THE LARGEST CLIQUE */
/*        FOUND DURING THE ORDERING. */

    if (mindeg + 1 == numord && *maxclq == 0) {
	*maxclq = numord;
    }

/*        TERMINATION TEST. */

    --numord;
    if (numord == 0) {
	goto L80;
    }

/*        DELETE VERTEX J FROM THE MINDEG LIST. */

    iwa1[mindeg] = iwa3[j];
    if (iwa3[j] > 0) {
	iwa2[iwa3[j]] = 0;
    }

/*        DETERMINE ALL THE NEIGHBORS OF VERTEX J WHICH PRECEDE J */
/*        IN THE SUBGRAPH SPANNED BY THE UN-ORDERED VERTICES. */

    i__1 = npntrp[j + 1] - 1;
    for (k = npntrp[j]; k <= i__1; ++k) {
	i__ = nghbrp[k];

/*           UPDATE THE POINTERS TO THE CURRENT DEGREE LISTS. */

	numdeg = -listp[i__];
	if (numdeg >= 0) {
	    ++listp[i__];
/* Computing MIN */
	    i__2 = mindeg, i__3 = -listp[i__];
	    mindeg = GMIN(i__2,i__3);

/*              DELETE VERTEX I FROM THE NUMDEG LIST. */

	    if (iwa2[i__] == 0) {
		iwa1[numdeg] = iwa3[i__];
	    } else {
		iwa3[iwa2[i__]] = iwa3[i__];
	    }
	    if (iwa3[i__] > 0) {
		iwa2[iwa3[i__]] = iwa2[i__];
	    }

/*              ADD VERTEX I TO THE NUMDEG-1 LIST. */

	    iwa2[i__] = 0;
	    iwa3[i__] = iwa1[numdeg - 1];
	    if (iwa1[numdeg - 1] > 0) {
		iwa2[iwa1[numdeg - 1]] = i__;
	    }
	    iwa1[numdeg - 1] = i__;
	}
/* L60: */
    }

/*        DETERMINE ALL THE NEIGHBORS OF VERTEX J WHICH SUCCEED J */
/*        IN THE SUBGRAPH SPANNED BY THE UN-ORDERED VERTICES. */

    i__1 = npntrs[j + 1] - 1;
    for (k = npntrs[j]; k <= i__1; ++k) {
	i__ = nghbrs[k];

/*           UPDATE THE POINTERS TO THE CURRENT DEGREE LISTS. */

	numdeg = -listp[i__];
	if (numdeg >= 0) {
	    ++listp[i__];
/* Computing MIN */
	    i__2 = mindeg, i__3 = -listp[i__];
	    mindeg = GMIN(i__2,i__3);

/*              DELETE VERTEX I FROM THE NUMDEG LIST. */

	    if (iwa2[i__] == 0) {
		iwa1[numdeg] = iwa3[i__];
	    } else {
		iwa3[iwa2[i__]] = iwa3[i__];
	    }
	    if (iwa3[i__] > 0) {
		iwa2[iwa3[i__]] = iwa2[i__];
	    }

/*              ADD VERTEX I TO THE NUMDEG-1 LIST. */

	    iwa2[i__] = 0;
	    iwa3[i__] = iwa1[numdeg - 1];
	    if (iwa1[numdeg - 1] > 0) {
		iwa2[iwa1[numdeg - 1]] = i__;
	    }
	    iwa1[numdeg - 1] = i__;
	}
/* L70: */
    }

/*        END OF ITERATION LOOP. */

    goto L30;
L80:
    return 0;

/*     LAST CARD OF SUBROUTINE SLOG. */

} /* slog_ */

/* Subroutine */
#ifdef __cplusplus
extern "C"
#endif
int fdhs_(int *n, int *indrow, int *jpntr, 
	int *indcol, int *ipntr, int *listp, int *ngrp, 
	int *maxgrp, int *numgrp, double *eta, double *fhesd, 
	double *fhes, int *iwa)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int i__, j, k, l, ip, jp;
    static double sum;
    static int numg, numl, irow;

/*     ********** */

/*     SUBROUTINE FDHS */

/*     THIS SUBROUTINE COMPUTES AN APPROXIMATION TO THE (SYMMETRIC) */
/*     HESSIAN MATRIX OF A FUNCTION BY A SUBSTITUTION METHOD. */
/*     THE LOWER TRIANGULAR PART OF THE APPROXIMATION IS STORED */
/*     WITH A COLUMN-ORIENTED DEFINITION OF THE SPARSITY PATTERN. */

/*     THIS SUBROUTINE REQUIRES A SYMMETRIC PERMUTATION OF THE */
/*     HESSIAN MATRIX AND A PARTITION OF THE COLUMNS OF THE HESSIAN */
/*     MATRIX CONSISTENT WITH THE DETERMINATION OF THE HESSIAN */
/*     MATRIX BY A LOWER TRIANGULAR SUBSTITUTION METHOD. */
/*     THIS INFORMATION CAN BE PROVIDED BY SUBROUTINE DSSM. */

/*     THE SYMMETRIC PERMUTATION OF THE HESSIAN MATRIX IS DEFINED */
/*     BY THE ARRAY LISTP. THIS ARRAY IS ONLY USED INTERNALLY. */

/*     THE PARTITION OF THE HESSIAN MATRIX IS DEFINED BY THE ARRAY */
/*     NGRP BY SETTING NGRP(J) TO THE GROUP NUMBER OF COLUMN J. */
/*     THE USER MUST PROVIDE AN APPROXIMATION TO THE COLUMNS OF */
/*     THE HESSIAN MATRIX IN EACH GROUP BY SPECIFYING A DIFFERENCE */
/*     PARAMETER VECTOR ETA AND AN APPROXIMATION TO H*D WHERE H IS */
/*     THE HESSIAN MATRIX AND THE VECTOR D IS DEFINED BY THE */
/*     FOLLOWING SECTION OF CODE. */

/*           DO 10 J = 1, N */
/*              D(J) = 0.0 */
/*              IF (NGRP(J) .EQ. NUMGRP) D(J) = ETA(J) */
/*        10    CONTINUE */

/*     IN THE ABOVE CODE NUMGRP IS A GROUP NUMBER AND ETA(J) IS THE */
/*     DIFFERENCE PARAMETER USED TO APPROXIMATE COLUMN J OF THE */
/*     HESSIAN MATRIX. SUITABLE VALUES FOR ETA(J) MUST BE PROVIDED. */

/*     AS MENTIONED ABOVE, AN APPROXIMATION TO H*D MUST BE PROVIDED. */
/*     FOR EXAMPLE, IF GRAD(X) IS THE GRADIENT OF THE FUNCTION AT X, */
/*     THEN */

/*           GRAD(X+D) - GRAD(X) */

/*     CORRESPONDS TO THE FORWARD DIFFERENCE APPROXIMATION. */

/*     THE LOWER TRIANGULAR SUBSTITUTION METHOD REQUIRES THAT THE */
/*     APPROXIMATIONS TO H*D FOR ALL THE GROUPS BE STORED IN SPECIAL */
/*     LOCATIONS OF THE ARRAY FHES. THIS IS DONE BY CALLING FDHS */
/*     SUCCESSIVELY WITH NUMGRP = 1,2,...,MAXGRP. ON THE CALL WITH */
/*     NUMGRP = MAXGRP, THE SUBROUTINE THEN PROCEEDS TO OVERWRITE */
/*     FHES WITH THE APPROXIMATION TO THE LOWER TRIANGULAR PART OF */
/*     THE HESSIAN MATRIX. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE FDHS(N,INDROW,JPNTR,INDCOL,IPNTR,LISTP,NGRP, */
/*                       MAXGRP,NUMGRP,ETA,FHESD,FHES,IWA) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE ORDER */
/*         OF THE HESSIAN MATRIX. */

/*       INDROW IS AN INT INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE LOWER TRIANGULAR PART */
/*         OF THE HESSIAN MATRIX. */

/*       JPNTR IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZEROES */
/*         IN THE LOWER TRIANGULAR PART OF THE HESSIAN MATRIX. */

/*       INDCOL IS AN INT INPUT ARRAY WHICH CONTAINS THE COLUMN */
/*         INDICES FOR THE NON-ZEROES IN THE LOWER TRIANGULAR PART */
/*         OF THE HESSIAN MATRIX. */

/*       IPNTR IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZEROES */
/*         IN THE LOWER TRIANGULAR PART OF THE HESSIAN MATRIX. */

/*       LISTP IS AN INT INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE SYMMETRIC PERMUTATION OF THE HESSIAN MATRIX. ELEMENT */
/*         (I,J) OF THE HESSIAN MATRIX IS THE (LISTP(I),LISTP(J)) */
/*         ELEMENT OF THE PERMUTED HESSIAN. */

/*       NGRP IS AN INT INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE PARTITION OF THE COLUMNS OF THE HESSIAN MATRIX. */
/*         COLUMN J BELONGS TO GROUP NGRP(J). */

/*       MAXGRP IS A POSITIVE INT INPUT VARIABLE WHICH SPECIFIES */
/*         THE NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF */
/*         THE HESSIAN MATRIX. */

/*       NUMGRP IS A POSITIVE INT INPUT VARIABLE SET TO A GROUP */
/*         NUMBER IN THE PARTITION. */

/*       ETA IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS THE */
/*         DIFFERENCE PARAMETER VECTOR. */

/*       FHESD IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS AN */
/*         APPROXIMATION TO H*D, WHERE H IS THE HESSIAN MATRIX */
/*         AND D IS THE DIFFERENCE VECTOR FOR GROUP NUMGRP. */

/*       FHES IS AN OUTPUT ARRAY OF LENGTH NNZ, WHERE NNZ IS THE */
/*         NUMBER OF NON-ZERO ELEMENTS IN THE LOWER TRIANGULAR PART */
/*         OF THE HESSIAN MATRIX. ON OUTPUT WITH NUMGRP LESS THAN */
/*         MAXGRP, THE FHESD ARRAY FOR GROUP NUMGRP HAS BEEN STORED */
/*         IN FHES. WHEN NUMGRP = MAXGRP THE SUBROUTINE OVERWRITES */
/*         FHES WITH AN APPROXIMATION TO THE LOWER TRIANGULAR PART */
/*         OF THE HESSIAN MATRIX. THE APPROXIMATION IS STORED IN */
/*         FHES WITH A COLUMN-ORIENTED DEFINITION OF THE SPARSITY */
/*         PATTERN. THUS THE ELEMENTS IN COLUMN J OF THE LOWER */
/*         TRIANGULAR PART OF THE HESSIAN MATRIX ARE */

/*               FHES(K), K = JPNTR(J),...,JPNTR(J+1)-1, */

/*         AND THE ROW INDICES FOR THESE ELEMENTS ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*       IWA IS AN INT WORK ARRAY OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... ABS */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. DECEMBER 1984. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     STORE THE I-TH ELEMENT OF GRADIENT DIFFERENCE FHESD */
/*     CORRESPONDING TO GROUP NUMGRP IF THERE IS A POSITION */
/*     (I,J) SUCH THAT NGRP(J) = NUMGRP AND (I,J) IS MAPPED */
/*     ONTO THE LOWER TRIANGULAR PART OF THE PERMUTED MATRIX. */

    /* Parameter adjustments */
    --iwa;
    --fhesd;
    --eta;
    --ngrp;
    --listp;
    --ipntr;
    --jpntr;
    --indrow;
    --indcol;
    --fhes;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (ngrp[j] == *numgrp) {
	    numl = listp[j];
	    i__2 = ipntr[j + 1] - 1;
	    for (ip = ipntr[j]; ip <= i__2; ++ip) {
		i__ = indcol[ip];
		if (listp[i__] > numl) {
		    i__3 = jpntr[i__ + 1] - 1;
		    for (jp = jpntr[i__]; jp <= i__3; ++jp) {
			if (indrow[jp] == j) {
			    fhes[jp] = fhesd[i__];
			    goto L20;
			}
/* L10: */
		    }
L20:
		    ;
		}
/* L30: */
	    }
	    i__2 = jpntr[j + 1] - 1;
	    for (jp = jpntr[j]; jp <= i__2; ++jp) {
		i__ = indrow[jp];
		if (listp[i__] >= numl) {
		    fhes[jp] = fhesd[i__];
		}
/* L40: */
	    }
	}
/* L50: */
    }

/*     EXIT IF THIS IS NOT THE LAST GROUP. */

    if (*numgrp < *maxgrp) {
	return 0;
    }

/*     MARK ALL COLUMN INDICES J SUCH THAT (I,J) IS MAPPED ONTO */
/*     THE LOWER TRIANGULAR PART OF THE PERMUTED MATRIX. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	numl = listp[i__];
	i__2 = ipntr[i__ + 1] - 1;
	for (ip = ipntr[i__]; ip <= i__2; ++ip) {
	    j = indcol[ip];
	    if (numl >= listp[j]) {
		indcol[ip] = -indcol[ip];
	    }
/* L60: */
	}
	i__2 = jpntr[i__ + 1] - 1;
	for (jp = jpntr[i__]; jp <= i__2; ++jp) {
	    j = indrow[jp];
	    if (numl > listp[j]) {
		indrow[jp] = -indrow[jp];
	    }
/* L70: */
	}
/* L80: */
    }

/*     INVERT THE ARRAY LISTP. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	iwa[listp[j]] = j;
/* L90: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	listp[j] = iwa[j];
/* L100: */
    }

/*     DETERMINE THE LOWER TRIANGULAR PART OF THE ORIGINAL MATRIX. */

    for (irow = *n; irow >= 1; --irow) {
	i__ = listp[irow];

/*        FIND THE POSITIONS OF THE ELEMENTS IN THE I-TH ROW OF THE */
/*        LOWER TRIANGULAR PART OF THE ORIGINAL MATRIX THAT HAVE */
/*        ALREADY BEEN DETERMINED. */

	i__1 = ipntr[i__ + 1] - 1;
	for (ip = ipntr[i__]; ip <= i__1; ++ip) {
	    j = indcol[ip];
	    if (j > 0) {
		i__2 = jpntr[j + 1] - 1;
		for (jp = jpntr[j]; jp <= i__2; ++jp) {
		    if (indrow[jp] == i__) {
			iwa[j] = jp;
			goto L120;
		    }
/* L110: */
		}
L120:
		;
	    }
/* L130: */
	}

/*        DETERMINE THE ELEMENTS IN THE I-TH ROW OF THE LOWER */
/*        TRIANGULAR PART OF THE ORIGINAL MATRIX WHICH GET MAPPED */
/*        ONTO THE LOWER TRIANGULAR PART OF THE PERMUTED MATRIX. */

	i__1 = ipntr[i__ + 1] - 1;
	for (k = ipntr[i__]; k <= i__1; ++k) {
	    j = -indcol[k];
	    if (j > 0) {
		indcol[k] = j;

/*              DETERMINE THE (I,J) ELEMENT. */

		numg = ngrp[j];
		sum = 0.f;
		i__2 = ipntr[i__ + 1] - 1;
		for (ip = ipntr[i__]; ip <= i__2; ++ip) {
		    l = (i__3 = indcol[ip], abs(i__3));
		    if (ngrp[l] == numg && l != j) {
			sum += fhes[iwa[l]] * eta[l];
		    }
/* L140: */
		}
		i__2 = jpntr[i__ + 1] - 1;
		for (jp = jpntr[i__]; jp <= i__2; ++jp) {
		    l = (i__3 = indrow[jp], abs(i__3));
		    if (ngrp[l] == numg && l != j) {
			sum += fhes[jp] * eta[l];
		    }
/* L150: */
		}

/*              STORE THE (I,J) ELEMENT. */

		i__2 = jpntr[j + 1] - 1;
		for (jp = jpntr[j]; jp <= i__2; ++jp) {
		    if (indrow[jp] == i__) {
			fhes[jp] = (fhes[jp] - sum) / eta[j];
			goto L170;
		    }
/* L160: */
		}
L170:
		;
	    }
/* L180: */
	}

/*        DETERMINE THE ELEMENTS IN THE I-TH ROW OF THE STRICT UPPER */
/*        TRIANGULAR PART OF THE ORIGINAL MATRIX WHICH GET MAPPED */
/*        ONTO THE LOWER TRIANGULAR PART OF THE PERMUTED MATRIX. */

	i__1 = jpntr[i__ + 1] - 1;
	for (k = jpntr[i__]; k <= i__1; ++k) {
	    j = -indrow[k];
	    if (j > 0) {
		indrow[k] = j;

/*              DETERMINE THE (I,J) ELEMENT. */

		numg = ngrp[j];
		sum = 0.f;
		i__2 = ipntr[i__ + 1] - 1;
		for (ip = ipntr[i__]; ip <= i__2; ++ip) {
		    l = (i__3 = indcol[ip], abs(i__3));
		    if (ngrp[l] == numg) {
			sum += fhes[iwa[l]] * eta[l];
		    }
/* L190: */
		}
		i__2 = jpntr[i__ + 1] - 1;
		for (jp = jpntr[i__]; jp <= i__2; ++jp) {
		    l = (i__3 = indrow[jp], abs(i__3));
		    if (ngrp[l] == numg && l != j) {
			sum += fhes[jp] * eta[l];
		    }
/* L200: */
		}

/*              STORE THE (I,J) ELEMENT. */

		fhes[k] = (fhes[k] - sum) / eta[j];
	    }
/* L210: */
	}
/* L220: */
    }

/*     RE-INVERT THE ARRAY LISTP. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	iwa[listp[j]] = j;
/* L230: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	listp[j] = iwa[j];
/* L240: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE FDHS. */

} /* fdhs_ */

/* Subroutine */ int degr_(int *n, int *indrow, int *jpntr, 
	int *indcol, int *ipntr, int *ndeg, int *iwa)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int ic, ip, jp, ir, jcol;

/*     ********** */

/*     SUBROUTINE DEGR */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, */
/*     THIS SUBROUTINE DETERMINES THE DEGREE SEQUENCE FOR */
/*     THE INTERSECTION GRAPH OF THE COLUMNS OF A. */

/*     IN GRAPH-THEORY TERMINOLOGY, THE INTERSECTION GRAPH OF */
/*     THE COLUMNS OF A IS THE LOOPLESS GRAPH G WITH VERTICES */
/*     A(J), J = 1,2,...,N WHERE A(J) IS THE J-TH COLUMN OF A */
/*     AND WITH EDGE (A(I),A(J)) IF AND ONLY IF COLUMNS I AND J */
/*     HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY DEGR AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE DEGR(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,IWA) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INT INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INT INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       NDEG IS AN INT OUTPUT ARRAY OF LENGTH N WHICH */
/*         SPECIFIES THE DEGREE SEQUENCE. THE DEGREE OF THE */
/*         J-TH COLUMN OF A IS NDEG(J). */

/*       IWA IS AN INT WORK ARRAY OF LENGTH N. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa;
    --ndeg;
    --jpntr;
    --indrow;
    --indcol;
    --ipntr;

    /* Function Body */
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	ndeg[jp] = 0;
	iwa[jp] = 0;
/* L10: */
    }

/*     COMPUTE THE DEGREE SEQUENCE BY DETERMINING THE CONTRIBUTIONS */
/*     TO THE DEGREES FROM THE CURRENT(JCOL) COLUMN AND FURTHER */
/*     COLUMNS WHICH HAVE NOT YET BEEN CONSIDERED. */

    i__1 = *n;
    for (jcol = 2; jcol <= i__1; ++jcol) {
	iwa[jcol] = *n;

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

	i__2 = jpntr[jcol + 1] - 1;
	for (jp = jpntr[jcol]; jp <= i__2; ++jp) {
	    ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	    i__3 = ipntr[ir + 1] - 1;
	    for (ip = ipntr[ir]; ip <= i__3; ++ip) {
		ic = indcol[ip];

/*              ARRAY IWA MARKS COLUMNS WHICH HAVE CONTRIBUTED TO */
/*              THE DEGREE COUNT OF COLUMN JCOL. UPDATE THE DEGREE */
/*              COUNTS OF THESE COLUMNS AS WELL AS COLUMN JCOL. */

		if (iwa[ic] < jcol) {
		    iwa[ic] = jcol;
		    ++ndeg[ic];
		    ++ndeg[jcol];
		}
/* L20: */
	    }
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE DEGR. */

} /* degr_ */

/* Subroutine */ int ido_(int *m, int *n, int *indrow, int *
	jpntr, int *indcol, int *ipntr, int *ndeg, int *list, 
	int *maxclq, int *iwa1, int *iwa2, int *iwa3, int 
	*iwa4)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Local variables */
    static int ic, ip, jp, ir, jcol, ncomp, maxinc, numinc, numord, 
	    maxlst, numwgt, numlst;
    extern /* Subroutine */ int numsrt_(int *, int *, int *, 
	    int *, int *, int *, int *);

/*     ********** */

/*     SUBROUTINE IDO */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS */
/*     SUBROUTINE DETERMINES AN INCIDENCE-DEGREE ORDERING OF THE */
/*     COLUMNS OF A. */

/*     THE INCIDENCE-DEGREE ORDERING IS DEFINED FOR THE LOOPLESS */
/*     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE */
/*     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF */
/*     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     THE INCIDENCE-DEGREE ORDERING IS DETERMINED RECURSIVELY BY */
/*     LETTING LIST(K), K = 1,...,N BE A COLUMN WITH MAXIMAL */
/*     INCIDENCE TO THE SUBGRAPH SPANNED BY THE ORDERED COLUMNS. */
/*     AMONG ALL THE COLUMNS OF MAXIMAL INCIDENCE, IDO CHOOSES A */
/*     COLUMN OF MAXIMAL DEGREE. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE IDO(M,N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST, */
/*                      MAXCLQ,IWA1,IWA2,IWA3,IWA4) */

/*     WHERE */

/*       M IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF ROWS OF A. */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INT INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INT INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       NDEG IS AN INT INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN */
/*         OF A IS NDEG(J). */

/*       LIST IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE INCIDENCE-DEGREE ORDERING OF THE COLUMNS OF A. THE J-TH */
/*         COLUMN IN THIS ORDER IS LIST(J). */

/*       MAXCLQ IS AN INT OUTPUT VARIABLE SET TO THE SIZE */
/*         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING. */

/*       IWA1,IWA2,IWA3, AND IWA4 ARE INT WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       MINPACK-SUPPLIED ... NUMSRT */

/*       FORTRAN-SUPPLIED ... MAX */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. AUGUST 1984. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     SORT THE DEGREE SEQUENCE. */

    /* Parameter adjustments */
    --ipntr;
    --iwa4;
    --iwa3;
    --iwa2;
    --list;
    --ndeg;
    --jpntr;
    --indrow;
    --indcol;

    /* Function Body */
    i__1 = *n - 1;
    numsrt_(n, &i__1, &ndeg[1], &c_n1, &iwa4[1], &iwa2[1], &iwa3[1]);

/*     INITIALIZATION BLOCK. */

/*     CREATE A DOUBLY-LINKED LIST TO ACCESS THE INCIDENCES OF THE */
/*     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS. */

/*     EACH UN-ORDERED COLUMN IC IS IN A LIST (THE INCIDENCE LIST) */
/*     OF COLUMNS WITH THE SAME INCIDENCE. */

/*     IWA1(NUMINC) IS THE FIRST COLUMN IN THE NUMINC LIST */
/*     UNLESS IWA1(NUMINC) = 0. IN THIS CASE THERE ARE */
/*     NO COLUMNS IN THE NUMINC LIST. */

/*     IWA2(IC) IS THE COLUMN BEFORE IC IN THE INCIDENCE LIST */
/*     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST */
/*     COLUMN IN THIS INCIDENCE LIST. */

/*     IWA3(IC) IS THE COLUMN AFTER IC IN THE INCIDENCE LIST */
/*     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST */
/*     COLUMN IN THIS INCIDENCE LIST. */

/*     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE */
/*     INCIDENCE OF IC TO THE GRAPH INDUCED BY THE ORDERED */
/*     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL) */
/*     IS THE INCIDENCE-DEGREE ORDER OF COLUMN JCOL. */

    maxinc = 0;
    for (jp = *n; jp >= 1; --jp) {
	ic = iwa4[jp];
	iwa1[*n - jp] = 0;
	iwa2[ic] = 0;
	iwa3[ic] = iwa1[0];
	if (iwa1[0] > 0) {
	    iwa2[iwa1[0]] = ic;
	}
	iwa1[0] = ic;
	iwa4[jp] = 0;
	list[jp] = 0;
/* L10: */
    }

/*     DETERMINE THE MAXIMAL SEARCH LENGTH FOR THE LIST */
/*     OF COLUMNS OF MAXIMAL INCIDENCE. */

    maxlst = 0;
    i__1 = *m;
    for (ir = 1; ir <= i__1; ++ir) {
/* Computing 2nd power */
	i__2 = ipntr[ir + 1] - ipntr[ir];
	maxlst += i__2 * i__2;
/* L20: */
    }
    maxlst /= *n;
    *maxclq = 0;
    numord = 1;

/*     BEGINNING OF ITERATION LOOP. */

L30:

/*        CHOOSE A COLUMN JCOL OF MAXIMAL DEGREE AMONG THE */
/*        COLUMNS OF MAXIMAL INCIDENCE MAXINC. */

L40:
    jp = iwa1[maxinc];
    if (jp > 0) {
	goto L50;
    }
    --maxinc;
    goto L40;
L50:
    numwgt = -1;
    i__1 = maxlst;
    for (numlst = 1; numlst <= i__1; ++numlst) {
	if (ndeg[jp] > numwgt) {
	    numwgt = ndeg[jp];
	    jcol = jp;
	}
	jp = iwa3[jp];
	if (jp <= 0) {
	    goto L70;
	}
/* L60: */
    }
L70:
    list[jcol] = numord;

/*        UPDATE THE SIZE OF THE LARGEST CLIQUE */
/*        FOUND DURING THE ORDERING. */

    if (maxinc == 0) {
	ncomp = 0;
    }
    ++ncomp;
    if (maxinc + 1 == ncomp) {
	*maxclq = GMAX(*maxclq,ncomp);
    }

/*        TERMINATION TEST. */

    ++numord;
    if (numord > *n) {
	goto L100;
    }

/*        DELETE COLUMN JCOL FROM THE MAXINC LIST. */

    if (iwa2[jcol] == 0) {
	iwa1[maxinc] = iwa3[jcol];
    } else {
	iwa3[iwa2[jcol]] = iwa3[jcol];
    }
    if (iwa3[jcol] > 0) {
	iwa2[iwa3[jcol]] = iwa2[jcol];
    }

/*        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL. */

    iwa4[jcol] = *n;

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

    i__1 = jpntr[jcol + 1] - 1;
    for (jp = jpntr[jcol]; jp <= i__1; ++jp) {
	ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	i__2 = ipntr[ir + 1] - 1;
	for (ip = ipntr[ir]; ip <= i__2; ++ip) {
	    ic = indcol[ip];

/*              ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO */
/*              COLUMN JCOL. */

	    if (iwa4[ic] < numord) {
		iwa4[ic] = numord;

/*                 UPDATE THE POINTERS TO THE CURRENT INCIDENCE LISTS. */

		numinc = list[ic];
		++list[ic];
/* Computing MAX */
		i__3 = maxinc, i__4 = list[ic];
		maxinc = GMAX(i__3,i__4);

/*                 DELETE COLUMN IC FROM THE NUMINC LIST. */

		if (iwa2[ic] == 0) {
		    iwa1[numinc] = iwa3[ic];
		} else {
		    iwa3[iwa2[ic]] = iwa3[ic];
		}
		if (iwa3[ic] > 0) {
		    iwa2[iwa3[ic]] = iwa2[ic];
		}

/*                 ADD COLUMN IC TO THE NUMINC+1 LIST. */

		iwa2[ic] = 0;
		iwa3[ic] = iwa1[numinc + 1];
		if (iwa1[numinc + 1] > 0) {
		    iwa2[iwa1[numinc + 1]] = ic;
		}
		iwa1[numinc + 1] = ic;
	    }
/* L80: */
	}
/* L90: */
    }

/*        END OF ITERATION LOOP. */

    goto L30;
L100:

/*     INVERT THE ARRAY LIST. */

    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	iwa2[list[jcol]] = jcol;
/* L110: */
    }
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	list[jp] = iwa2[jp];
/* L120: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE IDO. */

} /* ido_ */

/* Subroutine */ int numsrt_(int *n, int *nmax, int *num, int 
	*mode, int *index, int *last, int *next)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, j, k, l, jl, ju, jinc;

/*     **********. */

/*     SUBROUTINE NUMSRT */

/*     GIVEN A SEQUENCE OF INTS, THIS SUBROUTINE GROUPS */
/*     TOGETHER THOSE INDICES WITH THE SAME SEQUENCE VALUE */
/*     AND, OPTIONALLY, SORTS THE SEQUENCE INTO EITHER */
/*     ASCENDING OR DESCENDING ORDER. */

/*     THE SEQUENCE OF INTS IS DEFINED BY THE ARRAY NUM, */
/*     AND IT IS ASSUMED THAT THE INTS ARE EACH FROM THE SET */
/*     0,1,...,NMAX. ON OUTPUT THE INDICES K SUCH THAT NUM(K) = L */
/*     FOR ANY L = 0,1,...,NMAX CAN BE OBTAINED FROM THE ARRAYS */
/*     LAST AND NEXT AS FOLLOWS. */

/*           K = LAST(L) */
/*           WHILE (K .NE. 0) K = NEXT(K) */

/*     OPTIONALLY, THE SUBROUTINE PRODUCES AN ARRAY INDEX SO THAT */
/*     THE SEQUENCE NUM(INDEX(I)), I = 1,2,...,N IS SORTED. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE NUMSRT(N,NMAX,NUM,MODE,INDEX,LAST,NEXT) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE. */

/*       NMAX IS A POSITIVE INT INPUT VARIABLE. */

/*       NUM IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS THE */
/*         SEQUENCE OF INTS TO BE GROUPED AND SORTED. IT */
/*         IS ASSUMED THAT THE INTS ARE EACH FROM THE SET */
/*         0,1,...,NMAX. */

/*       MODE IS AN INT INPUT VARIABLE. THE SEQUENCE NUM IS */
/*         SORTED IN ASCENDING ORDER IF MODE IS POSITIVE AND IN */
/*         DESCENDING ORDER IF MODE IS NEGATIVE. IF MODE IS 0, */
/*         NO SORTING IS DONE. */

/*       INDEX IS AN INT OUTPUT ARRAY OF LENGTH N SET SO */
/*         THAT THE SEQUENCE */

/*               NUM(INDEX(I)), I = 1,2,...,N */

/*         IS SORTED ACCORDING TO THE SETTING OF MODE. IF MODE */
/*         IS 0, INDEX IS NOT REFERENCED. */

/*       LAST IS AN INT OUTPUT ARRAY OF LENGTH NMAX + 1. THE */
/*         INDEX OF NUM FOR THE LAST OCCURRENCE OF L IS LAST(L) */
/*         FOR ANY L = 0,1,...,NMAX UNLESS LAST(L) = 0. IN */
/*         THIS CASE L DOES NOT APPEAR IN NUM. */

/*       NEXT IS AN INT OUTPUT ARRAY OF LENGTH N. IF */
/*         NUM(K) = L, THEN THE INDEX OF NUM FOR THE PREVIOUS */
/*         OCCURRENCE OF L IS NEXT(K) FOR ANY L = 0,1,...,NMAX */
/*         UNLESS NEXT(K) = 0. IN THIS CASE THERE IS NO PREVIOUS */
/*         OCCURRENCE OF L IN NUM. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     DETERMINE THE ARRAYS NEXT AND LAST. */

    /* Parameter adjustments */
    --next;
    --index;
    --num;

    /* Function Body */
    i__1 = *nmax;
    for (i__ = 0; i__ <= i__1; ++i__) {
	last[i__] = 0;
/* L10: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = num[k];
	next[k] = last[l];
	last[l] = k;
/* L20: */
    }
    if (*mode == 0) {
	return 0;
    }

/*     STORE THE POINTERS TO THE SORTED ARRAY IN INDEX. */

    i__ = 1;
    if (*mode > 0) {
	jl = 0;
	ju = *nmax;
	jinc = 1;
    } else {
	jl = *nmax;
	ju = 0;
	jinc = -1;
    }
    i__1 = ju;
    i__2 = jinc;
    for (j = jl; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	k = last[j];
L30:
	if (k == 0) {
	    goto L40;
	}
	index[i__] = k;
	++i__;
	k = next[k];
	goto L30;
L40:
/* L50: */
	;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE NUMSRT. */

} /* numsrt_ */

/* Subroutine */ int seq_(int *n, int *indrow, int *jpntr, 
	int *indcol, int *ipntr, int *list, int *ngrp, 
	int *maxgrp, int *iwa)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int j, ic, ip, jp, ir, jcol;

/*     ********** */

/*     SUBROUTINE SEQ */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS */
/*     SUBROUTINE DETERMINES A CONSISTENT PARTITION OF THE */
/*     COLUMNS OF A BY A SEQUENTIAL ALGORITHM. */

/*     A CONSISTENT PARTITION IS DEFINED IN TERMS OF THE LOOPLESS */
/*     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE */
/*     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF */
/*     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     A PARTITION OF THE COLUMNS OF A INTO GROUPS IS CONSISTENT */
/*     IF THE COLUMNS IN ANY GROUP ARE NOT ADJACENT IN THE GRAPH G. */
/*     IN GRAPH-THEORY TERMINOLOGY, A CONSISTENT PARTITION OF THE */
/*     COLUMNS OF A CORRESPONDS TO A COLORING OF THE GRAPH G. */

/*     THE SUBROUTINE EXAMINES THE COLUMNS IN THE ORDER SPECIFIED */
/*     BY THE ARRAY LIST, AND ASSIGNS THE CURRENT COLUMN TO THE */
/*     GROUP WITH THE SMALLEST POSSIBLE NUMBER. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SEQ AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,LIST,NGRP,MAXGRP, */
/*                      IWA) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INT INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INT INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       LIST IS AN INT INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE ORDER TO BE USED BY THE SEQUENTIAL ALGORITHM. */
/*         THE J-TH COLUMN IN THIS ORDER IS LIST(J). */

/*       NGRP IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE PARTITION OF THE COLUMNS OF A. COLUMN JCOL BELONGS */
/*         TO GROUP NGRP(JCOL). */

/*       MAXGRP IS AN INT OUTPUT VARIABLE WHICH SPECIFIES THE */
/*         NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF A. */

/*       IWA IS AN INT WORK ARRAY OF LENGTH N. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa;
    --ngrp;
    --list;
    --jpntr;
    --indrow;
    --indcol;
    --ipntr;

    /* Function Body */
    *maxgrp = 0;
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	ngrp[jp] = *n;
	iwa[jp] = 0;
/* L10: */
    }

/*     BEGINNING OF ITERATION LOOP. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jcol = list[j];

/*        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL. */

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

	i__2 = jpntr[jcol + 1] - 1;
	for (jp = jpntr[jcol]; jp <= i__2; ++jp) {
	    ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	    i__3 = ipntr[ir + 1] - 1;
	    for (ip = ipntr[ir]; ip <= i__3; ++ip) {
		ic = indcol[ip];

/*              ARRAY IWA MARKS THE GROUP NUMBERS OF THE */
/*              COLUMNS WHICH ARE ADJACENT TO COLUMN JCOL. */

		iwa[ngrp[ic]] = j;
/* L20: */
	    }
/* L30: */
	}

/*        ASSIGN THE SMALLEST UN-MARKED GROUP NUMBER TO JCOL. */

	i__2 = *maxgrp;
	for (jp = 1; jp <= i__2; ++jp) {
	    if (iwa[jp] != j) {
		goto L50;
	    }
/* L40: */
	}
	++(*maxgrp);
L50:
	ngrp[jcol] = jp;
/* L60: */
    }

/*        END OF ITERATION LOOP. */

    return 0;

/*     LAST CARD OF SUBROUTINE SEQ. */

} /* seq_ */

/* Subroutine */ int setr_(int *m, int *n, int *indrow, int *
	jpntr, int *indcol, int *ipntr, int *iwa)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int jp, ir, jcol;

/*     ********** */

/*     SUBROUTINE SETR */

/*     GIVEN A COLUMN-ORIENTED DEFINITION OF THE SPARSITY PATTERN */
/*     OF AN M BY N MATRIX A, THIS SUBROUTINE DETERMINES A */
/*     ROW-ORIENTED DEFINITION OF THE SPARSITY PATTERN OF A. */

/*     ON INPUT THE COLUMN-ORIENTED DEFINITION IS SPECIFIED BY */
/*     THE ARRAYS INDROW AND JPNTR. ON OUTPUT THE ROW-ORIENTED */
/*     DEFINITION IS SPECIFIED BY THE ARRAYS INDCOL AND IPNTR. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SETR(M,N,INDROW,JPNTR,INDCOL,IPNTR,IWA) */

/*     WHERE */

/*       M IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF ROWS OF A. */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INT INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INT OUTPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INT OUTPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(1) IS SET TO 1 AND THAT IPNTR(M+1)-1 IS */
/*         THEN THE NUMBER OF NON-ZERO ELEMENTS OF THE MATRIX A. */

/*       IWA IS AN INT WORK ARRAY OF LENGTH M. */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE ROWS. */

    /* Parameter adjustments */
    --iwa;
    --ipntr;
    --jpntr;
    --indrow;
    --indcol;

    /* Function Body */
    i__1 = *m;
    for (ir = 1; ir <= i__1; ++ir) {
	iwa[ir] = 0;
/* L10: */
    }
    i__1 = jpntr[*n + 1] - 1;
    for (jp = 1; jp <= i__1; ++jp) {
	++iwa[indrow[jp]];
/* L20: */
    }

/*     SET POINTERS TO THE START OF THE ROWS IN INDCOL. */

    ipntr[1] = 1;
    i__1 = *m;
    for (ir = 1; ir <= i__1; ++ir) {
	ipntr[ir + 1] = ipntr[ir] + iwa[ir];
	iwa[ir] = ipntr[ir];
/* L30: */
    }

/*     FILL INDCOL. */

    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	i__2 = jpntr[jcol + 1] - 1;
	for (jp = jpntr[jcol]; jp <= i__2; ++jp) {
	    ir = indrow[jp];
	    indcol[iwa[ir]] = jcol;
	    ++iwa[ir];
/* L40: */
	}
/* L50: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SETR. */

} /* setr_ */

/* Subroutine */ int slo_(int *n, int *indrow, int *jpntr, 
	int *indcol, int *ipntr, int *ndeg, int *list, 
	int *maxclq, int *iwa1, int *iwa2, int *iwa3, int 
	*iwa4)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Local variables */
    static int ic, ip, jp, ir, jcol, mindeg, numdeg, numord;

/*     ********** */

/*     SUBROUTINE SLO */

/*     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS */
/*     SUBROUTINE DETERMINES THE SMALLEST-LAST ORDERING OF THE */
/*     COLUMNS OF A. */

/*     THE SMALLEST-LAST ORDERING IS DEFINED FOR THE LOOPLESS */
/*     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE */
/*     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF */
/*     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION. */

/*     THE SMALLEST-LAST ORDERING IS DETERMINED RECURSIVELY BY */
/*     LETTING LIST(K), K = N,...,1 BE A COLUMN WITH LEAST DEGREE */
/*     IN THE SUBGRAPH SPANNED BY THE UN-ORDERED COLUMNS. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SLO AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SLO(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST, */
/*                      MAXCLQ,IWA1,IWA2,IWA3,IWA4) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       INDROW IS AN INT INPUT ARRAY WHICH CONTAINS THE ROW */
/*         INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       JPNTR IS AN INT INPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW. */
/*         THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       INDCOL IS AN INT INPUT ARRAY WHICH CONTAINS THE */
/*         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A. */

/*       IPNTR IS AN INT INPUT ARRAY OF LENGTH M + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL. */
/*         THE COLUMN INDICES FOR ROW I ARE */

/*               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1. */

/*         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO */
/*         ELEMENTS OF THE MATRIX A. */

/*       NDEG IS AN INT INPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN */
/*         OF A IS NDEG(J). */

/*       LIST IS AN INT OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES */
/*         THE SMALLEST-LAST ORDERING OF THE COLUMNS OF A. THE J-TH */
/*         COLUMN IN THIS ORDER IS LIST(J). */

/*       MAXCLQ IS AN INT OUTPUT VARIABLE SET TO THE SIZE */
/*         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING. */

/*       IWA1,IWA2,IWA3, AND IWA4 ARE INT WORK ARRAYS OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... MIN */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. AUGUST 1984. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     INITIALIZATION BLOCK. */

    /* Parameter adjustments */
    --iwa4;
    --iwa3;
    --iwa2;
    --list;
    --ndeg;
    --jpntr;
    --indrow;
    --indcol;
    --ipntr;

    /* Function Body */
    mindeg = *n;
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	iwa1[jp - 1] = 0;
	iwa4[jp] = *n;
	list[jp] = ndeg[jp];
/* Computing MIN */
	i__2 = mindeg, i__3 = ndeg[jp];
	mindeg = GMIN(i__2,i__3);
/* L10: */
    }

/*     CREATE A DOUBLY-LINKED LIST TO ACCESS THE DEGREES OF THE */
/*     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS. */

/*     EACH UN-ORDERED COLUMN IC IS IN A LIST (THE DEGREE LIST) */
/*     OF COLUMNS WITH THE SAME DEGREE. */

/*     IWA1(NUMDEG) IS THE FIRST COLUMN IN THE NUMDEG LIST */
/*     UNLESS IWA1(NUMDEG) = 0. IN THIS CASE THERE ARE */
/*     NO COLUMNS IN THE NUMDEG LIST. */

/*     IWA2(IC) IS THE COLUMN BEFORE IC IN THE DEGREE LIST */
/*     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST */
/*     COLUMN IN THIS DEGREE LIST. */

/*     IWA3(IC) IS THE COLUMN AFTER IC IN THE DEGREE LIST */
/*     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST */
/*     COLUMN IN THIS DEGREE LIST. */

/*     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE */
/*     DEGREE OF IC IN THE GRAPH INDUCED BY THE UN-ORDERED */
/*     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL) */
/*     IS THE SMALLEST-LAST ORDER OF COLUMN JCOL. */

    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	numdeg = ndeg[jp];
	iwa2[jp] = 0;
	iwa3[jp] = iwa1[numdeg];
	if (iwa1[numdeg] > 0) {
	    iwa2[iwa1[numdeg]] = jp;
	}
	iwa1[numdeg] = jp;
/* L20: */
    }
    *maxclq = 0;
    numord = *n;

/*     BEGINNING OF ITERATION LOOP. */

L30:

/*        CHOOSE A COLUMN JCOL OF MINIMAL DEGREE MINDEG. */

L40:
    jcol = iwa1[mindeg];
    if (jcol > 0) {
	goto L50;
    }
    ++mindeg;
    goto L40;
L50:
    list[jcol] = numord;

/*        MARK THE SIZE OF THE LARGEST CLIQUE */
/*        FOUND DURING THE ORDERING. */

    if (mindeg + 1 == numord && *maxclq == 0) {
	*maxclq = numord;
    }

/*        TERMINATION TEST. */

    --numord;
    if (numord == 0) {
	goto L80;
    }

/*        DELETE COLUMN JCOL FROM THE MINDEG LIST. */

    iwa1[mindeg] = iwa3[jcol];
    if (iwa3[jcol] > 0) {
	iwa2[iwa3[jcol]] = 0;
    }

/*        FIND ALL COLUMNS ADJACENT TO COLUMN JCOL. */

    iwa4[jcol] = 0;

/*        DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND */
/*        TO NON-ZEROES IN THE MATRIX. */

    i__1 = jpntr[jcol + 1] - 1;
    for (jp = jpntr[jcol]; jp <= i__1; ++jp) {
	ir = indrow[jp];

/*           FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC) */
/*           WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX. */

	i__2 = ipntr[ir + 1] - 1;
	for (ip = ipntr[ir]; ip <= i__2; ++ip) {
	    ic = indcol[ip];

/*              ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO */
/*              COLUMN JCOL. */

	    if (iwa4[ic] > numord) {
		iwa4[ic] = numord;

/*                 UPDATE THE POINTERS TO THE CURRENT DEGREE LISTS. */

		numdeg = list[ic];
		--list[ic];
/* Computing MIN */
		i__3 = mindeg, i__4 = list[ic];
		mindeg = GMIN(i__3,i__4);

/*                 DELETE COLUMN IC FROM THE NUMDEG LIST. */

		if (iwa2[ic] == 0) {
		    iwa1[numdeg] = iwa3[ic];
		} else {
		    iwa3[iwa2[ic]] = iwa3[ic];
		}
		if (iwa3[ic] > 0) {
		    iwa2[iwa3[ic]] = iwa2[ic];
		}

/*                 ADD COLUMN IC TO THE NUMDEG-1 LIST. */

		iwa2[ic] = 0;
		iwa3[ic] = iwa1[numdeg - 1];
		if (iwa1[numdeg - 1] > 0) {
		    iwa2[iwa1[numdeg - 1]] = ic;
		}
		iwa1[numdeg - 1] = ic;
	    }
/* L60: */
	}
/* L70: */
    }

/*        END OF ITERATION LOOP. */

    goto L30;
L80:

/*     INVERT THE ARRAY LIST. */

    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	iwa2[list[jcol]] = jcol;
/* L90: */
    }
    i__1 = *n;
    for (jp = 1; jp <= i__1; ++jp) {
	list[jp] = iwa2[jp];
/* L100: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SLO. */

} /* slo_ */

/* Subroutine */ int srtdat_(int *n, int *nnz, int *indrow, 
	int *indcol, int *jpntr, int *iwa)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, j, k, l;

/*     ********** */

/*     SUBROUTINE SRTDAT */

/*     GIVEN THE NON-ZERO ELEMENTS OF AN M BY N MATRIX A IN */
/*     ARBITRARY ORDER AS SPECIFIED BY THEIR ROW AND COLUMN */
/*     INDICES, THIS SUBROUTINE PERMUTES THESE ELEMENTS SO */
/*     THAT THEIR COLUMN INDICES ARE IN NON-DECREASING ORDER. */

/*     ON INPUT IT IS ASSUMED THAT THE ELEMENTS ARE SPECIFIED IN */

/*           INDROW(K),INDCOL(K), K = 1,...,NNZ. */

/*     ON OUTPUT THE ELEMENTS ARE PERMUTED SO THAT INDCOL IS */
/*     IN NON-DECREASING ORDER. IN ADDITION, THE ARRAY JPNTR */
/*     IS SET SO THAT THE ROW INDICES FOR COLUMN J ARE */

/*           INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SRTDAT AND IS */
/*     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT. */

/*     THE SUBROUTINE STATEMENT IS */

/*       SUBROUTINE SRTDAT(N,NNZ,INDROW,INDCOL,JPNTR,IWA) */

/*     WHERE */

/*       N IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF COLUMNS OF A. */

/*       NNZ IS A POSITIVE INT INPUT VARIABLE SET TO THE NUMBER */
/*         OF NON-ZERO ELEMENTS OF A. */

/*       INDROW IS AN INT ARRAY OF LENGTH NNZ. ON INPUT INDROW */
/*         MUST CONTAIN THE ROW INDICES OF THE NON-ZERO ELEMENTS OF A. */
/*         ON OUTPUT INDROW IS PERMUTED SO THAT THE CORRESPONDING */
/*         COLUMN INDICES OF INDCOL ARE IN NON-DECREASING ORDER. */

/*       INDCOL IS AN INT ARRAY OF LENGTH NNZ. ON INPUT INDCOL */
/*         MUST CONTAIN THE COLUMN INDICES OF THE NON-ZERO ELEMENTS */
/*         OF A. ON OUTPUT INDCOL IS PERMUTED SO THAT THESE INDICES */
/*         ARE IN NON-DECREASING ORDER. */

/*       JPNTR IS AN INT OUTPUT ARRAY OF LENGTH N + 1 WHICH */
/*         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN THE OUTPUT */
/*         INDROW. THE ROW INDICES FOR COLUMN J ARE */

/*               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1. */

/*         NOTE THAT JPNTR(1) IS SET TO 1 AND THAT JPNTR(N+1)-1 */
/*         IS THEN NNZ. */

/*       IWA IS AN INT WORK ARRAY OF LENGTH N. */

/*     SUBPROGRAMS CALLED */

/*       FORTRAN-SUPPLIED ... MAX */

/*     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983. */
/*     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE' */

/*     ********** */

/*     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE COLUMNS. */

    /* Parameter adjustments */
    --iwa;
    --jpntr;
    --indcol;
    --indrow;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	iwa[j] = 0;
/* L10: */
    }
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	++iwa[indcol[k]];
/* L20: */
    }

/*     SET POINTERS TO THE START OF THE COLUMNS IN INDROW. */

    jpntr[1] = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jpntr[j + 1] = jpntr[j] + iwa[j];
	iwa[j] = jpntr[j];
/* L30: */
    }
    k = 1;

/*     BEGIN IN-PLACE SORT. */

L40:
    j = indcol[k];
    if (k >= jpntr[j]) {

/*           CURRENT ELEMENT IS IN POSITION. NOW EXAMINE THE */
/*           NEXT ELEMENT OR THE FIRST UN-SORTED ELEMENT IN */
/*           THE J-TH GROUP. */

/* Computing MAX */
	i__1 = k + 1, i__2 = iwa[j];
	k = GMAX(i__1,i__2);
    } else {

/*           CURRENT ELEMENT IS NOT IN POSITION. PLACE ELEMENT */
/*           IN POSITION AND MAKE THE DISPLACED ELEMENT THE */
/*           CURRENT ELEMENT. */

	l = iwa[j];
	++iwa[j];
	i__ = indrow[k];
	indrow[k] = indrow[l];
	indcol[k] = indcol[l];
	indrow[l] = i__;
	indcol[l] = j;
    }
    if (k <= *nnz) {
	goto L40;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SRTDAT. */

} /* srtdat_ */


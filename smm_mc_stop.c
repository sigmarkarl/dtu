/* M. Nielsen Mar 2004 mniel@cbs.dtu.dk */

/* 
Copyright (C) 2008-2015 Danish Technical University

This suite of programs and library routine is free software. You can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

In other words, you are free to modify, copy, or redistribute this
source code and its documentation in any way you like, but you must
distribute all derivative versions as free software under the same
terms that I've provided my code to you (i.e. the GNU General Public
License). This precludes any use of the code in proprietary or
commercial software unless your source code is made freely available.

If you wish to use the code under a different Open Source license
that's not compatible with the GPL (like the Artistic License, BSD
license, or the Netscape Public License), please contact me
(Morten Nielsen, mniel@cbs.dtu.dk) for permission.

Incorporation into commercial software under non-GPL terms is possible
by obtaining a specially licensed version from The Danish Technical University.
Contact Morten Nielsen (mniel@cbs.dtu.dk) to arrange licensing terms.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

int	p_seed;
float	p_lambda;
float	p_tmin;
float	p_tmax;
int	p_nt;
int	p_niter;
int	p_verbose;
float	p_dr;
float	p_scale;

PARAM   smm_mc_param[] = {
	"-s", VINT      p_seed, "Seed [-1] Default, [0] Time [>0] Specific seed", "-1",
	"-l", VFLOAT	p_lambda, "Lambda for minization", "0.05",
	"-ts", VFLOAT   p_tmax, "Start temperature", "0.001",
        "-te", VFLOAT   p_tmin, "Final temperature", "0.000001",
        "-nt", VINT     p_nt, "Number of temperature steps", "5",
	"-i", VINT      p_niter, "Number of iterations per T cycle", "2000",
	"-v", VSWITCH	p_verbose, "Verbose mode", "0",
	"-r", VFLOAT	p_dr, "Step size for MC move", "0.1",
	"-scale", VFLOAT p_scale, "1/(Energy scale) for MC move", "100.0",
	0
};

float   fvector_xyerror_nonorm2( int n, float *v1, float *v2 ) {
        int     i;
        float   err, tmp;

        err = 0.0;

        for ( i=0;i<n;i++ ) {
                tmp = v1[i] - v2[i];
                err += tmp*tmp;
        }

        return( err );
}

float	cal_error( float *target, float **inp, float *w, int n, float lambda, float *p, int nin )

/* Calculate error between target and predictec values */

/* 

	E = 0.5 * sum_i ( p_i - t_i )^2 + lambda` * sum_k w_k^2

	E = sum_i E_i

	E_i = 0.5 * ( ( p_i - t_i )^2 + lambda * sum_k w_k^2 )

	where lambda = lambda'/n is the per-target lambda value used in the smm gradient
	decent algorirthm

*/

{

	float	e;
	int	i;
	float	lambda_prime;

	e = 0.0;
	lambda_prime = lambda * n;

	/* Calculate prediction values */
	for ( i=0; i<n; i++ )
		p[i] = fvector_dot( w, inp[i], nin );

	/* Calculate error */
	e = 0.5 * fvector_xyerror_nonorm2( n, p, target ) + lambda_prime * fvector_dot( w, w, nin);

	return( e );
}

void print_matrix( float *w, int alen, int len )

{
	int	k,j;

	printf( "\nLast position-specific scoring matrix computed       \n");

	printf( "      " );
        for (k=0;k<alen;k++)
                printf( "%7c ", PROFILE_ORDER[k] );
        printf( "\n");

        for (j=0;j<len;j++){
                printf( "%3d A ",j+1);
                for (k=0;k<alen;k++)
                        printf( "%7.3f ", w[k+j*alen]);
                printf( "\n");
        }

}

int smm_mc_main( int argc, char *argv[] ) {
	PEPLIST	*peplist, *peplist2, *pl;
	int		len, alen;
	int	i,j, k;
	float	dt;
	int	nn, nacc, iter;
	float	e, e2, ne, ne2, de, de2, fe, fe2;
	float	old_mat1, old_mat2;
	int	ip1, ip2, ix, it;
	float	dr;
	float	*w, **inp, **inp2, *target, *target2, *p, *p2;
	float	t, dmat;
	float	best_err;
	int	n,n2,nin, c;

	pparse( &argc, &argv, smm_mc_param, 1, "trainfile" );

	/* Initialize seed for random number generator */
        if ( p_seed >= 0 )
                setseed( p_seed );

	/* Read peptides */
	peplist = peplist_read( argv[1] );
	peplist2 = peplist_read( argv[2] );

	if ( peplist == NULL ) {
		printf( "Error. Not elements in peplist %s\n", argv[1] );
		exit( 1 );
	}

	len = peplist->len;
	alen = 20;

	/* Find number of peptides */
	for ( n=0,pl = peplist; pl; pl=pl->next, n++ );
	for ( n2=0,pl = peplist2; pl; pl=pl->next, n2++ );

	nin = peplist->len * alen;

	printf( "# NIN %i N %i\n", nin, n );

	/* Allocate vectors for smm matrix, target, and prediction */
	/* Allocate SMM matrix. The matrix is allocated as a vector of size len*20 where len is the peptide length */

	w = fvector( 0, nin-1 );
	target = fvector( 0, n-1 );
	target2 = fvector( 0, n2-1 );
	p = fvector( 0, n-1 );
	p2 = fvector( 0, n2-1 );

	/* Allocate matrix to store input. Each row is a peptide */
        /* For 9mer peptides each row contain 180 values. Each amino acids is
        sparse encode. That is each amino acids is encoded as
                A 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                R 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                N 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                .....
                V 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1

        */

	inp = fmatrix( 0, n-1, 0, nin-1 );
	inp2 = fmatrix( 0, n2-1, 0, nin-1 );

	for ( i=0, pl = peplist; pl; pl=pl->next, i++ ) {

		if ( pl->len != len ) {
			printf( "Error. %s pl->len %i != len %i\n", pl->pep, pl->len, len );
			exit( 1 );
		}

		ix = 0;
		for ( j=0; j<pl->len; j++ ) {
			c = strpos( PROFILE_ORDER, pl->pep[j] );

			if ( c < 0 ) {
				printf( "Error. Unknown character %c %s\n", pl->pep[j], PROFILE_ORDER );
				exit( 1 );
			}

			for ( k=0; k<20; k++ ) {
                                inp[i][ix] = ( c == k ? 1.0 : 0 );
                                ix++;
                        } 
		}

                target[i] = pl->score;

	}

	for ( i=0, pl = peplist2; pl; pl=pl->next, i++ ) {
		if ( pl->len != len ) {
			printf( "Error. %s pl->len %i != len %i\n", pl->pep, pl->len, len );
			exit( 1 );
		}

		ix = 0;
		for ( j=0; j<pl->len; j++ ) {
			c = strpos( PROFILE_ORDER, pl->pep[j] );

			if ( c < 0 ) {
				printf( "Error. Unknown character %c %s\n", pl->pep[j], PROFILE_ORDER );
				exit( 1 );
			}

			for ( k=0; k<20; k++ ) {
				inp2[i][ix] = ( c == k ? 1.0 : 0 );
				ix++;
			}
		}

		target2[i] = pl->score;
	}

	/* Assign random initial values between -0.1 and 0.1 for the w vector */
	ix = 0;	
	for ( i=0; i<len; i++ )
	for ( j=0; j<alen; j++ ) 
		w[ix++] = 0.1* ( 2*drand48() - 1 );

	/* Define size of MC move */
	dr = p_dr; 

	/* Calculate error */
	e = cal_error( target, inp, w, n, p_lambda, p, nin );
	e2 = cal_error( target2, inp2, w, n2, p_lambda, p2, nin );

	printf( "# Initial E %f PCC %f\n", e, fvector_xycorr( nin, target, p ));

	/* Define temperature step for cooling */
	dt = ( p_tmax - p_tmin)/(p_nt > 0 ? p_nt : 1 );

	/* Define best (lowest) error */
	best_err = 9999.9;

	/* Do temperature loop */
	for ( it=0, t=p_tmax; it <= p_nt; it++, t-=dt ) {

		/* Do iteration loop for each temperature sted */
		for ( iter=0;iter<p_niter;iter++) {

			/* Select random weight */
			ip1 = (int)( drand48()*nin );

			/* Save old matrix value */
			old_mat1 = w[ip1];

			/* Find random weight position not equal to old weight position */
			ip2 = (int)( drand48()*nin );	
			while ( ip2 == ip1 )
				ip2 = (int)( drand48()*nin );

			/* Save old matrix value */
			old_mat2 = w[ip2];
	
			/* Define size of MC move */
			dmat = (2*drand48() - 1)*dr;

			/* Update matrix values. */
			w[ip1] += dmat;
			w[ip2] -= dmat;

			/* Calculate error value after MC move */
			ne = cal_error( target, inp, w, n, p_lambda, p, nin );
			ne2 = cal_error( target2, inp2, w, n2, p_lambda, p2, nin );

			/* Scale is a factor multiplied to the de values to place dE and temperature
                        and a comparable scale, like 1/kB the Boltzman factor of statistical mechanics */

			de = (ne - e)*p_scale;
			de2 = (ne2 - e2)*p_scale;

			/* if de < 0 or boltzmann > random number, accept move */ 
			/* Make sure you get the sign of the energy right */

			if( de < 0 || ( exp( -de/t ) > drand48() ) ) {
				nacc++;
				e = ne;
				e2 = ne2;

				fe = fvector_xyerror( n, p, target );
				fe2 = fvector_xyerror( n2, p2, target2 );
			}
			else {
				/* Reject move and reset w matrix to old values */
				w[ip1] = old_mat1;
				w[ip2] = old_mat2;
			}

			nn++;

			/* Update MC step size */
			if ( nn > 100 ) {
				if ( nacc*1.0/nn > 0.5 )
					dr*=1.1;
				else
					dr*=0.9;
				nn = 0;
				nacc = 0;

				if ( p_verbose )
					printf( "# DR %f\n", dr );
			}

			if ( p_verbose ) 
				printf( "# NT %i ITR %i E %f E2 %f PCC %f\n", it, iter, fe, fe2, fvector_xycorr( n, target, p ) );

		}
	}

	printf( "# Best energy %f PCC %f SME %f\n", e, fvector_xycorr( n, target, p ), fvector_xyerror( n, target, p ) );

	/* Print syn matrix */
	print_matrix( w, alen, len );

	exit( 0 );
}

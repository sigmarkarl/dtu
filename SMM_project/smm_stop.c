/* M.Nielsen October 2009 mniel@cbs.dtu.dk */

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

int	p_verbose;
int     p_seed;
int	p_ncycles;
float	p_eta;
float	p_lambda, p_eta2lambda;
float   p_bplim;

PARAM   smm_param[] = {
	"-v", VSWITCH	p_verbose,	"Verbose mode", "0",
        "-s", VINT      p_seed, "Seed [-1] Default, [0] Time [>0] Specific seed", "-1",
	"-nc", VINT     p_ncycles, "Number of iterations", "500",
	"-eta", VFLOAT	p_eta, "Eta for weight update", "0.01",
	"-bp", VFLOAT	p_bplim, "Limit for backpropagation", "0.00001",
	"-l",  VFLOAT   p_lambda, "Lambda", "0.05",
	0
};

float   fvector_xyerror_nonorm( int n, float *v1, float *v2 )

{
        int     i;
        float   err, tmp;

        err = 0.0;

        for ( i=0;i<n;i++ ) {
                tmp = v1[i] - v2[i];
                err += tmp*tmp;
        }

        return( err );
}

inline void	update_weight( float *w, float *inp, float d_o, int n ) {
	int	j;

	/* Update weights */

	/* w = w - eta * dw
	
	dw[j] = dE/dw_j

	E = 1/2 * (p - t)^2 + lambda * sum_j w[j]*w[j] = 1/2 * d_o^2 + lambda * sum_j w[j]*w[j]

	*/

	float p_etad_o = p_eta * d_o;
	for ( j=0; j<n; j++ ) {
		w[j] -= p_etad_o * inp[j] + p_eta2lambda * w[j];
	}
}

int smm_main( int argc, char *argv[] ) {

	float	*w;
	int	i, j;
	float	**inp, **inp2;
	float	*target, *target2, *p, *p2;
	int	 nin, n, n2;
	int	ncy;
    float   err, pcc, tpcc, gerr;
	int	*order, ix;
	float	d_o;
	PEPLIST	*peplist, *peplist2, *pl;
	int	c, k;
	int	alen, plen;

	pparse( &argc, &argv, smm_param, 1, "inputfile" );

	/* Define seed for random number generator */
	if ( p_seed >= 0 )
		setseed( p_seed );

	/* Read peptide input */
	peplist = peplist_read( argv[1] );
	peplist2 = peplist_read( argv[2] );

	if ( ! peplist ) {
		printf( "Error. No peptides read from file %s. Exit\n", argv[1] );
	}

	/* Count number of peptides read */
	for ( n=0,pl=peplist; pl; pl=pl->next, n++ );
	for ( n2=0,pl=peplist2; pl; pl=pl->next, n2++ );

	alen = 20; 
	plen = peplist->len;

	/* nin is equal to the number of input values */
	nin = plen * alen;

	/* Allocate matrix to store input. Each row is a peptide */
	/* For 9mer peptides each row contains 180 values. Each amino acid is
	sparse encode. That is each amino acids is encoded as
		A 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
		R 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
		N 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
		.....
		V 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1

	*/
			
	inp = fmatrix( 0, n-1, 0, nin - 1 );
	inp2 = fmatrix( 0, n2-1, 0, nin - 1 );

	/* Allocate vector for target values */
	target = fvector( 0, n-1 );
	target2 = fvector( 0, n2-1 );

	for ( i=0,pl=peplist; pl; pl=pl->next, i++ ) {

		ix = 0;

		for ( j=0; j<plen; j++ ) {

			c = strpos( PROFILE_ORDER, pl->pep[j] );

			if ( c < 0 ) {
				printf( "Error. Sequence letter unknown %c %s\n",
                                        	pl->pep[j],  PROFILE_ORDER);
                               	exit( 1 );
                       	}

			/* Do the sparse encoding */

			for ( k=0; k<20; k++ ) {
                inp[i][ix] = ( c == k ? 1.0 : 0 );
				ix++;
			}

		}

		target[i] = pl->score;
	}

	for ( i=0,pl=peplist2; pl; pl=pl->next, i++ ) {

			ix = 0;

			for ( j=0; j<plen; j++ ) {

				c = strpos( PROFILE_ORDER, pl->pep[j] );

				if ( c < 0 ) {
					printf( "Error. Sequence letter unknown %c %s\n",
	                                        	pl->pep[j],  PROFILE_ORDER);
	                               	exit( 1 );
	                       	}

				/* Do the sparse encoding */

				for ( k=0; k<20; k++ ) {
	                inp2[i][ix] = ( c == k ? 1.0 : 0 );
					ix++;
				}

			}

			target2[i] = pl->score;
		}

	printf( "# Lambda %f\n", p_lambda );
	p_eta2lambda = p_eta * 2 * p_lambda;

	printf( "# NIN %i N %i\n", nin, n );

	/* Allocate SMM matrix. The matrix is allocated as a vector of size len*20 where len is the peptide length
	This format makes the gradient decent more easy */

	w = fvector( 0, nin-1 );

	if ( w == NULL) {
                printf( "Error. Cannot allocate SYNAPS\n");
                exit(1);
        }

	/* Assign random initial values between -0.1 and 0.1 for the w vector */
	for ( i=0; i<nin; i++ ) {
		w[i] = 0.2 * ( 0.5 - drand48());

		printf("%f\n",w[i]);
	}

	/* Allocate vector to store prediction values for each peptide */
	p = fvector( 0, n-1 );
	p2 = fvector( 0, n-1 );

	/* Allocate vector to store order in which peptides are used for training
	This order is changed before each training cycle. The ivector_ramp routine
	returns a vector with values order[0] = 0, order[1] = 1, order[2] = 2, ... 
	*/

	order = ivector_ramp( 0, n-1 );

	for ( ncy=0;ncy<p_ncycles; ncy++ ) {

		/* Change the order of the peptides */
		ivector_rerandomize( order, 0, n-1 );

		for ( i = 0; i < n; i++ ) {

			ix = order[i];

			// Calculate the predicted binding
			p[ix] = fvector_dot( inp[ix], w, nin );

			// Calculate the error of the prediction d_o = p - t
			d_o = p[ix] - target[ix];

			// Only update weights if error is greather than p_bplim
			if ( SQR( d_o ) > p_bplim )
				update_weight( w, inp[ix], d_o, nin );

		}

		if ( p_verbose ) {
			for ( i = 0; i < n2; i++ ) {
				// Calculate the predicted binding
				p2[i] = fvector_dot( inp2[i], w, nin );

				// Calculate the error of the prediction d_o = p - t
				//d_o = p2[ix] - target2[ix];

				// Only update weights if error is greather than p_bplim
				//if ( SQR( d_o ) > p_bplim )
				//	update_weight( w, inp2[ix], d_o, nin );
			}

			pcc = fvector_xyerror(n, p, target);
			tpcc = fvector_xyerror( n2, p2, target2 );
			//pcc = fvector_xycorr( n2, p, target );
			//tpcc = fvector_xycorr( n2, p2, target2 );
            err = fvector_xyerror( n, p, target );

			/* Note that the lambda used here is the per_target normalized lambda value
			To get the global lambda value you must multiply with n */

			gerr = 0.5 * fvector_xyerror_nonorm( n, p, target ) + 
				n * p_lambda * fvector_dot( w, w, nin );

			printf( "Ncycle %i PCC %f TPCC %f Err %f TERR %f\n", ncy, pcc, tpcc, err, gerr );
		}
	}

	for ( i=0; i<n; i++ )
		p[i] = fvector_dot( inp[i], w, nin );

	pcc = fvector_xycorr( n, p, target );
	err = fvector_xyerror( n, p, target );
	gerr = 0.5 * fvector_xyerror_nonorm( n, p, target ) + 
		n * p_lambda * fvector_dot( w, w, nin );

	/* Write out SMM vector in matrix form */

	printf( "# Ncyle %i PCC %f Err %f MSE %f\n", p_ncycles, pcc, gerr, err );
	printf( "\nLast position-specific scoring matrix computed       \n");

        printf( "      " );
        for (i=0;i<20;i++)
                printf( "%7c ", PROFILE_ORDER[i]);
        printf( "\n");

	ix = 0;
	for ( i=0; i<plen; i++ ) {

                printf( "%3d A ",i+1);

                for (j=0; j<20; j++ ) 
			printf( "%7.4f ", w[i*20+j] );

                printf( "\n");
        }

	exit( 0 );
}

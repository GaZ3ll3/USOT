/*
 * BCAssem.c
 *
 *  Created on: Oct 6, 2014
 *      Author: lurker
 */
#include "mex.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

/*
 *
 *  Boundary Condition Assembly.
 *
 *  Do integral of test function and trial function on boundary.
 *
 *  Input : Edges
 *
 *  Output: [I, J, V] for building sparse matrix
 *
 */



void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	REAL *np = mxGetPr(prhs[0]);
	REAL *ne = mxGetPr(prhs[1]);
	REAL *nodes = mxGetPr(prhs[2]);
	int *edges = (int *)mxGetData(prhs[3]);
	/* each boundary node takes 2 edges.*/
	int nzmax = 4*(int)*ne;

	plhs[0] = mxCreateNumericMatrix(nzmax, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(nzmax, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(nzmax, 1,  mxDOUBLE_CLASS, mxREAL);

	REAL *pI = mxGetPr(plhs[0]);
	REAL *pJ = mxGetPr(plhs[1]);
	REAL *pV = mxGetPr(plhs[2]);

	mwSize i, j, k;
	mwSize node_1, node_2;
	REAL length;

#pragma omp parallel for
	for (i = 0; i < *ne; i++)
	{
		node_1 = edges[2*i];
		node_2 = edges[2*i + 1];

		length = sqrt(pow(nodes[2*node_1] - nodes[2*node_2],2) + pow(nodes[2*node_1 + 1] - nodes[2*node_2 + 1], 2));


		for (j = 0; j < 2; j++)
		{
			for (k = 0; k < 2; k++)
			{
				*pI++ = edges[2*i + j] + 1;
				*pJ++ = edges[2*i + k] + 1;

				if (j != k)
				{
					*pV++ = length/3.0;
				}
				else
				{
					*pV++ = length/6.0;
				}

			}
		}
	}
}

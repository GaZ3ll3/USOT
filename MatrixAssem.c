/*
 * MatrixAssem.c
 *
 *  Created on: Oct 1, 2014
 *      Author: lurker
 */

#include "mex.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*
 *  Generate COO format sparse matrix vectors
 *
 *  [I, J, V, U] = MatrixAssem(Mesh Property)
 *
 *  I : rows
 *  J : cols
 *  V : vals
 *
 *  nzmax as length of V.
 */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	/*
	 *  Number of nodes
	 *  Number of triangles
	 *  All nodes
	 *  All triangles, use int32.
	 */
	REAL *np = mxGetPr(prhs[0]);
	REAL *nt = mxGetPr(prhs[1]);
	REAL *nodes = mxGetPr(prhs[2]);
	int *tri = (int *)mxGetData(prhs[3]);
	int nzmax = 9*(int)*nt;

	plhs[0] = mxCreateNumericMatrix(nzmax, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(nzmax, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(nzmax, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[3] = mxCreateNumericMatrix(nzmax, 1,  mxDOUBLE_CLASS, mxREAL);

	REAL *pI = mxGetPr(plhs[0]);
	REAL *pJ = mxGetPr(plhs[1]);
	REAL *pV = mxGetPr(plhs[2]);
	REAL *pU = mxGetPr(plhs[3]);


	mwSize i, j, k;
	mwSize node_1, node_2, node_3;
	double grad[2][3];
	REAL area;
	REAL det;

#pragma omp parallel for
	for (i = 0; i < *nt; i++)
	{
		node_1 = tri[3*i];
		node_2 = tri[3*i + 1];
		node_3 = tri[3*i + 2];

		/*
		 * already ordered
		 */
		det = (nodes[node_2*2] - nodes[node_1*2])*(nodes[node_3*2 + 1] - nodes[node_1*2 + 1]) -
				(nodes[node_2*2 + 1] - nodes[node_1*2 + 1])*(nodes[node_3*2] - nodes[node_1*2]);
		area = 0.5*fabs(det);

		/*
		 * scaled gradient
		 */
		grad[0][0] = (nodes[node_2*2 + 1] - nodes[node_3*2 + 1]);
		grad[1][0] = (nodes[node_3*2] - nodes[node_2*2]);

		grad[0][1] = (nodes[node_3*2 + 1] - nodes[node_1*2 + 1]);
		grad[1][1] = (nodes[node_1*2] - nodes[node_3*2]);

		grad[0][2] = (nodes[node_1*2 + 1] - nodes[node_2*2 + 1]);
		grad[1][2] = (nodes[node_2*2] - nodes[node_1*2]);

		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				*pI++ = tri[3*i + j] + 1 ;

				*pJ++ = tri[3*i + k] + 1;

				*pU++ = (grad[0][j]*grad[0][k] + grad[1][j]*grad[1][k])/4.0/area;

				if (j != k)
				{
					*pV++ = area*1.0/12.0;
				}
				else
				{
					*pV++ = area*1.0/6.0;
				}
			}

		}
	}

	/*
	 * Output as three vectors, length is approximated (number of triangles)*3
	 */






}




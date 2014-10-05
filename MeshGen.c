/*****************************************************************************/
/* Modified from tricall.(Example code from Triangle)
/*****************************************************************************/

/* If SINGLE is defined when triangle.o is compiled, it should also be       */
/*   defined here.  If not, it should not be defined here.                   */

/* #define SINGLE */

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include "mex.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "triangle/triangle.h"

/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/

void report(io, markers, reporttriangles, reportneighbors, reportsegments,
            reportedges, reportnorms)
struct triangulateio *io;
int markers;
int reporttriangles;
int reportneighbors;
int reportsegments;
int reportedges;
int reportnorms;
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}

/*****************************************************************************/
/*                                                                           */
/*  main()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/

void  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  struct triangulateio in, mid, out;

  /* Define input points. */

  char *options = mxArrayToString(prhs[0]);
  REAL *Nodes = mxGetPr(prhs[1]);


  



  in.numberofpoints = (size_t) mxGetN(prhs[1]);
  in.numberofpointattributes = 0;
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

  mwSize i;

  for (i = 0; i < in.numberofpoints; i++)
  {
  	in.pointlist[2*i] = Nodes[2*i];
  	in.pointlist[2*i + 1] = Nodes[2*i + 1];
  }


/*
 * input arguments
 *
 *
 */ 
  in.pointmarkerlist = NULL;
  in.numberofsegments = 0;
  in.numberofholes = 0;
  in.numberofregions = 0;
  in.regionlist = NULL;

/*
 *
 * coarse triangluation
 * 
 * output argument setup
 */
  mid.pointlist = (REAL *) NULL;            
  mid.pointmarkerlist = (int *) NULL; 
  mid.trianglelist = (int *) NULL;          
  mid.triangleattributelist = (REAL *) NULL;
  mid.neighborlist = (int *) NULL;         
  mid.segmentlist = (int *) NULL;
  mid.segmentmarkerlist = (int *) NULL;
  mid.edgelist = (int *) NULL;             
  mid.edgemarkerlist = (int *) NULL;   

 /*
  * p: polygon
  * c: hull
  * z: from zero
  */
  triangulate("pcz", &in, &mid, (struct triangulateio *) NULL);

/*
 *
 * Refinement setup
 *
 */
  mid.trianglearealist = NULL;
  out.pointlist = (REAL *) NULL;            
  out.pointattributelist = (REAL *) NULL;
  out.trianglelist = (int *) NULL;          
  out.triangleattributelist = (REAL *) NULL;
  out.edgelist = (int *) NULL;             
  out.edgemarkerlist = (int *) NULL;
  out.segmentlist = (int *) NULL;

  /*
   *  option example as "prq30.0a0.1zBCV"
   */
  triangulate(options, &mid, &out, (struct triangulateio *) NULL);

  /*printf("Refined triangulation:\n\n");
  report(&out, 0, 1, 0, 0, 0, 0);*/

  /*
   * export to MATLAB
   *
   */
   plhs[0] = mxCreateDoubleScalar(mxREAL);
   *mxGetPr(plhs[0]) = out.numberofpoints;
   plhs[1] = mxCreateDoubleScalar(mxREAL);
   *mxGetPr(plhs[1]) = out.numberoftriangles;
   plhs[2] = mxCreateDoubleScalar(mxREAL);
   *mxGetPr(plhs[2]) = out.numberofsegments;;
   plhs[3] = mxCreateDoubleMatrix(2, out.numberofpoints, mxREAL);
   memcpy(mxGetPr(plhs[3]), out.pointlist,out.numberofpoints*2*sizeof(REAL));
   plhs[4] = mxCreateNumericMatrix(3, out.numberoftriangles, mxINT32_CLASS, mxREAL);
   memcpy(mxGetPr(plhs[4]), out.trianglelist ,out.numberoftriangles*3*sizeof(int));
   plhs[5] = mxCreateNumericMatrix(2, out.numberofsegments,mxINT32_CLASS, mxREAL);
   memcpy(mxGetPr(plhs[5]), out.segmentlist ,out.numberofsegments*2*sizeof(int));




  /* Free all allocated arrays, including those allocated by Triangle. */

  free(in.pointlist);
  free(in.pointmarkerlist);
  free(in.regionlist);
  free(mid.pointlist);
  free(mid.pointmarkerlist);
  free(mid.trianglelist);
  free(mid.triangleattributelist);
  free(mid.trianglearealist);
  free(mid.neighborlist);
  free(mid.segmentlist);
  free(mid.segmentmarkerlist);
  free(mid.edgelist);
  free(mid.edgemarkerlist);
  free(out.pointlist);
  free(out.pointattributelist);
  free(out.trianglelist);
  free(out.triangleattributelist);
  free(out.segmentlist);
}

/*
 * Mesh.cc
 *
 *  Created on: Oct 8, 2014
 *      Author: lurker
 */
#include <mexplus.h>
#include <iostream>
#include <pprint.h>
#include "Mesh.h"


#define ANSI_DECLARATORS 1

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#define VOID int

extern "C"
{
#include "triangle/triangle.h"
}


using namespace std;
using namespace mexplus;


class Mesh {
public:
	Mesh(mxArray* boundary);
	virtual ~Mesh();
	/*
	 * Members
	 */
	triangulateio MeshData;
	/*
	 * Methods
	 */
	void Refine();
private:
	void clear();
};


Mesh::Mesh(mxArray* boundary){
	size_t coordinate_length = mxGetM(boundary)/2;
	auto boundary_ptr = mxGetPr(boundary);
	/*
	 * Two step meshing
	 */
	struct triangulateio  input, mid;
	input.numberofpoints = coordinate_length;
	input.numberofpointattributes = 0;
	input.pointlist = (REAL *) malloc(coordinate_length * 2 * sizeof(REAL));

	for (size_t i = 0; i < coordinate_length; i++)
	{
		input.pointlist[2*i] = *(boundary_ptr++);
		input.pointlist[2*i + 1] = *(boundary_ptr++);
	}

	input.pointmarkerlist = nullptr;
	input.numberofsegments = 0;
	input.numberofholes = 0;
	input.numberofregions = 0;
	input.regionlist = nullptr;

	mid.pointlist = (REAL *) nullptr;
	mid.pointmarkerlist = (int *) nullptr;
	mid.trianglelist = (int *) nullptr;
	mid.triangleattributelist = (REAL *) nullptr;
	mid.neighborlist = (int *) nullptr;
	mid.segmentlist = (int *) nullptr;
	mid.segmentmarkerlist = (int *) nullptr;
	mid.edgelist = (int *) nullptr;
	mid.edgemarkerlist = (int *) nullptr;


/*
 * Coarse meshing
 */
	triangulate("pcz", &input, &mid, (struct triangulateio *) nullptr);
/*
 * Fine meshing
 */
	mid.trianglearealist = nullptr;
	MeshData.pointlist = (REAL *) nullptr;
	MeshData.pointattributelist = (REAL *) nullptr;
	MeshData.trianglelist = (int *) nullptr;
	MeshData.triangleattributelist = (REAL *) nullptr;
	MeshData.edgelist = (int *) nullptr;
	MeshData.edgemarkerlist = (int *) nullptr;
	MeshData.segmentlist = (int *) nullptr;
	triangulate("prq30.0a0.00001zBC", &mid, &MeshData, (struct triangulateio *) nullptr);
/*
 * Clean up
 */
	free(input.pointlist);
	free(input.pointmarkerlist);
	free(input.regionlist);
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
} // Mesh


Mesh::~Mesh()
{
	clear();
}
void Mesh::Refine(){

} // Refine

void Mesh::clear(){
	free(MeshData.pointlist);
	MeshData.pointlist = nullptr;
	free(MeshData.pointattributelist);
	MeshData.pointattributelist = nullptr;
	free(MeshData.trianglelist);
	MeshData.trianglelist = nullptr;
	free(MeshData.triangleattributelist);
	MeshData.triangleattributelist = nullptr;
	free(MeshData.segmentlist);
	MeshData.segmentlist = nullptr;
	free(MeshData.edgemarkerlist);
	MeshData.edgemarkerlist = nullptr;
	free(MeshData.edgelist);
	MeshData.edgelist = nullptr;
} // clear

template class mexplus::Session<Mesh>;

namespace {

// Create a new instance of Mesh and return its session id.
MEX_DEFINE(new) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Mesh>::create(new Mesh(const_cast<mxArray*>(input.get(0)))));
}

// Delete the instance by id
MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Mesh>::destroy(input.get(0));
}

MEX_DEFINE(report)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	const Mesh& mesh = Session<Mesh>::getConst(input.get(0));
	mexPrintf("Mesh Generated:\n");
	mexPrintf("\tPoints  : %15d\n", mesh.MeshData.numberofpoints);
	mexPrintf("\tElements: %15d\n", mesh.MeshData.numberoftriangles);
	mexPrintf("\tEdges   : %15d\n", mesh.MeshData.numberofedges);
	mexPrintf("\tSegments: %15d\n", mesh.MeshData.numberofsegments);
	mexPrintf("\tRegions : %15d\n", mesh.MeshData.numberofregions);
	mexPrintf("\tHoles   : %15d\n", mesh.MeshData.numberofholes);
	mexPrintf("\tCorners : %15d\n", mesh.MeshData.numberofcorners);

}

} // namespace

MEX_DISPATCH



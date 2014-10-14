/*
 * Mesh.cc
 *
 *  Created on: Oct 8, 2014
 *      Author: lurker
 */

#include "Mesh.h"


Mesh::Mesh(mxArray* boundary, mxArray* min_area){
	/*
	 * Import use-specified Minimal Area
	 */
	Min_Area = *(double *)(mxGetData(min_area));
	topology = new Topology;


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

	triangulate(const_cast<char*>(("prq30.0a" + std::to_string(Min_Area) + "ezBC").c_str()), &mid, &MeshData, (struct triangulateio *) nullptr);
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


Mesh::~Mesh(){
	clear();
}

void Mesh::Refine() {
/*
 * mid point refinement:
 *
 * each edges will produce a new node to insert
 */

}

/* TODO:
 *
 * Need to promote all boundary edges.
 */

void Mesh::Promote(int deg){

	/*
	 *  Clear topology for next use
	 */
	topology->nodes.clear();
	topology->elems.clear();
	topology->edges.clear();
	topology->boundary.clear();

	double P1_Coord_X, P1_Coord_Y, P2_Coord_X,P2_Coord_Y;
	double Delta_X, Delta_Y;
	if (deg <= 0) {mexWarnMsgTxt("Degree has to be positive. Use first order instead.\n");}
	else if (deg >= 20) {mexPrintf("Degree Too Large\n");}
	else {
		int id = (deg + 1) * (deg + 2)/2;
		/*
		 *  Allocate memory for topology
		 */
		topology->nodes.resize(2*( MeshData.numberofpoints +
				(deg - 1)*MeshData.numberofedges +
				(id - 3*deg)*MeshData.numberoftriangles));

		topology->elems.resize(((deg + 1) * (deg + 2)/2) * MeshData.numberoftriangles);

		topology->boundary.resize((deg + 1) * MeshData.numberofsegments);

		/*
		 * Fill nodes
		 */
		for (size_t i = 0; i < MeshData.numberofpoints; i++){
			topology->nodes[2*i] = *(MeshData.pointlist + 2*i);
			topology->nodes[2*i + 1] = *(MeshData.pointlist + 2*i + 1);
		}

		for (size_t i = 0; i < MeshData.numberofedges ; i++){
			P1_Coord_X = *(MeshData.pointlist + 2*(*(MeshData.edgelist + 2*i)));
			P1_Coord_Y = *(MeshData.pointlist + 2*(*(MeshData.edgelist + 2*i)) + 1);
			P2_Coord_X = *(MeshData.pointlist + 2*(*(MeshData.edgelist + 2*i + 1)));
			P2_Coord_Y = *(MeshData.pointlist + 2*(*(MeshData.edgelist + 2*i + 1)) + 1);

			Delta_X = (P2_Coord_X - P1_Coord_X)/(double)deg;
			Delta_Y = (P2_Coord_Y - P1_Coord_Y)/(double)deg;

			vector<int> EdgePoints(deg - 1);

			for (size_t j = 0; j < deg - 1; j++){
				topology->nodes[2*(MeshData.numberofpoints + i*(deg - 1) + j)] =
						P1_Coord_X + Delta_X*(j + 1);
				topology->nodes[2*(MeshData.numberofpoints + i*(deg - 1) + j) + 1] =
						P1_Coord_Y + Delta_Y*(j + 1);

				EdgePoints[j] = MeshData.numberofpoints + i*(deg - 1) + j;
			}

			topology->edges.insert(make_pair(std::to_string(*(MeshData.edgelist + 2*i))
			+ "_" + std::to_string(*(MeshData.edgelist + 2*i + 1)), EdgePoints));
		}

		size_t counter = 2*MeshData.numberofpoints + 2*(deg - 1)*MeshData.numberofedges;

		std::unordered_map<std::string, std::vector<int>>::const_iterator Iterator;

		for (size_t index = 0; index < MeshData.numberoftriangles; index++){
			/*
			 * Insert vertex
			 */
			topology->elems[(deg + 1)*(deg + 2)/ 2 * index ] = *(MeshData.trianglelist + 3*index);
			topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 1 ] = *(MeshData.trianglelist + 3*index + 1);
			topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 2 ] = *(MeshData.trianglelist + 3*index + 2);



			Iterator = topology->edges.find(std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index ]) + "_" +
							std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 1 ]));

			if (Iterator != topology->edges.end()){
				std::copy(Iterator->second.begin(), Iterator->second.end(), topology->elems.begin() + (deg + 1)*(deg + 2)/ 2 * index + 3);
			}

			if (Iterator == topology->edges.end()){
				Iterator = topology->edges.find(std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 1]) + "_" +
								std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index  ]));
				std::reverse_copy(Iterator->second.begin(), Iterator->second.end(), topology->elems.begin() + (deg + 1)*(deg + 2)/ 2 * index + 3);
			}



			Iterator = topology->edges.find(std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 1]) + "_" +
							std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 2 ]));
			if (Iterator != topology->edges.end()){
				std::copy(Iterator->second.begin(), Iterator->second.end(), topology->elems.begin() + (deg + 1)*(deg + 2)/ 2 * index + 3 + (deg - 1));
			}


			if (Iterator == topology->edges.end()){
				Iterator = topology->edges.find(std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 2]) + "_" +
								std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 1 ]));
				std::reverse_copy(Iterator->second.begin(), Iterator->second.end(), topology->elems.begin() + (deg + 1)*(deg + 2)/ 2 * index + 3 + (deg - 1));
			}

			Iterator = topology->edges.find(std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 2]) + "_" +
							std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index ]));

			if (Iterator != topology->edges.end()){
				std::copy(Iterator->second.begin(), Iterator->second.end(), topology->elems.begin() + (deg + 1)*(deg + 2)/ 2 * index + 3 + 2*(deg - 1));
			}

			if (Iterator == topology->edges.end()){
				Iterator = topology->edges.find(std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index ]) + "_" +
								std::to_string(topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 2 ]));
				std::reverse_copy(Iterator->second.begin(), Iterator->second.end(), topology->elems.begin() + (deg + 1)*(deg + 2)/ 2 * index + 3 + 2*(deg - 1));

			}




			size_t internal_counter = 0;
			for (size_t i  = 1; i < deg - 1; i++){
				for (size_t j = 1; j < deg - i; j++){
					size_t k = deg - i - j;

					topology->elems[(deg + 1)*(deg + 2)/ 2 * index + 3*deg + internal_counter] = counter/2;
					internal_counter += 1;

					topology->nodes[counter] =
							(double)k/(double)deg * (*(MeshData.pointlist + 2*(*(MeshData.trianglelist + 3*index)))) +
							(double)i/(double)deg * (*(MeshData.pointlist + 2*(*(MeshData.trianglelist + 3*index + 1)))) +
							(double)j/(double)deg * (*(MeshData.pointlist + 2*(*(MeshData.trianglelist + 3*index + 2))));

					topology->nodes[counter + 1] =
							(double)k/(double)deg * (*(MeshData.pointlist + 2*(*(MeshData.trianglelist + 3*index)) + 1)) +
							(double)i/(double)deg * (*(MeshData.pointlist + 2*(*(MeshData.trianglelist + 3*index + 1)) + 1)) +
							(double)j/(double)deg * (*(MeshData.pointlist + 2*(*(MeshData.trianglelist + 3*index + 2)) + 1));
					counter += 2;
				}
			}
		}

		for (size_t index = 0; index < MeshData.numberofsegments; index++){
			/*
			 * Orientation reverse
			 */
			topology->boundary[(deg + 1)*index + 1] = *(MeshData.segmentlist + 2*index);
			topology->boundary[(deg + 1)*index    ] = *(MeshData.segmentlist + 2*index + 1);
			/*
			 * In order insert the interpolated nodes.
			 */
			Iterator = topology->edges.find(std::to_string(topology->boundary[(deg + 1)*index]) + "_" +
										std::to_string(topology->boundary[(deg + 1)*index + 1]));

			if (Iterator != topology->edges.end()){
				std::copy(Iterator->second.begin(), Iterator->second.end(), topology->boundary.begin() + (deg + 1)*index + 2);
			}

			if (Iterator == topology->edges.end()){
				Iterator = topology->edges.find(std::to_string(topology->boundary[(deg + 1)*index + 1]) + "_" +
								std::to_string(topology->boundary[(deg + 1)*index]));
				std::reverse_copy(Iterator->second.begin(), Iterator->second.end(), topology->boundary.begin() + (deg + 1)*index + 2);
			}
		}// End of for
	}
}

void Mesh::clear(){
	/*
	 *  Should not free same pointer twice.
	 */
	if (MeshData.pointlist != nullptr) {free(MeshData.pointlist); MeshData.pointlist = nullptr;}
	if (MeshData.pointattributelist != nullptr) {free(MeshData.pointattributelist); MeshData.pointattributelist = nullptr;}
	if (MeshData.trianglelist != nullptr) {free(MeshData.trianglelist); MeshData.trianglelist = nullptr;}
	if (MeshData.triangleattributelist != nullptr) {free(MeshData.triangleattributelist); MeshData.triangleattributelist = nullptr;}
	if (MeshData.segmentlist != nullptr) {free(MeshData.segmentlist); MeshData.segmentlist = nullptr;}
	if (MeshData.edgemarkerlist != nullptr) {free(MeshData.edgemarkerlist); MeshData.edgemarkerlist = nullptr;}
	if (MeshData.edgelist != nullptr) {free(MeshData.edgelist); MeshData.edgelist = nullptr;}
} // clear

template class mexplus::Session<Mesh>;


namespace {

// Create a new instance of Mesh and return its session id.
MEX_DEFINE(new) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Mesh>::create(new Mesh(const_cast<mxArray*>(input.get(0)),
			const_cast<mxArray*>(input.get(1)))));
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
	Mesh* mesh = Session<Mesh>::get(input.get(0));

	bool verbose = false;

	mexPrintf("Mesh Generated:\n");
	mexPrintf("\tPoints  : %15d\n", (*mesh).MeshData.numberofpoints);
	mexPrintf("\tElements: %15d\n", (*mesh).MeshData.numberoftriangles);
	mexPrintf("\tEdges   : %15d\n", (*mesh).MeshData.numberofedges);
	mexPrintf("\tSegments: %15d\n", (*mesh).MeshData.numberofsegments);
	mexPrintf("\tRegions : %15d\n", (*mesh).MeshData.numberofregions);
	mexPrintf("\tHoles   : %15d\n", (*mesh).MeshData.numberofholes);
	mexPrintf("\tCorners : %15d\n", (*mesh).MeshData.numberofcorners);

	if (verbose){

	/*
	 * Check triangles, if they are oriented.
	 */


		for (size_t index = 0; index < mesh->MeshData.numberoftriangles; index++){
		   if ((mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index + 2] + 1] -
				   mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index] + 1])*
		   (mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index + 1]    ] -
				   mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index]   ]) -
		   (mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index + 2]    ] -
				   mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index]   ])*
		   (mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index + 1] + 1] -
				   mesh->MeshData.pointlist[2*mesh->MeshData.trianglelist[3*index] + 1]) > 0){

			   mexPrintf("Oriented\n");
		   }
		   else{
			   mexPrintf("UnOriented\n");
		   }
		}

	/*
	 * Check segments, if they are oriented.
	 */


	}


}

// TODO
MEX_DEFINE(promote)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 3);
	Mesh* mesh = Session<Mesh>::get(input.get(0));
	int deg = input.get<int>(1);
	mesh->Promote(deg);
//	std::cout << (*mesh).topology->elems << std::endl;
	plhs[0] = mxCreateNumericMatrix(2, (*mesh).topology->nodes.size()/2, mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[0]), &(*mesh).topology->nodes[0], (*mesh).topology->nodes.size()*sizeof(REAL));
	plhs[1] = mxCreateNumericMatrix((deg + 1)*(deg + 2)/2, (*mesh).topology->elems.size()/((deg + 1)*(deg + 2)/2), mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[1]), &(*mesh).topology->elems[0], (*mesh).topology->elems.size()*sizeof(int));


	plhs[2] = mxCreateNumericMatrix((deg + 1), (*mesh).topology->boundary.size()/((deg + 1)), mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[2]), &(*mesh).topology->boundary[0], (*mesh).topology->boundary.size()*sizeof(int));
}


MEX_DEFINE(export)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 7);
	Mesh* mesh = Session<Mesh>::get(input.get(0));
	plhs[0] = MxArray::from((*mesh).MeshData.numberofpoints);
	plhs[1] = MxArray::from((*mesh).MeshData.numberoftriangles);
	plhs[2] = MxArray::from((*mesh).MeshData.numberofsegments);
	plhs[3] = mxCreateNumericMatrix(2, (*mesh).MeshData.numberofpoints, mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[3]), (*mesh).MeshData.pointlist,(*mesh).MeshData.numberofpoints*2*sizeof(REAL));
	plhs[4] = mxCreateNumericMatrix(3, (*mesh).MeshData.numberoftriangles, mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[4]), (*mesh).MeshData.trianglelist ,(*mesh).MeshData.numberoftriangles*3*sizeof(int));
	plhs[5] = mxCreateNumericMatrix(2, (*mesh).MeshData.numberofsegments, mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[5]), (*mesh).MeshData.segmentlist ,(*mesh).MeshData.numberofsegments*2*sizeof(int));
	plhs[6] = mxCreateNumericMatrix(2, (*mesh).MeshData.numberofedges, mxINT32_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[6]), (*mesh).MeshData.edgelist ,(*mesh).MeshData.numberofedges*2*sizeof(int));
}
} // namespace


MEX_DISPATCH



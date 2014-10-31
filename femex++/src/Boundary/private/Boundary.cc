/*
 * Boundary.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: lurker
 */

#include "Boundary.h"

namespace MEX {

Boundary::Boundary(MatlabPtr _type, MatlabPtr _edges) {
	// TODO Auto-generated constructor stub


	b_type = static_cast<Boundary_t>((int32_t)*Matlab_Cast<Real_t>(_type));
	b_edges.resize(mxGetM(_edges)*mxGetN(_edges));
	memcpy(&b_edges[0], mxGetPr(_edges), mxGetM(_edges)*mxGetN(_edges)*sizeof(int32_t));
}

Boundary::~Boundary() {
	// TODO Auto-generated destructor stub
}

} /* namespace MEX */



using namespace MEX;

template class mexplus::Session<Boundary>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Boundary>::create(new Boundary(
			const_cast<MatlabPtr>(input.get(0)),
						const_cast<MatlabPtr>(input.get(1))
						)));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Boundary>::destroy(input.get(0));
}

MEX_DEFINE(export) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto boundary = Session<Boundary>::get(input.get(0));

	mexPrintf("%d\n", boundary->b_type);
	for (size_t i = 0; i < boundary->b_edges.size(); i++){
		mexPrintf("%d\n", boundary->b_edges[i]);
	}
}
}

MEX_DISPATCH

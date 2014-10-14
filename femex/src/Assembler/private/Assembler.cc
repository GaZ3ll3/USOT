/*
 * Assembler.cc
 *
 *  Created on: Oct 9, 2014
 *      Author: lurker
 */
#include "Assembler.h"


#define MAKE_VALUE(pointer) \
  shared_ptr<mxArray>(pointer, mxDestroyArray)
// Initialize vector<const mxArray*> rhs.
#define MAKE_RHS(rhs, ...) \
  const shared_ptr<mxArray> kFixture[] = { __VA_ARGS__ }; \
  const size_t kFixtureSize = sizeof(kFixture) / sizeof(shared_ptr<mxArray>); \
  vector<const mxArray*> rhs(kFixtureSize); \
  for (int i = 0; i < kFixtureSize; ++i) \
    rhs[i] = kFixture[i].get();
// Initialize vector<mxArray*> lhs.
#define MAKE_LHS(lhs, ...) \
  const shared_ptr<mxArray> kFixture[] = { __VA_ARGS__ }; \
  const size_t kFixtureSize = sizeof(kFixture) / sizeof(shared_ptr<mxArray>); \
  vector<mxArray*> lhs(kFixtureSize); \
  for (int i = 0; i < kFixtureSize; ++i) \
    lhs[i] = kFixture[i].get();

Assembler::Assembler(){


}

Assembler::~Assembler() {

}

template class mexplus::Session<Assembler>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Assembler>::create(new Assembler()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Assembler>::destroy(input.get(0));
}

MEX_DEFINE(evalNodalInfo)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofpoints = mxGetN(prhs[1]);
	size_t numberofqpoints = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);


	mxArray* Vander = mxCreateNumericMatrix(numberofpoints,numberofpoints,mxDOUBLE_CLASS, mxREAL);
	mxArray* VanderF = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	mxArray* VanderX = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	mxArray* VanderY = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);

	int deg = round((sqrt(8*numberofpoints + 1) - 3)/2);

	if ((deg + 1) * (deg + 2)/2 != numberofpoints){
		mexErrMsgTxt("Invalid length of input nodes\n");
	}

	// extract all nodes from promoted nodes.
	double* Vander_ptr = mxGetPr(Vander);
	double* VanderF_ptr = mxGetPr(VanderF);
	double* VanderX_ptr = mxGetPr(VanderX);
	double* VanderY_ptr = mxGetPr(VanderY);

	double* nodes_ptr = mxGetPr(prhs[1]);
	double* qnodes_ptr = mxGetPr(prhs[2]);

	/*
	 * nodes is 2 x N matrix
	 *
	 * qnodes is 2 x M matrix, with column majority
	 */

	/*
	 * Assembler use x decreasing order
	 */
	for (size_t col = 0; col < numberofpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*Vander_ptr++ = pow(*(nodes_ptr + 2*col), i - j)*pow(*(nodes_ptr + 2*col + 1), j);
			}
		}
	}

	for (size_t col = 0; col < numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*VanderF_ptr++ = pow(*(qnodes_ptr + 2*col), i - j)*pow(*(qnodes_ptr + 2*col + 1), j);
			}
		}
	}

	for (size_t col = 0; col < numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == i){
					*VanderX_ptr++ = 0.;
				}
				else{
					*VanderX_ptr++ = (i - j) * pow(*(qnodes_ptr + 2*col), i - j - 1)*pow(*(qnodes_ptr + 2*col + 1), j);
				}
			}
		}
	}

	for (size_t col = 0; col < numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == 0){
					*VanderY_ptr++ = 0.;
				}
				else{
					*VanderY_ptr++ = j* pow(*(qnodes_ptr + 2*col), i - j)*pow(*(qnodes_ptr + 2*col + 1), j - 1);
				}
			}
		}
	}

	mxArray* RHS_f[2];
	RHS_f[0] = Vander;
	RHS_f[1] = VanderF;

	mxArray** LHS_f = &plhs[0];
	mexCallMATLAB(1, LHS_f, 2, RHS_f, "mldivide");

	mxArray* RHS_x[2];
	RHS_x[0] = Vander;
	RHS_x[1] = VanderX;

	mxArray** LHS_x = &plhs[1];
	mexCallMATLAB(1, LHS_x, 2, RHS_x, "mldivide");

	mxArray** RHS_y;
	RHS_y =  (mxArray **)mxCalloc(2, sizeof(mxArray*));
	RHS_y[0] = Vander;
	RHS_y[1] = VanderY;

	mxArray** LHS_y = &plhs[2];
	mexCallMATLAB(1, LHS_y, 2, RHS_y, "mldivide");

	mxDestroyArray(Vander);
	mxDestroyArray(VanderF);
	mxDestroyArray(VanderX);
	mxDestroyArray(VanderY);
}

MEX_DEFINE(assems)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	/*
	 * all nodes
	 */
	double* pnodes_ptr = mxGetPr(prhs[1]);
	/*
	 * all elems
	 */
	int* pelem_ptr = (int*)mxGetPr(prhs[2]);
	size_t numberofelem = mxGetN(prhs[2]);
	size_t numberofnodesperelem = mxGetM(prhs[2]);

	/*
	 * nodes in a reference elem
	 * evaluated at qnodes
	 *
	 * better take it as Num of qnodes * Num of basis
	 */
	double* reference_x = mxGetPr(prhs[3]);
	double* reference_y = mxGetPr(prhs[4]);
	size_t numberofqnodes = mxGetM(prhs[3]);

	if (numberofqnodes != mxGetM(prhs[4])){
		mexErrMsgTxt("Invalid Input: Reference element gradient size does not match.\n");
	}
	/*
	 * Num of qnodes
	 */
	double* weights   = mxGetPr(prhs[5]);

	/*
	 * Num of qnodes * Num of elems
	 */
//	double* ExternalFunction = mxGetPr(prhs[5]);
//
//

	/*
	 * I
	 */
	plhs[0] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);
	/*
	 * J
	 */
	plhs[1] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);
	/*
	 * V
	 */
	plhs[2] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);


	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area;
	double Jacobian[2][2];

	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i];
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1];
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2];

		Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
		Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
		Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
		Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		for (size_t j = 0; j < numberofnodesperelem; j++){
			for (size_t k = 0; k < numberofnodesperelem; k++){
				*pI = *(pelem_ptr + i*numberofnodesperelem + j) + 1;
				*pJ = *(pelem_ptr + i*numberofnodesperelem + k) + 1;
				*pV = 0.;
				for (size_t l = 0; l < numberofqnodes; l++){
					*pV = *pV + (
							(Jacobian[0][0]*reference_x[j*numberofqnodes + l] + Jacobian[0][1]*reference_y[j*numberofqnodes + l])*
							(Jacobian[0][0]*reference_x[k*numberofqnodes + l] + Jacobian[0][1]*reference_y[k*numberofqnodes + l]) +

							(Jacobian[1][0]*reference_x[j*numberofqnodes + l] + Jacobian[1][1]*reference_y[j*numberofqnodes + l])*
							(Jacobian[1][0]*reference_x[k*numberofqnodes + l] + Jacobian[1][1]*reference_y[k*numberofqnodes + l])
							)*weights[l]/4.0/area;
				}
				pI++; pJ++; pV++;
			}
		}
	}
}

MEX_DEFINE(assema)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 5);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	/*
	 * all nodes
	 */
	double* pnodes_ptr = mxGetPr(prhs[1]);
	/*
	 * all elems
	 */
	int* pelem_ptr = (int*)mxGetPr(prhs[2]);
	size_t numberofelem = mxGetN(prhs[2]);
	size_t numberofnodesperelem = mxGetM(prhs[2]);

	/*
	 * nodes in a reference elem
	 * evaluated at qnodes
	 *
	 * better take it as Num of qnodes * Num of basis
	 */
	double* reference = mxGetPr(prhs[3]);
	size_t numberofqnodes = mxGetM(prhs[3]);
	/*
	 * Num of qnodes
	 */
	double* weights   = mxGetPr(prhs[4]);

	/*
	 * Num of qnodes * Num of elems
	 */
//	double* ExternalFunction = mxGetPr(prhs[5]);
//
//

	/*
	 * I
	 */
	plhs[0] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);
	/*
	 * J
	 */
	plhs[1] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);
	/*
	 * V
	 */
	plhs[2] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);


	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area;

	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i];
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1];
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2];

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		for (size_t j = 0; j < numberofnodesperelem; j++){
			for (size_t k = 0; k < numberofnodesperelem; k++){
				*pI = *(pelem_ptr + i*numberofnodesperelem + j) + 1;
				*pJ = *(pelem_ptr + i*numberofnodesperelem + k) + 1;
				*pV = 0.;
				for (size_t l = 0; l < numberofqnodes; l++){
					*pV = *pV + reference[j*numberofqnodes + l]*reference[k*numberofqnodes + l]*weights[l]*area;
				}
				pI++; pJ++; pV++;
			}
		}
	}


}

MEX_DEFINE(assemsa)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));
	/*
	 * do both at same time, if it is needed.(always do).
	 */
}



}


MEX_DISPATCH



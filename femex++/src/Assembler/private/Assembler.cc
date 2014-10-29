/*
 * Assembler.cc
 *
 *  Created on: Oct 9, 2014
 *      Author: lurker
 */
#include "Assembler.h"

Assembler::Assembler(){


}

Assembler::~Assembler() {

}


/*
 * Extract information of basis on Reference Triangle
 */
void Assembler::Reference(MatlabPtr &F, MatlabPtr &DX, MatlabPtr &DY,
		MatlabPtr Points, MatlabPtr QPoints){

	std::size_t _numberofpoints = mxGetN(Points);
	std::size_t _numberofqpoints = mxGetN(QPoints);

	/*
	 *  Temporary Arrays, will be destroyed later.
	 */
	mxArray* Vander = mxCreateNumericMatrix(_numberofpoints,_numberofpoints,mxDOUBLE_CLASS, mxREAL);
	mxArray* VanderF = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	mxArray* VanderX = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	mxArray* VanderY = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);


	int deg = round((sqrt(8*_numberofpoints + 1) - 3)/2);

	if ((deg + 1) * (deg + 2)/2 != _numberofpoints){
		mexErrMsgTxt("Invalid length of input nodes\n");
	}

	// extract all nodes from promoted nodes.
	Real_t* Vander_ptr = mxGetPr(Vander);
	Real_t* VanderF_ptr = mxGetPr(VanderF);
	Real_t* VanderX_ptr = mxGetPr(VanderX);
	Real_t* VanderY_ptr = mxGetPr(VanderY);

	Real_t* nodes_ptr = mxGetPr(Points);
	Real_t* qnodes_ptr = mxGetPr(QPoints);

	// basis on all nodes

	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*Vander_ptr++ = pow(nodes_ptr[2*col], i - j)*pow(nodes_ptr[2*col + 1], j);
			}
		}
	}



	// basis on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*VanderF_ptr++ = pow(qnodes_ptr[2*col], i - j)*pow(qnodes_ptr[2*col + 1], j);
			}
		}
	}

	// partial x on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == i){
					*VanderX_ptr++ = 0.;
				}
				else{
					*VanderX_ptr++ = (i - j) * pow(qnodes_ptr[2*col], i - j - 1)*pow(qnodes_ptr[2*col + 1], j);
				}
			}
		}
	}

	// partial y on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == 0){
					*VanderY_ptr++ = 0.;
				}
				else{
					*VanderY_ptr++ = j* pow(qnodes_ptr[2*col], i - j)*pow(qnodes_ptr[2*col + 1], j - 1);
				}
			}
		}
	}



	mxArray* RHS_f[] = {Vander, VanderF};
	mexCallMATLAB(1, &F, 2, RHS_f, "mldivide");

	mxArray* RHS_x[] = {Vander, VanderX};
	mexCallMATLAB(1, &DX, 2, RHS_x, "mldivide");

	mxArray* RHS_y[] = {Vander, VanderY};
	mexCallMATLAB(1, &DY, 2, RHS_y, "mldivide");


	mxDestroyArray(Vander);
	mxDestroyArray(VanderF);
	mxDestroyArray(VanderX);
	mxDestroyArray(VanderY);


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

MEX_DEFINE(reference)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofpoints = mxGetN(prhs[1]);
	size_t numberofqpoints = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);

	assembler->Reference(plhs[0], plhs[1], plhs[2], const_cast<MatlabPtr>(prhs[1]), const_cast<MatlabPtr>(prhs[2]));
}

MEX_DEFINE(assems)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	/*
	 * all nodes
	 */
	Real_t* pnodes_ptr = mxGetPr(prhs[1]);
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
	Real_t* reference_x = mxGetPr(prhs[3]);
	Real_t* reference_y = mxGetPr(prhs[4]);
	size_t numberofqnodes = mxGetM(prhs[3]);

	if (numberofqnodes != mxGetM(prhs[4])){
		mexErrMsgTxt("Invalid Input: Reference element gradient size does not match.\n");
	}
	/*
	 * Num of qnodes
	 */
	Real_t* weights   = mxGetPr(prhs[5]);

	/*
	 * Num of qnodes * Num of elems
	 */

	/*
	 * Pointer not initialized
	 */
	Real_t* Interp = (Real_t*) nullptr;

	if (mxIsFunctionHandle(prhs[6])){
		mexErrMsgTxt("Handle not acceptable yet.\n");
	}
	else if(mxIsDouble(prhs[6])){
		size_t InterpM = mxGetM(prhs[6]);
		size_t InterpN = mxGetN(prhs[6]);

		if (InterpM != numberofqnodes){
			mexErrMsgTxt("Invalid Input: fifth argument' M not sufficient.\n");
		}
		if (InterpN != numberofelem){
			mexErrMsgTxt("Invalid Input: fifth argument' N not sufficient.\n");
		}

		Interp = mxGetPr(prhs[6]);
	}
	else{
		mexErrMsgTxt("The fifth Argument has to be a function handle or a vector");
	}

	/*
	 * I
	 */
	plhs[0] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);
	/*
	 * J
	 */
	plhs[1] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);
	/*
	 * V
	 */
	plhs[2] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];

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
					*pV = *pV + Interp[i*numberofqnodes + l]*(
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
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	/*
	 * all nodes
	 */
	Real_t* pnodes_ptr = mxGetPr(prhs[1]);
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
	Real_t* reference = mxGetPr(prhs[3]);
	size_t numberofqnodes = mxGetM(prhs[3]);
	/*
	 * Num of qnodes
	 */
	Real_t* weights   = mxGetPr(prhs[4]);

	/*
	 * Num of qnodes * Num of elems
	 */

	/*
	 * Pointer not initialized
	 */
	Real_t* Interp = (Real_t*) nullptr;

	if (mxIsFunctionHandle(prhs[5])){
		mexErrMsgTxt("Handle not acceptable yet.\n");
	}
	else if(mxIsDouble(prhs[5])){
		size_t InterpM = mxGetM(prhs[5]);
		size_t InterpN = mxGetN(prhs[5]);

		if (InterpM != numberofqnodes){
			mexErrMsgTxt("Invalid Input: fifth argument' M not sufficient.\n");
		}
		if (InterpN != numberofelem){
			mexErrMsgTxt("Invalid Input: fifth argument' N not sufficient.\n");
		}

		Interp = mxGetPr(prhs[5]);
	}
	else{
		mexErrMsgTxt("The fifth Argument has to be a function handle or a vector");
	}


	/*
	 * I
	 */
	plhs[0] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);
	/*
	 * J
	 */
	plhs[1] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);
	/*
	 * V
	 */
	plhs[2] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;

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
					*pV = *pV + Interp[i*numberofqnodes + l] * reference[j*numberofqnodes + l]*reference[k*numberofqnodes + l]*weights[l]*area;
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



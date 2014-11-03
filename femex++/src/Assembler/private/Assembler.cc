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

	auto _numberofpoints = mxGetN(Points);
	auto _numberofqpoints = mxGetN(QPoints);

	/*
	 *  Temporary Arrays, will be destroyed later.
	 */
	auto Vander = mxCreateNumericMatrix(_numberofpoints,_numberofpoints,mxDOUBLE_CLASS, mxREAL);
	auto VanderF = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderX = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderY = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);


	int deg = round((sqrt(8*_numberofpoints + 1) - 3)/2);

	if ((deg + 1) * (deg + 2)/2 != _numberofpoints){
		mexErrMsgTxt("Invalid length of input nodes\n");
	}

	// extract all nodes from promoted nodes.
	auto Vander_ptr  = mxGetPr(Vander);
	auto VanderF_ptr = mxGetPr(VanderF);
	auto VanderX_ptr = mxGetPr(VanderX);
	auto VanderY_ptr = mxGetPr(VanderY);

	auto nodes_ptr   = mxGetPr(Points);
	auto qnodes_ptr  = mxGetPr(QPoints);

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

void Assembler::AssembleMass(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr Fcn){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;


	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		// Due to symmetric property, only need half of the work load.
		for (size_t j = 0; j < numberofnodesperelem; j++){
			for (size_t k = 0; k < j + 1; k++){
				*pI = pelem_ptr[i*numberofnodesperelem + j];
				*pJ = pelem_ptr[i*numberofnodesperelem + k];
				*pV = 0.;
				for (size_t l = 0; l < numberofqnodes; l++){
					*pV = *pV + Interp[i*numberofqnodes + l]*
							reference[j+ l*numberofnodesperelem]*
							reference[k+ l*numberofnodesperelem]*
							weights[l];
				}
				*pV  = (*pV)*area;

				pI++; pJ++; pV++;
				if (k != j) {
					*pI = *(pJ - 1);
					*pJ = *(pI - 1);
					*pV = *(pV - 1);
					pI++; pJ++; pV++;
				}

			}
		}
	}
}

void Assembler::Qnodes(Real_t*& Coords, MatlabPtr Nodes, MatlabPtr QNodes, MatlabPtr Elems){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(QNodes);


	mwSize vertex_1, vertex_2 , vertex_3;

	for (size_t i = 0; i < numberofelem; i++) {

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		for (size_t l = 0; l < numberofqnodes; l++) {
			Coords[2*(i*numberofqnodes + l)] =
					pnodes_ptr[2 * vertex_1] *(1 - qnodes_ptr[2*l] - qnodes_ptr[2*l + 1]) +
					pnodes_ptr[2 * vertex_2]*qnodes_ptr[2*l] +
					pnodes_ptr[2 * vertex_3]*qnodes_ptr[2*l + 1];
			Coords[2*(i*numberofqnodes + l) + 1] =
					pnodes_ptr[2 * vertex_1 + 1] *(1 - qnodes_ptr[2*l] - qnodes_ptr[2*l + 1]) +
					pnodes_ptr[2 * vertex_2 + 1]*qnodes_ptr[2*l] +
					pnodes_ptr[2 * vertex_3 + 1]*qnodes_ptr[2*l + 1];
		}
	}
}


// calculate integral on boundary
void Assembler::AssembleLoad(Real_t*& pLoad, MatlabPtr Nodes,
		MatlabPtr QNodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr Fcn) {

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area, tmp;


	// Fcn cannot be a function handle, too slow
	auto Fcn_ptr = Matlab_Cast<Real_t>(Fcn);
	// linear interpolation
	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		if (mxGetNumberOfElements(Fcn) != mxGetN(Nodes)) {
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp +=  (Fcn_ptr[vertex_1] *(1 - qnodes_ptr[2*l] - qnodes_ptr[2*l + 1]) +
							 Fcn_ptr[vertex_2]*qnodes_ptr[2*l] +
							 Fcn_ptr[vertex_3]*qnodes_ptr[2*l + 1])*reference[j+ l*numberofnodesperelem]*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += tmp*area;
			}
		}
		else {
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp +=  Fcn_ptr[i*numberofqnodes + l]*reference[j+ l*numberofnodesperelem]*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += tmp*area;
			}
		}
	}
}


void Assembler::AssembleStiff(Real_t* &pI, Real_t* &pJ, Real_t*&pV,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr RefX,
		MatlabPtr RefY, MatlabPtr Weights, MatlabPtr Fcn) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];

	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
		Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
		Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
		Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

		// Orientation corrected.
		det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
		area = 0.5*fabs(det);

		// Due to symmetric property, half of work load can be reduced
		for (size_t j = 0; j < numberofnodesperelem; j++){
			for (size_t k = 0; k < j + 1; k++){
				*pI = pelem_ptr[i*numberofnodesperelem + j];
				*pJ = pelem_ptr[i*numberofnodesperelem + k];
				*pV = 0.;
				for (size_t l = 0; l < numberofqnodes; l++){
					*pV = *pV + Interp[i*numberofqnodes + l] *(
							(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
							(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
							+
							(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
							(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
							)*weights[l];
				}
				*pV = (*pV)/4.0/area;
				pI++; pJ++; pV++;
				if (j != k){
					*pI = *(pJ - 1);
					*pJ = *(pI - 1);
					*pV = *(pV - 1);
					pI++; pJ++; pV++;
				}
			}
		}
	}
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


	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler->AssembleStiff(pI, pJ, pV, const_cast<MatlabPtr>(prhs[1]),
			const_cast<MatlabPtr>(prhs[2]), const_cast<MatlabPtr>(prhs[3]),
			const_cast<MatlabPtr>(prhs[4]), const_cast<MatlabPtr>(prhs[5]),
			const_cast<MatlabPtr>(prhs[6]));



}

MEX_DEFINE(assema)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(1, numberofnodesperelem * numberofnodesperelem * numberofelem, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler->AssembleMass(pI, pJ, pV, const_cast<MatlabPtr>(prhs[1]),
			const_cast<MatlabPtr>(prhs[2]), const_cast<MatlabPtr>(prhs[3]),
			const_cast<MatlabPtr>(prhs[4]), const_cast<MatlabPtr>(prhs[5]));


}

MEX_DEFINE(assemsa)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 9);
	OutputArguments output(nlhs, plhs, 4);
	Assembler* assembler = Session<Assembler>::get(input.get(0));
	/*
	 * do both at same time, if it is needed.(always do).
	 */
}

MEX_DEFINE(asseml) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   = mxGetN(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes,1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pLoad = mxGetPr(plhs[0]);
	assembler->AssembleLoad(pLoad, const_cast<MatlabPtr>(prhs[1]),
			const_cast<MatlabPtr>(prhs[2]), const_cast<MatlabPtr>(prhs[3]),
			const_cast<MatlabPtr>(prhs[4]), const_cast<MatlabPtr>(prhs[5]),
			const_cast<MatlabPtr>(prhs[6]));

}

MEX_DEFINE(qnodes)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));


	size_t numberofelem           = mxGetN(prhs[3]);
	size_t numberofqnodes         = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(2, numberofelem * numberofqnodes,  mxDOUBLE_CLASS, mxREAL);
	Real_t* Coords = mxGetPr(plhs[0]);
	assembler->Qnodes(Coords,const_cast<MatlabPtr>(prhs[1]), const_cast<MatlabPtr>(prhs[2]),
			const_cast<MatlabPtr>(prhs[3]));
}
}


MEX_DISPATCH



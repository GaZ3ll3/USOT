/*
 * Integrator.h
 *
 *  Created on: Oct 10, 2014
 *      Author: lurker
 */

#ifndef INTEGRATOR_PRIVATE_INTEGRATOR_CC_
#define INTEGRATOR_PRIVATE_INTEGRATOR_CC_


#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */


#include <iostream>
#include <iterator>
#include <string.h>

#include <mexplus.h>
#include <pprint.h>

#include <vector>
#include <cstdlib>

#include "utils.h"

using namespace std;
using namespace mexplus;

class Integrator {
public:
	Integrator(int dim, int Degree) {
		prec = Degree;
		if (dim == 2){
			_dim = dim;
			QuadratureData();
		}
		else if (dim == 1){
			_dim = dim;
			GaussData();
		}
		else {
			_dim = 0;
			mexErrMsgTxt("Dimension exceeds maximum 2. \n");
		}
	}
	virtual ~Integrator() {if (_dim != 0) {clear();}}

	/*
	 * public members
	 */
	std::size_t _dim;
	std::vector<Real_t> qwts;
	std::vector<Real_t> qpts;
	std::size_t prec;

	/*
	 * Evaluate nodal basis function on nodes
	 */
private:
	void QuadratureData();
	void GaussData();
	void clear();

};

#endif /* INTEGRATOR_PRIVATE_INTEGRATOR_CC_ */

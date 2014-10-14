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


using namespace std;
using namespace mexplus;

class Integrator {
public:
	Integrator(int Degree) {
		prec = Degree;
		QuadratureData();
	}
	virtual ~Integrator() {}

	/*
	 * public members
	 */
	std::vector<double> qwts;
	std::vector<double> qpts;
	std::size_t prec;

	/*
	 * Evaluate nodal basis function on nodes
	 */
	void QuadratureData();

};

#endif /* INTEGRATOR_PRIVATE_INTEGRATOR_CC_ */

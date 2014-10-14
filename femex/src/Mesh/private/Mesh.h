/*
 * Mesh.h
 *
 *  Created on: Oct 10, 2014
 *      Author: lurker
 */

#ifndef MESH_PRIVATE_MESH_H_
#define MESH_PRIVATE_MESH_H_

#include <iostream>
#include <iterator>
#include <string.h>

#include <unordered_map>
#include <unordered_set>

#include <mexplus.h>
#include <pprint.h>


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


#define FEMEX_EXPECT(condition) if (!(condition)) \
    mexErrMsgTxt(#condition " not true.")


typedef	struct Topology {
		vector<double> nodes;
		unordered_map<std::string, vector<int>> edges;
		vector<int>    elems;
} Topology;

class Mesh {
public:
	Mesh(mxArray* boundary, mxArray* min_area);
	virtual ~Mesh();
/*
 * Members
 */
	struct triangulateio MeshData;
	REAL Min_Area;
	Topology* topology;
/*
 * Methods
 *
 * Refine: generate linear Lagrange element mesh
 *
 *
 */
	void Refine();
	void Promote(int);

private:
	void clear();
};


#endif /* MESH_PRIVATE_MESH_H_ */

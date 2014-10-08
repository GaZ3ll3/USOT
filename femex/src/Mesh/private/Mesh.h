/*
 * Mesh.h
 *
 *  Created on: Oct 8, 2014
 *      Author: lurker
 */

#ifndef MESH_H_
#define MESH_H_


using namespace std;

template <class T> class Matrix
{
public:
	T* data;
	size_t m;
	size_t n;

	T& operator()(size_t i, size_t j) {return data[(i - 1)*m + (j - 1)];};
};


#endif /* MESH_H_ */

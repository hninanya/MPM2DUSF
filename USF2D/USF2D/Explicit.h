#ifndef EXPLICIT_H
#define EXPLICIT_H

#include "Readfile.h"
#include "mesh.h"

const double pi = 3.1415926535897932384626433832795028841971694;


class SolutionPhase:public ClaBoundaryC, public BoundaryConditions, public MeshMPM2D//, public ClaElements
{
public:
 	void solution();
	void dampupdate(int , double* , double* , int, MeshMPM2D*, double);
		
	//void particlesInside(const int&,const int&, MeshMPM2D*, MeshMPM2D*, ClaElements*);
	//bool check_inside_point(int, ClaElements&, MeshMPM2D&);
};


#endif // !EXPLICIT_H

#ifndef MESH_H
#define MESH_H

#include "Readfile.h"
#include "BasisFunc.h"
#include "Newton.h"

class MeshMPM2D;
class BoundaryConditions
{
public:
	void BCVelocitySolid(ClaBoundaryC*, int, MeshMPM2D*);
};

class MeshMPM2D: public ClaParticles, public Shape2D, public ClaElements
{
public:
	///MPM2D - USF

	///--------------------
	double b;										//for solving the nonlinear/linear os system equation
	///--------------------

	int	   particle_id;
	///particle global coordinates
	double particle_x_i;
	double particle_x_j;
	double particle_x_k;

	int	   pid;
	double pcoordx;
	double pcoordy;
	double pcoordz;
	///Natural coordinates
	double part_coord_xi;
	double part_coord_eta;
	double part_coord_zeta;
	///Nodal velocity
	double nvelx;
	double nvely;
	double nvelz;
	///Nodal acceleration
	double naccx;
	double naccy;
	double naccz;
	///Nodal mass
	double nmass;
	///Nodal mommentum
	double nmomentumx;
	double nmomentumy;
	///Damping nodal forces
	double ndforcex;
	double ndforcey;
	///Internal nodal forces
	double niforcex;
	double niforcey;
	///External nodal forces
	double neforcex;
	double neforcey;
	///Total nodal forces
	double ntforcex;
	double ntforcey;
	//External nodal forces (gravity)
	double pfgx;
	double pfgy;


	int dependences[100];							//Stablish the max number of particle per element
	int inside_particles;	

	int belong;

	double ele_area_ref;

	double pvol;
	double pvol0;
	double pmass;
	double pvelx;
	double pvely;
	double pvelz;
	//stresses
	double psigxx;
	double psigxy;
	double psigyx;
	double psigyy;

	double pVSxx;
	double pVSxy;
	double pVSyx;
	double pVSyy;
	//particle velocity gradient 
	double pVG11;
	double pVG12;
	double pVG21;
	double pVG22;
	//particle increment gradient deformation
	double pIGD11;
	double pIGD12;
	double pIGD21;
	double pIGD22;
	//particle gradient deformation
	double pGD11;
	double pGD12;
	double pGD21;
	double pGD22;
	//particle delta deformation d\epsilon
	double de11;
	double de12;
	double de21;
	double de22;
	//particle material
	double young;
	double poisson;
	double pden;
	//particle strain
	double e11;
	double e12;
	double e21;
	double e22;
	//particle momemntum
	double pmox;
	double pmoy;
	//detF
	double detF;

		
	double	CreatModelMPM(const int&, const int&, const int&, int);
	void    CreateParticles(int,int, MeshMPM2D*, ClaParticles*);
	void	VelocityParticles(int, MeshMPM2D*, ClaPartVel*);
	bool	check_inside_point(int, ClaElements&, MeshMPM2D&);
	void	particlesInside(const int&, const int&, MeshMPM2D*, MeshMPM2D*, ClaElements*);
	void	MatrixAVectorb(int, int, MeshMPM2D&);
	void	NaturalCoordinates(int, int);
	void	NodalMassMom(const int&, MeshMPM2D*, MeshMPM2D*, MeshMPM2D*, MeshMPM2D*, ClaElements*);
	void	ConstitutiveModel(const int&, MeshMPM2D*, MeshMPM2D*, double , MeshMPM2D*, ClaElements*);
	void	NodesToParticles2(const int&, MeshMPM2D*, MeshMPM2D*, double, MeshMPM2D*, ClaElements*);
	void	NodalAccVel(MeshMPM2D*, int, int, double);
	void	GridNodalMomentumUpdate(MeshMPM2D*, double, const int&);
	void	NodalExternalforces(const int &, MeshMPM2D*, MeshMPM2D*, MeshMPM2D*, MeshMPM2D*, ClaElements*);
	void	NodalInternalforces(const int &, MeshMPM2D*, MeshMPM2D*, MeshMPM2D*, MeshMPM2D*, ClaElements*);

	void	ResetMesh(MeshMPM2D*, int);
	int		sign(double);

};

extern Solver		*pointer_unknonwn;
extern MeshMPM2D	*particlePointer;
extern MeshMPM2D	*meshParticles;
extern MeshMPM2D	*elemParticles;		  //Is is recommended that thi pointer belongs to Readfile.cpp
extern Shape2D		*GShape;

extern MeshMPM2D	*Vector;
extern double		**MatrixA;
extern double		**Matrix;

extern MeshMPM2D	*nodePointer;


#endif	//!MESH_H    

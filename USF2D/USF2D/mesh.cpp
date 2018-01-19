#include "mesh.h"
#include "Readfile.h"
#include <iomanip>
#include <iostream>
#include "BasisFunc.h"
#include "Newton.h"
#include "material.h"
#include "Explicit.h" 
#include <utility>


using namespace std;

Solver		 *pointer_unknonwn = NULL;
MeshMPM2D	 *particlePointer  = NULL;
MeshMPM2D	 *meshParticles    = NULL;
MeshMPM2D	 *elemParticles    = NULL;	
MeshMPM2D	 *partElements	   = NULL;
MeshMPM2D	 *Vector		   = NULL;
MeshMPM2D	 *nodePointer      = NULL;
Shape2D		 *GShape		  = NULL;


double		 **MatrixA;
double		 **Matrix;
double		error;
double		emax = pow(10, -10);



void MeshMPM2D::VelocityParticles(int nparvel, MeshMPM2D*part, ClaPartVel*preadvel)
{
	int pid = 0;
	for (int i = 0; i < nparvel; i++)
	{
		pid = preadvel[i].id;
		part[pid - 1].pvelx = preadvel[i].vx;
		part[pid - 1].pvely = preadvel[i].vy;
		part[pid - 1].pvelz = preadvel[i].vz;
	}
}

void MeshMPM2D::CreateParticles(int npar,int nele, MeshMPM2D*ppointer, ClaParticles*pread)
{
	if (strcmp(readParticles, "GENERATE") == 0)
	{
		//----------------------------------------------------------------------------------------------------------------
		//CALCULATING ************** 1 **************** MATERIAL POINT WITHIN EACH ELEMENT (AS CENTER)
		//----------------------------------------------------------------------------------------------------------------
		for (int i = 0; i < nele; i++)
		{
			meshParticles[i].particle_x_i = (vectorNodes[elementVector[i].n[0] - 1].x + vectorNodes[elementVector[i].n[1] - 1].x +
				vectorNodes[elementVector[i].n[2] - 1].x + vectorNodes[elementVector[i].n[3] - 1].x) / 4;
			meshParticles[i].particle_x_j = (vectorNodes[elementVector[i].n[0] - 1].y + vectorNodes[elementVector[i].n[1] - 1].y +
				vectorNodes[elementVector[i].n[2] - 1].y + vectorNodes[elementVector[i].n[3] - 1].y) / 4;

			meshParticles[i].particle_id = i + 1;
			ppointer[i].pid     = meshParticles[i].particle_id;
			ppointer[i].pcoordx = meshParticles[i].particle_x_i;
			ppointer[i].pcoordy = meshParticles[i].particle_x_j;
		}
	}
	else if (strcmp(readParticles, "READ") == 0)
	{
		for (int k = 0; k < npar; k++)
		{
			ppointer[k].pid     = pread[k].part_id;
			ppointer[k].pcoordx = pread[k].part_x;
			ppointer[k].pcoordy = pread[k].part_y;
			ppointer[k].pcoordz = pread[k].part_z;
		}
	}
}

double MeshMPM2D::CreatModelMPM(const int& numNodes, const int& numElements, const int& numParticles, int dimension)
{
	particlePointer = (MeshMPM2D*)calloc(numParticles, sizeof(MeshMPM2D));
	if (particlePointer == NULL)	exit(1);

	nodePointer = (MeshMPM2D*)calloc(numNodes, sizeof(MeshMPM2D));
	if (nodePointer == NULL) exit(1);

	meshParticles = (MeshMPM2D*)calloc(numParticles, sizeof(MeshMPM2D));
	if (meshParticles == NULL) exit(1);

	elemParticles = (MeshMPM2D*)calloc(numElements, sizeof(MeshMPM2D));
	if (elemParticles == NULL)exit(1);

	///For Newton-Raphson---------------------------------------------------
	Vector = (MeshMPM2D*)calloc(dimension, sizeof(MeshMPM2D));
	if (Vector == NULL) exit(1);

	pointer_unknonwn = (Solver*)calloc(dimension, sizeof(Solver));
	if (pointer_unknonwn == NULL) exit(1);
	///For Newton-Raphson---------------------------------------------------
	GShape = (Shape2D*)calloc(4, sizeof(Shape2D));		///4 Corresponds to a element of four nodes									  
	if (GShape == NULL) exit(1);

	double total_volume = 0;

	CreateParticles(numParticles,numElements, particlePointer, partVector);
	particlesInside(numElements,numParticles, elemParticles, particlePointer, elementVector);

	for (int i = 0; i < numElements; i++)
	{
		//----------------------------------------------------------------------------------------------------------------
		//CALCULATING THE AREA/VOLUME OF EACH ELEMENT 
		//----------------------------------------------------------------------------------------------------------------

		elemParticles[i].ele_area_ref = (vectorNodes[elementVector[i].n[0] - 1].x*vectorNodes[elementVector[i].n[1] - 1].y +
			vectorNodes[elementVector[i].n[1] - 1].x*vectorNodes[elementVector[i].n[2] - 1].y +
			vectorNodes[elementVector[i].n[2] - 1].x*vectorNodes[elementVector[i].n[3] - 1].y +
			vectorNodes[elementVector[i].n[3] - 1].x*vectorNodes[elementVector[i].n[0] - 1].y -
			vectorNodes[elementVector[i].n[1] - 1].x*vectorNodes[elementVector[i].n[0] - 1].y -
			vectorNodes[elementVector[i].n[2] - 1].x*vectorNodes[elementVector[i].n[1] - 1].y -
			vectorNodes[elementVector[i].n[3] - 1].x*vectorNodes[elementVector[i].n[2] - 1].y -
			vectorNodes[elementVector[i].n[0] - 1].x*vectorNodes[elementVector[i].n[3] - 1].y) / 2;

		total_volume = total_volume + elemParticles[i].ele_area_ref * 1;
		//cout << particlePointer[i].ele_area_ref << endl;

	}

	for (int i = 0; i < numParticles; i++)
	{
		//----------------------------------------------------------------------------------------------------------------
		//CALCULATING THE MASS/VOLUME OF EACH PARTICLE AS FUNCTION OF ELEMENTAL AREA
		//----------------------------------------------------------------------------------------------------------------
		particlePointer[i].pmass = elemParticles[particlePointer[i].belong - 1].ele_area_ref*matVector[elementVector[particlePointer[i].belong -
			1].matid - 1].density / 4;
		//Density
		particlePointer[i].pden = matVector[elementVector[particlePointer[i].belong - 1].matid - 1].density;
		//Volume
		particlePointer[i].pvol = particlePointer[i].pmass / particlePointer[i].pden;
		particlePointer[i].pvol0 = particlePointer[i].pvol;
		//Initial Stress
		particlePointer[i].psigxx = model.iniSxx;
		particlePointer[i].psigxy = model.iniSxy;
		particlePointer[i].psigyy = model.iniSyy;
		particlePointer[i].psigyx = model.iniSxy;


		//----------> (Here^) Divided by one particle per element


		//----------------------------------------------------------------------------------------------------------------
		//velocity of each particle
		//----------------------------------------------------------------------------------------------------------------
		particlePointer[i].pvelx = 0; //Initilize velocity for each particle
		particlePointer[i].pvely = 0;

		//----------------------------------------------------------------------------------------------------------------
		//Particle gradient deformation -> identity beacause  there is nt deformation so F=I
		//----------------------------------------------------------------------------------------------------------------
		particlePointer[i].pGD11 = 1;
		particlePointer[i].pGD12 = 0;
		particlePointer[i].pGD21 = 0;
		particlePointer[i].pGD22 = 1;

		//----------------------------------------------------------------------------------------------------------------
		//Particle material
		//----------------------------------------------------------------------------------------------------------------
		particlePointer[i].young = matVector[elementVector[particlePointer[i].belong - 1].matid - 1].young;
		particlePointer[i].poisson = matVector[elementVector[particlePointer[i].belong - 1].matid - 1].poisson;

	}

	return 1;
}


bool MeshMPM2D::check_inside_point(int polyCorners, ClaElements &eleVector, MeshMPM2D &coordPart)
{
	int   i, j = polyCorners - 1;
	bool  oddNodes = false;

	for (i = 0; i<polyCorners; i++) {
		if ((vectorNodes[eleVector.n[i] - 1].y< coordPart.pcoordy && vectorNodes[eleVector.n[j] - 1].y >= coordPart.pcoordy
			|| vectorNodes[eleVector.n[j] - 1].y< coordPart.pcoordy && vectorNodes[eleVector.n[i] - 1].y >= coordPart.pcoordy)
			&& (vectorNodes[eleVector.n[i] - 1].x <= coordPart.pcoordx || vectorNodes[eleVector.n[j] - 1].x <= coordPart.pcoordx)) {
			if (vectorNodes[eleVector.n[i] - 1].x + (coordPart.pcoordy - vectorNodes[eleVector.n[i] - 1].y) / (vectorNodes[eleVector.n[j] - 1].y -
				vectorNodes[eleVector.n[i] - 1].y)*(vectorNodes[eleVector.n[j] - 1].x - vectorNodes[eleVector.n[i] - 1].x)<coordPart.pcoordx) {
				oddNodes = !oddNodes;
			}
		}
		j = i;
	}
	return oddNodes;
}

void MeshMPM2D::particlesInside(const int& numElements, const  int& numParticles, MeshMPM2D*elepart, MeshMPM2D*pPointer, ClaElements*elPointer)
{
	//----------------------------------------------------------------------------------------------------------------
	//PARTICLES INSIDE ELEMENT ith                (This loop determines how many particles are inside of each element)
	//----------------------------------------------------------------------------------------------------------------	
	for (int iele = 0; iele < numElements; iele++)											  //Loop over all elements
	{
		int counter = 0;
		for (int ipar = 0; ipar < numParticles; ipar++)
		{
			int inside;
			inside = check_inside_point(elPointer[iele].nnodes, elPointer[iele], pPointer[ipar]);
			if (inside == 1)
			{
				counter++;
				elepart[iele].inside_particles = counter;//It determines how many particles are located in each element
				elepart[iele].dependences[counter - 1] = pPointer[ipar].pid;//It determines which part is located in each elem
				pPointer[ipar].belong = elPointer[iele].id; //It determines in which elem each part is located
			}
		}
	}
}

void MeshMPM2D::ResetMesh(MeshMPM2D*Node, int numNodes)
{
	for (int i = 0; i < numNodes; i++)
	{
		Node[i].niforcex		= 0;
		Node[i].niforcey		= 0;

		Node[i].neforcex		= 0;
		Node[i].neforcey		= 0;

		Node[i].ntforcex		= 0;
		Node[i].ntforcey		= 0;
		
		Node[i].nmass			= 0;

		Node[i].nmomentumx		= 0;
		Node[i].nmomentumy		= 0;

		Node[i].naccx			= 0;
		Node[i].naccy			= 0;

		Node[i].nvelx			= 0;
		Node[i].nvely			= 0;
	}	
}

int MeshMPM2D::sign(double num)
{
	if (num < 0) return -1;
	else if (num == 0) return 0;
	else if (num > 0) return 1;

	return 0;
}

void BoundaryConditions::BCVelocitySolid(ClaBoundaryC*BCvel, int numBCvel, MeshMPM2D*Node)
{
		Node[0].nvelx = 0;
		Node[0].nvely = 0;

		Node[0].nmomentumx = 0;
		Node[0].nmomentumy = 0;

		Node[0].ntforcex = 0;
		Node[0].ntforcey = 0;

		Node[0].naccx = 0;
		Node[0].naccy = 0;

		Node[1].nvelx = 0;
		Node[1].nvely = 0;

		Node[1].nmomentumx = 0;
		Node[1].nmomentumy = 0;

		Node[1].ntforcex = 0;
		Node[1].ntforcey = 0;

		Node[1].naccx = 0;
		Node[1].naccy = 0;


	//for (int i = 0; i < numBCvel; i++)
	//{
	//	Node[BCvel[i].idBC - 1].nvelx = BCvel[i].vx;
	//	Node[BCvel[i].idBC - 1].nvely = BCvel[i].vy;

	//	Node[BCvel[i].idBC - 1].nmomentumx = BCvel[i].vx;
	//	Node[BCvel[i].idBC - 1].nmomentumy = BCvel[i].vy;

	//	Node[BCvel[i].idBC - 1].ntforcex = 0;
	//	Node[BCvel[i].idBC - 1].ntforcey = 0;

	//	Node[BCvel[i].idBC - 1].naccx = 0;
	//	Node[BCvel[i].idBC - 1].naccy = 0;

	//}
}



void MeshMPM2D::MatrixAVectorb(int ipar, int dimension, MeshMPM2D&pVector)
{
	double Matrix[2][2] = {};
	Phi_4(pointer_unknonwn[0].C, pointer_unknonwn[1].C);
	GradPhi_4(pointer_unknonwn[0].C, pointer_unknonwn[1].C);

	Vector[0].b = -pVector.pcoordx;
	Vector[1].b = -pVector.pcoordy;


	for (int j = 0; j < 4; j++)												  //4 corresponds to type of element
	{
		Vector[0].b += GShape[j].Sh*vectorNodes[elementVector[pVector.belong - 1].n[j] - 1].x;
		Vector[1].b += GShape[j].Sh*vectorNodes[elementVector[pVector.belong - 1].n[j] - 1].y;
	}

	for (int j = 0; j < 4; j++)
	{
		Matrix[0][0] += GShape[j].GSh_xi*vectorNodes[elementVector[pVector.belong - 1].n[j] - 1].x;
		Matrix[0][1] += GShape[j].GSh_eta*vectorNodes[elementVector[pVector.belong - 1].n[j] - 1].x;
		Matrix[1][0] += GShape[j].GSh_xi*vectorNodes[elementVector[pVector.belong - 1].n[j] - 1].y;
		Matrix[1][1] += GShape[j].GSh_eta*vectorNodes[elementVector[pVector.belong - 1].n[j] - 1].y;
	}
	//MatrixA = Matrix;
	error = normVector<MeshMPM2D*, int>(Vector, dimension);
	LUFactorization<MeshMPM2D*, int, double[2][2]>(Vector, dimension, Matrix);
}

void MeshMPM2D::NaturalCoordinates(int numParticles, int dimension)
{
	for (int ipar = 0; ipar < numParticles; ipar++)
	{
		//Newton-Raphson
		initialGuess(dimension,0);			 //0 = (0,0,0) Initial Guess at center of element [Numada,2013]
		do
		{
			MatrixAVectorb(ipar, dimension, particlePointer[ipar]);
		} while (error > emax);
		meshParticles[ipar].part_coord_xi = pointer_unknonwn[0].C;
		meshParticles[ipar].part_coord_eta = pointer_unknonwn[1].C;

		//std::cout << "[ipar= " << ipar + 1 << "]= " << meshParticles[ipar].part_coord_xi << "  " << meshParticles[ipar].part_coord_eta << std::endl;
	}
}
void MeshMPM2D::NodalMassMom(const int & numElements, MeshMPM2D*elePart, MeshMPM2D*meshPart, MeshMPM2D*Node, MeshMPM2D*Part, ClaElements*EleVector)
{
	for (int i = 0; i < numParticles; i++)
	{
		Node[i].nmomentumx = 0;
		Node[i].nmomentumy = 0;
		Node[i].nmass = 0;


	}

	for (int i = 0; i < numParticles; i++)
	{
		Part[i].pmox = Part[i].pvelx*Part[i].pmass*Part[i].detF;
		Part[i].pmoy = Part[i].pvely*Part[i].pmass*Part[i].detF;
	}

	for (int iele = 0; iele < numElements; iele++)
	{
		int mpts = elePart[iele].inside_particles;
		for (int ipar = 0; ipar < mpts; ipar++)
		{
			int idpart = elePart[iele].dependences[ipar] - 1;

			Phi_4(meshPart[idpart].part_coord_xi, meshPart[idpart].part_coord_eta);
			GradPhi_4(meshPart[idpart].part_coord_xi, meshPart[idpart].part_coord_eta);
			GradPhi_4xy(iele);

			for (int i = 0; i < 4; i++)	///											// 4 define the nodes of each element
			{
				int idnode = EleVector[iele].n[i] - 1;

				//Nodal mass
				Node[idnode].nmass += GShape[i].Sh *Part[idpart].pmass;

				//Nodal momentum
				Node[idnode].nmomentumx += GShape[i].Sh *Part[idpart].pmox;
				Node[idnode].nmomentumy += GShape[i].Sh *Part[idpart].pmoy;

			}
		}
	}

	for (int i = 0; i < numNodes; i++)
	{
		nodePointer[i].nvelx = nodePointer[i].nmomentumx / nodePointer[i].nmass;
		nodePointer[i].nvely = nodePointer[i].nmomentumy / nodePointer[i].nmass;
	}

}

void MeshMPM2D::NodalExternalforces(const int & numElements, MeshMPM2D*elePart, MeshMPM2D*meshPart, MeshMPM2D*Node, MeshMPM2D*Part, ClaElements*EleVector)
{
	for (int i = 0; i < numParticles; i++)
	{
		Node[i].neforcex = 0;
		Node[i].neforcey = 0;
	}

	for (int i = 0; i < numParticles; i++)
	{
		Part[i].pfgx = model.gravx*Part[i].pmass;
		Part[i].pfgy = model.gravy*Part[i].pmass;
	}


	for (int iele = 0; iele < numElements; iele++)
	{
		int mpts = elePart[iele].inside_particles;
		for (int ipar = 0; ipar < mpts; ipar++)
		{
			int idpart = elePart[iele].dependences[ipar] - 1;

			Phi_4(meshPart[idpart].part_coord_xi, meshPart[idpart].part_coord_eta);
			GradPhi_4(meshPart[idpart ].part_coord_xi, meshPart[idpart ].part_coord_eta);
			GradPhi_4xy(iele);

			for (int i = 0; i < 4; i++)	///											// 4 define the nodes of each element
			{
				int idnode = EleVector[iele].n[i] - 1;

				Node[idnode ].neforcex += GShape[i].Sh*Part[idpart].pfgx;
				Node[idnode ].neforcey += GShape[i].Sh*Part[idpart].pfgy;
						
			}
		}
	}

}

void MeshMPM2D::NodalInternalforces(const int & numElements, MeshMPM2D*elePart, MeshMPM2D*meshPart, MeshMPM2D*Node, MeshMPM2D*Part, ClaElements*EleVector)
{
	for (int i = 0; i < numNodes; i++)
	{
		Node[i].niforcex = 0;
		Node[i].niforcey = 0;
	}


	for (int iele = 0; iele < numElements; iele++)
	{
		int mpts = elePart[iele].inside_particles;
		for (int ipar = 0; ipar < mpts; ipar++)
		{
			int idpart = elePart[iele].dependences[ipar] - 1;

			Phi_4(meshPart[idpart].part_coord_xi, meshPart[idpart].part_coord_eta);
			GradPhi_4(meshPart[idpart].part_coord_xi, meshPart[idpart].part_coord_eta);
			GradPhi_4xy(iele);

			for (int i = 0; i < 4; i++)	///											// 4 define the nodes of each element
			{
				int idnode = EleVector[iele].n[i] - 1;

				//Nodal internal force
				Node[idnode].niforcex -= (	Part[idpart].pVSxx*GShape[i].GSh_x +Part[idpart].pVSxy*GShape[i].GSh_y);

				Node[idnode].niforcey -=(					Part[idpart].pVSyx*GShape[i].GSh_x +		Part[idpart].pVSyy*GShape[i].GSh_y);


			}
		}
	}
	///Total nodal forces
	for (int i = 0; i < numNodes; i++)
	{
		Node[i].ntforcex = Node[i].niforcex + Node[i].neforcex;
		Node[i].ntforcey = Node[i].niforcey + Node[i].neforcey;
	}
}

void MeshMPM2D::ConstitutiveModel
(const int&numElements, MeshMPM2D*Part, MeshMPM2D*Node, double dtime, MeshMPM2D*elePart, ClaElements*EleVector)
{
	//for (int i = 0; i < numNodes; i++)
	//{
	//	cout << Node[i].nvelx << " " << Node[i].nvely << endl;
	//}

	for (int iele = 0; iele < numElements; iele++)
	{
		int mpts = elePart[iele].inside_particles;
		for (int ipar = 0; ipar < mpts; ipar++)
		{
			int idpart = elePart[iele].dependences[ipar];
			Phi_4(meshParticles[idpart - 1].part_coord_xi, meshParticles[idpart - 1].part_coord_eta);
			GradPhi_4(meshParticles[idpart - 1].part_coord_xi, meshParticles[idpart - 1].part_coord_eta);
			GradPhi_4xy(iele);

			Part[idpart - 1].pVG11 = 0;
			Part[idpart - 1].pVG12 = 0;
			Part[idpart - 1].pVG21 = 0;
			Part[idpart - 1].pVG22 = 0;
			
			for (int i = 0; i < 4; i++)												// 4 define the nodes of each element
			{
				int idnode = EleVector[iele].n[i];

				//Particle velocity gradient
				Part[idpart - 1].pVG11 += GShape[i].GSh_x*Node[idnode - 1].nvelx;
				Part[idpart - 1].pVG12 += GShape[i].GSh_x*Node[idnode - 1].nvely;
				Part[idpart - 1].pVG21 += GShape[i].GSh_y*Node[idnode - 1].nvelx;
				Part[idpart - 1].pVG22 += GShape[i].GSh_y*Node[idnode - 1].nvely;

				//Part[idpart - 1].pVG11 += GShape[i].GSh_x*Node[idnode - 1].nvelx;
				//Part[idpart - 1].pVG12 += GShape[i].GSh_y*Node[idnode - 1].nvelx;
				//Part[idpart - 1].pVG21 += GShape[i].GSh_x*Node[idnode - 1].nvely;
				//Part[idpart - 1].pVG22 += GShape[i].GSh_y*Node[idnode - 1].nvely;
			}

			///Gradiente deformation F=(I+L.dt)F
			Part[idpart - 1].pIGD11 = (1 + Part[idpart - 1].pVG11*dtime);
			Part[idpart - 1].pIGD12 = (Part[idpart - 1].pVG12*dtime);
			Part[idpart - 1].pIGD21 = (Part[idpart - 1].pVG21*dtime);
			Part[idpart - 1].pIGD22 = (1 + Part[idpart - 1].pVG22*dtime);


			///Gradiente deformation F=(I+L.dt)F
			Part[idpart - 1].pGD11 = (1 + Part[idpart - 1].pVG11*dtime)*Part[idpart - 1].pGD11 + Part[idpart - 1].pVG12*dtime*Part[idpart - 1].pGD21;
			Part[idpart - 1].pGD12 = (1 + Part[idpart - 1].pVG11*dtime)*Part[idpart - 1].pGD12 + Part[idpart - 1].pVG12*dtime*Part[idpart - 1].pGD22;
			Part[idpart - 1].pGD21 = Part[idpart - 1].pVG21*dtime*Part[idpart - 1].pGD11 + (1 + Part[idpart - 1].pVG22*dtime)*Part[idpart - 1].pGD21;
			Part[idpart - 1].pGD22 = Part[idpart - 1].pVG21*dtime*Part[idpart - 1].pGD12 + (1 + Part[idpart - 1].pVG22*dtime)*Part[idpart - 1].pGD22;

			///Updating volume
			Part[idpart - 1].detF = (Part[idpart - 1].pGD11*Part[idpart - 1].pGD22 - Part[idpart - 1].pGD21*Part[idpart - 1].pGD12);
			Part[idpart - 1].pvol = Part[idpart - 1].detF*Part[idpart - 1].pvol0;

			///symmetric part of Gradient velocity for updating deformation de_t+dt=0.5(L+L')dt
			Part[idpart - 1].de11 = 0.5*(Part[idpart - 1].pVG11 + Part[idpart - 1].pVG11)*dtime;
			Part[idpart - 1].de12 = 0.5*(Part[idpart - 1].pVG12 + Part[idpart - 1].pVG21)*dtime;
			Part[idpart - 1].de21 = 0.5*(Part[idpart - 1].pVG21 + Part[idpart - 1].pVG12)*dtime;
			Part[idpart - 1].de22 = 0.5*(Part[idpart - 1].pVG22 + Part[idpart - 1].pVG22)*dtime;
			
			///Constitutive model
			planeStrainLinear(particlePointer[idpart - 1]);

			//Strain at particles
			Part[idpart - 1].e11 += Part[idpart - 1].de11;
			Part[idpart - 1].e12 += Part[idpart - 1].de12;
			Part[idpart - 1].e21 += Part[idpart - 1].de21;
			Part[idpart - 1].e22 += Part[idpart - 1].de22;
		}
	}

	for (int i = 0; i < numParticles; i++)
	{
		Part[i].pVSxx = Part[i].psigxx*Part[i].pvol;
		Part[i].pVSxy = Part[i].psigxy*Part[i].pvol;
		Part[i].pVSyx = Part[i].psigyx*Part[i].pvol;
		Part[i].pVSyy = Part[i].psigyy*Part[i].pvol;

	}

}

void MeshMPM2D::NodesToParticles2(const int&numElements, MeshMPM2D*Part, MeshMPM2D*Node, double dtime, MeshMPM2D*elePart, ClaElements*EleVector)
{
	//VelocityParticles(numPartVel, particlePointer, PartVelocityVector);

	for (int iele = 0; iele < numElements; iele++)
	{
		int mpts = elePart[iele].inside_particles;
		for (int ipar = 0; ipar < mpts; ipar++)
		{
			int idpart = elePart[iele].dependences[ipar];
			Phi_4(meshParticles[idpart - 1].part_coord_xi, meshParticles[idpart - 1].part_coord_eta);
			GradPhi_4(meshParticles[idpart - 1].part_coord_xi, meshParticles[idpart - 1].part_coord_eta);
			GradPhi_4xy(iele);

			for (int i = 0; i < 4; i++)												// 4 define the nodes of each element
			{
				int idnode = EleVector[iele].n[i];

				///Particle velocity
				Part[idpart - 1].pvelx += dtime*GShape[i].Sh*Node[idnode - 1].naccx;
				Part[idpart - 1].pvely += dtime*GShape[i].Sh*Node[idnode - 1].naccy;
				///Particle coordinate
				Part[idpart - 1].pcoordx += dtime*GShape[i].Sh*Node[idnode - 1].nvelx;
				Part[idpart - 1].pcoordy += dtime*GShape[i].Sh*Node[idnode - 1].nvely;

			}
		}
	}
}

void MeshMPM2D::NodalAccVel(MeshMPM2D* Node, int numNodes, int step, double dt)
{
	for (int i = 0; i < numNodes; i++)
	{
		Node[i].naccx = (Node[i].ntforcex - model.damp*Node[i].nmomentumx) / Node[i].nmass;
		Node[i].naccy = (Node[i].ntforcey-model.damp*Node[i].nmomentumy) / Node[i].nmass;
	}

	if (step == 0) 
	{
		for (int j = 0; j < numNodes; j++)
		{
			Node[j].naccx *= .5;
			Node[j].naccy *= .5;
		}
	}


}

void MeshMPM2D::GridNodalMomentumUpdate(MeshMPM2D*Node, double dtime, const int&numNodes)
{///Explicit approach - Fordward-euler Method
	
	for (int i = 0; i < numNodes; i++)
	{
		//Node[i].nmomentumx += Node[i].ntforcex * dtime;
		//Node[i].nmomentumy += Node[i].ntforcey * dtime;

		Node[i].nvelx += Node[i].naccx * dtime;
		Node[i].nvely += Node[i].naccy * dtime;

	}
}

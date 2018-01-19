#define _CRT_SECURE_NO_WARNINGS


#include <iostream>
#include <iomanip>
#include "Explicit.h"
#include "Readfile.h"
#include "mesh.h"
#include "Newton.h"	
#include "BasisFunc.h"
#include "tensor.h"
#include "material.h"
#include "InOuFiles.h"
#include <fstream>


int			dimension = 2;
Solver		nonLinear;
MeshMPM2D	objectMeshMPM2D;
Matrix2		ObjectTest(1,2,3,4);
double		last_ce = 0.0;
double		last_cew = 0.0;

void SolutionPhase::dampupdate(int step, double *last_ce, double *last_cew, int numParticles, MeshMPM2D* Part, double dt)
{

	if (step == 0)
	{
		model.damp =  2.0 * pi * model.dampfreq * model.dampfrac;
		model.dampw = 2.0 * pi * model.dampfreq * model.dampfrac;
		return;
	}

	double curr_ce = 0.0;
	double mass = 0.0;
	double damp_par = 0.0;
	double curr_cew = 0.0;
	double massw = 0.0;
	double damp_parw = 0.0;

	for (int i = 0; i < numParticles; i++)
	{
		mass = Part[i].pden * Part[i].pvol;
		curr_ce += 0.5 * mass * Part[i].pvelx * Part[i].pvelx;
		curr_ce += 0.5 * mass * Part[i].pvely * Part[i].pvely;

		massw = 0;// pch.pporos[i] * pch.mat[i].WDensity * pch.pV[i] * pch.pJ[i];
		curr_cew += 0;// 0.5 * massw * pch.pvw[i].x * pch.pvw[i].x;
		curr_cew += 0;// 0.5 * massw * pch.pvw[i].y * pch.pvw[i].y;
	}

	double d_ce = curr_ce - (*last_ce);
	double d_cew = curr_cew - (*last_cew);

	if (d_ce != 0.0)
	{
		damp_par = (2.0 * model.damp * curr_ce * dt) / fabs(d_ce);

		if (damp_par < model.dampadmul)
			model.damp *= model.dampadmul;
		else
			model.damp /= model.dampadmul;
	}
	if (d_cew != 0.0)
	{
		damp_parw = (2.0 * model.dampw * curr_cew * dt) / fabs(d_cew);

		if (damp_parw < model.dampadmul)
			model.dampw *= model.dampadmul;
		else
			model.dampw /= model.dampadmul;
	}

	(*last_ce) = curr_ce;
	(*last_cew) = curr_cew;
}




void SolutionPhase::solution( )
{
	int step = 0;
	double time = 0;
	
	//void timestep()
	//{
	double Ymod;
	double Poisson;
	double Density;
	double Lambda;
	double G;
	double mindt = 1.0e20;
	double timestep = 0.0;
	double dt;

	//for (int i = 1; i<pch.Npart(); i += 1)
	//{
	Ymod = particlePointer[0].young;
	Poisson = particlePointer[0].poisson;;
	Density = particlePointer[0].pden;
	Lambda = (Ymod * Poisson) / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
	G = Ymod / (2.0 * (1.0 + Poisson));
	timestep = 1.0 / sqrt((Lambda + 2.0 * G) / Density);
	if (timestep < mindt) mindt = timestep;
	//}

	//pch.dt = min(pch.dx,pch.dy) * timestep * mpmModel.dtfrac;
	dt = 0.25* timestep * model.dtfrac;
	//dt = 0.000025551942895096763;// *0 + 0.01;
	//dt = 0.01;

	//}

	//-------------------------------------------------------------------------------------------------------------
	//EXPORT
	//-------------------------------------------------------------------------------------------------------------
	ExternFiles export_1;
	std::cout << "__________________________________________________________________" << endl;

	std::cout << "\nEXPORT RESULTS" << "\n" << endl;

	char *Data[3] = { NULL, NULL,NULL };
	Data[1] = (char *)calloc(BUFSIZ, sizeof(char));

	Data[2] = (char *)calloc(BUFSIZ, sizeof(char));

	cout << "\nEnter a name for the exported file" << setw(4) << " : ";
	scanf("%s", Data[1]);

	cout << "Which extension (.vtk or .txt)?   " << setw(4) << " : ";
	scanf("%s", Data[2]);
	strcat(Data[1], Data[2]);
	ofstream write(Data[1]);
	//export_1.WriteFileHeading(numParticles,  write);
	//export_1.PrintHistory(numParticles, write, step,time);

	cout << "dt " << dt << endl;

		while (time <= totaltime)
	{
		//--------------------------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------------------------
		particlesInside(numElements, numParticles, elemParticles, particlePointer, elementVector);
		//cout << endl;
		//cout << "elemParticles[iele].inside_particles" << " " << "elemParticles[iele].dependences[counter - 1]" << endl;
		//for (int iele = 0; iele < numElements; iele++)
		//{
		//	cout << elemParticles[iele].inside_particles << " " << elemParticles[iele].dependences[0] << " " << elemParticles[iele].dependences[1]
		//		<< " " << elemParticles[iele].dependences[2] << " " << elemParticles[iele].dependences[3] << " " <<
		//		elementVector[iele].n[0] << " " <<
		//		elementVector[iele].n[1] << " " <<
		//		elementVector[iele].n[2] << " " <<
		//		elementVector[iele].n[3] << endl;
		//}
		//--------------------------------------------------------------------------------------------------------------------
		NaturalCoordinates(numParticles, dimension);
		//cout << endl;
		//cout << "particlePointer[i].part_coord_xi" << " " << "particlePointer[i].part_coord_eta" << endl;
		//for (int i = 0; i < numParticles; i++)
		//{
		//	cout << meshParticles[i].part_coord_xi << " " << meshParticles[i].part_coord_eta << endl;
		//}
		//--------------------------------------------------------------------------------------------------------------------
		//ResetMesh(nodePointer,numNodes);///OK
		//cout << endl;
		//cout << "Reset" << endl;
		//for (int i = 0; i < numNodes; i++)
		//{
		//	cout << nodePointer[i].niforcex << " " << nodePointer[i].niforcey <<
		//		" " << nodePointer[i].neforcex << " " << nodePointer[i].neforcey <<
		//		" " << nodePointer[i].ntforcex << " " << nodePointer[i].ntforcey <<
		//		" "<< nodePointer[i].nmass << " " << nodePointer[i].nmomentumx << " " << nodePointer[i].nmomentumy<<endl;
		//}
		//cout << model.damp << endl;
		//--------------------------------------------------------------------------------------------------------------------
		dampupdate(step, &last_ce, &last_cew, numParticles, particlePointer, dt);
		NodalExternalforces(numElements, elemParticles, meshParticles, nodePointer, particlePointer, elementVector);
		ConstitutiveModel(numElements, particlePointer, nodePointer, dt, elemParticles, elementVector);
		NodalMassMom(numElements, elemParticles, meshParticles, nodePointer, particlePointer, elementVector);
		//cout << endl;
		//cout << "particlePointer[i].pvelx" << " " << "particlePointer[i].pvely" << endl;
		//for (int i = 0; i < numParticles; i++)
		//{
		//	cout << particlePointer[i].pVSxx << " " << particlePointer[i].pVSyy << endl;
		//}
		//cout << endl;
		//cout << "particlePointer[i].pvelx" << " " << "particlePointer[i].pvely" << endl;
		//for (int i = 0; i < numParticles; i++)
		//{
		//	cout << particlePointer[i].detF <<" " << particlePointer[i].detF << endl;
		//}
		
		NodalInternalforces(numElements, elemParticles, meshParticles, nodePointer, particlePointer, elementVector);
		//cout << endl;
		//cout << "nodePointer[i].naccx" << " " << "nodePointer[i].naccy" << endl;
		//for (int i = 0; i < numNodes; i++)
		//{
		//	cout << nodePointer[i].niforcex << " " << nodePointer[i].niforcex << endl;
		//}
		
		NodalAccVel(nodePointer, numNodes, step, dt);
		BCVelocitySolid(vectorBCvelocity, numBCvel, nodePointer);


		GridNodalMomentumUpdate(nodePointer, dt, numNodes);


		NodesToParticles2(numElements, particlePointer, nodePointer, dt, elemParticles, elementVector);

	


		//cout << endl;
		//cout << "nodePointer[i].nmass"<< " " << "nodePointer[i].nmomentumx" << " " << "nodePointer[i].nmomentumy" << endl;
		//for (int i = 0; i < numNodes; i++)
		//{
		//	cout<< nodePointer[i].nmass << " " << nodePointer[i].nmomentumx << " " << nodePointer[i].nmomentumy << endl;
		//}
		//--------------------------------------------------------------------------------------------------------------------
	

	
		//--------------------------------------------------------------------------------------------------------------------
		//cout << endl;
		//cout << "nodePointer[i].nmass" << " " << "nodePointer[i].nmomentumx" << " " << "nodePointer[i].nmomentumy" << endl;
		//for (int i = 0; i < numNodes; i++)
		//{
		//	cout <<nodePointer[i].nmomentumx << " " << nodePointer[i].nmomentumy << " " << nodePointer[i].ntforcex << " " << nodePointer[i].ntforcey << endl;
		//}
	
		//--------------------------------------------------------------------------------------------------------------------




		////--------------------------------------------------------------------------------------------------------------------
		//cout << endl;
		//cout << "nodePointer[i].nvelx" << " " << "nodePointer[i].nvely" << endl;
		//for (int i = 0; i < numNodes; i++)
		//{
		//	cout << nodePointer[i].nvelx << " " << nodePointer[i].nvely << endl;
		//}

		//--------------------------------------------------------------------------------------------------------------------
		//cout << endl;
		//cout << "nodePointer[i].nmass" << " " << "nodePointer[i].nvelx" << " " << "nodePointer[i].nvely" << endl;
		//for (int i = 0; i < numNodes; i++)
		//{
		//	cout << nodePointer[i].nmass << " " << nodePointer[i].nvelx << " " << nodePointer[i].nvely << endl;
		//}

		//--------------------------------------------------------------------------------------------------------------------
		//cout << endl;
		//for (int i = 0; i < numParticles; i++)
		//{
		//	cout << particlePointer[i].pcoordx << " " << particlePointer[i].pcoordy << endl;
		//}

		//cout << "particlePointer[i].pe11" << " " << "particlePointer[i].pe12" << " " << "particlePointer[i].pe21" << " "
		//	<< "particlePointer[i].pe22" << endl;
		//for (int i = 0; i < numParticles; i++)
		//{
		//	cout << particlePointer[i].psigxx << " " << particlePointer[i].psigxy << " " << particlePointer[i].psigyx << " " 
		//		<< particlePointer[i].psigyy << endl;
		//}

		//export_1.PrintHistory(numParticles, write, step,time);
		step++;
		time +=  dt;
		
	}
	//
	//export_1.PrintHistory(numParticles, write, step,time);
	//export_1.OutputField(numParticles, write, step, time);
	export_1.vtkfile(numParticles, write, step, time, particlePointer);

	
	

	cout << "\n\nThe file haz been created           : " << Data[1] << endl;

	std::cout << "\n__________________________________________________________________" << endl;
}

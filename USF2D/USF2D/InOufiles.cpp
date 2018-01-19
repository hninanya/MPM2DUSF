#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include "InOuFiles.h"
#include "mesh.h"
#include "Newton.h"	
#include "Readfile.h"
 


using namespace std;

static char inputfile[BUFSIZ];
static FILE *_ioData = NULL;



int ExternFiles::OpenFile(int argc, char*argv[])
{
	if (argc>1)
	{
		strcpy(inputfile, argv[1]);
	}
	else
	{
		return 1;
	}

	if ((_ioData=fopen(inputfile,"r"))==NULL)
	{
		return 0;
	}	 

	return 1;
}

int ExternFiles::ReadFile()
{
	pointerFile(_ioData);
	subtitles();
	
	return 1;
}

void ExternFiles::vtkfile(int numParticles, std::ostream&os, int step, double time, MeshMPM2D*part)
{
	os << "# vtk DataFile Version 1.0" << endl;
	os << "3D triangulation data" << endl;
	os << "ASCII" << endl;
	os << " " << endl;
	os << "DATASET POLYDATA" << endl;
	os << "POINTS" << " " << numParticles << " " << "float" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].pcoordx << " " << part[i].pcoordy << " " << part[i].pcoordz << endl;
	}
	os << " " << endl;
	os << "POINT_DATA" << " " << numParticles << endl;
	os << "SCALARS" << " " << "uxx" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].pcoordx - partVector[i].part_x << endl;
	}
	os << " " << endl;
	os << "SCALARS" << " " << "uyy" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].pcoordy - partVector[i].part_y << endl;
	}
	os << " " << endl;
	os << "SCALARS" << " " << "psigxx" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].psigxx << endl;
	}
	os << " " << endl;
	os << "SCALARS" << " " << "psigxy" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].psigxy << endl;
	}
	os << " " << endl;
	os << "SCALARS" << " " << "psigyy" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].psigyy << endl;
	}
	os << " " << endl;
	os << "SCALARS" << " " << "e11" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].e11 << endl;
	}
	os << " " << endl;
	os << "SCALARS" << " " << "e12" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].e12 << endl;
	}
	os << " " << endl;
	os << "SCALARS" << " " << "e22" << " " << "float" << " " << 1 << endl;
	os << "LOOKUP_TABLE" << " " << "default" << endl;
	for (int i = 0; i < numParticles; i++)
	{
		os << part[i].e22 << endl;
	}
}
void ExternFiles::WriteFileHeading(int numPartciles, std::ostream&os)
{



	////---------------------------------------------------------------------------------------------------------------------------------
	//os << "\n" << endl;
	//os << "****************************************************************************************************" << endl;
	//os << "MESH" << endl;
	//os << "****************************************************************************************************" << endl;
	//os << "#element" << setw(10) << "node(1)" << setw(10) << "node(2)" << setw(10) << "node(3)" << setw(10) << "node(4)" << endl;
	//for (int i = 0; i < numElements; i++)
	//{
	//	os << i + 1;
	//	for (int j = 0; j < 4; j++)
	//	{
	//		os << setw(11) << elementVector[i].n[j];
	//	}
	//	os << endl;
	//}
}

void ExternFiles::OutputField(int numPartciles, std::ostream&os, int step, double time)
{
	//if (step == 0)
	//{
	//	os << "\n" << endl;
	//	os << "****************************************************************************************************" << endl;
	//	os << "PARTICLES" << "Step - " << step+1 << endl;
	//	os << "****************************************************************************************************" << endl;
	//	os << " id" << setw(11) << "coord_x" << setw(11) << "coord_y"  << setw(12) << "e11" << setw(21)
	//		<< "e22" << setw(21) << "e12" << setw(25) << "sigma11" << setw(21) << "sigma22" << setw(21) << "sigma12" << endl;
	//}
	//for (int i = 0; i < numPartciles; i++)
	//{
	//	//os << setw(3) << i + 1 << setw(11) << fixed << setprecision(0) << elemParticles[i].inside_particles << endl;
	//	os << setw(3) << particlePointer[i].pid << setw(11) << fixed << setprecision(5) << particlePointer[i].pcoordx << setw(11) <<
	//		fixed << setprecision(5) << particlePointer[i].pcoordy << setw(21) << fixed <<
	//		setprecision(10) << particlePointer[i].e11 << setw(21) << fixed << setprecision(10) << particlePointer[i].e22 <<
	//		setw(21) << fixed << setprecision(10) << particlePointer[i].e12 << setw(21) << fixed << setprecision(10) <<
	//		particlePointer[i].psigxx << setw(21) << fixed << setprecision(10) << particlePointer[i].psigyy << setw(21) <<
	//		fixed << setprecision(10) << particlePointer[i].psigxy << endl;

	//}

	os << time << " " << particlePointer[10].pcoordy << " " << particlePointer[10].e22 << " " << particlePointer[10].psigyy << endl;

	if (step == 0)
	{
		os << "\n" << endl;
		os << "****************************************************************************************************" << endl;
		os << "PARTICLES" << "Step - " << step + 1 << endl;
		os << "****************************************************************************************************" << endl;
		os << setw(11) << "coord_x" << setw(11) << "coord_y" << setw(12) << "e11" << setw(21)
			<< "e22" << setw(21) << "e12" << setw(25) << "sigma11" << setw(21) << "sigma22" << setw(21) << "sigma12" << endl;
	}
	os << setw(11) << "coord_x" << setw(11) << "coord_y" << setw(12) << "e11" << setw(21)
		<< "e22" << setw(21) << "e12" << setw(25) << "sigma11" << setw(21) << "sigma22" << setw(21) << "sigma12" << endl;
	for (int i = 0; i < numPartciles; i++)
	{
		os << setw(11) << fixed << setprecision(5) << particlePointer[i].pcoordx << setw(11) <<
			fixed << setprecision(5) << particlePointer[i].pcoordy << setw(21) << fixed <<
			setprecision(10) << particlePointer[i].e11 << setw(21) << fixed << setprecision(10) << particlePointer[i].e22 <<
			setw(21) << fixed << setprecision(10) << particlePointer[i].e12 << setw(21) << fixed << setprecision(10) <<
			particlePointer[i].psigxx << setw(21) << fixed << setprecision(10) << particlePointer[i].pVSyy << setw(21) <<
			fixed << setprecision(10) << particlePointer[i].psigxy << endl;

	}

	//for (int i = 0; i < numElements; i++)
	//{
	//	os << setw(3) << elementVector[i].id << setw(12) << fixed << setprecision(0) << elemParticles[i].inside_particles << setw(11);
	//	for (int ipar = 0; ipar < elemParticles[i].inside_particles; ipar++)
	//	{
	//		os << elemParticles[i].dependences[ipar] << setw(10);
	//	}

	//	os << endl;
	//}


	//---------------------------------------------------------------------------------------------------------------------------------

	//if (step == 0)
	//{
	//	os << "\n" << endl;
	//	os << "****************************************************************************************************" << endl;
	//	os << "ELEMENTS - PARTICLES" << endl;
	//	os << "****************************************************************************************************" << endl;
	//	os << " ElementID" << setw(12) << "dependences" << setw(11) << "ParticleID" << endl;
	//}
	//for (int i = 0; i < numElements; i++)
	//{
	//	os << setw(3) << elementVector[i].id << setw(12) << fixed << setprecision(0) << elemParticles[i].inside_particles << setw(11);
	//	for (int ipar = 0; ipar < elemParticles[i].inside_particles; ipar++)
	//	{
	//		os << elemParticles[i].dependences[ipar] << setw(10);
	//	}

	//	os << endl;
	//}



	//---------------------------------------------------------------------------------------------------------------------------------
	//if (step == 0)
	//{
	//	os << "\n" << endl;
	//	os << "****************************************************************************************************" << endl;
	//	os << "PARTICLES - ELEMENTS" << endl;
	//	os << "****************************************************************************************************" << endl;
	//	os << "ParticleID" << setw(12) << "ElememtID" << endl;
	//}
	//for (int i = 0; i < numElements; i++)
	//{
	//	os << setw(3) << meshParticles[i].particle_id << setw(12) << fixed << setprecision(0) << meshParticles[i].belong << endl;
	//}
}
void ExternFiles::PrintHistory(int numPartciles, std::ostream&os, int step, double time)
{


		os <<time<<" "<< particlePointer[10].pcoordy <<" "<<particlePointer[10].e22<<" "<<particlePointer[10].psigyy<< endl;


}


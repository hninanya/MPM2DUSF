#define _CRT_SECURE_NO_WARNINGS

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>


#include "BasisFunc.h"
#include "mesh.h"
#include "Newton.h"
#include "InOuFiles.h"
#include "Explicit.h"
#include <omp.h>

#define DEBUG

//MPM2DUSF - Testando dede XPS01131487641986

#include "Readfile.h"


using namespace std;
using namespace std::chrono;

ExternFiles openfile;
ExternFiles readfile;



int main(int argc, char *argv[]){

#ifdef _OPENMP
	int nProcessors = omp_get_max_threads();
	printf("Num_threads:%3d\n", nProcessors);
	omp_set_num_threads(nProcessors);
#endif



	std::cout << "Institution	: Pontifical Catholic University of Rio de Janeiro" << endl;	
	std::cout << "Date		:";
	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	cout <<" "<< (now->tm_year + 1900) << '-'
		<< (now->tm_mon + 1) << '-'
		<< now->tm_mday
		<< endl;
	std::cout << "__________________________________________________________________" << endl;
	std::cout << endl;

	

	//INPUT FILE
	char *Files_[2] = { NULL, NULL };
	Files_[1] = (char*)calloc(BUFSIZ, sizeof(char));

	if (argc < 2)
	{
		printf("Enter the input file with its extension : ");
		scanf("%s", Files_[1]);
	}
	else
	{
		strcpy(Files_[1], argv[1]);
	}

	if (!openfile.OpenFile(2,Files_))
	{
		return 0;
	}

	//READ FILE

	readfile.ReadFile();

//#ifdef DEBUG
	//INITIALIZE BACKGROUND GRID NODAL MASS AND MOMENTUM
	high_resolution_clock::time_point t1 = high_resolution_clock::now(); //<-----Run Time
	MeshMPM2D inibackground;
	SolutionPhase explicit_solution;

	inibackground.CreatModelMPM(numNodes, numElements, numParticles, 2); ///2 = dimension
	explicit_solution.solution();

	high_resolution_clock::time_point t2 = high_resolution_clock::now(); //<-----Run Time
	auto duration = duration_cast<microseconds>(t2 - t1).count();		 //<-----Run Time

	cout << "\nExecution Time = " << duration << " microseconds" << endl;//<-----Run Time
//#endif // !DEBUG





		   



	//-------------------------------------------------------------------------------------------------------------
	//EXIT
	//-------------------------------------------------------------------------------------------------------------
	cout << "\n\nDone!!! \n\n";
	printf("Start data : %s  -  Hora : %s", __DATE__, __TIME__);
	cout << "\n\n";

	system("pause");

	return 0;
}									  


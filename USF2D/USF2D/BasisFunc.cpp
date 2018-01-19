#include <iostream>
#include <iomanip>	
#include "BasisFunc.h"
#include "Readfile.h"
#include "mesh.h"
using namespace std;

					 

void Shape2D::Phi_4( double xi,  double eta)
{
	GShape[0].Sh = (xi - 1)*(eta - 1) / 4;
	GShape[1].Sh = -(xi + 1)*(eta - 1) / 4;
	GShape[2].Sh = (xi + 1)*(eta + 1) / 4;
	GShape[3].Sh = -(xi - 1)*(eta + 1) / 4;



	//0.06250
	//	0.18750
	//	0.56250
	//	0.18750


}

void Shape2D::GradPhi_4( double xi,  double eta)		 //Arrange according to element's numeration
{
	GShape[0].GSh_xi = (eta - 1) / 4;	
	GShape[1].GSh_xi = -(eta - 1) / 4;	
	GShape[2].GSh_xi = (eta + 1) / 4;	
	GShape[3].GSh_xi = -(eta + 1) / 4;	

	//-0.12500
	//	0.12500
	//	0.37500
	//	- 0.37500


	GShape[0].GSh_eta = (xi - 1) / 4;	
	GShape[1].GSh_eta = -(xi + 1) / 4;	
	GShape[2].GSh_eta = (xi + 1) / 4;	
	GShape[3].GSh_eta = -(xi - 1) / 4;	

	//-0.12500
	//	- 0.37500
	//	0.37500
	//	0.12500

}

void Shape2D::GradPhi_4xy(int iele)
{
	double J[2][2] = { {0} };
	double invJ[2][2] = { { 0 } };
	
	for (int i = 0; i < 4; i++)
	{
		J[0][0] += GShape[i].GSh_xi*vectorNodes[elementVector[iele].n[i] - 1].x;
		J[0][1] += GShape[i].GSh_xi*vectorNodes[elementVector[iele].n[i] - 1].y;

		J[1][0] += GShape[i].GSh_eta*vectorNodes[elementVector[iele].n[i] - 1].x;
		J[1][1] += GShape[i].GSh_eta*vectorNodes[elementVector[iele].n[i] - 1].y;
	}
	invMatrix2x2(*J, *invJ);
	
	GShape[0].GSh_x = invJ[0][0] * GShape[0].GSh_xi + invJ[0][1] * GShape[0].GSh_eta;
	GShape[1].GSh_x = invJ[0][0] * GShape[1].GSh_xi + invJ[0][1] * GShape[1].GSh_eta;
	GShape[2].GSh_x = invJ[0][0] * GShape[2].GSh_xi + invJ[0][1] * GShape[2].GSh_eta;
	GShape[3].GSh_x = invJ[0][0] * GShape[3].GSh_xi + invJ[0][1] * GShape[3].GSh_eta;

	GShape[0].GSh_y = invJ[1][0] * GShape[0].GSh_xi + invJ[1][1] * GShape[0].GSh_eta;
	GShape[1].GSh_y = invJ[1][0] * GShape[1].GSh_xi + invJ[1][1] * GShape[1].GSh_eta;
	GShape[2].GSh_y = invJ[1][0] * GShape[2].GSh_xi + invJ[1][1] * GShape[2].GSh_eta;
	GShape[3].GSh_y = invJ[1][0] * GShape[3].GSh_xi + invJ[1][1] * GShape[3].GSh_eta;



	//  gradx1 - 2.5
	//	grady1 - 2.5


	//	gradx2	2.5
	//	grady2 - 7.5


	//	gradx3	7.5
	//	grady3	7.5


	//	gradx4 - 7.5
	//	grady4	2.5


}








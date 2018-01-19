#include "material.h"
#include "mesh.h"

#ifndef PI
#define PI 3.141592654
#endif

void planeStrainLinear(MeshMPM2D&part)
{
	double Ymod = part.young;
	double nu	= part.poisson;

	double cm[3][3];

	cm[0][0] = cm[1][1] = (Ymod * (1.0 - nu)) / ((1.0 + nu) * (1.0 - (2.0 * nu)));
	cm[0][1] = cm[1][0] = (Ymod * nu) / ((1.0 + nu) * (1.0 - (2.0 * nu)));
	cm[2][2] = Ymod / (1.0 + nu);
	cm[0][2] = cm[1][2] = cm[2][0] = cm[2][1] = 0.0;
	//the Cauchy stress update
	part.psigxx += (cm[0][0] * (part.de11) + cm[0][1] * (part.de22) + cm[0][2] * (part.de12));
	part.psigxy += (cm[2][0] * (part.de11) + cm[2][1] * (part.de22) + cm[2][2] * (part.de12));
	part.psigyx += (cm[2][0] * (part.de11) + cm[2][1] * (part.de22) + cm[2][2] * (part.de21));
	part.psigyy += (cm[1][0] * (part.de11) + cm[1][1] * (part.de22) + cm[1][2] * (part.de21));
}
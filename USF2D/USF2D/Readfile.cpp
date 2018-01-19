#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <cstdio>
#include <iomanip>
#include "Readfile.h"

using namespace std;

ClaModel					 model;	   
ClaCoordNodes				*vectorNodes		= NULL;
ClaMatVector				*matVector			= NULL;
ClaElements					*elementVector		= NULL;
ClaParticles				*partVector			= NULL;	
ClaBoundaryC				*vectorBCvelocity	= NULL;
ClaPartVel					*PartVelocityVector = NULL;



char						 readParticles[BUFSIZ];

int		                     numPartVel			= 0;
int							 numNodes			= 0;
int							 numBCvel			= 0;
int							 numMat				= 0;
int							 numElements		= 0;
int							 numParticles		= 0;
double						 dtime				= 0;
double						 totaltime			= 0;

static size_t				 valuereader		= 0;
static FILE					*_readData			= NULL;									 


int Readsubtitles::subtitles()
{
	char subtitle[100];
	while (1)
	{
		if (findsubtitle(subtitle) == 0){ return 0; }
		else if (strcmp(subtitle,"HEADING")==0)
		{
		ReadHeading();
		}
		else if (strcmp(subtitle, "NCELLS") == 0)
		{
		ReadNCells();
		}
		else if (strcmp(subtitle, "LIMITS") == 0)
		{
		ReadLimits();
		}
		else if (strcmp(subtitle, "TOTAL.TIME") == 0)
		{
		ReadTotalTime();
		}
		else if (strcmp(subtitle, "PRINT.STEP") == 0)
		{
		//ReadPrintStep();
		}
		else if (strcmp(subtitle, "GRAVITY") == 0)
		{
		ReadGravity();
		}
		else if (strcmp(subtitle, "INITIAL.STRESS") == 0)
		{
		ReadHasInitStress();
		}
		else if (strcmp(subtitle, "POINT.RESULT") == 0)
		{
		//ReadPointResults();
		}
		else if (strcmp(subtitle, "DAMP") == 0)
		{
		ReadDamp();
		}
		else if (strcmp(subtitle, "TOTAL.NODES") == 0)
		{
		ReadTotalNodes();
		}
		else if (strcmp(subtitle, "TOTAL.PARTICLES") == 0)
		{
		ReadTotalParticles();
		}
		else if (strcmp(subtitle, "PARTICLES") == 0)
		{
		ReadParticlesCoords();
		}
		else if (strcmp(subtitle, "NODES") == 0)
		{
		ReadNodeCoord();
		}
		else if (strcmp(subtitle, "MATERIAL") == 0)
		{
		ReadNumMaterial();
		}
		else if (strcmp(subtitle, "MATERIAL.LINEAR") == 0)
		{
		ReadMatLinear();
		}
		else if (strcmp(subtitle, "MATERIAL.FLOW") == 0)
		{
		ReadMatFlow();
		}
		else if (strcmp(subtitle, "ELEMENTS") == 0)
		{
		ReadNumOfElements();
		}
		else if (strcmp(subtitle, "ELEMENT.Q4") == 0)
		{
		ReadElementQ4();
		}
		else if (strcmp(subtitle, "BC.VELOCITY") == 0)
		{
		ReadBCVelocity();
		}
		else if (strcmp(subtitle, "POINT.VELOCITY") == 0)
		{
		ReadPartVelocity();
		}
		else if (strcmp(subtitle, "ELEMENTS") == 0)
		{
			cout << "ELEMENTS Subtitle has been found" << endl;
		}
		else if (strcmp(subtitle, "TIME") == 0)
		{
			cout << "TIME Subtitle has been found" << endl;
		}
		else if (strcmp(subtitle, "MATERIAL") == 0)
		{
			cout << "MATERIAL Subtitle has been found" << endl;
		}
		else if (strcmp(subtitle, "GRAVITY") == 0)
		{
			cout << "GRAVITY Subtitle has been found" << endl;
		}
		//else if (strcmp(subtitle, "NAME") == 0){ cout << "NAME ..." << endl; }
		else if (strcmp(subtitle, "END") == 0) { break; }

	}

	return 1;
}

void Readsubtitles::pointerFile(FILE *file)
{
	_readData = file;
}

int Readsubtitles::findsubtitle(char *subtitle)
{
	int c;
	while ((c = fgetc(_readData)) != '%')
		if (c == EOF)
			return 0;
	return (fscanf(_readData, "%s", subtitle));
}

void Readsubtitles::ReadHeading()
{
	char heading[BUFSIZ];
	fscanf(_readData, "%s", heading);

	printf("\nHEADING : %s\n", heading);
}

void Readsubtitles::ReadNCells()
{
	int ncells_x, ncells_y;

	valuereader = fscanf(_readData, "%d%d", &ncells_x, &ncells_y);
	model.ncx = ncells_x;
	model.ncy = ncells_y;

}

void Readsubtitles::ReadBCVelocity()
{
	int id, totalBCvelocity;
	double vx, vy, vz;

	valuereader = fscanf(_readData, "%d", &totalBCvelocity);
	/*if (totalBCvelocity == 0)
	{
		return;
	}*/
	numBCvel = totalBCvelocity;
	vectorBCvelocity = (ClaBoundaryC*)calloc(totalBCvelocity, sizeof(ClaBoundaryC));

	for (int i = 0; i < totalBCvelocity; i++)
	{
		valuereader = fscanf(_readData, "%d%lf%lf%lf", &id, &vx, &vy, &vz);
		
		vectorBCvelocity[i].idBC	= id;
		vectorBCvelocity[i].vx      = vx;
		vectorBCvelocity[i].vy      = vy;
		vectorBCvelocity[i].vz		= vz;
	}
}
void Readsubtitles::ReadPartVelocity(void)
{
	int npveloc = 0;
	int id;
	double vx, vy, vz;

	valuereader = fscanf(_readData, "%d", &npveloc);

	/* Set mpmNumPointVelocity
	*/
	numPartVel = npveloc;

	/* allocate mpmNumPointVelocity */
	PartVelocityVector = (ClaPartVel*)calloc(npveloc, sizeof(ClaPartVel));

	for (int i = 0; i < npveloc; i++)
	{
		valuereader = fscanf(_readData, "%d%lf%lf%lf", &id, &vx, &vy, &vz);
		PartVelocityVector[i].id = id;
		PartVelocityVector[i].vx = vx;
		PartVelocityVector[i].vy = vy;
		PartVelocityVector[i].vz = vz;
	}

} /* End of _mpmReadPointVelocity */



void Readsubtitles::ReadLimits(void)
{
	double minx, miny, minz, maxx, maxy, maxz;

	valuereader = fscanf(_readData, "%lf%lf%lf%lf%lf%lf", &minx, &miny, &minz, &maxx, &maxy, &maxz);
	model.minx = minx;
	model.miny = miny;
	model.minz = minz;
	model.maxx = maxx;
	model.maxy = maxy;
	model.maxz = maxz;

}

void Readsubtitles::ReadTotalTime(void)
{
	//Total time of analysis
	double ttime;
	//double dt;

	valuereader = fscanf(_readData, "%lf", &ttime);// , &dt);
	totaltime = ttime;
	//dtime = dt;
}

//Read total particles
void Readsubtitles::ReadTotalParticles(void)
{
	int nparticles;
	char particles[BUFSIZ];
	valuereader = fscanf(_readData, "%s%d",particles, &nparticles);
	strcpy(readParticles, particles);
	numParticles = nparticles;
	partVector = (ClaParticles*)calloc(nparticles, sizeof(ClaParticles));
}

void Readsubtitles::ReadParticlesCoords(void)
{
	int totalparticles, pid;
	double pcoordx, pcoordy, pcoordz;
	valuereader = fscanf(_readData, "%d", &totalparticles);


	for (int i = 0; i < totalparticles; i++)
	{
		valuereader = fscanf(_readData, "%d%lf%lf%lf", &pid, &pcoordx, &pcoordy, &pcoordz);
		partVector[i].part_id = pid;
		partVector[i].part_x = pcoordx;
		partVector[i].part_y = pcoordy;
		partVector[i].part_z = pcoordz;
	}
}

void Readsubtitles::ReadGravity(void)
{
	double hasgravity_x;
	double hasgravity_y;


	valuereader = fscanf(_readData, "%lf%lf", &hasgravity_x, &hasgravity_y);
	model.gravx = hasgravity_x;
	model.gravy = hasgravity_y;

}

void Readsubtitles::ReadDamp(void)
{
	double dfreq, dfrac, dadmul, dtfrac;

	valuereader = fscanf(_readData, "%lf%lf%lf%lf", &dfreq, &dfrac, &dadmul, &dtfrac);
	model.dampfreq = dfreq;
	model.dampfrac = dfrac;
	model.dampadmul = dadmul;
	model.dtfrac = dtfrac;
}

void Readsubtitles::ReadTotalNodes(void)
{
	int nnode;

	valuereader = fscanf(_readData, "%d", &nnode);
	numNodes = nnode; // Set number of nodes
	if (nnode == 0)
	{
		return;
	}
	vectorNodes = (ClaCoordNodes*)calloc(nnode, sizeof(ClaCoordNodes));
}

void Readsubtitles::ReadNodeCoord(void)
{
	int totalnodes, id;
	double x, y, z;

	valuereader = fscanf(_readData, "%d", &totalnodes);
	for (int i = 0; i < totalnodes; i++)
	{
		valuereader = fscanf(_readData, "%d%lf%lf%lf", &id, &x, &y, &z);

		vectorNodes[i].id = id;
		vectorNodes[i].x = x;
		vectorNodes[i].y = y;
		vectorNodes[i].z = z;
	}
	//printf("File Grid Coordinate has been uploaded\n");
	//for (int i = 0; i < totalnodes; i++)
	//{
	//	cout <<setprecision(2)<<vectorNodes[i].id << " " << vectorNodes[i].x << " " << vectorNodes[i].y
	//		<< " " << vectorNodes[i].z << endl;
	//}
}

void Readsubtitles::ReadNumMaterial(void)
{
	int nmat;

	valuereader = fscanf(_readData, "%d", &nmat);
	numMat = nmat; //Set Number of materials
	matVector = (ClaMatVector*)calloc(nmat, sizeof(ClaMatVector));
}

void Readsubtitles::ReadMatLinear(void)
{
	int nmat, id;
	double e, nu, rho;

	valuereader = fscanf(_readData, "%d", &nmat);//Here, I am definning properties for each material
	for (int i = 0; i < nmat; i++)
	{
		valuereader = fscanf(_readData, "%d%lf%lf%lf", &id, &e, &nu, &rho);
		matVector[id - 1].id = id;
		matVector[id - 1].type = 5;
		matVector[id - 1].young = e;
		matVector[id - 1].poisson = nu;
		matVector[id - 1].density = rho;
	}
}
void Readsubtitles::ReadMatFlow(void)
{
	int nmat, id;
	double kxx, kyy, kw, poros, rho;

	valuereader = fscanf(_readData, "%d", &nmat);
	for (int i = 0; i < nmat; i++)
	{
		valuereader = fscanf(_readData, "%d%lf%lf%lf%lf%lf", &id, &kxx, &kyy, &kw, &poros, &rho);
		matVector[id - 1].id = id;
		matVector[id - 1].kxx = kxx;
		matVector[id - 1].kyy = kyy;
		matVector[id - 1].kw = kw;
		matVector[id - 1].poros = poros;
		matVector[id - 1].wDensity = rho;
	}
}

void Readsubtitles::ReadHasInitStress(void)
{
	int hasinitstress = 1;
	double istress_xx = 0.0;
	double istress_yy = 0.0;
	double istress_xy = 0.0;

	valuereader = fscanf(_readData, "%lf %lf %lf", &istress_xx, &istress_yy, &istress_xy);

	model.hasinitstress = hasinitstress;
	model.iniSxx = istress_xx;
	model.iniSxy = istress_yy;
	model.iniSyy = istress_xy;

} 


void Readsubtitles::ReadNumOfElements(void)
{
	int nelem;

	valuereader = fscanf(_readData, "%d", &nelem);
	numElements = nelem; //Set Number of elements;
	elementVector = (ClaElements*)calloc(nelem, sizeof(ClaElements));	

}

void Readsubtitles::ReadElementQ4(void)
{
	int nelem;
	int matid = 0;
	int ordid = 0;
	int t = 0;
	int id, n;

	valuereader = fscanf(_readData, "%d", &nelem);
	for (int i = 0; i < nelem; i++)
	{
		valuereader = fscanf(_readData, "%d%d%d%d", &id, &matid, &ordid, &t);
		elementVector[id - 1].id = id;
		elementVector[id - 1].matid = matid;
		elementVector[id - 1].type = 2;
		elementVector[id - 1].nnodes = 4;
		for (int j = 0; j < 4; j++)
		{
			valuereader = fscanf(_readData, "%d", &n);
			elementVector[id - 1].n[j] = n;
		}
	}

}

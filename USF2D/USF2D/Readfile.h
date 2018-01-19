#ifndef READFILE_H
#define READFILE_H

#include <cstdio>

class Readsubtitles
{
public:		
	void pointerFile(FILE *);
	int subtitles(void);
	int findsubtitle(char *);

	static void ReadHeading(void);
	static void	ReadNCells(void);
	static void ReadLimits();
	static void ReadTotalTime(void);
	static void ReadPrintStep(void);
	static void	ReadGravity(void);
	static void	ReadPointResults(void);
	static void ReadDamp(void);
	static void	ReadTotalNodes(void);
	static void ReadNodeCoord(void);
	static void ReadNumMaterial(void);
	static void ReadMatLinear(void);
	static void ReadMatFlow(void);
	static void ReadNumOfElements(void);
	static void ReadElementQ4(void);
	static void	ReadTotalParticles(void);
	static void ReadHasInitStress(void);
	static void ReadBCVelocity(void);
	static void ReadParticlesCoords(void);
	static void ReadPartVelocity(void);



};

class ClaModel
{
public:
	int ncx;
	int ncy;

	double minx;
	double miny;
	double minz;
	double maxx;
	double maxy;
	double maxz;


	int printstep;

	double gravx;
	double gravy;

	double dampfreq;
	double dampfrac;
	double dampadmul;
	double dtfrac;
	double damp;
	double dampw;

	int hasinitstress;
	double iniSxx;
	double iniSyy;
	double iniSxy;
};

class ClaCoordNodes
{
public:
	int id;
	double x, y, z;
};

class ClaMatVector
{
public:

	//Linear Material
	int id;
	int type;
	double young;
	double poisson;
	double density;

	//Material Flow
	double kxx;
	double kyy;
	double kw;
	double poros;
	double wDensity;

};

class ClaElements
{
public:
	int type;	// 2 = Q4; 1 = T3 
	int nnodes; // number of nodes to element
	int id;		// element id
	int n[4];	// element connectivity
	int matid;	// element material

	char*status;
};

class ClaParticles
{
public:
	int		nparticles; //number of particles
  //int		inside_particles;
	int		part_id;
	double	part_x;
	double	part_y;
	double	part_z;
};


class ClaPartVel
{
public:
	int           id;
	double        vx;
	double        vy;
	double        vz;
} ;

class ClaBoundaryC
{
public:
	int	   idBC;
	double vx;
	double vy;
	double vz;
};
extern int		             numPartVel;
extern int					 numBCvel;
extern int					 numNodes;
extern int					 numElements;
extern int					 numParticles;
extern double				 dtime;
extern double				 totaltime;


extern char					 readParticles[BUFSIZ];

extern ClaModel				 model;
extern ClaCoordNodes		*vectorNodes;
extern ClaMatVector			*matVector;
extern ClaElements			*elementVector;
extern ClaParticles			*partVector;
extern ClaPartVel			*PartVelocityVector;


extern ClaBoundaryC			*vectorBCvelocity;






#endif // !READFILE_H

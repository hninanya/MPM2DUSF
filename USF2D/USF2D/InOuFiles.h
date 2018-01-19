#ifndef INOUFILES_H
#define INOUFILES_H

#include <string>
#include "Readfile.h"
#include <fstream>
#include "mesh.h"
//public functions

class ExternFiles: public Readsubtitles
{
public:
	int OpenFile(int, char *[]);
	int ReadFile(void);
	int closeFile(void);
	void WriteFileHeading(int, std::ostream&);
	void PrintHistory(int, std::ostream&, int, double);
	void OutputField(int, std::ostream&, int, double);
	void vtkfile(int, std::ostream&, int, double, MeshMPM2D*);



};

#endif // !INOUFILES_H

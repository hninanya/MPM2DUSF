#include <iostream>
#include <iomanip>
#include "Newton.h"
#include "mesh.h"

using namespace std;

void initialGuess(int size_Vector, double initial_guess_value)
{
	
	for (int i = 0; i < size_Vector; i++)
	{
		pointer_unknonwn[i].C = initial_guess_value;
	}
}	
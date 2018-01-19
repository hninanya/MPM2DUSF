#ifndef NEWTON_H
#define NEWTON_H
#include <iostream>
using namespace std;
class Solver
{
public:
	//double	emax = pow(10,-10) ;
	double	iteration;
	int		counter_iteration;

	//Forward elimination
	double	bp;
	double	x_unknown;

	double	C;

	//Friend functions&Template functions
	template <class R, class S, class U>
	friend double LUFactorization(R& Vector, const S& size_Vector, U Matrix);	
};

//public functions
void	initialGuess(int, double);

template <class R, class S, class U> double LUFactorization(R& Vector, const S& size_Vector, U Matrix)
{
	double Lower[2][2] = {};
	double Upper[2][2] = {};

	double em;
	for (int i = 0; i < size_Vector; i++)
	{
		for (int j = 0; j < size_Vector; j++)
		{
			Lower[i][j] = Matrix[i][j];
		}
	}
	for (int k = 0; k < size_Vector - 1; k++)
	{
		for (int i = k + 1; i < size_Vector; i++)
		{
			em = Lower[i][k] / Lower[k][k];
			//cout << em << endl;
			Lower[i][k] = em;

			for (int j = k + 1; j < size_Vector; j++)
			{
				Lower[i][j] = Lower[i][j] - em*Lower[k][j];
			}
		}
	}
	for (int i = 0; i < size_Vector; i++)
	{
		for (int j = 0; j < size_Vector; j++)
		{
			if (j >= i)
			{
				Upper[i][j] = Lower[i][j];
			}
			else
			{
				Upper[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < size_Vector; i++)
	{
		for (int j = 0; j < size_Vector; j++)
		{
			if (j > i)
			{
				Lower[i][j] = 0;
			}
			if (i == j)
			{
				Lower[i][j] = 1;
			}
		}
	}
	//Forward elmination step to calculate b'=Vector_bp
	pointer_unknonwn[0].bp = Vector[0].b;
	for (int i = 1; i < size_Vector; i++)
	{
		pointer_unknonwn[i].bp = Vector[i].b;
		for (int j = 0; j < i; j++)
		{	
			pointer_unknonwn[i].bp = pointer_unknonwn[i].bp - Lower[i][j] * pointer_unknonwn[j].bp;
		}
	}
	//Back substitution step to calculate x_unknown = delta*C
	pointer_unknonwn[size_Vector - 1].x_unknown = pointer_unknonwn[size_Vector -
		1].bp / Upper[size_Vector - 1][size_Vector - 1];
	for (int i = size_Vector - 2; i >= 0; i--)
	{
		pointer_unknonwn[i].x_unknown = pointer_unknonwn[i].bp;
		for (int j = size_Vector - 1; j >= i + 1; j--)
		{
			pointer_unknonwn[i].x_unknown = pointer_unknonwn[i].x_unknown -
				Upper[i][j] * pointer_unknonwn[j].x_unknown;
		}
		pointer_unknonwn[i].x_unknown = pointer_unknonwn[i].x_unknown / Upper[i][i];
	}

	//Updating solution C = C + delta*C for Newton-Raphson Method

	for (int i = 0; i < size_Vector; i++)
	{
		pointer_unknonwn[i].C = pointer_unknonwn[i].C - pointer_unknonwn[i].x_unknown;
		
	}
	return 1;
}
template <class R, class S> double normVector(R& Vector, const S& size_Vector)
{
	double sumsqrt = 0;
	double norm = 0;
	
	for (int i = 0; i < size_Vector; i++)
	{
		sumsqrt += pow(Vector[i].b, 2);
	}
	norm = sqrt(sumsqrt);
	return norm;
}




#endif // !NEWTON_H

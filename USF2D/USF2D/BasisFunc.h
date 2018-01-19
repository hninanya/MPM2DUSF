#ifndef BASISFUNC_H
#define BASISFUNC_H

class Shape2D
{
public:


	double Sh;
	double GSh_xi;
	double GSh_eta;
	double GSh_x;
	double GSh_y;




	void Phi_4( double,  double);
	void GradPhi_4( double, double);
	void GradPhi_4xy(int );

	//void* Phi_9nodes(const double&, const double&);
	//void* GradPhi_9nodes(const double&, const double&);
};

inline double detMatrix2x2(double *matrix)
{
	double det2=0;
	det2 = *(matrix) **(matrix + 3) - *(matrix + 2) **(matrix + 1);
	return det2;
}

inline void invMatrix2x2(double* matrix, double *inv_matrix)
{
	double d = detMatrix2x2(matrix);
	*(inv_matrix) = *(matrix + 3)/d;
	*(inv_matrix+1) = -*(matrix + 1)/d;
	*(inv_matrix + 2) = -*(matrix + 2)/d;
	*(inv_matrix + 3) = *(matrix)/d;

	//std::cout << "Jacobian determinant = " << d << std::endl;
}

#endif // !BASISFUNC_H



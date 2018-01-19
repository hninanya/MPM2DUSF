#ifndef TENSOR_H
#define TENSOR_H

class Matrix2
{
public:
	double xx, xy, yx, yy;
	Matrix2() { xx = 0.; xy = 0; yx = 0; yy = 0; }
	Matrix2(double a) { xx = a; xy = a; yx = a; yy = a; }
	Matrix2(double a, double b, double c, double d) { xx = a; xy = b; yx = c; yy = d; }
	Matrix2(const Matrix2&m) { xx = m.xx; xy = m.xy; yx = m.yx; yy = m.yy; }
};

inline double det(const Matrix2&m) { return m.xx*m.yy - m.xy*m.yx; }
//jj


#endif // !TENSOR_H


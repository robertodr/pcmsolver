#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

#include "Getkw.h"
#include "taylor.hpp"
#include "GreensFunction.h"
#include "UniformDielectric.h"

template<class T>
UniformDielectric<T>::UniformDielectric(double dielConst) 
{
    epsilon = dielConst;
    this->uniformFlag = true;
}

template<class T>
UniformDielectric<T>::UniformDielectric(Section green) 
{
    epsilon = green.getDbl("Eps");
    this->uniformFlag = true;
}

template<class T>
double UniformDielectric<T>::evald(Vector3d &direction, Vector3d &p1, Vector3d &p2){
    double der =  direction.dot(this->gradientProbe(p1, p2))/direction.norm();
    return der * epsilon;
}

template<class T>
T UniformDielectric<T>::evalGreensFunction(T * sp, T * pp) {
	T distance = sqrt((sp[0] - pp[0]) * (sp[0] - pp[0]) +
						   (sp[1] - pp[1]) * (sp[1] - pp[1]) +
						   (sp[2] - pp[2]) * (sp[2] - pp[2]));
	return 1/(epsilon * distance);
}

template class UniformDielectric<double>;
template class UniformDielectric< taylor <double, 3, 1> >;
template class UniformDielectric< taylor <double, 3, 2> >;

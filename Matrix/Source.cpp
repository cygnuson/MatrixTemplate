

#include <iostream>

#include "Matrix.hpp"

int main()
{
	auto m = cg::Matrix<4, 1, float>::Create({
		{10},
		{10},
		{10},
		{0}
	});

	m = cg::Vector4<float>::RotationMatrixD(33, 33, 33) * m;

	std::cout << "\n" << m.ToString();


	int stop = 0;
	return 0;

}



#include <iostream>

#include "Matrix.hpp"

int main()
{
	auto x = cg::Vector4<int>::Create(2, 3, 4, 5);
	auto y = cg::Vector4<int>::Create(2, 3, 4, 5);
	auto z = x*y;


	auto m = cg::rmat::Create<3, 3, int>({
		1,2,3,
		3,4,5,
		5,6,7
	});
	auto m2 = cg::rmat::Create<3, 3, int>({
		1,2,3,
		3,4,5,
		5,6,7
	});
	auto i2 = cg::rmat::Create<3, 3, int>({
		2,0,0,
		0,2,0,
		0,0,2
	});


	std::cout << std::endl << cg::rmat::ToString(m);

	cg::rmat::ColumnSwap(m,0,2);

	std::cout << std::endl << cg::rmat::ToString(m);

	int stop = 0;
	return 0;

}

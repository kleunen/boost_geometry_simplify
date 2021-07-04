#include "simplify.hpp"
#include <iostream>
#include <chrono>

void benchmark()
{
	// Benchmark overlapping approach
	MultiPolygon poly;
	boost::geometry::read_wkt("MULTIPOLYGON (((10 50, 80 50, 80 70, 40 70, 40 30, 30 30, 30 80, 90 80, 90 40, 20 40, 20 20, 50 20, 50 90, 60 90, 60 10, 10 10, 10 50)))", poly);

  	auto const start = std::chrono::steady_clock::now();

	constexpr std::size_t count = 1000000;
	for(std::size_t i = 0; i < count; ++i) {
		MultiPolygon result;
    	boost::geometry::simplify(poly, result, 40);
	}

  	auto const end = std::chrono::steady_clock::now();
  	auto const ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  	std::cout << " time: " << (double)ms / count << std::endl;
}

int main()
{
	MultiPolygon mp;
	boost::geometry::read_wkt("MULTIPOLYGON(((0 10,10 10,10 0,0 0,0 10),(1.1 1.1,3.1 1.9,4.9 1.1,3.1 3.1,2.5 2.1,1.1 1.1)))", mp);

    // Simplify it, using distance of 0.5 units
    MultiPolygon simplified;
    boost::geometry::simplify(mp, simplified, 0.5);

    std::cout
        << "  original: " << boost::geometry::wkt(mp) << ", valid: " << std::boolalpha << boost::geometry::is_valid(mp) << std::endl
        << "simplified: " << boost::geometry::wkt(simplified) << ", valid: " << std::boolalpha << boost::geometry::is_valid(simplified) << std::endl; 
 
	benchmark();   
	return 0;
}

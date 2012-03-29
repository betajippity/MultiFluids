// Written by Peter Kutz.
// Used for math operations and constants that aren't included in the C++ standard library.

#pragma once

namespace BasicMath {

    extern const double PI;
	extern const double ONE_OVER_PI;
	extern const double TWO_PI;
	extern const double FOUR_PI;
	extern const double ONE_OVER_FOUR_PI;
	extern const double E;

	extern unsigned int mod(int x, int y); // Proper int mod with an always-positive result (unlike %).

	extern double mod(double x, double y); // Proper double mod with an always-positive result (unlike fmod).

	extern double radiansToDegrees(double radians);

	extern double degreesToRadians(double degrees);

	extern double average(double n1, double n2);
	
	extern double round(double n);

    extern double log2(double n); // std::log2 is not avaiable in all compilers.

    extern bool isNaN(double n); // std::isnan is not avaiable in all compilers.
	
}


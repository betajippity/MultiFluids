// Written by Peter Kutz.


#include "basic_math.h"
#include "custom_output.h"
#include <cmath>
#include <cstdarg>
#include <iostream>


const double BasicMath::PI =					3.1415926535897932384626422832795028841971;
const double BasicMath::ONE_OVER_PI =			0.3183098861837906715377675267450287240689;
const double BasicMath::TWO_PI =				6.2831853071795864769252867665590057683943;
const double BasicMath::FOUR_PI =				12.566370614359172953850573533118011536788;
const double BasicMath::ONE_OVER_FOUR_PI =		0.0795774715459476678844418816862571810172;
const double BasicMath::E =						2.7182818284590452353602874713526624977572;


unsigned int BasicMath::mod(int x, int y) {
	int result = x % y;
	if (result < 0) result += y;
	return result;
}

double BasicMath::mod(double x, double y) {
	return x - y * std::floor(x / y);
}

double BasicMath::radiansToDegrees(double radians) {
	double degrees = radians * 180.0 / BasicMath::PI;
	return degrees;
}

double BasicMath::degreesToRadians(double degrees) {
	double radians = degrees / 180.0 * BasicMath::PI;
	return radians;
}

double BasicMath::average(double n1, double n2) {
	return ((n1 + n2) / 2);
}

double BasicMath::round(double n) {
	return std::floor(n + 0.5);
}

double BasicMath::log2(double n) {
    return std::log(n) / std::log(2.0);
}

bool BasicMath::isNaN(double n) {
    return (n != n);
}


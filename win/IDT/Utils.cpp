#include "stdafx.h"
#include "Utils.h"


double str2num(string str) {
	istringstream iss(str);
	double num;
	iss >> num;
	return num;
}

string num2str(int num) {
	stringstream iss;
	iss << num;
	return iss.str();
}

double tanOfHalfAngle(double a, double b, double c) {
	return sqrt((a - b + c) * (a + b - c) / \
		((a + b + c) * (-a + b + c)));
}

#pragma once
#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <cmath>
using namespace std;

#define SQUARE(x) ((x)*(x))

double str2num(string str);

string num2str(int num);

double tanOfHalfAngle(double a, double b, double c);

#endif
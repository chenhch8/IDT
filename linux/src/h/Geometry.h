#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>
#include "Utils.h"
using namespace std;

typedef struct POINT {
	double x, y, z;
	vector<int> f;  // 顶点关联的面
	POINT() {}
	POINT(double a[]): x(a[0]), y(a[1]), z(a[2]) {}
	double distance(POINT &a, bool type = true) {
		double result = SQUARE(x - a.x) + SQUARE(y - a.y) + SQUARE(z - a.z);
		return type ? sqrt(result) : result;
	}
	POINT midPoint(POINT &a) {
		double _p[3];
		_p[0] = (x + a.x) / 2;
		_p[1] = (y + a.y) / 2;
		_p[2] = (z + a.z) / 2;
		return POINT(_p);
	}
	bool operator == (POINT &a) {
		return x == a.x && y == a.y && z == a.z;
	}
} Point;

typedef struct EDGE {
	int p[2];       // 边的两个端点
	int f[2];       // 边对应的两个三角面
	double length;  // 边长
	EDGE() {}
	EDGE(int _p1, int _p2, double _len, int _f1 = -1, int _f2 = -1) {
		set(_p1, _p2, _len, _f1, _f2);
	}
	EDGE(int _p1, int _p2, double _len, vector<int> &_f) {
		int _f_[2] = {-1, -1};
		for (short i = 0; i < _f.size(); ++i)
			_f_[i] = _f[i];
		set(_p1, _p2, _len, _f_[0], _f_[1]);
	}
	void set(int _p1, int _p2, double _len, int _f1 = -1, int _f2 = -1) {
		p[0] = _p1; p[1] = _p2;
		length = _len;
		f[0] = _f1; f[1] = _f2;
	}
} Edge;

typedef struct FACE {
	vector<int> e;  // 面对应的边
	vector<int> p;  // 面对应的点
	// 初始化数据，默认初始化 p
	FACE() {}
	FACE(vector<int> &d, bool type = true) {
		vector<int> &t = type ? p : e;
		t = d;
	}
	FACE(vector<int> &_e, vector<int> &_p) {
		e = _e;
		p = _p;
	}
} Face;

typedef struct EDGESLINK {
    int index;      // 边的另一个顶点
    double length;  // 边长
    vector<int> f;  // 边所属的三角面
    EDGESLINK(int _index, int _f, double _len) {
        index = _index;
        length = _len;
        f.push_back(_f);
    }
} ELink;

typedef vector<Point> PointsList;
typedef vector<EDGE>  EdgesList;
typedef vector<FACE>  FacesList;

#endif
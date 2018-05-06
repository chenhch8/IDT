#pragma once
#ifndef BASEMESH_H
#define BASEMESH_H

#include "Geometry.h"
#include "Utils.h"
#include <fstream>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <regex>
using namespace std;

class BaseMesh {
public:
	BaseMesh() {}
	~BaseMesh() {}
	// �����ļ�
	bool Load(const char *, string &);
	// д���ļ�
	bool Save(const char *, string &);
	// ��ȡ����
	void ParseData(PointsList &, FacesList &, string &);
	// ȥ���ر�
	vector<ELink> * EdgesUnique(PointsList &, FacesList &);
	// ɾ����ʱ��������
	void DeleteNodes(vector<ELink> *, int);
	// �������������ݺͱ�����
	void GetEdgesAndFaces(EdgesList &, FacesList &, vector<ELink> *, int);
};

#endif

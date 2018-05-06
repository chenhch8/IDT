#pragma once
#ifndef MESH_H
#define MESH_H

#include "Geometry.h"
#include "BaseMesh.h"
#include "IDT.h"
#include <iostream>
#include <stdio.h>
using namespace std;

class CMesh : BaseMesh {
public:
	CMesh() {}
	~CMesh() {}
	bool Load(const char *);  // ���������ļ�
	void Save(const char *);
	void iDT();               // ����iDT
private:
	PointsList points;  // ��������
	FacesList  faces;   // ����������
	EdgesList  edges;   // ������
};

#endif
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
 	bool Load(const char *);  // 解析数据文件
 	void Save(const char *);
 	void iDT();               // 启动iDT
 private:
 	PointsList points;  // 顶点链表
    FacesList  faces;   // 三角面链表
 	EdgesList  edges;   // 边链表
};

#endif
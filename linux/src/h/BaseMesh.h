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
    // 读入文件
    bool Load(const char *, string &);
    // 写出文件
    bool Save(const char *, string &);
    // 提取数据
    void ParseData(PointsList &, FacesList &, string &);
    // 去除重边
    vector<ELink> * EdgesUnique(PointsList &, FacesList &);
    // 删除临时辅助数据
    void DeleteNodes(vector<ELink> *, int);
    // 解析保存面数据和边数据
    void GetEdgesAndFaces(EdgesList &, FacesList &, vector<ELink> *, int);
};

#endif
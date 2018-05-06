#ifndef IDT_H
#define IDT_H

#include "Geometry.h"
#include <iostream>
#include <stack>
#include <queue>
#include <vector>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <cmath>
using namespace std;

typedef struct SPLITFACE {
    int fid[2];
    Face face[2];
    SPLITFACE(int fid1, Face &f1, int fid2, Face &f2) {
        fid[0] = fid1;
        fid[1] = fid2;
        face[0] = f1;
        face[1] = f2;
    }
} SplitFace;

typedef struct SPLITEDGE {
    int eid[3];
    int len;
    Edge edge[3];
    SPLITEDGE() {}
    SPLITEDGE(int id1, Edge &e1, int id2, Edge &e2, int id3, Edge &e3) {
        eid[0] = id1;
        edge[0] = e1;
        eid[1] = id2;
        edge[1] = e2;
        eid[2] = id3;
        edge[2] = e3;
        len = 3;
    }
    SPLITEDGE(int id, Edge &e) {
        eid[0] = id;
        edge[0] = e;
        len = 1;
    }
} SplitEdge;

typedef struct CROSSEDGES {
    vector<int> edgPath;
    set<int> hasPass;
    CROSSEDGES() {}
} CrossEdges;

class IDT {
 public:
 	IDT();
 	~IDT();
 	vector<int> Flip(const int edgeId);
 	bool Delaunay(const int edgeId);
 	bool IsEmpty();
 	void PushInStack();
 	void PushInStack(int id);
 	int PopStack();
 	bool GetMark(int id);
 	void SetMark(bool value);
 	void SetMark(int id, bool value);
    void initMesh(PointsList &points, EdgesList &es, FacesList &fs);
    double Melt(PointsList &ps, EdgesList &es, FacesList &fs);
 private:
    stack<int> cache;
    vector<bool> *isFliped;
    vector<bool> *mark;
    PointsList *points;
    EdgesList *edges;
    FacesList *faces;
    void clearMark();
    void clearMesh();
    int findCrossEdgesAndCrossPoints(int edgId, EdgesList &oglEgs, FacesList &oglFcs, vector<Point> &_cosPts, vector<int> &_cosEgs);
    vector<Point> findAllCrossPoints(int start, int end, const int edgId, const vector<int> &es, const EdgesList &oglEgs);
    double findCrossPoint(Point left, Point right, Point &p0, Point &p1, Point &mid);
    bool updateMesh(vector<Point> &crossPoints, vector<int> &egsId,
                    PointsList &ps, EdgesList &es, FacesList &fs,
                    vector<vector<int>> &newInOldEgs, int start, int end);
    void splitFace(bool type, int fId, Face &face, EdgesList &es, int pKey, int key1,
                   int key2, vector<int> &newPtsId, vector<int> &newEgsId,
                   vector<int> &newFcsId, int ppt, int ept, int fpt,
                   vector<SplitFace> &newFaces, bool reverse = false);
    void splitEdge(int eId, Edge &edge, int pKey, vector<int> &newPtsId,
                   vector<int> &newEgsId, vector<int> &newFcsId, int ppt,
                   int ept, int fpt, vector<SplitEdge> &newEdges);
};

#endif
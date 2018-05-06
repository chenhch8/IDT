#include "stdafx.h"
#include "BaseMesh.h"


bool BaseMesh::Load(const char *filename, string &result) {
	fstream fp;
	fp.open(filename, ios::in);
	if (!fp.is_open()) {
		cout << "Error: File No Exits!" << endl;
		return false;
	}
	// 获取文件长度
	fp.seekg(0, fp.end);
	unsigned fileLen = fp.tellg();
	fp.seekg(0, fp.beg);
	char *buffer = new char[fileLen];
	// 读取文件内容
	cout << "Reading File: " << filename << "..." << flush;
	fp.read(buffer, fileLen);
	fp.close();
	cout << "done" << endl;
	result = string(buffer);
	delete buffer;
	buffer = NULL;
	return true;
}

bool BaseMesh::Save(const char *filename, string &result) {
	fstream fp;
	fp.open(filename, ios::out);
	if (!fp.is_open()) {
		cout << "Error: Writing Fail!" << endl;
		return false;
	}
	cout << "Writing File: " << filename << "..." << flush;
	fp << result;
	fp.close();
	cout << "done" << endl;
	return true;
}

void BaseMesh::ParseData(PointsList &points, FacesList &faces, string &content) {
	regex reg("((v(\\s-?(\\d*\\.)?\\d+(e-\\d+)?){3})|(f(\\s\\d+){3}))", regex::icase),
		reg2("-?(\\d*\\.)?\\d+(e-\\d+)?", regex::icase);
	string ss;
	double nums[3];
	short i = 0, j;
	for (sregex_iterator it(content.begin(), content.end(), reg), end_it; it != end_it; it++) {
		ss = it->str();
		for (sregex_iterator it2(ss.begin(), ss.end(), reg2); it2 != end_it; it2++, i = ++i % 3)
			nums[i] = str2num(it2->str());
		if (ss[0] == 'v') {
			points.push_back(Point(nums));
		}
		else {
			vector<int> nums2(3);
			for (j = 0; j < 3; ++j) {
				nums2[j] = nums[j] - 1;
				points[nums2[j]].f.push_back(faces.size());
			}
			Face tmp = Face(nums2, true);
			sort(nums2.begin(), nums2.end());
			tmp.e = nums2;
			faces.push_back(tmp);
		}
	}
}

vector<ELink> * BaseMesh::EdgesUnique(PointsList &points, FacesList &faces) {
	vector<ELink> *node = new vector<ELink>[points.size()];
	int i, j, k, u, v1, v2;
	double len;
	for (i = 0; i < faces.size(); ++i) {
		for (j = 0; j < 2; ++j) {
			v1 = faces[i].e[j];
			for (k = j + 1; k < 3; ++k) {
				v2 = faces[i].e[k];
				len = points[v1].distance(points[v2]);
				if (node[v1].empty()) {
					node[v1].push_back(ELink(v2, i, len));
				}
				else {
					u = 0;
					while (u < node[v1].size() && node[v1][u].index != v2) {
						u++;
					}
					if (u == node[v1].size()) {
						node[v1].push_back(ELink(v2, i, len));
					}
					else {
						node[v1][u].f.push_back(i);
					}
				}
			}
		}
	}
	return node;
}

void BaseMesh::DeleteNodes(vector<ELink> *node, int len) {
	if (!node) return;
	delete[]node;
	node = NULL;
}

void BaseMesh::GetEdgesAndFaces(EdgesList &edges, FacesList &faces, vector<ELink> *node, int len) {
	vector<short> fe(faces.size(), 0);
	int i, j, k;
	for (i = 0; i < len; ++i) {
		for (j = 0; j < node[i].size(); ++j) {
			edges.push_back(Edge(i, node[i][j].index, node[i][j].length, node[i][j].f));
			for (k = 0; k < node[i][j].f.size(); ++k) {
				faces[node[i][j].f[k]].e[fe[node[i][j].f[k]]++] = edges.size() - 1;
			}
		}
	}
}

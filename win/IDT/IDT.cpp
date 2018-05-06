#include "stdafx.h"
#include "IDT.h"

IDT::IDT() {
	mark = NULL;
	points = NULL;
	edges = NULL;
	faces = NULL;
	isFliped = NULL;
}

IDT::~IDT() {
	clearMark();
	clearMesh();
}

void IDT::clearMesh() {
	if (faces) {
		delete points;
		delete faces;
		delete edges;
		delete isFliped;
		points = NULL;
		faces = NULL;
		edges = NULL;
		isFliped = NULL;
	}
}

void IDT::clearMark() {
	if (mark) {
		delete mark;
		mark = NULL;
	}
}

vector<int> IDT::Flip(const int edgId) {
	(*isFliped)[edgId] = true;
	// [1] 找出三角形的三条边 a, b, c 和边 edgId 对应的两个对顶角
	const int *fId = (*edges)[edgId].f,
		*pId = (*edges)[edgId].p;
	int i, j;
	int a[2] = { -1, -1 }, c[2] = { -1, -1 };
	int p[2];  // 四边形的非边 edgId 上的两个顶点编号
			   // 遍历两个三角形
	for (i = 0; i < 2; ++i) {
		vector<int> &e = (*faces)[fId[i]].e;
		// 遍历三角形的三条边
		for (j = 0; j < 3; ++j) {
			if (e[j] == edgId) continue;
			if ((*edges)[e[j]].p[0] == pId[0] || (*edges)[e[j]].p[1] == pId[0]) {
				c[i] = e[j];
				p[i] = (*edges)[e[j]].p[0] == pId[0] ? (*edges)[e[j]].p[1] : (*edges)[e[j]].p[0];
			}
			else {
				a[i] = e[j];
			}
		}
	}
	// [2] 计算新边长度
	double tanA[2], tanB, cosB, len;
	for (i = 0; i < 2; ++i)
		tanA[i] = tanOfHalfAngle((*edges)[a[i]].length, (*edges)[c[i]].length, (*edges)[edgId].length);
	tanB = (tanA[0] + tanA[1]) / (1 - tanA[0] * tanA[1]);
	cosB = (1 - SQUARE(tanB)) / (1 + SQUARE(tanB));
	len = sqrt(SQUARE((*edges)[c[0]].length) + SQUARE((*edges)[c[1]].length) - \
		2 * (*edges)[c[0]].length * (*edges)[c[1]].length * cosB);
	// [3] 边翻转——边面信息更新
	// 更新面
	for (i = 0; i < 3; ++i) {  // 更新面的三顶点
		if ((*faces)[fId[0]].p[i] == pId[1])
			(*faces)[fId[0]].p[i] = p[1];
		if ((*faces)[fId[1]].p[i] == pId[0])
			(*faces)[fId[1]].p[i] = p[0];
	}
	int tmp1[] = { edgId, c[0], c[1] },
		tmp2[] = { edgId, a[0], a[1] };
	(*faces)[fId[0]].e = vector<int>(tmp1, tmp1 + 3);  // 更新面的三条边
	(*faces)[fId[1]].e = vector<int>(tmp2, tmp2 + 3);
	// 更新边
	for (i = 0; i < 2; ++i) {
		if ((*edges)[c[1]].f[i] == fId[1])
			(*edges)[c[1]].f[i] = fId[0];
		if ((*edges)[a[0]].f[i] == fId[0])
			(*edges)[a[0]].f[i] = fId[1];
	}
	(*edges)[edgId].set(p[0], p[1], len, fId[0], fId[1]);
	// [4] 返回边
	int tmp[4] = { a[0], a[1], c[0], c[1] };
	return vector<int>(tmp, tmp + 4);
}

bool IDT::Delaunay(const int edgId) {
	const int *fId = (*edges)[edgId].f;
	if (fId[1] == -1) {  // boundary edge
		return true;
	}
	else {  // interior edge
		double cots[2], half_angle,
			a, b, c;
		int i, j;
		for (i = 0; i < 2; ++i) {
			vector<int> &e = (*faces)[fId[i]].e;
			a = (*edges)[edgId].length;
			b = c = -1;
			for (j = 0; j < 3; ++j) {
				if (e[j] == edgId) continue;
				if (b == -1) {
					b = (*edges)[e[j]].length;
				}
				else {
					c = (*edges)[e[j]].length;
				}
			}
			half_angle = tanOfHalfAngle(a, b, c);
			cots[i] = (1 - SQUARE(half_angle)) / (2 * half_angle);
		}
		return cots[0] + cots[1] >= 0 ? true : false;
	}
}

double IDT::Melt(PointsList &ps, EdgesList &es, FacesList &fs) {
	EdgesList oglEgs(es.begin(), es.end());
	FacesList oglFcs(fs.begin(), fs.end());
	vector<vector<int>> newInOldEgs(es.size());
	vector<int> cEgs;
	vector<Point> cPts;
	int start, end;
	int count = 0, sum = 0;
	for (int i = 0; i < (*isFliped).size(); ++i)
		if ((*isFliped)[i]) {
			sum++;
			cPts.clear();
			cEgs.clear();
			start = findCrossEdgesAndCrossPoints(i, oglEgs, oglFcs, cPts, cEgs);
			end = (*edges)[i].p[(*edges)[i].p[0] == start ? 1 : 0];
			if (cPts.empty()) {
				count++;
			}
			else if (!updateMesh(cPts, cEgs, ps, es, fs, newInOldEgs, start, end)) {
				count++;
			}
		}
	return (double)count / sum;
}

void IDT::PushInStack() {
	for (int i = 0; i < (*edges).size(); ++i)
		cache.push(i);
}

void IDT::PushInStack(int id) {
	cache.push(id);
}

int IDT::PopStack() {
	int id = cache.top();
	cache.pop();
	return id;
}

bool IDT::IsEmpty() {
	return cache.empty();
}

bool IDT::GetMark(int id) {
	return (*mark)[id];
}

void IDT::SetMark(int id, bool value = true) {
	(*mark)[id] = value;
}

void IDT::SetMark(bool value = true) {
	clearMark();
	mark = new vector<bool>((*edges).size(), value);
}

void IDT::initMesh(PointsList &ps, EdgesList &es, FacesList &fs) {
	clearMesh();
	points = new PointsList(ps.begin(), ps.end());
	edges = new EdgesList(es.begin(), es.end());
	faces = new FacesList(fs.begin(), fs.end());
	isFliped = new vector<bool>(es.size(), false);
}

int IDT::findCrossEdgesAndCrossPoints(int edgId, EdgesList &oglEgs, FacesList &oglFcs,
	vector<Point> &_cosPts, vector<int> &_cosEgs) {
	_cosEgs.clear();
	_cosPts.clear();
	int start = (*edges)[edgId].p[0],
		end = (*edges)[edgId].p[1];
	int i, j, countI, countJ;
	stack<CrossEdges> que;
	const int *pt, *ft;
	vector<int> &f = (*points)[start].f;
	countI = 0;
	vector<bool> iflag(f.size(), false);
	for (i = rand() % f.size(); countI < f.size(); i = rand() % f.size()) {
		if (iflag[i]) continue;
		iflag[i] = true;
		countI++;
		vector<int> &e = oglFcs[f[i]].e;
		CrossEdges tmp;
		tmp.hasPass.insert(f[i]);
		countJ = 0;
		vector<bool> jflag(f.size(), false);
		for (j = rand() % e.size(); countJ < e.size(); j = ++j % e.size()) {
			if (jflag[j]) continue;
			jflag[j] = true;
			countJ++;
			pt = oglEgs[e[j]].p;
			if (pt[0] != start && pt[1] != start) {
				tmp.edgPath.push_back(e[j]);
				break;
			}
		}
		que.push(tmp);
	}
	int last;
	// vector<vector<Point>> crossPoints;
	// vector<vector<int>> crossEdges;
	bool stop = false;
	while (!que.empty() && !stop) {
		// cout << que.size() << endl;
		CrossEdges ce = que.top();
		que.pop();
		if (ce.hasPass.size() > 5) continue;
		last = ce.edgPath.size() - 1;
		ft = oglEgs[ce.edgPath[last]].f;
		countI = 0;
		for (i = rand() % 2; countI < 2; i = ++i % 2) {
			countI++;
			if (ft[i] != -1 && ce.hasPass.find(ft[i]) == ce.hasPass.end()) {
				vector<int> &p = oglFcs[ft[i]].p;
				for (j = 0; j < p.size(); ++j) {
					if (p[j] == end) {
						vector<Point> cspts = findAllCrossPoints(start, end, edgId, ce.edgPath, oglEgs);
						if (!cspts.empty()) {
							// crossEdges.push_back(ce.edgPath);
							// crossPoints.push_back(cspts);
							_cosEgs = ce.edgPath;
							_cosPts = cspts;
							stop = true;
						}
						break;
					}
				}
				if (j == p.size()) {
					vector<int> &e = oglFcs[ft[i]].e;
					vector<bool> flag(e.size(), false);
					countJ = 0;
					for (j = rand() % e.size(); countJ < e.size() - 1; j = rand() % e.size()) {
						if (e[j] == ce.edgPath[last] || flag[j]) continue;
						flag[j] = true;
						countJ++;
						CrossEdges next;
						next.hasPass = ce.hasPass;
						next.edgPath = ce.edgPath;
						next.hasPass.insert(ft[i]);
						next.edgPath.push_back(e[j]);
						que.push(next);
					}
				}
			}
		}
	}
	// if (crossPoints.size() > 1 || !crossPoints.size())
	//     cout << crossPoints.size() << endl;
	// if (crossPoints.size() == 1) {
	//     _cosPts = crossPoints[0];
	//     _cosEgs = crossEdges[0];
	// }
	return start;
}

vector<Point> IDT::findAllCrossPoints(int start, int end, const int edgId, const vector<int> &es,
	const EdgesList &oglEgs) {
	vector<Point> result(es.size() + 2);
	int i;
	const int *p;
	result[0] = (*points)[start];
	result[result.size() - 1] = (*points)[end];
	end = result.size() - 1;
	for (i = 1; i < end; ++i) {
		p = oglEgs[es[i - 1]].p;
		result[i] = (*points)[p[0]].midPoint((*points)[p[1]]);
	}
	double lenMin = 1024000, len, gap;
	int thread = 1000;
	do {
		len = 0;
		for (i = 1; i < end; ++i) {
			p = oglEgs[es[i - 1]].p;
			len += findCrossPoint((*points)[p[0]], (*points)[p[1]],
				result[i - 1], result[i + 1], result[i]);
		}
		len += result[end - 1].distance(result[end]);
		if (len <= lenMin) {
			gap = lenMin - len;
			lenMin = len;
		}
	} while (gap && thread--);
	if (gap != 0 || fabs(lenMin - (*edges)[edgId].length) > LENTHREAD) result.clear();
	// #if DEBUG == ON
	//     if (!gap)
	//         cout << lenMin << " " << (*edges)[edgId].length << (fabs(lenMin - (*edges)[edgId].length) > 1e-10) << endl;
	// #endif
	return result;
}

double IDT::findCrossPoint(Point left, Point right, Point &p0,
	Point &p1, Point &mid) {
	mid = left.midPoint(right);
	double lenL, lenR, lenMin = 1024000, gap;
	do {
		Point m1 = left.midPoint(mid),
			m2 = right.midPoint(mid);
		lenL = m1.distance(p0) + m1.distance(p1);
		lenR = m2.distance(p0) + m2.distance(p1);
		if (lenL < lenR) {
			right = mid;
			mid = m1;
			gap = fabs(lenL - lenMin);
			lenMin = lenL;
		}
		else if (lenL > lenR) {
			left = mid;
			mid = m2;
			gap = fabs(lenR - lenMin);
			lenMin = lenR;
		}
		else {
			break;
		}
	} while (gap);
	return mid.distance(p0);
}

bool IDT::updateMesh(vector<Point> &crossPoints, vector<int> &egsId,
	PointsList &ps, EdgesList &es, FacesList &fs,
	vector<vector<int>> &newInOldEgs, const int start, const int end) {
	int i, j, k, eId;
	// [1] 找出旧边上与目标边相交的子边
	double disA, disB, disC;
	vector<int> sonEgsId(egsId.size(), -1);
	for (i = 0; i < egsId.size(); ++i) {
		eId = egsId[i];
		vector<int> &carry = newInOldEgs[eId];
		if (carry.empty())
			carry.push_back(eId);
		for (j = 0; j < carry.size(); ++j) {
			const int *p = es[carry[j]].p;
			disA = ps[p[0]].distance(ps[p[1]], false);
			disB = ps[p[0]].distance(crossPoints[i + 1], false);
			disC = ps[p[1]].distance(crossPoints[i + 1], false);
			if (disA > disB && disA > disC) {
				sonEgsId[i] = carry[j];
				break;
			}
		}
	}
	for (i = 0; i < sonEgsId.size(); ++i) {
		if (sonEgsId[i] == -1) {
			// cout << "Error here: 256" << endl;
			return false;
		}
	}
	// [2] 设置新边、面、点 id
	vector<int> newEgsId(2 * crossPoints.size() - 3),
		newFcsId(crossPoints.size() - 1),
		newPtsId(crossPoints.size() - 2);
	for (i = 0; i < newEgsId.size(); ++i)
		newEgsId[i] = es.size() + i;
	for (i = 0; i < newFcsId.size(); ++i)
		newFcsId[i] = fs.size() + i;
	for (i = 0; i < newPtsId.size(); ++i)
		newPtsId[i] = ps.size() + i;
	// [3] 添加新点
	for (i = 1; i < crossPoints.size() - 1; ++i)
		ps.push_back(crossPoints[i]);
	// [4] 更新边、面信息
	const int *f = NULL;
	int ppt = 0, ept = 0, fpt = 0, fId;
	int sonEgId, polygon, pKey, pKey2;
	vector<SplitFace> newFaces;
	vector<SplitEdge> newEdges;
	// 对每条分裂边
	for (i = 0; i < sonEgsId.size(); ++i) {
		sonEgId = sonEgsId[i];
		f = es[sonEgId].f;  // f 在后面的指向不能变！
							// 4.1 确定要分裂的面
							// 对面f[0]的每条边/点
		vector<int> &tmpE = fs[f[0]].e;
		vector<int> &tmpP = fs[f[0]].p;
		polygon = tmpE.size();
		for (j = 0; j < polygon; ++j) {
			if (tmpP[j] == start || i > 0 && tmpE[j] == sonEgsId[i - 1])
				break;
		}
		fId = f[j == polygon ? 1 : 0];  // fId 即为目标分裂面
		vector<int> &e = fs[fId].e;
		vector<int> &p = fs[fId].p;
		polygon = e.size();
#if DEBUG == ON
		for (j = 0; j < polygon; ++j) {
			if (p[j] == start || i > 0 && e[j] == sonEgsId[i - 1])
				break;
		}
		if (j == polygon) {
			cout << "Error here: 304" << endl;
			exit(1);
		}
#endif
		if (i == 0) {  // 点边分割
					   // 更新 pKey
			pKey = es[sonEgId].p[0];
			pKey2 = es[sonEgId].p[1];
			splitFace(true, fId, fs[fId], es, pKey, start, sonEgId,
				newPtsId, newEgsId, newFcsId, ppt, ept, fpt,
				newFaces);
		}
		else {  // 边边分割
				// 更新 pKey
			for (j = 0;; j = ++j % polygon) {
				if (e[j] == sonEgId || e[j] == sonEgsId[i - 1])
					continue;
				if (pKey == es[sonEgId].p[0] || pKey == es[sonEgId].p[1])
					break;
				if (es[e[j]].p[0] == pKey && es[e[j]].p[1] != pKey2)
					pKey = es[e[j]].p[1];
				else if (es[e[j]].p[1] == pKey && es[e[j]].p[0] != pKey2)
					pKey = es[e[j]].p[0];

			}
			pKey2 = es[sonEgId].p[es[sonEgId].p[0] == pKey ? 1 : 0];
			splitFace(false, fId, fs[fId], es, pKey, sonEgsId[i - 1], sonEgId,
				newPtsId, newEgsId, newFcsId, ppt, ept, fpt, newFaces);
		}
		splitEdge(sonEgId, es[sonEgId], pKey, newPtsId, newEgsId, newFcsId,
			ppt, ept, fpt, newEdges);
		Edge e1;
		if (i == 0)
			e1.set(start, newPtsId[ppt], -1, fId, newFcsId[fpt]);
		else
			e1.set(newPtsId[ppt - 1], newPtsId[ppt], -1, fId, newFcsId[fpt]);
		SplitEdge &spe = newEdges[newEdges.size() - 1];
		spe.eid[0] = newEgsId[ept];
		spe.edge[0] = e1;
		ppt += 1;
		fpt += 1;
		ept += 2;
	}
	fId = f[fId == f[0] ? 1 : 0];
	splitFace(true, fId, fs[fId], es, pKey, end, sonEgId,
		newPtsId, newEgsId, newFcsId, ppt - 1, ept, fpt,
		newFaces, true);
	Edge e1(newPtsId[ppt - 1], end, -1, fId, newFcsId[fpt]);
	newEdges.push_back(SplitEdge(newEgsId[ept], e1));
	// [5] 更新 newInOldEgs 和网格
	int last = newFaces.size() - 1;
	for (i = 0; i < last; ++i) {
		SplitFace &spf = newFaces[i];
		SplitEdge &spe = newEdges[i];

		fs[spf.fid[0]] = spf.face[0];
		fs.push_back(spf.face[1]);
		es[spe.eid[1]] = spe.edge[1];
		es.push_back(spe.edge[0]);
		es.push_back(spe.edge[2]);
		newInOldEgs[egsId[i]].push_back(spe.eid[2]);
#if DEBUG == ON
		bool error = false;
		if (es.size() - 1 != spe.eid[2]) {
			cout << "Error here: 365" << endl;
			error = true;
		}
		if (fs.size() - 1 != spf.fid[1]) {
			cout << "Error here: 369" << endl;
			error = true;
		}
		if (error) exit(1);
#endif
	}
	SplitFace &spf = newFaces[last];
	SplitEdge &spe = newEdges[last];
	fs[spf.fid[0]] = spf.face[0];
	fs.push_back(spf.face[1]);
	es.push_back(spe.edge[0]);
	return true;
}

// type=true:  点边分割—— key1 为点 id，key2 为边 id
// type=false: 边边分割—— key1 和 key2 均为边 id
void IDT::splitFace(bool type, int fId, Face &face, EdgesList &es, int pKey, int key1,
	int key2, vector<int> &newPtsId, vector<int> &newEgsId,
	vector<int> &newFcsId, int ppt, int ept, int fpt,
	vector<SplitFace> &newFaces, bool reverse) {
	Face f1, f2;
	int polygon = face.p.size();
	// 更新f1、f2
	bool done1 = false, done2 = false, key2In = false;
	const int *p;
	int *f;
	int pkey = pKey,
		pkey2 = es[key2].p[es[key2].p[0] == pKey ? 1 : 0];
	for (int i = 0; !done1 || !done2; i = ++i % polygon) {
		if (face.e[i] == key2) {
			if (!done1 && !key2In) {
				f1.e.push_back(face.e[i]);
				key2In = true;
			}
			continue;
		}
		p = es[face.e[i]].p;
		if ((p[0] == pkey || p[1] == pkey) && !done1) {
			f1.p.push_back(pkey);
			f1.e.push_back(face.e[i]);
			pkey = p[p[0] == pkey ? 1 : 0];
			if (key1 == (type ? pkey : face.e[i])) {
				if (!key2In) {
					f1.e.push_back(key2);
					key2In = true;
				}
				done1 = true;
			}
		}
		else if ((p[0] == pkey2 || p[1] == pkey2) && !done2) {
			f2.p.push_back(pkey2);
			if (type || !type && face.e[i] != key1) {
				f2.e.push_back(face.e[i]);
				f = es[face.e[i]].f;
				if (f[0] == fId)
					f[0] = newFcsId[fpt];
#if DEBUG == OFF
				else
					f[1] = newFcsId[fpt];
#elif DEBUG == ON
				else if (f[1] == fId)
					f[1] = newFcsId[fpt];
				else {
					cout << "f 错误: 426" << endl;
					exit(1);
				}
#endif
			}
			pkey2 = p[p[0] == pkey2 ? 1 : 0];
			if (key1 == (type ? pkey2 : face.e[i]))
				done2 = true;
		}
	}
	f1.e.push_back(newEgsId[ept]);
	f2.e.push_back(newEgsId[ept]);
	if (!reverse)
		f2.e.push_back(newEgsId[ept + 1]);
	else
		f2.e.push_back(newEgsId[ept - 1]);
	int tmpP;
	if (!type) {
		f2.e.push_back(newEgsId[ept - 1]);
		tmpP = newPtsId[ppt - 1];
	}
	else {
		tmpP = key1;
	}
	f1.p.push_back(tmpP);
	f2.p.push_back(tmpP);
	f1.p.push_back(newPtsId[ppt]);
	f2.p.push_back(newPtsId[ppt]);
	// 顶点逆时针存储
	vector<int> *ptr = &f1.p;
	if (f2.p.size() > 3)
		ptr = &f2.p;
	for (int i = 0;; ++i)
		if (face.p[i] == (*ptr)[0]) {
			if (face.p[(i + 1) % polygon] == (*ptr)[1])
				ptr = &(ptr == &f1.p ? f2 : f1).p;
			break;
		}
	for (int i = 0; i < (*ptr).size() / 2; ++i) {
		(*ptr)[i] += (*ptr)[(*ptr).size() - 1 - i];
		(*ptr)[(*ptr).size() - 1 - i] = (*ptr)[i] - (*ptr)[(*ptr).size() - 1 - i];
		(*ptr)[i] -= (*ptr)[(*ptr).size() - 1 - i];
	}
	set<int> help1, help2;
	vector<int> pp, ee;
	pp = f1.p; f1.p.clear();
	ee = f1.e; f1.e.clear();
	for (int i = 0; i < pp.size(); ++i) {
		if (help1.find(pp[i]) == help1.end()) {
			help1.insert(pp[i]);
			f1.p.push_back(pp[i]);
		}
		if (help2.find(ee[i]) == help2.end()) {
			help2.insert(ee[i]);
			f1.e.push_back(ee[i]);
		}
	}
	help1.clear();
	help2.clear();
	pp = f2.p; f2.p.clear();
	ee = f2.e; f2.e.clear();
	for (int i = 0; i < pp.size(); ++i) {
		if (help1.find(pp[i]) == help1.end()) {
			help1.insert(pp[i]);
			f2.p.push_back(pp[i]);
		}
		if (help2.find(ee[i]) == help2.end()) {
			help2.insert(ee[i]);
			f2.e.push_back(ee[i]);
		}
	}
	newFaces.push_back(SplitFace(fId, f1, newFcsId[fpt], f2));
}

// type=false: 更新三条边
void IDT::splitEdge(int eId, Edge &edge, int pKey, vector<int> &newPtsId,
	vector<int> &newEgsId, vector<int> &newFcsId, int ppt,
	int ept, int fpt, vector<SplitEdge> &newEdges) {
	SplitEdge spe;
	Edge e2 = edge,
		e3 = edge;
	if (edge.p[0] == pKey) {
		e2.p[1] = newPtsId[ppt];
		e3.p[0] = newPtsId[ppt];
	}
	else {
		e2.p[0] = newPtsId[ppt];
		e3.p[1] = newPtsId[ppt];
	}
	e3.f[0] = newFcsId[fpt];
	e3.f[1] = newFcsId[fpt + 1];
	e2.length = e3.length = -1;
	spe.eid[1] = eId;
	spe.edge[1] = e2;
	spe.eid[2] = newEgsId[ept + 1];
	spe.edge[2] = e3;
	newEdges.push_back(spe);
}

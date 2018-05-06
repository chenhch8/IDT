#include "../h/CMesh.h"

bool CMesh::Load(const char *filename) {
    string filecontent;
    // [1] 读取 obj 文件
    if (!BaseMesh::Load(filename, filecontent)) {
        return false;
    }
    // [2] 解析 obj 文件，将结果保存入 points、edges 和 faces 链表中
    cout << "Parsingn Obj file..." << flush;
    // [2.1] 解析数据
    ParseData(points, faces, filecontent);
    // [2.2] 数据归类
    // [2.2.1] 边查重，建立并返回临时链表
    vector<ELink> *head = EdgesUnique(points, faces);
    // [2.2.2] 存入边、面数据
    GetEdgesAndFaces(edges, faces, head, points.size());
    // [2.2.3] 删除临时链表
    DeleteNodes(head, points.size());
    cout << "done" << endl;
    cout << "\tPoints: " << points.size() << "\n\tEdges: " << edges.size() << "\n\tFaces: " << faces.size() << endl;
    return true;
}

void CMesh::Save(const char *filename) {
    string filecontent = "####\n"\
                         "#\n"\
                         "# Vertices: " + num2str(points.size()) + "\n"\
                         "# Faces: " + num2str(faces.size()) + "\n"\
                         "#\n"\
                         "####\n";
    unsigned i, j;
    #define SIZE 10240
    char buffer[SIZE];
    for (i = 0; i < points.size(); ++i) {
        snprintf(buffer, SIZE, "v %lf %lf %lf\n", points[i].x, points[i].y, points[i].z);
        filecontent += string(buffer);
    }
    filecontent += "# " + num2str(points.size()) + " vertices, 0 vertices normals\n\n";
    for (i = 0; i < faces.size(); ++i) {
        string str = "f";
        for (j = 0; j < faces[i].p.size(); ++j) {
            str += " " + num2str(faces[i].p[j] + 1);
        }
        str += "\n";
        filecontent += str;
    }
    #undef SIZE
    filecontent += "# " + num2str(faces.size()) + " faces, 0 coords texture\n\n" + \
                   "# End of file";
    BaseMesh::Save(filename, filecontent);
}

void CMesh::iDT() {
    int e, i;
    cout << "start IDT..." << flush;
    IDT *idt = new IDT();
    idt->initMesh(points, edges, faces);
    idt->SetMark(true);
    idt->PushInStack();
    while (!idt->IsEmpty()) {
        e = idt->PopStack();
        idt->SetMark(e, false);
        if (!idt->Delaunay(e)) {
            vector<int> es = idt->Flip(e);
            for (i = 0; i < es.size(); ++i) {
                if (!idt->GetMark(es[i])) {
                    idt->PushInStack(es[i]);
                    idt->SetMark(es[i], true);
                }
            }
        }
    }
    int p = points.size(), f = faces.size();
    e = edges.size();
    cout << "done\nproduct new mesh..." << flush;
    double miss = idt->Melt(points, edges, faces);
    cout << "done" << endl;
    cout << "\tPoints: " << points.size() << " (+" << (points.size() - p) << ")" \
         << "\n\tEdges: " << edges.size() << " (+" << (edges.size() - e) << ")" \
         << "\n\tFaces: " << faces.size() << " (+" << (faces.size() - f) << ")" << endl;
    cout << "新边丢失率：" << miss << endl;
    delete idt;
    idt = NULL;
}
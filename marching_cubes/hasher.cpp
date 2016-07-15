
#include "hasher.h"

using namespace std;
using namespace boost;

typedef unsigned long hash_id;

//construct hashmap and points list
pointHashmap::pointHashmap(vector<TRIANGLE> triangles, gridPtr g){
    dims = g->dims;
    //go through list of triangles and generate unordered map and points list
    numPoints=0;
    for(int i=0; i<triangles.size(); i++){
        for(int j=0; j<3; j++){
            hash_id id = hashXYZ(triangles[i].p[j]);
            if(map.find(id)==map.end()){
                map[id]=numPoints;
                points.push_back(triangles[i].p[j]);
                numPoints++;
            }
        }
    }
}

//get index in points list of point
//returns -1 if not in list
int pointHashmap::index(const XYZ& pnt){
    hash_id id = hashXYZ(pnt);
    if(map.find(id)==map.end()) return -1;
    return map.at(id);
}

#ifndef HASHER_H
#define HASHER_H

#include <Eigen/Eigen>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <boost/unordered_map.hpp>

#include "grid.h"

using namespace std;
using namespace boost;

typedef unsigned long hash_id;

class XYZ{
public:
    float x;
    float y;
    float z;

    XYZ(){}
    XYZ(float x_, float y_, float z_){
        x=x_; y=y_; z=z_;
    }

    static float dist(const XYZ& lhs, const XYZ& rhs){
        float d = (rhs.x-lhs.x)*(rhs.x-lhs.x);
        d+=(rhs.y-lhs.y)*(rhs.y-lhs.y);
        d+=(rhs.z-lhs.z)*(rhs.z-lhs.z);
        return sqrt(d);
    }

    friend bool operator==(const XYZ& lhs, const XYZ& rhs){
        return dist(lhs, rhs)<0.000001;
    }
    friend bool operator!=(const XYZ& lhs, const XYZ& rhs){
        return !(lhs==rhs);
    }
    friend bool operator<(const XYZ& lhs, const XYZ& rhs){
        if (lhs.x < rhs.x)
            return true;
        else if (lhs.x > rhs.x)
            return false;

        if (lhs.y < rhs.y)
            return true;
        else if (lhs.y > rhs.y)
            return false;

        if (lhs.z < rhs.z)
            return true;
        else if (lhs.z > rhs.z)
            return false;

        return false;
    }
    friend bool operator<=(const XYZ& lhs, const XYZ& rhs){
        return ((lhs<rhs)||(lhs==rhs));
    }
    friend bool operator>(const XYZ& lhs, const XYZ& rhs){
        return !(lhs<=rhs);
    }
    friend bool operator>=(const XYZ& lhs, const XYZ& rhs){
        return !(lhs<rhs);
    }

    XYZ& operator=(const XYZ& other){
        if (this != &other) {
            this->x=other.x;
            this->y=other.y;
            this->z=other.z;
        }
        return *this;
    }
};

typedef struct {
   XYZ p[3];
} TRIANGLE;


class pointHashmap{
public:
    vector<XYZ> points;
    int numPoints;

    //construct hashmap and points list
    pointHashmap(vector<TRIANGLE> triangles, gridPtr g);

    //get index in points list of point
    //returns -1 if not in list
    int index(const XYZ& pnt);

private:
    Eigen::Vector3i dims;
    unordered_map<hash_id, int> map;
    hash_id hashXYZ(const XYZ& s){
        hash_id val = (hash_id)((double)s.x*1000000000.0*(double)(dims[1]*dims[2]));
        val+=(hash_id)((double)s.y*1000000.0*(double)dims[2]);
        val+=(hash_id)((double)s.z*1000.0);
        return val;
    }

};

#endif

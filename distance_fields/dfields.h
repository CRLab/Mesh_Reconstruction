#ifndef GETSQDIST_H
#define GETSQDIST_H

#include <Eigen/Eigenvalues>

#include "grid.h"

using namespace std;

typedef boost::shared_ptr<grid> gridPtr;

//extract perimeter of binary volume
gridPtr fastPerim(gridPtr volume_grid);

//implementation of Saito and Toriwaki distance field
gridPtr getsqdist(gridPtr volume_grid);

//returns linear indexes of closest 0-value voxel
gridPtr getsqdist_index(gridPtr volume_grid);

//get square root of grid
gridPtr getsqrt(gridPtr g);

//get linear indices of neighboring voxels
vector<int> getNeighbors(gridPtr g, const Eigen::Vector3i &pnt);

//check if voxel is on surface
bool isSurface(gridPtr g, const Eigen::Vector3i &pnt);

//get linear indices of neighboring voxels on the surface
vector<int> getSurfaceNeighbors(gridPtr g, const Eigen::Vector3i &pnt);

//get centroid of a list of points
Eigen::Vector3f getCentroid(gridPtr g, const vector<int> &pnts);

//get covariance matrix of surface point
Eigen::Matrix3f getCovariance(gridPtr g, const Eigen::Vector3i &pnt);

//get normal vector from covariance matrix
Eigen::Vector3f getNormalVector(const Eigen::Matrix3f &covariance);

//orient normal vector with respect to local centroid
void orientNormal(gridPtr g, Eigen::Vector3f &normal, const Eigen::Vector3i &pnt);

//get indexes of surface points of volume grid
vector<int> getSurface(gridPtr volume);

//get surface normals
vector<Eigen::Vector3f> getSurfaceNormals(gridPtr volume, const vector<int> &surface);

#endif

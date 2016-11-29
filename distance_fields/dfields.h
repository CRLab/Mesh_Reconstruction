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

//*****************************************************************************************************
//functions for computing normals of surface

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

//returns index in vector of value, -1 if not contained
int contains(vector<int> vec, int val);

//create normals and visualize in pcl viewer
void visualizeNormals(gridPtr volume);

//*****************************************************************************************************
//feature detection using normals

//input is binary volume pre-smoothing
//binary grid: 1=feature, 0=no feature
//features can be ignored during smoothing
//reasonable threshold ~0.9
gridPtr getFeatureMap(gridPtr volume, gridPtr surfaceMap, const vector<Eigen::Vector3f> &normals, float threshold);

gridPtr dfield_gpu(gridPtr volume_grid);

#endif

#ifndef GRID_H
#define GRID_H

#include <pcl/point_types.h>
#include <boost/thread/thread.hpp>
#include <pcl/filters/voxel_grid.h>

#include <stdlib.h>

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;

class grid{
private:
    float*** voxels;

    void allocGrid();
    void deallocGrid();

public:
    Eigen::Vector3i t_;
    Eigen::Vector3i dims;
    //proxy class for float**
    class y_z{
    private:
        float** yz_voxels;
        int dims_y;
        int dims_z;
    public:
        //proxy class for single float*
        class z{
        private:
            float* z_voxels;
            int dims_z;
        public:
            z(float* z_voxels_in, int dimz_in);

            //index operator for returning
            float& operator[](int index);
        }; //end of z proxy class

        y_z(float** yz_voxels_in, int dimy_in, int dimz_in);

        //index operator for returning proxy float*
        z operator[](int index);
    }; //end of y_z proxy class

    //index operator for returning proxy float**
    y_z operator[](int index);

    //construct unfilled grid with dimensions and translation
    grid(Eigen::Vector3i dims_in, Eigen::Vector3i t_in);
    //copy constructor
    grid(grid& other);

    //create grid from point cloud
    grid(pcl::PointCloud<pcl::InterestPoint>::Ptr grid_cloud, VoxelGridPtr vox);

    ~grid();

    //copy assignment
    grid& operator=(const grid& other);

    //linear index = z*num_x*num_y + y*num_x + x;
    //convert linear index to vector subscript
    Eigen::Vector3i ind2sub(int linear_index);
    //convert vector subscript to linear index
    int sub2ind(Eigen::Vector3i subs);
    //get value using linear indexing
    float& operator()(int index);

    grid operator+(const grid& rhs);

    void fillGrid(int res_factor);


};

#endif

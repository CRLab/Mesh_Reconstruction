#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <boost/thread/thread.hpp>
#include <pcl/filters/voxel_grid.h>

#include <stdlib.h>

#include "grid.h"

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;

grid::y_z::z::z(float* z_voxels_in, int dimz_in){
    z_voxels=z_voxels_in;
    dims_z=dimz_in;
}

float& grid::y_z::z::operator[](int index){
    if(index<0 || index>=dims_z){
        cerr<<"z index of grid out of range [0:"<<dims_z<<"): "<<index<<endl;
        exit(1);
    }
    return z_voxels[index];
}

grid::y_z::y_z(float** yz_voxels_in, int dimy_in, int dimz_in){
    yz_voxels=yz_voxels_in;
    dims_y = dimy_in;
    dims_z = dimz_in;
}

grid::y_z::z grid::y_z::operator[](int index){
    if(index<0 || index>=dims_y){
        cerr<<"y index of grid out of range [0:"<<dims_y<<"): "<<index<<endl;
        exit(1);
    }
    return z(yz_voxels[index], dims_z);
}

void grid::allocGrid(){
    voxels = new float**[dims[0]];
    for(int i=0; i<dims[0]; i++){
        voxels[i] = new float*[dims[1]];
        for(int j=0; j<dims[1]; j++){
            voxels[i][j] = new float[dims[2]];
        }
    }
}
void grid::deallocGrid(){
    for(int i=0; i<dims[0]; i++){
        for(int j=0; j<dims[1]; j++){
            delete [] voxels[i][j];
        }
        delete [] voxels[i];
    }
    delete [] voxels;
}

grid::grid(const Eigen::Vector3i &dims_in, const Eigen::Vector3f &scale_in, const Eigen::Vector3f &shift_in, const int pad_in){
    dims=dims_in;
    scale = scale_in;
    shift = shift_in;
    pad=pad_in;
    allocGrid();
}

grid::grid(grid& other){
    dims=other.dims;
    scale=other.scale;
    shift=other.shift;
    pad=other.pad;
    allocGrid();
    for(int i=0; i<dims[0]; i++){
        for(int j=0; j<dims[1]; j++){
            for(int k=0; k<dims[2]; k++){
                voxels[i][j][k]=other[i][j][k];
            }
        }
    }
}

//create grid from point cloud
grid::grid(pcl::PointCloud<pcl::InterestPoint>::Ptr grid_cloud, VoxelGridPtr vox){
    Eigen::Vector3i min_box = vox->getMinBoxCoordinates();
    Eigen::Vector3i num_divisions=vox->getNrDivisions();
    dims=num_divisions;
    Eigen::Vector3f high_vals;
    high_vals[0]=-10000.0; high_vals[1]=-10000.0; high_vals[2]=-10000.0;
    //pad the grid with empty voxels on each side
    pad= 6;
    shift[0]=10000.0; shift[1]=10000.0; shift[2]=10000.0;
    dims[0]+=pad*2; dims[1]+=pad*2; dims[2]+=pad*2;
    allocGrid();
    for(int i=0; i<dims[0]; i++){
        for(int j=0; j<dims[1]; j++){
            for(int k=0; k<dims[2]; k++){
                Eigen::Vector3i pnt;
                pnt[0]=i+min_box[0]-pad; pnt[1]=j+min_box[1]-pad; pnt[2]=k+min_box[2]-pad;
                int index = vox->getCentroidIndexAt(pnt);
                if(index!=-1){
                    if(grid_cloud->points[index].x<shift[0]) shift[0]=grid_cloud->points[index].x;
                    if(grid_cloud->points[index].y<shift[1]) shift[1]=grid_cloud->points[index].y;
                    if(grid_cloud->points[index].z<shift[2]) shift[2]=grid_cloud->points[index].z;
                    if(grid_cloud->points[index].x>high_vals[0]) high_vals[0]=grid_cloud->points[index].x;
                    if(grid_cloud->points[index].y>high_vals[1]) high_vals[1]=grid_cloud->points[index].y;
                    if(grid_cloud->points[index].z>high_vals[2]) high_vals[2]=grid_cloud->points[index].z;
                    voxels[i][j][k] = grid_cloud->points[index].strength;
                }
                else voxels[i][j][k]=-1.0;

                if(i<pad||i>=dims[0]-pad||j<pad||j>=dims[1]-pad||k<pad||k>=dims[2]-pad){
                    voxels[i][j][k]=-1.0;
                }
            }
        }
    }
    scale[0]=(high_vals[0]-shift[0])/((float)num_divisions[0]);
    scale[1]=(high_vals[1]-shift[1])/((float)num_divisions[1]);
    scale[2]=(high_vals[2]-shift[2])/((float)num_divisions[2]);
}

grid::~grid(){
    deallocGrid();
}

//copy assignment
grid& grid::operator=(const grid& other){
    if(this!=&other){
        this->pad=other.pad;
        if(this->dims[0]!=other.dims[0]||this->dims[1]!=other.dims[1]||this->dims[2]!=other.dims[2]){
            deallocGrid();
            this->dims = other.dims;
            allocGrid();
        }
        if(this->scale[0]!=other.scale[0]||this->scale[1]!=other.scale[1]||this->scale[2]!=other.scale[2]){
            this->scale = other.scale;
        }
        if(this->shift[0]!=other.shift[0]||this->shift[1]!=other.shift[1]||this->shift[2]!=other.shift[2]){
            this->shift = other.shift;
        }
        for(int i=0; i<dims[0]; i++){
            for(int j=0; j<dims[1]; j++){
                for(int k=0; k<dims[2]; k++){
                    voxels[i][j][k]=other.voxels[i][j][k];
                }
            }
        }
    }
    return *this;
}

grid::y_z grid::operator[](int index){
    if(index<0 || index>=dims[0]){
        cerr<<"x index of grid out of range [0:"<<dims[0]<<"): "<<index<<endl;
        exit(1);
    }
    return y_z(voxels[index], dims[1], dims[2]);
}

//linear index = z*num_x*num_y + y*num_x + x;
//convert linear index to vector subscript
Eigen::Vector3i grid::ind2sub(int linear_index){
    Eigen::Vector3i subs;
    subs[2] = linear_index/(dims[0]*dims[1]);
    linear_index = linear_index%(dims[0]*dims[1]);
    subs[1] = linear_index/dims[0];
    subs[0] = linear_index%dims[0];
    return subs;
}
//convert vector subscript to linear index
int grid::sub2ind(const Eigen::Vector3i &subs){
    return subs[2]*dims[0]*dims[1]+subs[1]*dims[0]+subs[0];
}
//get value using linear indexing
float& grid::operator()(int index){
    if(index<0 || index>=dims[0]*dims[1]*dims[2]){
        cerr<<"linear index out of range [0:"<<dims[0]*dims[1]*dims[2]<<"): "<<index<<endl;
        exit(1);
    }
    Eigen::Vector3i subs=ind2sub(index);
    return voxels[subs[0]][subs[1]][subs[2]];
}

grid grid::operator+(const grid& rhs){
    if(rhs.dims[0]!=dims[0] || rhs.dims[1]!=dims[1] || rhs.dims[2]!=dims[2]
            || rhs.scale[0]!=scale[0] || rhs.scale[1]!=scale[1] || rhs.scale[2]!=scale[2]
            || rhs.shift[0]!=shift[0] || rhs.shift[1]!=shift[1] || rhs.shift[2]!=shift[2]
            || rhs.pad!=pad){
        cerr<<"can't add grids with different dimensions";
    }
    grid out(dims, scale, shift, pad);
    for(int i=0; i<out.dims[0]; i++){
        for(int j=0; j<out.dims[1]; j++){
            for(int k=0; k<out.dims[2]; k++){
                out.voxels[i][j][k] = voxels[i][j][k]+rhs.voxels[i][j][k];
            }
        }
    }
    return out;
}

//get point in cloud corresponding to center of voxel pnt
Eigen::Vector3f grid::getCloudPoint(const Eigen::Vector3f &pnt){
    Eigen::Vector3f p;
    p[0]=((pnt[0]-(float)pad)*scale[0])+shift[0];
    p[1]=((pnt[1]-(float)pad)*scale[1])+shift[1];
    p[2]=((pnt[2]-(float)pad)*scale[2])+shift[2];
    return p;
}

//run this on the conf_grid to fill in any gaps between observed and completed clouds
void grid::fillGrid(int res_factor){
    int max_width = res_factor+1;
    for(int i=0; i<dims[0]; i++){
        for(int j=0; j<dims[1]; j++){
            int num_gap=0;
            bool first_hit=false;
            bool first_gap=false;
            bool end_hit=false;
            for(int k=0; k<dims[2]; k++){
                if(!end_hit){
                    //check if hit first occupied voxel
                    if(!first_hit){
                        if(voxels[i][j][k]!=-1.0){
                            first_hit=true;
                        }
                    }
                    //first occupied voxel has been hit
                    else{
                        //check if first gap has been hit
                        if(!first_gap){
                            if(voxels[i][j][k]==-1.0){
                                first_gap=true;
                                num_gap++;
                            }
                        }
                        //first gap has been hit
                        else{
                            //check if voxel is unoccupied
                            if(voxels[i][j][k]==-1.0){
                                num_gap++;
                            }
                            //next
                            else{
                                end_hit=true;
                                //check if within max width
                                if(num_gap<max_width){
                                    //fill in previous voxels
                                    for(int n=num_gap; n>0; n--){
                                        //compute confidence for new voxel from gaussian
                                        int dist = num_gap-n+1;
                                        float exp = -dist*dist/(4*res_factor);
                                        float conf = pow(2.71828f,exp);
                                        voxels[i][j][k-n]=conf;
                                    }
                                }
                            }
                        }
                    }
                }
            }//end of z column
        }
    }//end if iteration
}

//function for pcl visualization of grid
void grid::visualize(){
    //convert imbedding function to point cloud for visualization
    //create point clouds from bands for visualization
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pcl_grid (new pcl::PointCloud<pcl::PointXYZRGB>());
    for(int i=0; i<dims[0]; i++){
        for(int j=0; j<dims[1]; j++){
            for(int k=0; k<dims[2]; k++){
                if(voxels[i][j][k]>0){
                    pcl::PointXYZRGB pnt;
                    pnt.x=(float)i; pnt.y=(float)j; pnt.z=(float)k;
                    pnt.r=0;pnt.g=0;pnt.b=255;
                    pcl_grid->push_back(pnt);
                }
                else if(voxels[i][j][k]<0){
                    pcl::PointXYZRGB pnt;
                    pnt.x=(float)i; pnt.y=(float)j; pnt.z=(float)k;
                    pnt.r=0;pnt.g=0;pnt.b=0;
                    pcl_grid->push_back(pnt);
                }
                else{
                    pcl::PointXYZRGB pnt;
                    pnt.x=(float)i; pnt.y=(float)j; pnt.z=(float)k;
                    pnt.r=0;pnt.g=255;pnt.b=0;
                    pcl_grid->push_back(pnt);
                }
            }
        }
    }

    //visualize
    //display in visualizor
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0, 0, 0);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(pcl_grid);
    viewer->addPointCloud<pcl::PointXYZRGB> (pcl_grid, rgb, "cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "cloud");
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
}

//get create grid with confidences
gridPtr createGrid(pcl::PointCloud<pcl::InterestPoint>::Ptr grid_cloud, VoxelGridPtr vox, int res_factor){
    gridPtr gp(new grid(grid_cloud, vox));

    //fill in any holes between observed and completed clouds
    gp->fillGrid(res_factor);
    return gp;
}

//create binary volume grid from confidence grid
gridPtr getBinaryVolume(gridPtr grid_cloud){
    gridPtr g(new grid(grid_cloud->dims, grid_cloud->scale, grid_cloud->shift, grid_cloud->pad));
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                if((*grid_cloud)[i][j][k]>=0) (*g)[i][j][k]=1.0;
                else (*g)[i][j][k]=0.0;
            }
        }
    }
    return g;
}

//copy grid
gridPtr copyGrid(gridPtr in){
    gridPtr g(new grid(*in));
    return g;
}

//add 2 grids voxel by voxel
gridPtr addGrids(gridPtr in1, gridPtr in2){
    gridPtr g(new grid(in1->dims, in1->scale, in1->shift, in1->pad));
    *g = ((*in1)+(*in2));
    return g;
}

//get linear indices of all non-zero voxels in grid
vector<int> findIndexes(gridPtr band){
    vector<int> indexes_copy (band->dims[0]*band->dims[1]*band->dims[2], 0);
    int num_ind = 0;
    for(int i=0; i<band->dims[0]; i++){
        for(int j=0; j<band->dims[1]; j++){
            for(int k=0; k<band->dims[2]; k++){
                if((*band)[i][j][k]!=0.0){
                    Eigen::Vector3i pnt;
                    pnt[0]=i; pnt[1]=j; pnt[2]=k;
                    indexes_copy[num_ind] = band->sub2ind(pnt);
                    num_ind++;
                }
            }
        }
    }
    //resize indexes
    vector<int> indexes(num_ind, 0);
    for(int i=0; i<num_ind; i++){
        indexes[i] = indexes_copy[i];
    }
    return indexes;
}

//create index map
gridPtr getIndexMap(gridPtr band, const vector<int>& indexes){
    gridPtr map_(new grid(band->dims, band->scale, band->shift, band->pad));
    //set all values to -1
    for(int i=0; i<map_->dims[0]; i++){
        for(int j=0; j<map_->dims[1]; j++){
            for(int k=0; k<map_->dims[2]; k++){
                (*map_)[i][j][k]=-1.0;
            }
        }
    }
    //give i.d.'s to band location
    for(int i=0; i<indexes.size(); i++){
        (*map_)(indexes[i])=(float)i;
    }
    return map_;
}

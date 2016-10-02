
#include "dfields.h"

using namespace std;

#define GPU_CHECKERROR( err ) (gpuCheckError( err, __FILE__, __LINE__ ))
static void gpuCheckError( cudaError_t err,
                           const char *file,
                           int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}

#define getFlatIndex(x,y,z,xdim,ydim,zdim) z*ydim*xdim+y*xdim+x

#define get3DIndex(i,j,k,index,xdim,ydim,zdim) k=index/(xdim*ydim),index-=k*xdim*ydim,j=index/xdim,index-=j*xdim,i=index/1

float* flatten(gridPtr g) {
    float* g_flat = new float[g->dims[0] * g->dims[1] * g->dims[2]];
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                g_flat[getFlatIndex(i,j,k,g->dims[0],g->dims[1],g->dims[2])] = (*g)[i][j][k];
            }
        }
    }
    return g_flat;
}

gridPtr unflatten(float* g_flat, gridPtr volume_grid, int xdim, int ydim, int zdim) {
    gridPtr g(new grid(volume_grid->dims, volume_grid->scale, volume_grid->shift, volume_grid->pad));
    for(int i=0; i<xdim; i++) {
        for (int j = 0; j < ydim; j++) {
            for (int k = 0; k < zdim; k++) {
                (*g)[i][j][k] = g_flat[getFlatIndex(i,j,k,xdim,ydim,zdim)];
            }
        }
    }
    return g;
}

int numDiff(float *g, gridPtr real) {
    int numDiff = 0;
    for(int i=0; i<real->dims[0]; i++) {
        for (int j = 0; j < real->dims[1]; j++) {
            for (int k = 0; k < real->dims[2]; k++) {
                if((*real)[i][j][k] != g[getFlatIndex(i,j,k,real->dims[0],real->dims[1],real->dims[2])]) {
                    numDiff++;
                    cout << (*real)[i][j][k] << "," << g[getFlatIndex(i,j,k,real->dims[0],real->dims[1],real->dims[2])] << " | ";
                }
            }
        }
    }
    cout << endl;
    return numDiff;
}

__global__ void fastPerimGPU(float* src_grid, float* dest_grid, int xdim, int ydim, int zdim) {
    int center = blockIdx.x*blockDim.x+threadIdx.x;

    if(center < xdim * ydim * zdim) {
        dest_grid[center] = 1.0;
        int i,j,k,index=center;
        get3DIndex(i,j,k,index,xdim,ydim,zdim);

        int neg_i=i>0;
        int neg_j=j>0;
        int neg_k=k>0;
        int pos_i=i<xdim-1;
        int pos_j=j<ydim-1;
        int pos_k=k<zdim-1;

        float value = src_grid[center];
        for(int i_=i-neg_i; i_<=i+pos_i; i_++){
            for(int j_=j-neg_j; j_<=j+pos_j; j_++){
                for(int k_=k-neg_k; k_<=k+pos_k; k_++){
                    if(value!=src_grid[getFlatIndex(i_,j_,k_, xdim,ydim,zdim)]) {
                        dest_grid[center]=0.0;
                        return;
                    }
                }
            }
        }
    }
}

__global__ void sqdistpt1(float* dist_grid, float* dest_grid, int xdim, int ydim, int zdim) {
    int center = blockIdx.x*blockDim.x+threadIdx.x;

    if(center < xdim * ydim * zdim) {

        int i,j,k,var=center;
        get3DIndex(i,j,k,var,xdim,ydim,zdim);

        //get squared dist to closest 0 in row (i's)
        float min = 100000;
        for(int index=0; index<xdim; index++){
            float dist = 100000;
            if(dist_grid[getFlatIndex(index, j, k, xdim, ydim, zdim)]==0.0){
                dist = (index-i)*(index-i);
            }
            if(dist<min){
                min=dist;
            }
        }
        //set value of voxel to min q dist
        dest_grid[center]=(float)min;
    }
}

__global__ void sqdistpt2(float* dist_grid, float* dest_grid, int xdim, int ydim, int zdim) {
    int center = blockIdx.x*blockDim.x+threadIdx.x;

    if(center < xdim * ydim * zdim) {

        int i,j,k,var=center;
        get3DIndex(i,j,k,var,xdim,ydim,zdim);

        //get squared dist to closest 0 in row (i's)
        int min = 100000;
        for(int index=0; index<ydim; index++){
            int dist = dist_grid[getFlatIndex(i, index, k, xdim, ydim, zdim)] + (index - j) * (index - j);
            if(dist<min) min=dist;
        }
        //set value of voxel to min q dist
        dest_grid[center]=(float)min;
    }
}

__global__ void sqdistpt3(float* dist_grid, float* dest_grid, int xdim, int ydim, int zdim) {
    int center = blockIdx.x*blockDim.x+threadIdx.x;

    if(center < xdim * ydim * zdim) {

        int i,j,k,var=center;
        get3DIndex(i,j,k,var,xdim,ydim,zdim);

        //get squared dist to closest 0 in row (i's)
        int min = 100000;
        for(int index=0; index<zdim; index++){
            int dist = dist_grid[getFlatIndex(i, j, index, xdim, ydim, zdim)] + (index - k) * (index - k);
            if(dist<min) min=dist;
        }
        //set value of voxel to min q dist
        dest_grid[center]=(float)min;
    }
}

__global__ void sqrtgpu(float* dist_grid, float* dest_grid, int xdim, int ydim, int zdim) {
    int center = blockIdx.x*blockDim.x+threadIdx.x;

    if(center < xdim * ydim * zdim) {
        dest_grid[center] = sqrt(dist_grid[center]);
    }
}

gridPtr dfield_gpu(gridPtr volume_grid) {
    float* flat_g = flatten(volume_grid);
    int xdim = volume_grid->dims[0], ydim = volume_grid->dims[1], zdim = volume_grid->dims[2];
    int g_size = xdim * ydim * zdim;

    float* dest_grid, *src_grid;
    cout << "Initializing dfields arrays" << endl;
    GPU_CHECKERROR(cudaMalloc((void**) &dest_grid, g_size * sizeof(float)));
    GPU_CHECKERROR(cudaMalloc((void**) &src_grid, g_size * sizeof(float)));

    GPU_CHECKERROR(cudaMemcpy(src_grid, flat_g, g_size * sizeof(float), cudaMemcpyHostToDevice));

    int BLOCK_SIZE = 256;
    int NUM_BLOCKS = (int)ceil(double(g_size)/BLOCK_SIZE);

    //Alternate buffers so we don't have to move data around on the device
    //Buffers effectively alternate between input and output
    cout << "dfield Fast perim gpu" << endl;
    fastPerimGPU<<<NUM_BLOCKS, BLOCK_SIZE>>>(src_grid, dest_grid, xdim, ydim, zdim);

    cout << "dfield sqdist" << endl;
    sqdistpt1<<<NUM_BLOCKS, BLOCK_SIZE>>>(dest_grid, src_grid, xdim, ydim, zdim);
    sqdistpt2<<<NUM_BLOCKS, BLOCK_SIZE>>>(src_grid, dest_grid, xdim, ydim, zdim);
    sqdistpt3<<<NUM_BLOCKS, BLOCK_SIZE>>>(dest_grid, src_grid, xdim, ydim, zdim);

    cout << "dfield sqrt" << endl;
    sqrtgpu<<<NUM_BLOCKS, BLOCK_SIZE>>>(src_grid, dest_grid, xdim, ydim, zdim);

    //Copy data back to device
    GPU_CHECKERROR(cudaMemcpy(flat_g, dest_grid, g_size * sizeof(float), cudaMemcpyDeviceToHost));

    cudaThreadSynchronize();

    cudaFree(dest_grid);
    cudaFree(src_grid);

    cudaThreadSynchronize();

    return unflatten(flat_g, volume_grid, xdim, ydim, zdim);
}

//extract perimeter of binary volume
gridPtr fastPerim(gridPtr volume_grid){
    gridPtr g(new grid(volume_grid->dims, volume_grid->scale, volume_grid->shift, volume_grid->pad));
    for(int i=0; i<g->dims[0]; i++){
        for(int j=0; j<g->dims[1]; j++){
            for(int k=0; k<g->dims[2]; k++){
                (*g)[i][j][k]=1.0;
                int neg_i=1;int neg_j=1;int neg_k=1;
                int pos_i=1; int pos_j=1; int pos_k=1;
                if(i==0)neg_i=0;
                else if(i==g->dims[0]-1)pos_i=0;
                if(j==0)neg_j=0;
                else if(j==g->dims[1]-1)pos_j=0;
                if(k==0)neg_k=0;
                else if(k==g->dims[2]-1)pos_k=0;
                //search neighboring voxels for different value
                for(int i_=i-neg_i; i_<=i+pos_i; i_++){
                    for(int j_=j-neg_j; j_<=j+pos_j; j_++){
                        for(int k_=k-neg_k; k_<=k+pos_k; k_++){
                            if((*volume_grid)[i][j][k]!=(*volume_grid)[i_][j_][k_]){
                                (*g)[i][j][k]=0.0;
                            }
                        }
                    }
                }
                //end of neighbor search
            }
        }
    }

    return g;
}

//implementation of Saito and Toriwaki distance field
//squared euclidean distance to closest 0-value voxel
gridPtr getsqdist(gridPtr volume_grid){
    //make copies
    gridPtr dist_grid = copyGrid(volume_grid);
    gridPtr dist_grid_copy = copyGrid(dist_grid);
    //first transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                //get squared dist to closest 0 in row (i's)
                float min = 100000;
                for(int index=0; index<dist_grid->dims[0]; index++){
                    float dist = 100000;
                    if((*dist_grid)[index][j][k]==0.0){
                        dist = (index-i)*(index-i);
                    }
                    if(dist<min){
                        min=dist;
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=min;
            }
        }
    }//end of first transformation

    dist_grid=copyGrid(dist_grid_copy);
    //second transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                for(int index=0; index<dist_grid->dims[1]; index++){
                    int dist = (*dist_grid)[i][index][k] + (index-j)*(index-j);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
            }
        }
    }//end of second transformation

    dist_grid=copyGrid(dist_grid_copy);
    //third transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                for(int index=0; index<dist_grid->dims[2]; index++){
                    int dist = (*dist_grid)[i][j][index] + (index-k)*(index-k);
                    if(dist<min) min=dist;
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
            }
        }
    }//end of third transformation

    return dist_grid_copy;
}

//returns linear indexes of closest 0-value voxel
gridPtr getsqdist_index(gridPtr volume_grid){
    //make copies
    gridPtr index_grid;
    gridPtr index_grid_copy = copyGrid(volume_grid);
    gridPtr dist_grid = copyGrid(volume_grid);
    gridPtr dist_grid_copy = copyGrid(dist_grid);
    //first transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                //get squared dist to closest 0 in row (i's)
                float min = 100000;
                float x_ind = -1.0;
                for(int index=0; index<dist_grid->dims[0]; index++){
                    float dist = 100000;
                    if((*dist_grid)[index][j][k]==0.0){
                        dist = (index-i)*(index-i);
                    }
                    if(dist<min){
                        min=dist;
                        x_ind = index;
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=min;
                //set index value to x index
                (*index_grid_copy)[i][j][k]=x_ind;
            }
        }
    }//end of first transformation

    index_grid=copyGrid(index_grid_copy);
    dist_grid=copyGrid(dist_grid_copy);
    //second transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                float xy_ind = -1.0;
                for(int index=0; index<dist_grid->dims[1]; index++){
                    int dist = (*dist_grid)[i][index][k] + (index-j)*(index-j);
                    if(dist<min){
                        min=dist;
                        xy_ind = (float)(index*index_grid->dims[0])+(*index_grid)[i][index][k];
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
                (*index_grid_copy)[i][j][k] = xy_ind;
            }
        }
    }//end of second transformation

    index_grid=copyGrid(index_grid_copy);
    dist_grid=copyGrid(dist_grid_copy);
    //third transformation
    for(int i=0; i<dist_grid->dims[0]; i++){
        for(int j=0; j<dist_grid->dims[1]; j++){
            for(int k=0; k<dist_grid->dims[2]; k++){
                int min = 100000;
                float xyz_ind=-1.0;
                for(int index=0; index<dist_grid->dims[2]; index++){
                    int dist = (*dist_grid)[i][j][index] + (index-k)*(index-k);
                    if(dist<min){
                        min=dist;
                        xyz_ind = (float)(index*index_grid->dims[0]*index_grid->dims[1])+(*index_grid)[i][j][index];
                    }
                }
                //set value of voxel to min q dist
                (*dist_grid_copy)[i][j][k]=(float)min;
                (*index_grid_copy)[i][j][k] = xyz_ind;
            }
        }
    }//end of third transformation

    return index_grid_copy;

}

//get square root of grid
gridPtr getsqrt(gridPtr g){
    gridPtr s(new grid(g->dims, g->scale, g->shift, g->pad));
    for(int i=0; i<s->dims[0]; i++){
        for(int j=0; j<s->dims[1]; j++){
            for(int k=0; k<s->dims[2]; k++){
                (*s)[i][j][k]= (float) sqrt((double)(*g)[i][j][k]);
            }
        }
    }
    return s;
}


//*****************************************************************************************************
//functions for computing normals of surface

//get linear indices of neighboring voxels
vector<int> getNeighbors(gridPtr g, const Eigen::Vector3i &pnt){
    vector<int> indexes;
    for(int i=pnt[0]-1; i<=pnt[0]+1; i++){
        for(int j=pnt[1]-1; j<=pnt[1]+1; j++){
            for(int k=pnt[2]-1; k<=pnt[2]+1; k++){
                if((*g)[i][j][k]>0){
                    Eigen::Vector3i p;
                    p[0]=i; p[1]=j; p[2]=k;
                    indexes.push_back(g->sub2ind(p));
                }
            }
        }
    }
    return indexes;
}

//check if voxel is on surface
bool isSurface(gridPtr g, const Eigen::Vector3i &pnt){
    //search neighboring voxels for different value
    int i=pnt[0];
    int j=pnt[1];
    int k=pnt[2];
    if((*g)[i][j][k]!=1.0) return false;
    for(int i_=i-1; i_<=i+1; i_++){
        for(int j_=j-1; j_<=j+1; j_++){
            for(int k_=k-1; k_<=k+1; k_++){
                if((*g)[i_][j_][k_]==0.0){
                    return true;
                }
            }
        }
    }
    return false;
}

//get linear indices of neighboring voxels on the surface
vector<int> getSurfaceNeighbors(gridPtr g, const Eigen::Vector3i &pnt){
    vector<int> indexes;
    for(int i=pnt[0]-1; i<=pnt[0]+1; i++){
        for(int j=pnt[1]-1; j<=pnt[1]+1; j++){
            for(int k=pnt[2]-1; k<=pnt[2]+1; k++){
                if((*g)[i][j][k]>0){
                    Eigen::Vector3i p;
                    p[0]=i; p[1]=j; p[2]=k;
                    if(isSurface(g,p)){
                        indexes.push_back(g->sub2ind(p));
                    }
                }
            }
        }
    }
    return indexes;
}

//get centroid of a list of points
Eigen::Vector3f getCentroid(gridPtr g, const vector<int> &pnts){
    Eigen::Vector3f centroid;
    centroid[0]=0.0; centroid[1]=0.0; centroid[2]=0.0;
    for(int i=0; i<pnts.size(); i++){
        Eigen::Vector3i p = g->ind2sub(pnts[i]);
        centroid[0]+=(float)p[0]; centroid[1]+=(float)p[1]; centroid[2]+=(float)p[2];
    }
    centroid[0]=centroid[0]/(float)pnts.size();
    centroid[1]=centroid[1]/(float)pnts.size();
    centroid[2]=centroid[2]/(float)pnts.size();
    return centroid;
}

//get covariance matrix of surface point
Eigen::Matrix3f getCovariance(gridPtr g, const Eigen::Vector3i &pnt){
    vector<int> indexes = getSurfaceNeighbors(g, pnt);
    Eigen::Vector3f centroid = getCentroid(g, indexes);
    Eigen::Matrix3f covariance;
    covariance<<0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    for(int i=0; i<indexes.size(); i++){
        Eigen::Vector3f vec;
        vec[0]=(float)(g->ind2sub(indexes[i]))[0] - centroid[0];
        vec[1]=(float)(g->ind2sub(indexes[i]))[1] - centroid[1];
        vec[2]=(float)(g->ind2sub(indexes[i]))[2] - centroid[2];
        for(int row=0; row<3; row++){
            for(int col=0; col<3; col++){
                covariance(row,col) = covariance(row,col)+vec[row]*vec[col];
            }
        }
    }
    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            covariance(row,col) = covariance(row,col)/(float)indexes.size();
        }
    }
    return covariance;
}

//get normal vector from covariance matrix
Eigen::Vector3f getNormalVector(const Eigen::Matrix3f &covariance){
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es;
    es.compute(covariance);
    float min=10000;
    int index=0;
    for(int i=0; i<3; i++){
        int eigenvalue=es.eigenvalues()[i];
        if(eigenvalue<min){
            min=eigenvalue;
            index=i;
        }
    }
    Eigen::Vector3f vec = es.eigenvectors().col(index);
    float len = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
    vec[0]=vec[0]/len; vec[1]=vec[1]/len; vec[2]=vec[2]/len;
    return vec;
}

//orient normal vector with respect to local centroid
void orientNormal(gridPtr g, Eigen::Vector3f &normal, const Eigen::Vector3i &pnt){
    Eigen::Vector3f centroid = getCentroid(g, getNeighbors(g, pnt));
    centroid[0]-=(float)pnt[0]; centroid[1]-=(float)pnt[1]; centroid[2]-=(float)pnt[2];
    float prod = normal[0]*centroid[0]+normal[1]*centroid[1]+normal[2]*centroid[2];
    if(prod>0.0){
        normal[0]=-normal[0]; normal[1]=-normal[1]; normal[2]=-normal[2];
    }
}

//get indexes of surface points of volume grid
vector<int> getSurface(gridPtr volume){
    vector<int> indexes;
    gridPtr surface(new grid(volume->dims, volume->scale, volume->shift, volume->pad));
    for(int i=0; i<surface->dims[0]; i++){
        for(int j=0; j<surface->dims[1]; j++){
            for(int k=0; k<surface->dims[2]; k++){
                Eigen::Vector3i pnt;
                pnt[0]=i; pnt[1]=j; pnt[2]=k;
                if(isSurface(volume, pnt)){
                    indexes.push_back(volume->sub2ind(pnt));
                }
            }
        }
    }
    return indexes;
}

//get surface normals
vector<Eigen::Vector3f> getSurfaceNormals(gridPtr volume, const vector<int> &surface){
    vector<Eigen::Vector3f> normals;
    for(int i=0; i<surface.size(); i++){
        Eigen::Matrix3f covariance = getCovariance(volume, volume->ind2sub(surface[i]));
        Eigen::Vector3f norm = getNormalVector(covariance);
        normals.push_back(norm);
        orientNormal(volume, normals[i], volume->ind2sub(surface[i]));
    }
    return normals;
}

//returns index in vector of value, -1 if not contained
int contains(vector<int> vec, int val){
    for(int i=0; i<vec.size(); i++){
        if(vec[i]==val) return i;
    }
    return -1;
}


//create normals and visualize in pcl viewer
void visualizeNormals(gridPtr volume){
    vector<int> surface = getSurface(volume);
    gridPtr surfaceMap = getIndexMap(volume, surface);
    vector<Eigen::Vector3f> normals = getSurfaceNormals(volume, surface);

    //create point cloud from volume
    //convert imbedding function to point cloud for visualization
    //create point clouds from bands for visualization
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr pcl_grid (new pcl::PointCloud<pcl::PointXYZRGB>());
    for(int i=0; i<volume->dims[0]; i++){
        for(int j=0; j<volume->dims[1]; j++){
            for(int k=0; k<volume->dims[2]; k++){
                if((*volume)[i][j][k]>0){
                    Eigen::Vector3i pnt;
                    pnt[0]=i; pnt[1]=j; pnt[2]=k;
                    int index = volume->sub2ind(pnt);
                    int ind = contains(surface, index);
                    if(ind==-1){
                        pcl::PointXYZRGB p;
                        p.x=(float)i; p.y=(float)j; p.z=(float)k;
                        p.r=0;p.g=0;p.b=255;
                        pcl_grid->push_back(p);
                    }
                    else{
                        pcl::PointXYZRGB p;
                        p.x=(float)i; p.y=(float)j; p.z=(float)k;
                        p.r=255;p.g=0;p.b=0;
                        pcl_grid->push_back(p);
                    }
                }
                else{
                    pcl::PointXYZRGB p;
                    p.x=(float)i; p.y=(float)j; p.z=(float)k;
                    p.r=0;p.g=0;p.b=0;
                    pcl_grid->push_back(p);
                }
            }
        }
    }

    //create point cloud of normals
    pcl::PointCloud<pcl::Normal>::Ptr normal_grid (new pcl::PointCloud<pcl::Normal>());
    for(int i=0; i<volume->dims[0]; i++){
        for(int j=0; j<volume->dims[1]; j++){
            for(int k=0; k<volume->dims[2]; k++){
                int ind = (*surfaceMap)[i][j][k];
                if(ind==-1){
                    pcl::Normal norm;
                    norm.normal_x = norm.normal_y = norm.normal_z = norm.data_n[3] = 0.0f;
                    norm.curvature = 0;
                    normal_grid->push_back(norm);
                }
                else{
                    pcl::Normal norm;
                    norm.normal_x = normals[ind][0]; norm.normal_y = normals[ind][1]; norm.normal_z = normals[ind][2];
                    norm.curvature = 0;
                    norm.data_n[3]=0.0f;
                    normal_grid->push_back(norm);
                }
            }
        }
    }

    //visualize
    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (0, 0, 0);
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(pcl_grid);
    viewer->addPointCloud<pcl::PointXYZRGB> (pcl_grid, rgb, "sample cloud");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "sample cloud");
    viewer->addPointCloudNormals<pcl::PointXYZRGB, pcl::Normal> (pcl_grid, normal_grid, 1, 0.25, "normals");
    viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();
    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }
}

//********************************************************************************************************************
//feature detection using normals

//input is binary volume pre-smoothing
//binary grid: 1=feature, 0=no feature
//features can be ignored during smoothing
//reasonable threshold ~0.9
gridPtr getFeatureMap(gridPtr volume, gridPtr surfaceMap, const vector<Eigen::Vector3f> &normals, float threshold){
    gridPtr featureMap (new grid(volume->dims, volume->scale, volume->shift, volume->pad));
    for(int i=0; i<featureMap->dims[0]; i++){
        for(int j=0; j<featureMap->dims[1]; j++){
            for(int k=0; k<featureMap->dims[2]; k++){
                (*featureMap)[i][j][k]=0.0;
            }
        }
    }

    int edgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

    for(int i=0; i<volume->dims[0]-1; i++){
        for(int j=0; j<volume->dims[1]-1; j++){
            for(int k=0; k<volume->dims[2]-1; k++){
                //create gridcell
                int gridcell[8];
                for(int n=0; n<8; n++){
                    int i_ = ((n%4)==1 || (n%4)==2);
                    int j_ = (n%4)/2;
                    int k_ = n/4;
                    gridcell[n] = (int)(*volume)[i+i_][j+j_][k+k_];
                }
                int cubeindex=0;
                if (gridcell[0] < 1) cubeindex |= 1;
                if (gridcell[1] < 1) cubeindex |= 2;
                if (gridcell[2] < 1) cubeindex |= 4;
                if (gridcell[3] < 1) cubeindex |= 8;
                if (gridcell[4] < 1) cubeindex |= 16;
                if (gridcell[5] < 1) cubeindex |= 32;
                if (gridcell[6] < 1) cubeindex |= 64;
                if (gridcell[7] < 1) cubeindex |= 128;

                //check if cube is on the edge
                if(edgeTable[cubeindex]==0) continue;

                //check for feature in cube
                //get opening angles of normals
                float min=5.0;
                for(int n=0; n<7; n++){
                    for(int m=n+1; m<8; m++){
                        int i_1 = ((n%4)==1 || (n%4)==2);
                        int j_1 = (n%4)/2;
                        int k_1 = n/4;
                        int i_2 = ((m%4)==1 || (m%4)==2);
                        int j_2 = (m%4)/2;
                        int k_2 = m/4;
                        int index1 = (*surfaceMap)[i+i_1][j+j_1][k+k_1];
                        int index2 = (*surfaceMap)[i+i_2][j+j_2][k+k_2];
                        if(index1==-1 || index2==-1) continue;
                        Eigen::Vector3f norm1 = normals[index1];
                        Eigen::Vector3f norm2 = normals[index2];
                        float angle = norm1[0]*norm2[0]+norm1[1]*norm2[1]+norm1[2]*norm2[2];
                        if(angle<min) min=angle;
                    }
                }

                //check if within threshold for feature detection
                if(min>=threshold) continue;

                //set feature points to 1.0
                for(int i_=0; i_<2; i_++){
                    for(int j_=0; j_<2; j_++){
                        for(int k_=0; k_<2; k_++){
                            if((*surfaceMap)[i+i_][j+j_][k+k_]>-1.0){
                                (*featureMap)[i+i_][j+j_][k+k_]=1.0;
                            }
                        }
                    }
                }

            }
        }
    }//end of iteration through grid

    return featureMap;
}



#include "binvoxToPcl.h"

using namespace std;

typedef unsigned char byte;


int get_index(binvox vox, int x, int y, int z){
    int wxh = vox.width*vox.height;
    int index = x*wxh + z*vox.width + y;
    return index;
}

int read_binvox(binvox* vox, string filespec)
{

  ifstream *input = new ifstream(filespec.c_str(), ios::in | ios::binary);

  //
  // read header
  //
  string line;
  *input >> line;  // #binvox
  if (line.compare("#binvox") != 0) {
    cout << "Error: first line reads [" << line << "] instead of [#binvox]" << endl;
    delete input;
    return 0;
  }
  *input >> vox->version;
  cout << "reading binvox version " << vox->version << endl;

  vox->depth = -1;
  int done = 0;
  while(input->good() && !done) {
    *input >> line;
    if (line.compare("data") == 0) done = 1;
    else if (line.compare("dim") == 0) {
      *input >> vox->depth >> vox->height >> vox->width;
    }
    else if (line.compare("translate") == 0) {
      *input >> vox->tx >> vox->ty >> vox->tz;
    }
    else if (line.compare("scale") == 0) {
      *input >> vox->scale;
    }
    else {
      cout << "  unrecognized keyword [" << line << "], skipping" << endl;
      char c;
      do {  // skip until end of line
        c = input->get();
      } while(input->good() && (c != '\n'));

    }
  }
  if (!done) {
    cout << "  error reading header" << endl;
    return 0;
  }
  if (vox->depth == -1) {
    cout << "  missing dimensions in header" << endl;
    return 0;
  }

  vox->size = vox->width * vox->height * vox->depth;
  vox->voxels = new byte[vox->size];
  if (!vox->voxels) {
    cout << "  error allocating memory" << endl;
    return 0;
  }

  //
  // read voxel data
  //
  byte value;
  byte count;
  int index = 0;
  int end_index = 0;
  int nr_voxels = 0;

  input->unsetf(ios::skipws);  // need to read every byte now (!)
  *input >> value;  // read the linefeed char

  while((end_index < vox->size) && input->good()) {
    *input >> value >> count;

    if (input->good()) {
      end_index = index + count;
      if (end_index > vox->size) return 0;
      for(int i=index; i < end_index; i++) vox->voxels[i] = value;

      if (value) nr_voxels += count;
      index = end_index;
    }  // if file still ok

  }  // while

  input->close();
  cout << "  read " << nr_voxels << " voxels" << endl;

  return 1;

}

pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL(string filespec, int res_factor){
    if(res_factor<1){
        cout<<"Resolution factor must be at least 1"<<endl<<endl;
        exit(1);
    }
    binvox vox;
    if (!read_binvox(&vox, filespec)){
        cout << "Error reading [" << filespec << "]" <<endl <<endl;
        exit(1);
    }

    //convert voxels to point cloud
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ> ());
    for(int x=0; x<vox.depth; x++){
        for(int y=0; y<vox.width; y++){
            for(int z=0; z<vox.height; z++){
                int index = get_index(vox, x,y,z);
                if(vox.voxels[index]){
                    //check if voxel is on boundary of grid or object
                    bool x_low = (x==0); bool y_low = (y==0); bool z_low = (z==0);
                    bool x_high = (x==vox.depth-1); bool y_high = (y==vox.width-1); bool z_high = (z==vox.height-1);
                    if(!x_low) x_low = !vox.voxels[get_index(vox, x-1,y,z)];
                    if(!y_low) y_low = !vox.voxels[get_index(vox, x,y-1,z)];
                    if(!z_low) z_low = !vox.voxels[get_index(vox, x,y,z-1)];
                    if(!x_high) x_high = !vox.voxels[get_index(vox, x+1,y,z)];
                    if(!y_high) y_high = !vox.voxels[get_index(vox, x,y+1,z)];
                    if(!z_high) z_high = !vox.voxels[get_index(vox, x,y,z+1)];

                    //get bounds for edges
                    int x_lb = 0; int y_lb = 0; int z_lb = 0;
                    int x_ub = res_factor; int y_ub = res_factor; int z_ub = res_factor;

                    if(x_low) x_lb=res_factor/2; if(y_low) y_lb=res_factor/2; if(z_low) z_lb=res_factor/2;
                    if(x_high) x_ub=(res_factor+1)/2; if(y_high) y_ub=(res_factor+1)/2; if(z_high) z_ub=(res_factor+1)/2;

                    //create points evenly distributed inside of voxel
                    pcl::PointXYZ pnt[(x_ub-x_lb)*(y_ub-y_lb)*(z_ub-z_lb)];
                    for(int i=x_lb; i<x_ub; i++){
                        for(int j=y_lb; j<y_ub; j++){
                            for(int k=z_lb; k<z_ub; k++){
                                int ind=(i-x_lb)*(y_ub-y_lb)*(z_ub-z_lb)+(j-y_lb)*(z_ub-z_lb)+(k-z_lb);
                                pnt[ind].x = ((float)x+(0.5/(float)res_factor)+(float)i*(1/(float)res_factor))*vox.scale/((float)vox.depth)+vox.tx;
                                pnt[ind].y = ((float)y+(0.5/(float)res_factor)+(float)j*(1/(float)res_factor))*vox.scale/((float)vox.width)+vox.ty;
                                pnt[ind].z = ((float)z+(0.5/(float)res_factor)+(float)k*(1/(float)res_factor))*vox.scale/((float)vox.height)+vox.tz;
                                cloud->push_back(pnt[ind]);
                            }
                        }
                    }//finished generate points for voxel x,y,z

                }
            }
        }
    }//finished generating point cloud

    return cloud;
}

//function returns density by finding average distance between points in cloud
float getDensity(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud){
    //create new cloud missing every ten points
    pcl::PointCloud<pcl::PointXYZ>::Ptr new_cloud (new pcl::PointCloud<pcl::PointXYZ>());
    for(int i=0; i<cloud->points.size(); i++){
        if((i%10)>0){
            new_cloud->push_back(cloud->points[i]);
        }
    }


    //create kd tree
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(new_cloud);

    float distTotal = 0.0;
    int count = 0;
    for(int i=0; i<cloud->points.size(); i+=10){
        //perform nearest neighbor search
        vector<int> pointIdxNKNSearch(1);
        vector<float> pointNKNSquaredDistance(1);
        kdtree.nearestKSearch(cloud->points[i], 1, pointIdxNKNSearch, pointNKNSquaredDistance);
        //get distance
        distTotal += (float) sqrt((double)pointNKNSquaredDistance[0]);
        count++;
    }

    return (float)count/distTotal;
}

int getResolutionFactor(string filespec, pcl::PointCloud<pcl::PointXYZ>::Ptr otherCloud){
    float otherDensity = 3.0*getDensity(otherCloud)/4.0;
    pcl::PointCloud<pcl::PointXYZ>::Ptr initCloud = binvoxToPCL(filespec);
    float density = getDensity(initCloud);
    cout<<"Observed density: "<<otherDensity<<endl;
    cout<<"Completed density: "<<density<<endl;
    if(density>otherDensity){
        cout<<"Resolution Factor: 1"<<endl;
        return 1;
    }
    else{
        int res_factor = (int)ceil(otherDensity/density);
        cout<<"Resolution Factor: "<<res_factor<<endl;
        return res_factor;
    }
}

//automatically adjust point cloud resolutions based on density ratio:
pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL_autores(string filespec, pcl::PointCloud<pcl::PointXYZ>::Ptr otherCloud){
    int res_factor = getResolutionFactor(filespec, otherCloud);
    return binvoxToPCL(filespec, res_factor);
}


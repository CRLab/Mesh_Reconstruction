
#include "binvoxToPcl.h"

using namespace std;

typedef unsigned char byte;


int get_index(const binvox &vox, int x, int y, int z){
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
    for(int x=0; x<vox.depth-1; x++){
        for(int y=0; y<vox.width-1; y++){
            for(int z=0; z<vox.height-1; z++){
                //get indexes of cube
                vector<int> indexes;
                for(int i=0; i<=1; i++){
                    for(int j=0; j<=1; j++){
                        for(int k=0; k<=1; k++){
                            indexes.push_back(get_index(vox, x+i, y+j, z+k));
                        }
                    }
                }

                //create points evenly distributed inside of cube
                for(int i=0; i<res_factor; i++){
                    for(int j=0; j<res_factor; j++){
                        for(int k=0; k<res_factor; k++){
                            Eigen::Vector3f pnt;
                            pnt[0] = (float)x+((float)i)/((float)res_factor);
                            pnt[1] = (float)y+((float)j)/((float)res_factor);
                            pnt[2] = (float)z+((float)k)/((float)res_factor);
                            //check if pnt is closer to surface then outside
                            float val=0.0;
                            for(int n=0; n<8; n++){
                                //get distance to point (l1 norm)
                                int a=n/4;
                                int b=(n%4)/2;
                                int c=(n%4)%2;
                                float dist = fabs(pnt[0]-(float)(x+a));
                                dist+= fabs(pnt[1]-(float)(y+b));
                                dist+= fabs(pnt[2]-(float)(z+c));
                                dist+= 0.00001;
                                float weight = (float)(2*vox.voxels[indexes[n]]-1);
                                val+=(weight/dist);
                            }
                            if(val>0.0){
                                pcl::PointXYZ p;
                                p.x = (pnt[0]+0.5)*vox.scale/((float)vox.depth)+vox.tx;
                                p.y = (pnt[1]+0.5)*vox.scale/((float)vox.width)+vox.ty;
                                p.z = (pnt[2]+0.5)*vox.scale/((float)vox.height)+vox.tz;
                                cloud->push_back(p);
                            }
                        }
                    }
                }//finished generate points for cube x,y,z

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


#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/surface/mls.h>

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

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

pcl::PointCloud<pcl::PointXYZ>::Ptr binvoxToPCL(string filespec){
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
                    pcl::PointXYZ pnt;
                    pnt.x = ((float)x + 0.5)*vox.scale/((float)vox.depth) + vox.tx;
                    pnt.y = ((float)y + 0.5)*vox.scale/((float)vox.width) + vox.ty;
                    pnt.z = ((float)z + 0.5)*vox.scale/((float)vox.height) + vox.tz;
                    cloud->push_back(pnt);
                }
            }
        }
    }

    return cloud;

}

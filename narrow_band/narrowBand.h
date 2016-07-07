#ifndef NARROWBAND_H
#define NARROWBAND_H

#include "dfields.h"
#include "grid.h"

using namespace std;

typedef boost::shared_ptr< pcl::VoxelGrid<pcl::InterestPoint> > VoxelGridPtr;

struct bands{
    gridPtr band;
    gridPtr tight_band;
};
typedef boost::shared_ptr<bands> bandsPtr;

//morphological erosion with mask: [[000;010;000],[010,111,010],[000;010;000]]
gridPtr erode_grid(gridPtr gr);

//morphological dilation with mask: [[000;010;000],[010,111,010],[000;010;000]]
gridPtr dilate_grid(gridPtr gr);

//generate band and tight band (eroded band) using dist field "margin" and band_size
bandsPtr createBands(gridPtr margin, float band_size);

#endif

# Mesh_Reconstruction
Code for combining Partial Object View with a CNN generated Completion

##Usage
```
Usage: mesh_reconstruction <binvox file> <pcd file> <output file> <options>
Allowed options:
  --help                  produce help message
  --feature-detection     Toggle feature handling
  --cuda                  Toggle CUDA option
  --feature-threshold arg Increasing raises feature sensitivity. Default: 0.75
  --corner-threshold arg  Decreasing raises feature sensitivity. Default: 0.8
```

<binvox file> is the output of running our CNN, normally a 40x40x40 occupancy map representing our hypothesis about which voxels are occupied by the object.
<pcd file> is the captured partial view of the target object, this only contains the points visible to the camera and tends to be higher resolution than the completion we are combining it with
<output file> where to place to merged partial view and completion

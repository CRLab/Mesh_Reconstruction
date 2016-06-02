This archive contains the MATLAB implementation of the algorithm, performing surface extraction from binary volumes, as described in:

V. Lempitsky, "Surface Extraction from Binary Volumes with Higher-Order Smoothness", MSR-TR-2009-31 (also in submission).

LICENSE

See license.txt for the terms of use.

ACKNOWLEDGEMENT

If you use the code for the publication, please acknowledge it by citing the above-mentioned technical report or the subsequent publication.

CODE

The algorithm is contained in "ExtractSeparatingSurfaceCPU.m" (CPU verison) and "ExtractSeparatingSurfaceGPU.m" (GPU verison, requires CUDA-enabled NVIDIA card). 
Calling the functions without arguments will test the algorithm on a binary volume corresponding to a ball. See the files for further comments.

"FastPerim.cpp" and "RunQP.cpp" are the auxiliary C++ code pieces, the latter containing the numeric scheme for solving the quadratic program. 
Mexifications for Win32 are provided.

"QPGPU.cu" is the CUDA code for the same scheme to be run on an NVIDIA GPU (the mexification for Win32 is provided).

Due to memory inefficiency, the MATLAB version will fail on isosurface extraction from large volumes 
(an equivalent C++ version of the code was used to obtain respective images for the publication).

"bwsqdist" - an optional missing file from the distribution that computes the euclidean distance
transform fast and would accelerate the code considerably for large volumes. This is available on request.


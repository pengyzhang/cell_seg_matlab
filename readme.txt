******************************************************
Readme for sparse shape prior
10/26/2015
******************************************************
To run the program you must include seed files as in '.\data\seed_detection_result\' ,shape prior data '.\data\prior\' and the sparse-solving toolbox '.\l1_ls_matlab\'.

Functions for updating LSF:
LSbatch.m: entrance function. 
lse.m: update the lsf in three steps: updateSR; updatef; updateLSF.
sparse_solver.m: solving sparse representation.
distance_map.m: generate distance maps from list of contour coordinates.

Functions for performance evaluation:
contour_generate.m(previously named as temp.m): eliminate empty cells in boundarycoordinates.mat; generate figures from boundarycoordinates.mat.
genImageIphotoDrawXml.m: generate xml files.
CVPR_Performance.m: compute metrics.
******************************************************
# Sparse shape prior based nuclei segmentation
Paper citation:

1. Pengyue Zhang, Fusheng Wang, George Teodoro, Yanhui Liang, Daniel J. Brat, Jun Kong:
Automated level set segmentation of histopathologic cells with sparse shape prior support and dynamic occlusion constraint. ISBI 2017: 718-722

2. Jun Kong, Pengyue Zhang, Yanhui Liang, George Teodoro, Daniel J. Brat, Fusheng Wang:
Robust cell segmentation for histological images of Glioblastoma. ISBI 2016: 1041-1045
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

evaluate_performance.m: compute metrics.

******************************************************
# Sparse shape prior based nuclei segmentation
******************************************************
Paper citation:

1. Pengyue Zhang, Fusheng Wang, George Teodoro, Yanhui Liang, Daniel J. Brat, Jun Kong:
Automated level set segmentation of histopathologic cells with sparse shape prior support and dynamic occlusion constraint. ISBI 2017: 718-722

2. Jun Kong, Pengyue Zhang, Yanhui Liang, George Teodoro, Daniel J. Brat, Fusheng Wang:
Robust cell segmentation for histological images of Glioblastoma. ISBI 2016: 1041-1045
******************************************************
# Datasets
******************************************************
Two datasets are used to evaluate the experiment performance: DBM40 dataset from Emory Hospital Achive and TCGA dataset from online public resource. 512x512 image patches are used as processing units due to memory limits. Human annotation are performed on both datasets.
******************************************************
# Shape prior generation
******************************************************
annotation.m: manually select N points on image and save x-y coordinates of the points as 2*N dimensional vector
train_prior.m: load the vectors of shape priors. Interpolate uniformly on the vectores and reorder the points. Save the shape prior library as a matrix in which each colume represents a shape prior.
******************************************************
# Stage 1. Seed detection
******************************************************
Run getSeeds.m to detect seeds. Seeds will be used for level set function initialization in stage 2.
******************************************************
# Stage 2. Contour deformation
******************************************************
To run the program you must include seed files as in '.\data\your_dataset\seed_detection_result\' ,shape prior data '.\data\trainingShape_v3.mat' and the sparse-solving toolbox '.\utilities\l1_ls_matlab\'.

Functions for updating LSF:

LSbatch.m: entrance function. 

lse.m: update the lsf in three steps: updateSR; updatef; updateLSF.

sparse_solver.m: solving sparse representation.

distance_map.m: generate distance maps from list of contour coordinates.

Functions for evaluating performance:

evaluate_performance.m: compute metrics.

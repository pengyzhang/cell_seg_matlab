******************************************************
# Sparse shape prior based nuclei segmentation
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
# Stage 1. Seed detection
******************************************************
Run getSeeds.m to detect seeds. Seeds will be used for level set function initialization in step 2.
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

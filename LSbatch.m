%% @Pengyue Zhang
%% Entrance function for sparse shape prior based levelset function
% optimization
function LSbatch(pathname, imgname, seedPath, resultPath)
if nargin == 0
   clc; close all;
   addpath('./utilities/');
   addpath('./utilities/l1_ls_matlab/');
   pathname = './data/GBM40/raw_image/';
   seedPath = './data/GBM40/seed_detection_result/';
   resultPath = './data/GBM40/result/';
   if ~isdir(resultPath)
       mkdir(resultPath);
   end
   files = dir([pathname '*.tif']);
   imgnames = {files.name}.';
   imgIdx = [1,3,4,7,8,11,13,14,16,19,21,23,26];
   imgnames = imgnames(imgIdx);
   N = length(imgnames);
   for i = 1:N
       fprintf('processing %s--%02d out of %02d.../n', imgnames{i}, i, N);
       LSbatch(pathname, imgnames{i}, seedPath, resultPath);
       close all;
   end
   
   return;
end
% load training shape prior library and initialize shapes
load('./data/trainingShape_v3.mat');
initialPhi = [seedPath imgname '.mat'];
colorI = imread([pathname imgname]);
rowRange = 1:size(colorI,1); colRange = rowRange;
colorI = colorI(rowRange, colRange,:);
I = double(rgb2gray(colorI));
% load seeds from file
load(initialPhi, 'peaks');
[nPeaks, nDim] = size(peaks);
[u,allBW]= initial_phi(rowRange, colRange, peaks);

% parameters
mu=0.08*255^2; % coefficient of arc length term % 0.001
timestep=.1;
xi=2;  % coefficient of penalize term. for distance regularization term (regularize the level set function)
omega=2000;  % coefficient of exclusive term
nu=3000; % coefficient of sparse term

sigma = 4; % scale parameter that specifies the size of the neighborhood 4
lambdaU = 1;
lambdaB = 1; %1
iter_outer=15; 
iter_inner=5;   % inner iteration for level set evolution
epsilon=1;   %for Heaviside function
c0=5;


allInitialLSF = -c0*ones(size(I));
allInitialLSF(allBW) = c0;
allU = allInitialLSF;
figure(1); imagesc(colorI,[0, 255]); axis off; axis equal;
hold on; contour(allInitialLSF,[0 0],'c');
% title('Initial contour','FontSize', 20);

K=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel

Kx = conv2(K, [1 0 -1], 'valid');
Ky = conv2(K, [1 0 -1]', 'valid');
fx = conv2(I,Kx,'same');
fy = conv2(I,Ky,'same');

r = 50; g = exp(-0.5*(fx.*fx+fy.*fy)/r^2); %1.5

for n=1:iter_outer
    if ~exist('transform','var')
        for iPeaks = 1:nPeaks
            transform{iPeaks}.c = ones(size(R,1)/2,1)*[peaks(iPeaks,2),peaks(iPeaks,1)];
            theta = 0;
            transform{iPeaks}.T = [cos(theta),sin(theta);-sin(theta),cos(theta)];
            transform{iPeaks}.b = 1;
        end
    end
    % level set function optimization
    [u, allU, c1, c2, transform]= lse(u, allU, I, g, R, peaks, transform, lambdaU, lambdaB, mu, xi, omega, nu, timestep, epsilon, iter_inner);
        
    fprintf('n = %d/n', n);
    if mod(n,1)==0
        pause(0.001);
        figure,imagesc(colorI,[0, 255]); axis off; axis equal;
        hold on;
        for iPhi=1:size(u,2)
            contour(u{iPhi},[0 0], 'c');
        end
        
        iterNum=[num2str(n), ' iterations'];
        resultName = strcat(resultPath,'color_', imgname,'_iteration_',num2str(n),'.tif');
        hold off;
    end
 
   
end
saveas(gcf,resultName,'tif');
nucleiColor = colorI;
tempR = nucleiColor(:,:,1);
tempG = nucleiColor(:,:,2);
tempB = nucleiColor(:,:,3);
for iPhi = 1:size(u,2)
    bw = u{iPhi}>0;
    edgebw = edge(bw);
    tempR(edgebw)=0;
    tempG(edgebw)=255;
    tempB(edgebw)=0;
    nucleiColor(:,:,1)=tempR;
    nucleiColor(:,:,2)=tempG;
    nucleiColor(:,:,3)=tempB;
    
    b = bwboundaries(bw);
    boundaryCoordinate{iPhi} = b;
end
nucleiColor = uint8(nucleiColor);
figure, imagesc(nucleiColor,[0 255]); axis off; axis equal;
save([resultPath,imgname,'_boundary_coordinates.mat'],'boundaryCoordinate');
end

% initialize shapes with seeds
function [initPhi,allBW] = initial_phi(rowRange, colRange, peaks)
initialRadius = 5;
[nPeaks, nDim] = size(peaks);
[colGrid,rowGrid] = meshgrid(colRange,rowRange);
allBW = zeros(rowRange(end), colRange(end));
for iPeaks = 1:nPeaks
    BW = sqrt((rowGrid-peaks(iPeaks,2)).^2+(colGrid-peaks(iPeaks,1)).^2)<= initialRadius;
    BW = logical(BW);
    allBW = allBW+BW;
    edgeBW = edge(BW);
    dm = -bwdist(edgeBW);
    dm(BW==1) = -1*dm(BW==1);

    initPhi{iPeaks} = dm;
end

allBW = logical(allBW);
    
end
%% @Pengyue Zhang
% 06/03/2015
%% Organizing training set
function train_prior(imgPath,resultPath)
load('data\annotate_v3.mat');
if nargin == 0
    imagedir = 'data\300image\';
    imgPath = dir(fullfile(imagedir, '*.tif'));
    resultPath = 'data\prior\';
end
nFile=length(imgPath);
% nFile=10;
nShape=10;
% nPoint=10;
trainShapeColumn = [];
trainShapeRow = [];
imgResolution = 300;
Img = zeros(imgResolution,imgResolution);
for iFile = 1:nFile
%     Img = imread([imagedir,files(iFile).name]);
    for jShape = 1:nShape
        if fixNumberOfPoint == 0
            corShapeX = annotateX{(iFile-1)*nShape+jShape};
            corShapeY = annotateY{(iFile-1)*nShape+jShape};
            nPoint = length(corShapeX);
        else
            idxShape = (iFile-1)*nShape*nPoint+(jShape-1)*nPoint+1:(iFile-1)*nShape*nPoint+jShape*nPoint; % index of elements corresponding to current shape
            corShapeX = annotateX(idxShape); 
            corShapeY = annotateY(idxShape);
        end
        BW = roipoly(Img, corShapeX, corShapeY);
        imshow(BW);
        hold on
        STATS = regionprops(BW,'Orientation','Centroid','MajorAxisLength','Perimeter');
        plot(STATS.Centroid(1),STATS.Centroid(2),'b+'); % Region Centroids in X-Y coordinate system
        majorAxisAngle = STATS.Orientation*pi/180;
        startpointColumn = STATS.Centroid(1)+STATS.MajorAxisLength*cos(majorAxisAngle)/2; % Approximation of the intersection point between contour and the major axis, which is the startpoint of the shape
        startpointRow = STATS.Centroid(2)-STATS.MajorAxisLength*sin(majorAxisAngle)/2;
        while ~BW(floor(startpointRow), floor(startpointColumn)) % Refine the intersection point
            startpointColumn = startpointColumn - cos(majorAxisAngle)/2;
            startpointRow = startpointRow + sin(majorAxisAngle)/2;
%             plot(startpointX,startpointY,'b+');
        end
        
%         perimeter = 0;
        distFromStartpoint = sqrt((corShapeX-startpointColumn*ones(nPoint,1)).^2+(corShapeY-startpointRow*ones(nPoint,1)).^2); % Sort shape points according to Euclidean dist to startpoint
        [unused,minPoint] = min(distFromStartpoint);
        corShapeX = [corShapeX;corShapeX]; % Duplicate corShapeX
        corShapeY = [corShapeY;corShapeY];
        corShapeX = corShapeX(minPoint:minPoint+nPoint-1); % Rearrange point on ellipse
        corShapeY = corShapeY(minPoint:minPoint+nPoint-1);
        segLength = 0;
        for kPoint = 1:nPoint
%             if kPoint < nPoint
%                 perimeter = perimeter + sqrt((corShapeX(kPoint)-corShapeX(kPoint+1))^2+(corShapeY(kPoint)-corShapeY(kPoint+1))^2);
%             else 
%                 perimeter = perimeter + sqrt((corShapeX(kPoint)-corShapeX(1))^2+(corShapeY(kPoint)-corShapeY(1))^2);
%             end
            
            if kPoint == 1
                segLength(kPoint) = sqrt((corShapeX(kPoint)-startpointColumn)^2+(corShapeY(kPoint)-startpointRow)^2);
            elseif kPoint <= nPoint
                segLength(kPoint) = sqrt((corShapeX(kPoint)-corShapeX(kPoint-1))^2+(corShapeY(kPoint)-corShapeY(kPoint-1))^2);
            end
            
        end
        segLength(nPoint+1) = sqrt((corShapeX(nPoint)-startpointColumn)^2+(corShapeY(nPoint)-startpointRow)^2);
        perimeter = sum(segLength);
        arcLength = cumsum(segLength,2)./perimeter; % Uniform arclength in [0,1]
        arcLength = [0,arcLength];
        arcGrid = 0:0.01:1;
        interpColumn = interp1(arcLength',[startpointColumn;corShapeX;startpointColumn],arcGrid','spline'); % Interpolation with 101 points on each shape
        interpRow = interp1(arcLength',[startpointRow;corShapeY;startpointRow],arcGrid','spline');
        plot(interpColumn,interpRow,'c.');
        interpColumn = interpColumn - STATS.Centroid(1); % Relative position
        interpRow = interpRow - STATS.Centroid(2);
        
        trainShapeColumn = [trainShapeColumn interpColumn]; % Row-Column coordinate of training shapes 
        trainShapeRow = [trainShapeRow interpRow];
        
    end
end
%% Align all shape to a reference shape.
[nFeature, allShape] = size(trainShapeColumn);
referenceShape = [trainShapeColumn(:,1),trainShapeRow(:,1)]; % First shape as reference shape
Dim = 2;
vectorLength = nFeature*Dim;
R = zeros(vectorLength,nShape);
R(:,1) = reshape(referenceShape',vectorLength,1);
for iShape = 2:allShape % Align all shapes to reference shape
    tempShape = [trainShapeColumn(:,iShape),trainShapeRow(:,iShape)];
    [d1,alignedShape,transform] = procrustes(referenceShape,tempShape);
%     for jPoint = 1:101
%         plot(alignedShape((jPoint-1)*2+1),alignedShape(jPoint*2),'r.');
%     end
    R(:,iShape)=reshape(alignedShape',vectorLength,1); % Training shape set of size 202*100
end

save([resultPath,'trainingShape_v3.mat'],'R');
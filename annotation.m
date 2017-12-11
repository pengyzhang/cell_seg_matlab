%% @Pengyue Zhang
%% 05/28/2015
%% Extract shape priors manually
function annotation(imgPath,resultPath)
% load('annotate.mat');
if nargin == 0
    imagedir = 'data\500image\';
    imgPath = dir(fullfile(imagedir, '*.tif'));
    resultPath = 'data\annotation\';
end
nFile=length(imgPath);
nShape=10;

if ~exist('numberofPoint','var')
    fixNumberOfPoint = 0;
    for iFile = 1:nFile
        Img = imread([imagedir,imgPath(iFile).name]);
        figure,imshow(Img,'InitialMagnification',200);
        hold on
        for jShape = 1:nShape
            [x,y] = ginput;
            annotateX{(iFile-1)*nShape+jShape} = x;
            annotateY{(iFile-1)*nShape+jShape} = y;
            plot(x,y,'b*');
            disp(['Annotated:Files ',num2str(iFile),' Shape: ', num2str(jShape)]);
    %         plot(annotateX((i-1)*numberofShape*numberofPoint+(j-1)*numberofPoint+1:(i-1)*numberofShape*numberofPoint+j*numberofPoint),annotateY((i-1)*numberofShape*numberofPoint+(j-1)*numberofPoint+1:(i-1)*numberofShape*numberofPoint+j*numberofPoint));
        end
    end
else
    annotateX=[];
    annotateY=[];
    fixNumberOfPoint = 1;
    for iFile = 1:nFile
        Img = imread([imagedir,imgPath(iFile).name]);
        figure,imshow(Img,'InitialMagnification',200);
        hold on
        for jShape = 1:nShape
            [x,y] = ginput(numberofPoint);
            annotateX = [annotateX;x];
            annotateY = [annotateY;y];
            disp(['Annotated:Files ',num2str(iFile),' Shape: ', num2str(jShape)]);
    %         plot(annotateX((i-1)*numberofShape*numberofPoint+(j-1)*numberofPoint+1:(i-1)*numberofShape*numberofPoint+j*numberofPoint),annotateY((i-1)*numberofShape*numberofPoint+(j-1)*numberofPoint+1:(i-1)*numberofShape*numberofPoint+j*numberofPoint));
        end
    end
end
save([resultPath,'annotate_v3.mat'],'annotateX','annotateY','fixNumberOfPoint');
imgPath = '.\data\512image\';
boundaryPath = '.\NucleiSegLS\nu\3000\';
imgIdx = [1,3,4,7,8,11,13,14,16,19,21,23,26];
for i = 1:1
    colorI = imread([imgPath sprintf('%02d',imgIdx(i)) '.tif']);
    load([boundaryPath sprintf('%02d',imgIdx(i)) '_boundary_coordinates.mat']);
    emptyIdx = [];
    for j = 1:length(boundaryCoordinate)
        
        cell = boundaryCoordinate{j};
        
        if isempty(cell)
            emptyIdx = [emptyIdx j];
        end
        
    end
    boundaryCoordinate(emptyIdx) = [];
    save([boundaryPath sprintf('%02d',imgIdx(i)) '_boundary_coordinates.mat'],'boundaryCoordinate');
end
% colorMap = [255,0,0;0,255,0;255,255,51;255,255,0;255,0,255;0,255,255;255,255,255];
colorMap = [1,0,0;0,1,0;1,1,0.2;1,1,0;1,0,1;0,1,1;1,1,1];
% colorMap = [0,255,255;0,255,255;0,255,255;0,255,255;0,255,255;0,255,255;0,255,255];

% for i = 1:1
%     colorI = imread([imgPath sprintf('%02d',imgIdx(i)) '.tif']);
%     load([boundaryPath sprintf('%02d',imgIdx(i)) '_boundary_coordinates.mat']);
%     for j = 1:length(boundaryCoordinate)
%         cell = boundaryCoordinate{j};
%         object = cell{1,1};
%         rem = mod(j,7);
%         if rem == 0
%             for k = 1:size(object,1)
%                 colorI(object(k,1),object(k,2),1) = colorMap(7,1); 
%                 colorI(object(k,1),object(k,2),2) = colorMap(7,2);
%                 colorI(object(k,1),object(k,2),3) = colorMap(7,3); 
%             end
%         else
%             for k = 1:size(object,1)
%                 colorI(object(k,1),object(k,2),1) = colorMap(rem,1); 
%                 colorI(object(k,1),object(k,2),2) = colorMap(rem,2);
%                 colorI(object(k,1),object(k,2),3) = colorMap(rem,3); 
%             end
%         end
%         
%     end
%     figure, imshow(colorI);
% %     imwrite(colorI, [boundaryPath sprintf('%02d',i) '_final_contour_color.tif'],'tif');
% end
for i = 1:1
    colorI = imread([imgPath sprintf('%02d',imgIdx(i)) '.tif']);
    load([boundaryPath sprintf('%02d',imgIdx(i)) '_boundary_coordinates.mat']);
    figure, imagesc(colorI,[0 255]), axis off, axis equal, hold on;
    for j = 1:length(boundaryCoordinate)
        cell = boundaryCoordinate{j};
        object = cell{1,1};
        X = object(:,2);
        Y = object(:,1);
        BW = roipoly(colorI, X, Y);
        edgeBW = edge(BW);
        dm = -bwdist(edgeBW);
        dm(BW==1) = -1*dm(BW==1);
        rem = mod(j,7);
        if rem == 0
            [c,h] = contour(dm,[0,0],'color',colorMap(7,:));
            set(h,'linewidth',5);
        else
            [c,h] = contour(dm,[0,0],'color',colorMap(rem,:));
            set(h,'linewidth',5);
        end
    end
end
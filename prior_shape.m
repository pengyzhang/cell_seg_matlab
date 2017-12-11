close all;
dataPath = '.\data\prior\';
load([dataPath 'trainingShape_v3.mat']);
[coordinates nShape] = size(R);
points = coordinates/2;
imgSize = 80;

for i = 1:nShape
    img = ones(imgSize,imgSize);
    shape = R(:,i);
    for j = 1:points
        X(j) = round(shape((j-1)*2+1)+imgSize/2);
        Y(j) = round(shape(j*2)+imgSize/2);
        img(X(j),Y(j)) = 0;
    end
    
    imagesc(img); colormap(gray); axis off;
    set(gcf,'color','none');
end
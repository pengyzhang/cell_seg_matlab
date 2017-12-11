imgPath = 'E:\Research\Code\sparse shape prior\data\512image\';
for i = 1:40
    colorI = imread([imgPath sprintf('%02d',i) '.tif']);
    I = rgb2gray(colorI);
    I2 = imcomplement(I);
    imwrite(I2,[imgPath sprintf('%02d',i) '_complement.tif']);
end
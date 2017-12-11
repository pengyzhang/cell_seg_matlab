function getSeeds_old(inputFileName)

%process all .bmp images
if nargin==0
    dataPath = 'E:\Research\Code\sparse shape prior\data\512image\';
    
    inputFileName = dir([dataPath '*.tif']);
    inputFileName = {inputFileName.name}.';
    
    %k = strfind(inputFileName, '.');
    %imageNumber = cellfun(@(x,y) str2double(y(x(1)+1:x(2)-1)), k, inputFileName);
    %[~, ind] = sort(imageNumber);
    %inputFileName = inputFileName(ind);
    
    for i = 1:length(inputFileName)
        fprintf('Image %d (%s) is processed... \n', i, inputFileName{i});
        getSeeds_old([dataPath  inputFileName{i}]);
    end
    return;
end


close all; 
imtool close all;
warning('off', 'all');
clc;

%load image
I=imread(inputFileName);


%define OD matrix; each column is associated with one stain (i.e. Hematoxylin, Eosin, and red_marker)
%        Hemat      Eosin    Null
stains =[0.554688  0.380814  0;...   %Red
         0.781334  0.87215   0;...   %Green
         0.286075  0.307141  0];     %Blue
% stains =[0.554688  0.380814  0.13193;...   %Red
%          0.781334  0.87215   0.900097;...   %Green
%          0.286075  0.307141  0.415234];     %Blue
% stains =[0.554688  0.380814  0.191667;...   %Red
%          0.781334  0.87215   0.862937;...   %Green
%          0.286075  0.307141  0.46755];     %Blue

%calculate stain intensities using color deconvolution
[Deconvolved, colorImage] = ColorDeconvolution_FullNewVer(I, stains, [true true true]);
Hemat = Deconvolved(:,:,1); 
Eosin = Deconvolved(:,:,2); 
%HematColor = colorImage{1};
%EosinColor = colorImage{2};

%figure;imshow(I,[]);
%figure;imshow(HematColor,[]);
%figure;imshow(EosinColor,[]);
%figure;imshow(Hemat,[]);
%figure;imshow(Eosin,[]);



%complement and smooth image
sigma=5; %important!
G = fspecial('gaussian',round(3*sigma)*2+1,sigma); % Gaussian kernel
h = imcomplement(Hemat);
f = conv2(h, G, 'same');
%[fx, fy] = gradient(f);
%edge_canny = edge(f,'canny');


%morphological reconstruction
uint8f = uint8(f);
marker = imopen(uint8f, strel('disk',10));
recon = imreconstruct(marker,uint8f, 8);
dif = uint8f - recon;
%imtool(dif,[]);


% peaksBW = imregionalmax(dif);

%[fx,fy] = gradient(double(dif));
%firstOrder = sqrt(fx.^2 + fy.^2) < eps;
%secondOrder =  (Lambda1<0);
%peaksBW = firstOrder & secondOrder;

% figure;imshow(dif,[]); hold on;
% [r,c] = ind2sub(size(peaksBW), find(peaksBW));
% for i = 1:length(r)
%     plot(c(i), r(i), 'g*');
% end
    

%voting with sign of eigenvalues from Hessian matrix
vote = zeros(size(dif));
for sigma = 3:0.3:10;
%for sigma = 3:0.3:15
%for sigma = 2:0.3:10
[fxx,fxy,fyy] = Hessian2D(double(dif), sigma);
fxx = (sigma^2)*fxx;
fxy = (sigma^2)*fxy;
fyy = (sigma^2)*fyy;

tmp = sqrt((fxx - fyy).^2 + 4*fxy.^2);
Lambda1 = 0.5*(fxx + fyy + tmp);
%Lambda2 = 0.5*(fxx + fyy - tmp);
TF = Lambda1 < 0;
vote(TF) = vote(TF) + 1;
end

%imshow(vote,[]); axis ij;
%surf(vote); axis ij; colormap jet;



%otusu thresholding
level = graythresh(uint8f);
otsuBW = im2bw(uint8f, level);
minOtsuA = 150;
otsuBW = bwareaopen(otsuBW, minOtsuA, 8);


%peak detection on voting map (merge by connection and distance)
minD = 18; 
minA = 20; 
alpha = 10;

v = sort(unique(vote(:)),'descend');

vnorm = (flipud(v) - v(end))/(v(1)-v(end)); %normalize [vmax, vmin] to [0, 1]
peaks = []; %peaks(x, y);
for i = 1:length(v)
    a = vote >= v(i);
    L = bwlabel(a, 8);
    c = regionprops(L, 'Area', 'Centroid');
    d = round(reshape([c.Centroid], 2,length(c))');
    d = max(1, d);
    d(:,1) = min(size(vote,2), d(:,1));
    d(:,2) = min(size(vote,1), d(:,2));
    otsuTF = otsuBW((d(:,1)-1)*size(vote,1)+d(:,2))';
    b = ismember(L, find( otsuTF & ([c.Area]>(minA+exp(alpha*vnorm(i)))) ) );
    %b = ismember(L, find([c.Area]> max(20, minA-2*v(i))) );
    
    L = bwlabel(b, 8);
    c = regionprops(L, 'Centroid');
    cand_length = length(c);
    fprintf('new %d center candidates at voting threshold: %d\n', cand_length, v(i));
    
    if ~isempty(peaks)  %merge by connection
        removeL = L((peaks(:,1)-1) * size(vote,1) + peaks(:,2));
        removeTF = ismember(1:max(L(:)), removeL);
        c(removeTF) = [];
        fprintf('reduce %d centers covered by connected component of upper level centers:\n', sum(removeTF));
    else
        removeTF = [];
    end
    

    c = reshape([c.Centroid], 2,length(c))';
    c = round(c);
    c = max(1, c);
    c(:,1) = min(size(vote,2), c(:,1));
    c(:,2) = min(size(vote,1), c(:,2));
    
    c = [peaks; c];
    dist_length = 0;
    while true %merge by distance metric
        
        D = pdist(c,'euclidean');
        if all(D>minD)
            break;
        end
        
        D = squareform(D);
        mergeTF = (D<=minD) ;
        ind = ([1:size(D,1)]-1)*size(D,1) + [1:size(D,1)];
        mergeTF(ind) = true; %diagonal entries of D == 0, thus needs to be set to true mannually
        
        
        TF = false(1,size(D,2));
        for j = 1:size(D,1)-1
            
            lineTF = [false(1, j) mergeTF(j, j+1:end)];
            if any(TF&lineTF)
               continue;
            end
            
            TF = TF | lineTF;
        end
        c(TF',:) = [];
        fprintf('reduce %d centers by distance metric: minD=%d\n', sum(TF), minD);
        dist_length = dist_length + sum(TF);
        
        %imshow(dif,[]); hold on;
        %scatter(c(:,1), c(:,2), 100, 'g+');
    end
    
    
    orig_length = size(peaks,1);
    peaks = c;
    fprintf(['After voting threshold %d: '...
            '%d(original)+%d(new candidate)-%d(connection)-%d(distance)=%d\n\n'],...
             v(i), orig_length, cand_length, sum(removeTF), dist_length, size(peaks,1));
    if (orig_length+cand_length-sum(removeTF)-dist_length)~=size(peaks,1)
            error('center number does not match: orig_length+cand_length-sum(removeTF)-dist_length=%d; size(peaks)=%d',...
                  orig_length+cand_length-sum(removeTF)-dist_length, size(peaks,1) );
    end
    
end

%peaks(bottom+1:end, :) = [];
figure, imshow(I,[]); hold on;
scatter(peaks(:,1), peaks(:,2), 100, 'g+');
resultPath = 'E:\Research\Code\sparse shape prior\data\seed_detection_result\';
resultFileName = inputFileName(end-5:end);
imwrite(I,[resultPath resultFileName],'tif');
% saveas(gcf,[resultPath resultFileName],'tif');
save([resultPath resultFileName '.mat'],'peaks');
% keyboard;

% r=I(:,:,1);
% g=I(:,:,2);
% b=I(:,:,3);
% g(otsuBW)=150;
% figure; h=imshow(cat(3,r,g,b),[]);hold on;scatter(peaks(:,1), peaks(:,2), 100, 'g+');

end



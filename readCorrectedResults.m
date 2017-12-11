clearvars; close all; clc;

code_dir = '/Users/kongj/Renal_Chapman/';
addpath(genpath(code_dir));
if exist('openSlide_c.mexmaci64','file')~=3
    %mex -I/opt/local/include/openslide/ -L/opt/local/lib -lopenslide openSlide_c.cpp
    mex -I/usr/local/include/openslide/ -L/usr/local/lib -lopenslide openSlide_c.cpp
end

roidir = '/Users/kongj/Astrocyoma_IDH_nuclei/data/extracted_ROIs/';
xmldir = '/Users/kongj/Astrocyoma_IDH_nuclei/data/corrected_result_xmls/';

files = dir([xmldir '*.xml']);
dirIdx = [files.isdir];
files = {files(~dirIdx).name}';
N = length(files);

dashind = strfind(files, '-');
outputfiles = cellfun(@(x,y) [x(1:y(3))], files, dashind, 'UniformOutput', false); 

red_sizes = [];
stats = cell(N, 1);

if ~exist([xmldir 'features.mat'], 'file')
    save([xmldir 'features.mat'], 'stats', 'files');
end

for i = 1:N
   
   xmlname = [xmldir files{i}];
   matname = [xmlname(1:end-4) '.mat'];
   if ~exist(matname, 'file')
        XMLAnno = xmlStr2struct(xmlname);
        save(matname, 'XMLAnno');
        fprintf('save ''XMLAnno'' as %s\n', matname);
        continue;
   end
        
   %continue;
   
   imagename = [xmlname(1:end-4) '.tif'];
   %if exist(imagename, 'file')
   %    continue;
   %end
   
   fprintf('process %s...\n', matname);
   temp = load(matname, 'XMLAnno');
   XMLAnno = temp.XMLAnno;
   clear temp
   
   
   inputImageName = [roidir files{i}(1:end-9) '.tif'];
   info = imfinfo(inputImageName);
   if length(info)>1
       I = imread(inputImageName, 1);
   else
       I = imread(inputImageName);
   end
   red = I(:,:,1); green = I(:,:,2); blue = I(:,:,3);
   L = zeros(size(I,1), size(I,2));
   label = 0;
   
   Shape = XMLAnno.Document{1}.Layers{1}.Layer{1}.Shapes{1}.Shape;
   S = length(Shape);
   for j = 1:S
       Current_Shape = Shape{j};
       Point = Current_Shape.Data{1}.Points{1}.Point;
       P = length(Point);
       Coord = NaN(P, 2);
       for k = 1:P
          Current_Point = Point{k};
          Coord(k,:) = [str2double(Current_Point.ATTR.X) str2double(Current_Point.ATTR.Y)]+1;
       end
       
       Coord = round(Coord);
       Coord(:,1) = min(max(Coord(:,1), 1), size(I,2)); %X or col
       Coord(:,2) = min(max(Coord(:,2), 1), size(I,1)); %Y or row
       %figure;plot(Coord(:,1), Coord(:,2), 'o--'); hold on;
       
       try
           Coord = smooth_interp(Coord);
           %plot(Coord(:,1), Coord(:,2), 'r+--');
       catch
           keyboard;
       end
       
       cell_bw = roipoly(I(:,:,1), Coord(:,1), Coord(:,2)); 

       if str2double(Current_Shape.Settings{1}.Line{1}.Color{1}.ATTR.R) > 50
           red_sizes = [red_sizes sum(cell_bw(:))]; 
           continue;
       end
       
       label = label+1;
       L(cell_bw) = label;
       
       tempB = bwperim(cell_bw);
       red(tempB) = 0;
       green(tempB) = 255;
       blue(tempB) = 0;
       
   end
   
   %save images overlaid with corrected results
   if exist(imagename, 'file')
       continue;
   else
        imwrite(cat(3, red, green, blue), imagename);
   end
   
end



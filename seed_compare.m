clear all; close all;

Index = [1:40];
n = length(Index);
imgPath = '.\data\512image\';
imageJpath = '.\parvin\seed_result\';
peakPath = '.\seed_detection_result\';
resultPath = '.\seed_detection_human\';
for i = 1:2
    imageJ_filename = [imageJpath sprintf('%02d',Index(i)) '_coordinates'];
    imageJ_seed = load(imageJ_filename);
    peak_filename = [peakPath sprintf('%02d',Index(i)) '.tif.mat'];
    peak = load(peak_filename);
    peak = peak.peaks;
    human_filename = [resultPath sprintf('%02d',Index(i)) '_data.xml'];
    [human_seed unused] = read_xml(human_filename);
    colorI = imread([imgPath sprintf('%02d',Index(i)),'.tif']);
    figure,imagesc(colorI,[0, 255]); axis off; axis equal;
    hold on;
    scatter(imageJ_seed(:,2),imageJ_seed(:,3),50,'y+');
    scatter(peak(:,1),peak(:,2),50,'co');
%     scatter(human_seed(:,1), human_seed(:,2),100,'r+');
%     result_name = strcat(resultPath,sprintf('%02d',Index(i)));
%     saveas(gcf,result_name,'tif');
    
end
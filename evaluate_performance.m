%% @Pengyue Zhang
% Given the final xmls from Carol, five quantities will be computed for
% each of randomly selected cells:
% Jaccard Coef(A,B) = |A intersect B| / |A union B|
% P(A,B) = |A intersect B| / B
% R(A,B) = |A intersect B| / A
% F1(A,B) = 2*P*R/(P+R)
% Dh(A,B): Haussdorf distance


clear all; close all;

humanPath = '.\metrics\LS_human\';
machiPath = '.\metrics\LS_result\';
SRow = 512; SCol = 512;


if exist('.\validation.mat', 'file') ~= 0
    load('.\validation.mat');
else

    humanFiles = dir([humanPath '*.xml']);
    machiFiles = dir([machiPath '*.mat']);
    
    humanFiles = {humanFiles.name}.';
    machiFiles = {machiFiles.name}.';

    N = length(humanFiles);
    humanPts = [];
    
    J = cell(1,N);
    P = J;
    R = J;
    F1 = J;
    HausdorfD = J;
    
    miss = 0;
    for i = 1:13
        fprintf('Analyzing file %s ...\n', humanFiles{i});
        humanPts = xml2coordinate([humanPath humanFiles{i}], SRow, SCol);
        
        H = size(humanPts,1);
        Jtemp = NaN(H,1);
        Ptemp = Jtemp;
        Rtemp = Jtemp;
        F1temp = Jtemp;
        HausdorfDtemp = Jtemp;
        
        
        load([machiPath humanFiles{i}(1:2) '_boundary_coordinates.mat'], 'boundaryCoordinate');
        %load([machiPath machiFiles{i}], 'boundaryCoordinate');
        
        M = length(boundaryCoordinate);
        machiPts = cell(1,M);
        for j = 1: M
            machiPts{j} = boundaryCoordinate{j}{1};
        end
        
        
        
        for h = 1: H
           hp = humanPts{h}; 
           hbw = roipoly(512,512,hp(:,1), hp(:,2)); 
           hp = bwboundaries(hbw,8); hp = hp{1};
           
           indicator = zeros(M,1);
           for m = 1:M
               mp = machiPts{m};
               mbw = roipoly(SRow,SCol, mp(:,2), mp(:,1));
               
               indicator(m) = sum(hbw(:)&mbw(:)) / (sum(hbw(:)|mbw(:))+eps);
           end
           
           [maxIndicator, ind] = max(indicator);
           if maxIndicator==0
               miss = miss+1;
               continue;
           end
           
           mp = machiPts{ind};
           mbw = roipoly(SRow,SCol, mp(:,2), mp(:,1));
           Jtemp(h) = sum(hbw(:)&mbw(:)) / (sum(hbw(:)|mbw(:))+eps);
%            if Jtemp(h)<0.5
%               figure; imshow(imread([humanPath humanFiles{i}(1:2) '.tif']),[]);
%               hold on;
%               plot(hp(:,2), hp(:,1), 'r'); plot(mp(:,2),mp(:,1),'g');
%            end
           Ptemp(h) = sum(hbw(:)&mbw(:)) / (sum(mbw(:))+eps);
           Rtemp(h) = sum(hbw(:)&mbw(:)) / (sum(hbw(:))+eps);
           F1temp(h) = 2*Ptemp(h)*Rtemp(h)/ (Ptemp(h)+Rtemp(h)+eps);
           [HausdorfDtemp(h), ~] = HausdorffDist([hp(:,2) hp(:,1)], [mp(:,2) mp(:,1)]);
           
        end
        
        
        J{i} = Jtemp;
        P{i} = Ptemp;
        R{i} = Rtemp;
        F1{i} = F1temp;
        HausdorfD{i} = HausdorfDtemp;
        
        fprintf('J: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n', nanmean(Jtemp), nanstd(Jtemp), max(Jtemp), quantile(Jtemp,0.75), nanmedian(Jtemp), quantile(Jtemp,0.25), min(Jtemp) );
        fprintf('P: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n', nanmean(Ptemp), nanstd(Ptemp), max(Ptemp), quantile(Ptemp,0.75), nanmedian(Ptemp), quantile(Ptemp,0.25), min(Ptemp) );
        fprintf('R: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n', nanmean(Rtemp), nanstd(Rtemp), max(Rtemp), quantile(Rtemp,0.75), nanmedian(Rtemp), quantile(Rtemp,0.25), min(Rtemp) );
        fprintf('F1: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n',nanmean(F1temp), nanstd(F1temp), max(F1temp), quantile(F1temp,0.75), nanmedian(F1temp), quantile(F1temp,0.25), min(F1temp) );
        
    end
    
    %save('./validation.mat', 'J', 'P', 'R', 'F1', 'HausdorfD');
end


Jvec = cell2mat(J'); Jvec(isnan(Jvec),:) = [];
Pvec = cell2mat(P'); Pvec(isnan(Pvec),:) = [];
Rvec = cell2mat(R'); Rvec(isnan(Rvec),:) = [];
F1vec = cell2mat(F1'); F1vec(isnan(F1vec),:) = [];
Hvec = cell2mat(HausdorfD'); Hvec(isnan(Hvec),:) = [];

fprintf('J: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n', mean(Jvec), std(Jvec), median(Jvec), max(Jvec), min(Jvec) );
fprintf('P: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n', mean(Pvec), std(Pvec), median(Pvec), max(Pvec), min(Pvec) );
fprintf('R: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n', mean(Rvec), std(Rvec), median(Rvec), max(Rvec), min(Rvec) );
fprintf('F1: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n',mean(F1vec), std(F1vec), median(F1vec), max(F1vec), min(F1vec) );
fprintf('dH: %3.2f \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n',mean(Hvec), std(Hvec), median(Hvec), max(Hvec), min(Hvec) );

fprintf('Total number of annotated cells: %d \n', length(Jvec));

% figure;
% subplot(5,1, 1); bar(1-J(ind_evaluation)); ax = axis; axis([0 length(ind_evaluation)+0.5 ax(3:4)]); set(gca,'FontSize', 25); h=ylabel('1-J'); set(h,'FontSize',45);
% subplot(5,1, 2); bar(1-P(ind_evaluation)); ax = axis; axis([0 length(ind_evaluation)+0.5 ax(3:4)]); set(gca,'FontSize', 25); h=ylabel('1-P'); set(h,'FontSize',45);
% subplot(5,1, 3); bar(1-R(ind_evaluation)); ax = axis; axis([0 length(ind_evaluation)+0.5 ax(3:4)]); set(gca,'FontSize', 25); h=ylabel('1-R'); set(h,'FontSize',45);
% subplot(5,1, 4); bar(1-F1(ind_evaluation)); ax = axis; axis([0 length(ind_evaluation)+0.5 ax(3:4)]); set(gca,'FontSize', 25); h=ylabel('1-F_1'); set(h,'FontSize',45);
% subplot(5,1, 5); bar(HausdorfD(ind_evaluation));ax = axis; axis([0 length(ind_evaluation)+0.5 ax(3:4)]); set(gca,'FontSize', 25); h=ylabel('d_H'); set(h,'FontSize',45);
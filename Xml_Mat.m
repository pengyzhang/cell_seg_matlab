clear all; clc;

filename = '.\metrics\LS_human\01_data.xml';

try
    fid = fopen(filename);
catch
   error('Failed to open XML annotation file %s.',filename);
   return;
end

RNUM = 512; 
CNUM = 512;
INITIAL = 500;
Cells = cell(INITIAL,1);
i = 0;
tline = fgetl(fid);
line = 1;

startStr = '<Points>';
endStr = '</Points>';
while ischar(tline)
    fprintf('line #:%d\n',line);
    tline = fgetl(fid);
    line = line + 1;
    
    if length(tline)<length(startStr) || ~strcmp(tline(end-length(startStr)+1:end), startStr)
        continue;
    end
    
    i = i + 1;
    temp = [];
    tline = fgetl(fid);
    line = line+1;
    while ~strcmp(tline(end-length(endStr)+1:end), endStr)
        k = strfind(tline, '"');
        if length(k) ~= 4
            error(['X-, y-coordinates are wrong on line ' sprintf('%d', line)]);
        end
        x = round(str2double(tline(k(1)+1:k(2)-1)));
        y = round(str2double(tline(k(3)+1:k(4)-1)));
        temp = [temp; y, x];
        
        tline = fgetl(fid);
        line = line+1;
    end
    
    
    temp(temp(:,1)<1,1) = 1;
    temp(temp(:,1)>RNUM,1) = RNUM;
    temp(temp(:,2)<1,2) = 1;
    temp(temp(:,2)>CNUM,2) = CNUM;
    
    Cells{i} = temp;
end

fclose(fid);
Cells(cellfun(@isempty, Cells)) = [];

%poly2mask
%stack 2D masks
%bwlabeln for n-dim array
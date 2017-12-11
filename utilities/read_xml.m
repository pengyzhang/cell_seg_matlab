%% Read xml annotation files
function [Pts, colors]  = read_xml(filename)

try
    fid = fopen(filename);
catch execption
   error('Failed to open XML annotation file %s.',filename);
   
end

INITIAL = 1000;
Pts = NaN(INITIAL,1);
colors = cell(INITIAL,1);
colors(:) = {'green'};

i = 1;
OVER = false;
tline = fgetl(fid);
line = 1;


LINE_STR_Ver_21 = '<Line Width="1" Dash="Solid" DashOffset="0" Join="Round" OutlineType="Color">';
LINE_STR_Ver_18 = '<Line Width="1" Dash="Solid" Join="Round" OutlineType="Color">';
LINE_COLOR_STR =  '<Color Alpha="255" R="0" G="255" B="0" />';
Data_STR = '<Data Rotation="0">';

while ischar(tline)
    
   while ~(strcmpi(LINE_STR_Ver_21, tline) || strcmpi(LINE_STR_Ver_18, tline))
      tline = fgetl(fid);
      if ~ischar(tline)
          OVER = true;
          break;
      end
      tline = strtrim(tline);
      line = line + 1;
   end
   
   if OVER
       break;
   end
   
   tline = strtrim(fgetl(fid));
   line = line + 1;
    
   if ~strcmpi(LINE_COLOR_STR, tline)
       colors{i} = 'red';
   end
   
   tline = strtrim(fgetl(fid));
   line = line + 1;
   
   while ~strcmpi(Data_STR, tline)
      tline = strtrim(fgetl(fid));
      line = line + 1;
   end
   
   tline = strtrim(fgetl(fid));
   line = line + 1;
   
   k = strfind(tline,'"');
   if length(k) ~= 8
       error(sprintf('%s contains a line with wrong format for <Extent X,Y,Width,Height/>',filename));
   end
   
   x = str2double(tline((k(1)+1):(k(2)-1)));
   y = str2double(tline((k(3)+1):(k(4)-1)));
   width = str2double(tline((k(5)+1):(k(6)-1)));
   height = str2double(tline((k(7)+1):(k(8)-1)));
   Pts(i,1) = x+1+width/2;
   Pts(i,2) = y+1+height/2;
   
   i = i + 1;
end

fclose(fid);


nanTF = any(isnan(Pts),2);
Pts(nanTF,:) = [];
colors(nanTF) = [];




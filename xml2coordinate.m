function Pts = xml2coordinate(filename, RNUM, CNUM)
try
    fid = fopen(filename);
catch execption
   error('Failed to open XML annotation file %s.',filename);
   
end

%RNUM = 512; 
%CNUM = 512;
INITIAL = 500;
Pts = cell(INITIAL,1);

i = 0;
tline = fgetl(fid);
line = 1;


startStr = '<Points>';
endStr = '</Points>';

SkipTF = true;

while ischar(tline)
    %fprintf('line #:%d\n',line);
    tline = fgetl(fid);
    line = line + 1;
    
    if ischar(tline)
        tline = strtrim(tline);
    else
        break;
    end
    
    if strcmpi(startStr, tline)
        i = i + 1;
        temp = [];
        while ~strcmpi(endStr, tline)
            %fprintf('line #:%d\n',line);
            tline = strtrim(fgetl(fid));
            line = line + 1;
            
            if strcmpi(endStr, tline)
                break;
            end
            k = strfind(tline,'"');
            x = round(str2double(tline( (k(1)+1) : (k(2)-1) )));
            y = round(str2double(tline( (k(3)+1) : (k(4)-1) )));
            temp = [temp; x+1 y+1];
            
        end
        
         temp(temp(:,1)<1,1) = 1;
         temp(temp(:,1)>CNUM,1) = CNUM;
         temp(temp(:,2)<1,2) = 1;
         temp(temp(:,2)>RNUM,2) = RNUM;
         
         Pts{i} = temp;
        
    else
        continue;
    end
  
end

fclose(fid);
emptyTF = cellfun(@isempty, Pts);
Pts(emptyTF) = [];
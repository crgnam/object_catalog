matlabrc; clc; close all;

fid = fopen('lookup.csv');
fid2 = fopen('lookup2.csv','wt');
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end 
    fprintf(fid2, [deblank(line) ,',\n']);
end
fclose(fid);
fclose(fid2);
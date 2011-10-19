function p = readLandmarks(filename)
%-------------------------------------------------------
% function p = readLandmarks(filename)
%     filename: the filename of the landmark file

% read in landmark points written out by the markTags program.
% 
% hari sundar 06.07.07
%-------------------------------------------------------

fp=fopen(filename,'r');
nt = fscanf(fp,'%d', 1);

for i=1:nt
    nn = fscanf(fp, '%d', 1);
    for j=1:nn
        p(:, j, i)  = fscanf(fp, '%f %f', 2);
    end
end

fclose(fp);
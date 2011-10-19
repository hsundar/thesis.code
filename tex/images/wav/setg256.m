function []=SetG256
% Set the colormap to G256
if exist('G256')~= 1
   G256 = [0:1/255:1;0:1/255:1;0:1/255:1]';
end 	
colormap(G256);
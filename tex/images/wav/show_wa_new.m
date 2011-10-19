clear all
close all

%img = imread('show_wa_01.bmp','bmp');
img = imread('mr.png');

xx = 90;
yy = 90;

d1 = img(xx-56:xx+55,yy-56:yy+55);
d2 = img(xx-28:xx+27,yy-28:yy+27);

img(xx-1:xx,yy-1:yy) = 255;
%img(xx-56:xx+55,yy-56) = 255;
%img(xx-56:xx+55,yy+55) = 255;
%img(xx-56,yy-56:yy+55) = 255;
%img(xx+55,yy-56:yy+55) = 255;

figure;
imagesc(img);
setg256;
axis off;
axis image;
hold on;
line([xx-56,xx-56],[yy-56,yy+55],'LineWidth',2,'Color','blue');
line([xx-56,xx+55],[yy-56,yy-56],'LineWidth',2,'Color','blue');
line([xx+55,xx-56],[yy+55,yy+55],'LineWidth',2,'Color','blue');
line([xx+55,xx+55],[yy-56,yy+55],'LineWidth',2,'Color','blue');

line([xx-28,xx-28],[yy-28,yy+28],'LineWidth',2,'Color','red');
line([xx-28,xx+28],[yy-28,yy-28],'LineWidth',2,'Color','red');
line([xx+28,xx-28],[yy+28,yy+28],'LineWidth',2,'Color','red');
line([xx+28,xx+28],[yy-28,yy+28],'LineWidth',2,'Color','red');
[X,Y] = size(d1);

d1 = double(d1);
for i=1:X/2
    for j=1:Y/2
        nd1(i,j) = (d1(i*2,j*2)+d1(i*2-1,j*2)+d1(i*2,j*2-1)+d1(i*2-1,j*2-1))/4;
    end;        
end;
clear d1;
d1 = nd1;

[X,Y] = size(d1);

[ca1, ch1, cv1, cd1] = dwt2(d1, 'haar');

ca1 = ca1.*255./(max(max(abs(ca1))));
ch1 = ch1.*255./(max(max(abs(ch1))));
cv1 = cv1.*255./(max(max(abs(cv1))));
cd1 = cd1.*255./(max(max(abs(cd1))));

nimg(1:X+1, 1:Y+1) = 255;

nimg(1:X/2,1:Y/2)=abs(ca1);
nimg(1:X/2,Y/2+2:Y+1)=abs(cv1);
nimg(X/2+2:X+1,1:Y/2)=abs(ch1);
nimg(X/2+2:X+1,Y/2+2:Y+1)=abs(cd1);

figure;
imagesc(d1);
setg256;
axis off;
axis image;
hold on;
line([.75,.75],[.75,56.25],'LineWidth',3,'Color','blue');
line([.75,56.25],[.75,.75],'LineWidth',3,'Color','blue');
line([.75,56.25],[56.25,56.25],'LineWidth',3,'Color','blue');
line([56.25,56.25],[.75,56.25],'LineWidth',3,'Color','blue');

figure;
imagesc(nimg);
setg256;
axis off;
axis image;

nimgd = nimg(1:X/2, 1:Y/2);
nimg1 = sqrt(cv1.^2 + ch1.^2);

nimgd = nimgd.*255./(max(max(nimgd)));
nimg1 = nimg1.*255./(max(max(nimg1)));

nnimg(1:X+1, 1:Y/2) = 255;

nnimg(1:X/2,1:Y/2) = nimgd;
nnimg(X/2+2:X+1, 1:Y/2) = nimg1;

figure;
imagesc(nnimg);
setg256;
axis off;
axis image;

figure;
imagesc(nnimg);
setg256;
axis off;
axis image;
hold on;

rectangle('Curvature', [1 1], 'Position', [X/4-4.5,Y/4-4.5,9,9],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-9,Y/4-9,18,18],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-14,Y/4-14,28,28],'EdgeColor','white','LineWidth',3);

rectangle('Curvature', [1 1], 'Position', [X/4-4.5,Y/2+1+Y/4-4.5,9,9],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-9,Y/2+1+Y/4-9,18,18],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-14,Y/2+1+Y/4-14,28,28],'EdgeColor','white','LineWidth',3);

figure(1)
set(gcf, 'Color', 'white');
print -dpng show_wa_11.png
figure(2)

set(gcf, 'Color', 'white');
print -dpng show_wa_12.png
figure(3)
set(gcf, 'Color', 'white');
print -dpng show_wa_13.png
figure(4)
set(gcf, 'Color', 'white');
print -dpng show_wa_14.png
figure(5)
set(gcf, 'Color', 'white');
%text(30,Y/4,'{\bf v}_D','FontAngle', 'italic','FontSize',24);
%text(30,Y/4*3,'{\bf v}_1','FontAngle', 'italic','FontSize',24);
text(30,Y/2,'{\bf v}_w^{(1)}','FontAngle', 'italic','FontSize',24);
print -dpng show_wa_15.png


[X,Y] = size(d2);

[ca1, ch1, cv1, cd1] = dwt2(d2, 'haar');

ca1 = ca1.*255./(max(max(abs(ca1))));
ch1 = ch1.*255./(max(max(abs(ch1))));
cv1 = cv1.*255./(max(max(abs(cv1))));
cd1 = cd1.*255./(max(max(abs(cd1))));

nimg(1:X+1, 1:Y+1) = 255;

nimg(1:X/2,1:Y/2)=abs(ca1);
nimg(1:X/2,Y/2+2:Y+1)=abs(cv1);
nimg(X/2+2:X+1,1:Y/2)=abs(ch1);
nimg(X/2+2:X+1,Y/2+2:Y+1)=abs(cd1);

figure;
imagesc(d2);
setg256;
axis off;
axis image;
hold on;
line([.75,.75],[.75,56.25],'LineWidth',3,'Color','red');
line([.75,56.25],[.75,.75],'LineWidth',3,'Color','red');
line([.75,56.25],[56.25,56.25],'LineWidth',3,'Color','red');
line([56.25,56.25],[.75,56.25],'LineWidth',3,'Color','red');

figure;
imagesc(nimg);
setg256;
axis off;
axis image;

nimgd = nimg(1:X/2, 1:Y/2);
nimg1 = sqrt(cv1.^2 + ch1.^2);

nimgd = nimgd.*255./(max(max(nimgd)));
nimg1 = nimg1.*255./(max(max(nimg1)));

nnimg(1:X+1, 1:Y/2) = 255;

nnimg(1:X/2,1:Y/2) = nimgd;
nnimg(X/2+2:X+1, 1:Y/2) = nimg1;

figure;
imagesc(nnimg);
setg256;
axis off;
axis image;

figure;
imagesc(nnimg);
setg256;
axis off;
axis image;
hold on;

rectangle('Curvature', [1 1], 'Position', [X/4-4.5,Y/4-4.5,9,9],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-9,Y/4-9,18,18],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-14,Y/4-14,28,28],'EdgeColor','white','LineWidth',3);

rectangle('Curvature', [1 1], 'Position', [X/4-4.5,Y/2+1+Y/4-4.5,9,9],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-9,Y/2+1+Y/4-9,18,18],'EdgeColor','white','LineWidth',3);
rectangle('Curvature', [1 1], 'Position', [X/4-14,Y/2+1+Y/4-14,28,28],'EdgeColor','white','LineWidth',3);

figure(6)
set(gcf, 'Color', 'white');
print -dpng show_wa_16.png
figure(7)
set(gcf, 'Color', 'white');
print -dpng show_wa_17.png
figure(8)
set(gcf, 'Color', 'white');
print -dpng show_wa_18.png
figure(9)
set(gcf, 'Color', 'white');
%text(30,Y/4,'{\bf v}_D','FontAngle', 'italic','FontSize',24);
%text(30,Y/4*3,'{\bf v}_1','FontAngle', 'italic','FontSize',24);
text(30,Y/2,'{\bf v}_w^{(0)}','FontAngle', 'italic','FontSize',24);
print  -dpng show_wa_19.png

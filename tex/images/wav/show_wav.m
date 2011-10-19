clear all
close all

%img = imread('show_wa_01.bmp','bmp');
img = imread('mr.png');

xx = 90;
yy = 90;

d = img(xx-56:xx+55,yy-56:yy+55);

img = double(img)*.8;

figure;
imagesc(img);
setg256;
axis off;
axis image;
hold on;

plot(xx,yy,'ow','MarkerSize',4);

delta = 16;

for xi=-3:4
    for yi=-3:4
        x=xi*delta+xx-delta/2;
        y=yi*delta+yy-delta/2;
        plot(x,y,'xw');
    end;
end;

delta = 4*delta;

line([xx+delta,xx+delta], [yy-delta,yy+delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx-delta], [yy-delta,yy+delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx+delta], [yy-delta,yy-delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx+delta], [yy+delta,yy+delta], 'LineWidth', 2, 'Color', 'White');

delta = 8;

for xi=-3:4
    for yi=-3:4
        x=xi*delta+xx-delta/2;
        y=yi*delta+yy-delta/2;
        plot(x,y,'ow','MarkerSize',4);
    end;
end;

delta = 4*delta;

line([xx+delta,xx+delta], [yy-delta,yy+delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx-delta], [yy-delta,yy+delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx+delta], [yy-delta,yy-delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx+delta], [yy+delta,yy+delta], 'LineWidth', 2, 'Color', 'White');

delta = 4;

for xi=-3:4
    for yi=-3:4
        x=xi*delta+xx-delta/2;
        y=yi*delta+yy-delta/2;
        plot(x,y,'.w');
    end;
end;

delta = 4*delta;

line([xx+delta,xx+delta], [yy-delta,yy+delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx-delta], [yy-delta,yy+delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx+delta], [yy-delta,yy-delta], 'LineWidth', 2, 'Color', 'White');
line([xx-delta,xx+delta], [yy+delta,yy+delta], 'LineWidth', 2, 'Color', 'White');

[X,Y] = size(d);

[ca1, ch1, cv1, cd1] = dwt2(d, 'haar');

[ca2, ch2, cv2, cd2] = dwt2(ca1, 'haar');

[ca3, ch3, cv3, cd3] = dwt2(ca2, 'haar');

ca1 = ca1.*255./(max(max(abs(ca1))));
ch1 = ch1.*255./(max(max(abs(ch1))));
cv1 = cv1.*255./(max(max(abs(cv1))));
cd1 = cd1.*255./(max(max(abs(cd1))));

ca2 = ca2.*255./(max(max(abs(ca2))));
ch2 = ch2.*255./(max(max(abs(ch2))));
cv2 = cv2.*255./(max(max(abs(cv2))));
cd2 = cd2.*255./(max(max(abs(cd2))));

ca3 = ca3.*255./(max(max(abs(ca3))));
ch3 = ch3.*255./(max(max(abs(ch3))));
cv3 = cv3.*255./(max(max(abs(cv3))));
cd3 = cd3.*255./(max(max(abs(cd3))));

nimg(1:X+3, 1:Y+3) = 255;

nimg(1:X/8,1:Y/8)=abs(ca3);
nimg(1:X/8,Y/8+2:Y/4+1)=abs(cv3);
nimg(X/8+2:X/4+1,1:Y/8)=abs(ch3);
nimg(X/8+2:X/4+1,Y/8+2:Y/4+1)=abs(cd3);

nimg(1:X/4,Y/4+3:Y/2+2)=abs(cv2);
nimg(X/4+3:X/2+2,1:Y/4)=abs(ch2);
nimg(X/4+3:X/2+2,Y/4+3:Y/2+2)=abs(cd2);

nimg(1:X/2,Y/2+4:Y+3)=abs(cv1);
nimg(X/2+4:X+3,1:Y/2)=abs(ch1);
nimg(X/2+4:X+3,Y/2+4:Y+3)=abs(cd1);

n1 = nimg(1:X/8,1:Y/8);
n2 = sqrt(cv3.^2+ch3.^2);

%d(55:57,55:57)=255;
%%figure;
%imagesc(d);
%setg256;
%axis off;
%axis image;

figure;
imagesc(nimg);
setg256;
axis off;
axis image;

xx = 1;
yy = 1;
delta = X/8 - 1;

line([xx,xx],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx+delta, xx+delta],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy,yy],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy+delta,yy+delta],'LineWidth', 2, 'Color', 'White');

xx = X/8+2;
yy = 1;
delta = X/8 - 1;

line([xx,xx],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx+delta, xx+delta],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy,yy],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy+delta,yy+delta],'LineWidth', 2, 'Color', 'White');

xx = 1;
yy = Y/8 + 2;
delta = X/8 - 1;

line([xx,xx],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx+delta, xx+delta],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy,yy],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy+delta,yy+delta],'LineWidth', 2, 'Color', 'White');

xx = X/4+3+X/16;
yy = 1+X/16;
delta = X/8 - 1;

line([xx,xx],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx+delta, xx+delta],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy,yy],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy+delta,yy+delta],'LineWidth', 2, 'Color', 'White');

n31 = nimg(xx:xx+delta, yy:yy+delta);

yy = Y/4+3+Y/16;
xx = 1+Y/16;
delta = X/8 - 1;

line([xx,xx],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx+delta, xx+delta],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy,yy],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy+delta,yy+delta],'LineWidth', 2, 'Color', 'White');

n32 = nimg(xx:xx+delta, yy:yy+delta);

n3 = sqrt(n31.^2+n32.^2);

xx = X/2+4+X/4-X/16;
yy = 1+X/4-X/16;
delta = X/8 - 1;

line([xx,xx],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx+delta, xx+delta],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy,yy],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy+delta,yy+delta],'LineWidth', 2, 'Color', 'White');

n41 = nimg(xx:xx+delta, yy:yy+delta);

yy = Y/2+4+Y/4-Y/16;
xx = 1+Y/4-Y/16;
delta = X/8 - 1;

line([xx,xx],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx+delta, xx+delta],[yy,yy+delta],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy,yy],'LineWidth', 2, 'Color', 'White');
line([xx,xx+delta],[yy+delta,yy+delta],'LineWidth', 2, 'Color', 'White');

n42 = nimg(xx:xx+delta, yy:yy+delta);
n4 = sqrt(n41.^2+n42.^2);

figure(1);
set(gcf, 'Color', 'white');
figure(2);
set(gcf, 'Color', 'white');

figure;
imagesc(n1);
setg256;
axis off;
axis image;

figure;
imagesc(n2);
setg256;
axis off;
axis image;

figure;
imagesc(n3);
setg256;
axis off;
axis image;

figure;
imagesc(n4);
setg256;
axis off;
axis image;

figure(1)
print -dpng show_wav_fig1.png
figure(2)
print -dpng show_wav_fig2.png
figure(3)
rectangle('Curvature', [1 1], 'Position', [X/16-2.5,Y/16-2.5,6,6],'EdgeColor','white','LineWidth',6);
rectangle('Curvature', [1 1], 'Position', [X/16-5.5,Y/16-5.5,12,12],'EdgeColor','white','LineWidth',6);
set(gcf, 'Color', 'white');
print -dpng show_wav_fig3.png
figure(4)
rectangle('Curvature', [1 1], 'Position', [X/16-2.5,Y/16-2.5,6,6],'EdgeColor','white','LineWidth',6);
rectangle('Curvature', [1 1], 'Position', [X/16-5.5,Y/16-5.5,12,12],'EdgeColor','white','LineWidth',6);
set(gcf, 'Color', 'white');
print -dpng show_wav_fig4.png
figure(5)
rectangle('Curvature', [1 1], 'Position', [X/16-2.5,Y/16-2.5,6,6],'EdgeColor','white','LineWidth',6);
rectangle('Curvature', [1 1], 'Position', [X/16-5.5,Y/16-5.5,12,12],'EdgeColor','white','LineWidth',6);
set(gcf, 'Color', 'white');
print -dpng show_wav_fig5.png
figure(6)
rectangle('Curvature', [1 1], 'Position', [X/16-2.5,Y/16-2.5,6,6],'EdgeColor','white','LineWidth',6);
rectangle('Curvature', [1 1], 'Position', [X/16-5.5,Y/16-5.5,12,12],'EdgeColor','white','LineWidth',6);
set(gcf, 'Color', 'white');
print -dpng show_wav_fig6.png
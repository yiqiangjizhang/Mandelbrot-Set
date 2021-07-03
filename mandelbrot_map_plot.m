clc
clear all
close all

% Read csv file space separated
data = dlmread('output.csv');
x_size = size(data,2);
y_size = size(data,1);

x_vector = linspace(-2,1,x_size);
y_vector = linspace(-1.5,1.5,y_size);
img = pcolor(x_vector,y_vector,data);
set(img, 'EdgeColor', 'none');
colormap jet
colorbar;

% saveas(img,'image.png');
% imwrite(data, 'Mandelbrot_BW.png')
% print(gcf,'Mandelbrot_JET.png','-dpng','-r2000');         % -r dpi
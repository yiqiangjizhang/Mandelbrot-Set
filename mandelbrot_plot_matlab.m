clc
clear all
close all

% Read csv file space separated
data = dlmread('output.csv');
x = data(:,1);
y = data(:,2);


h = figure(1)
plot(x,y,'.')

colormap jet

print(gcf,'Mandelbrot_1D.png','-dpng','-r2000');         % -r dpi
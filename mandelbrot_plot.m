clc
clear all
close all

% Read csv file space separated
data = dlmread('output.csv');
x = data(:,1);
y = data(:,2);

plot(x,y,'.')
% plot(x, y, 'Marker', '.', 'MarkerSize', 1, ...
%      'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'1);
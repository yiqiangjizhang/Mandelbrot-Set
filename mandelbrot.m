%% Mandelbrot set 

%-------------------------------------------------------------------------%
% Mandelbrot set
%-------------------------------------------------------------------------%

% Date: 09/03/2021
% Author/s: Group 1
% Subject: High Performance Computing
% Professor: Manel Soria & Arnau Miro

% Clear workspace, command window and close windows
clear;
close all;
clc;
% Set interpreter to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%%

% X plane (real(c))
x0 = -2; % x inferior limit
xf = 1; % x superior limit
dx = 0.001; % Step size

% X plane (real(c))
y0 = -1.5;
yf = 1.5;
dy = 0.001; % Step size

% Vectors
x = x0:dx:xf;
y = y0:dy:yf;

% Create a 2D mesh
[X,Y] = meshgrid(x,y);

% Parameters
n_it = 100; % Number of iterations
R=10; % If abs(z)>R then z is considered as infty

% Complex number of Mandelbrot set
c = X + 1i*Y;

% Z function of Mandelbrot set
z = zeros(size(c));

% Mandelbrot equation
mandelbrot = @(z,c) z.^2 + c;


I=zeros(size(c));
% Iterate Mandelbrot set
for i=1:n_it
    z = mandelbrot(z,c);
    % bw=abs(z)<R;
    bw = abs(z)>2 & I==0;
    I(bw)=n_it - i;
end


% I=abs(z)<R; % logical image
imagesc(x,y,I);
%set(gca,'YDir','normal');
xlabel('Re(c)');
ylabel('Im(c)');
colormap jet
axis square






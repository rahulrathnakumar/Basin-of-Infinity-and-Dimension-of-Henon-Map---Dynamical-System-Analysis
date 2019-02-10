%% CODE TO ILLUSTRATE THE BASIN BOUNDARIES OF THE HENON MAP
%NAME: RAHUL RATHNAKUMAR
%EMAIL: rrathnak@asu.edu
%% Clearing all previous values
close all;
clear all;
clc;
%% Setting up the grid for computation
tic; %start timer
digits=64; % Setting 64 digits of precision
param = [2.12,-0.3,100,100]; % a, b, maxiter, L
%Setting the bounds and resolution for the grid
xmin = -3;
ymin = -3;
xmax = 3;
ymax = 3;
n = 1000; %divisions along X axis
m= 1000; %divisions along Y axis
%X is initially a (n*m) x 2 zero vector
x=zeros(n*m,2);
%Temporary column vector to store the values of yn while iterating
t=zeros(n*m,1); 
%Create the column 2 (Y-coordinates) of x by repeating the first
%element n times, and then proceed to the next subdivision (repeat this m
%times)
x(:,2)=repelem(linspace(ymin,ymax,m),n); 
%Create the column 1 (X-coordinates) of x by sequentially stacking the
%elements from xmin to xmax (n subdivisions) and repeat the sequence m
%times
x(:,1)=repmat(linspace(xmin,xmax,n),1,m);
%% Applying the Henon map formulation
    f=zeros(n*m,2); %Matrix to be used for iteration
    for k = 1:param(3) 
            %Check for elements that exceed the bounding value of L: if
            %such a set of row(s) is found, then they are made empty parallelly, thereby reducing
            %the number of rows left in x and f.
            bounds=find(abs(f(:,1))>100 | abs(f(:,2))>100);
            % map k to the corresponding (x0,y0)
            x(bounds,:)=0;
            f(bounds,:)=0;
      if(k == 1) %first iteration operates on x
            t=x(:,2); %store the initial value yn in t
            f(:,2)=x(:,1); %computation of y1
            f(:,1)=param(1)-x(:,1).^2+(param(2).*t(:,1)); %computation of y2
      else %for higher iterations, the operation happens recursively on f
            t=f(:,2);
            f(:,2)=f(:,1);
            f(:,1)=param(1)-f(:,1).^2+(param(2).*t(:,1));
      end
    end
toc; %timer ends
%% Plotting the results
plot(x(:,1),x(:,2),"."); %plots the values of X that remain bounded
%% References:
% [1] MATLAB Documentation
% [2] M. Hénon, Comm. Math. Phys. 50, 69 (1976)
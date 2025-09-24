% Part C: 2D Tensor Product Bases and Elevation Plots
% 2D Lagrange Bases

clear; clc; close all;
a = 0; b = 10;
nElem = 10;
p = 2;

endpoints = linspace(a,b,nElem+1);
midpoints = (endpoints(1:end-1)+endpoints(2:end))/2;
nodes1D  = sort([endpoints midpoints]);
nBasis1D = length(nodes1D)

% evaluation grid
xx = linspace(a,b,1000);
yy = linspace(a,b,1000);
[X,Y] = meshgrid(xx,yy);

% basis functions
phiX = zeros(nBasis1D,length(xx));
phiY = zeros(nBasis1D,length(yy));

for i = 1:nBasis1D
    for j = 1:length(xx)
        phiX(i,j) = lagrange_basis(i,xx(j),nodes1D);
    end
    for j = 1:length(yy)
        phiY(i,j) = lagrange_basis(i,yy(j),nodes1D);
    end
end

% Corner basis: NOde at (0,0) which is i=1, j=1
Phi_corner = phiX(1,:)'*phiY(1,:);

% Edge basis: node at (5,0)
[~, iMid] = min(abs(nodes1D-5));
Phi_edge = phiX(iMid,:)'*phiY(1,:);

% Interior basis: node at (5,5)
[~,jMid] = min(abs(nodes1D-5));
Phi_interior = phiX(iMid,:)'*phiY(jMid,:);




% Plot!!
figure;
surf(X,Y,Phi_corner,EdgeColor="none");
title("Corner Basis (0,0)"); 
xlabel("x"); 
ylabel("y");
zlabel("\Phi");
view(45,45);

figure;
surf(X,Y,Phi_edge,EdgeColor="none");
title("Edge Basis (5,0)"); 
xlabel("x"); 
ylabel("y");
zlabel("\Phi");
view(45,45);

figure;
surf(X,Y,Phi_interior,EdgeColor="none");
title("Interior Basis (5,5)"); 
xlabel("x"); 
ylabel("y");
zlabel("\Phi");
view(45,45);

% Lagrange Basis evaluator
function val = lagrange_basis(i,x,nodes)
    % evaluate the ith 1D Lagrange basis function at point x
    xi = nodes(i);
    others = nodes([1:i-1,i+1:end]);
    val = prod((x-others)./(xi-others));
end
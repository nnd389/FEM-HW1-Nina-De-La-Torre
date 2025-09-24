clear; clc; close all; 

a = 0; b = 10;
nElem = 10;
p = 2;

% Build open uniform knot vector
Xi = [zeros(1,p+1), 1:nElem-1, nElem, nElem*ones(1,p+1)];
nBasis1D = length(Xi)-p-1; % number of basis functions

% Evaluation grid
xx = linspace(a,b,80);
yy = linspace(a,b,80);
[X,Y] = meshgrid(xx,yy);

%Compute 1D B-spline basis values
phiX = zeros(nBasis1D, length(xx));
phiY = zeros(nBasis1D, length(yy));
for i = 1:nBasis1D
    for j = 1:length(xx)
        phiX(i,j) = CoxDeBoor(i,p,xx(j),Xi);
    end
    for j = 1:length(yy)
        phiY(i,j) = CoxDeBoor(i,p,yy(j),Xi);
    end
end

% Corner basis: support touches (0,0)
Phi_corner = phiX(1,:)'*phiY(1,:);

% Edge basis: centered near (5,0)
[~,iMid] = min(abs((0:nBasis1D-1)-nBasis1D/2)); % index near middle
Phi_edge = phiX(iMid,:)'*phiY(1,:);

% Interior basis: centered near (5,5)
Phi_interior = phiX(iMid,:)'*phiY(iMid,:);

% Plot!
figure;
surf(X,Y,Phi_corner,EdgeColor="none");
title("Corner B-spline basis touching (0,0)");
xlabel("x"); ylabel("y"); zlabel("B");

figure;
surf(X,Y,Phi_edge,EdgeColor="none");
title("Edge B-spline basis along a boundary edge");
xlabel("x"); ylabel("y"); zlabel("B");

figure;
surf(X,Y,Phi_interior,EdgeColor="none");
title("Interior B-spline basis near (5,5)");
xlabel("x"); ylabel("y"); zlabel("B");

% Cox-de Boor recursion
function val = CoxDeBoor(i,p,x,Xi)
    if p==0
        if (Xi(i) <= x && x < Xi(i+1)) || (x==Xi(end) && Xi(i) <= x && x <= Xi(i+1))
            val = 1;
        else
            val = 0;
        end
    else
        denom1 = Xi(i+p)-Xi(i);
        denom2 = Xi(i+p+1)-Xi(i+1);
        term1 = 0;
        term2 = 0;
        if denom1 > 0
            term1 = (x-Xi(i))/denom1*CoxDeBoor(i,p-1,x,Xi);
        end
        if denom2 > 0
            term2 = (Xi(i+p+1)-x)/denom2*CoxDeBoor(i+1,p-1,x,Xi);
        end
        val = term1 + term2;
    end
end



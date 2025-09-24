% Part B: Maximally SMooth B-Spline
% Plot all N_{i,p}(x) over [0.10] and verify partition of unity

clear; close all; clc; 

a = 0; b = 10; 
numElems = 10; 
intKnots = 1:(numElems - 1); % interior knots 1...9
xplot = linspace(a,b,2000); % Consider changing this to 1000?

for p = 1:3
    % build uniform knot vector captial Xi
    Xi = [zeros(1,p+1), intKnots, b*ones(1,p+1)];
    m = length(Xi)-1; % m = number of knot intervals
    nBasis = m-p; % number of basis functions
    fprintf("p = %d: length(Xi) = %d, m = %d, nBasis = %d \n", p, length(Xi), m, nBasis);

    % Cox-de Boor!!! (iterative)
    nx = length(xplot);
    N = zeros(m, nx); % initialize basis functions
    % use m instead of nBasis here^
    for i = 1:m
        N(i,:) = (Xi(i) <= xplot) & (xplot < Xi(i+1));
    end
    N(m, xplot==Xi(end)) = 1; % right endpoint

    for d = 1:p % build up to degree p
        for i = 1:(m-d)
            % left term
            denom1 = Xi(i+d)-Xi(i);
            if denom1 == 0
                term1 = zeros(1,nx);
            else
                term1 = ((xplot - Xi(i))/denom1).*N(i,:);
            end

            %right term
            denom2 = Xi(i+d+1) - Xi(i+1);
            if denom2 == 0
                term2 = zeros(1,nx);
            else
                term2 = ((Xi(i+d+1) - xplot) / denom2) .* N(i+1,:);
            end
            N(i,:) = term1 + term2; % update the basis function
        end
    end
    
    % Final B-spline basis functions are stores in rows 1...nBasis of N
    B = N(1:nBasis, :);

    % Plot!
    figure(Name=sprintf("B-Spline basis p=%d",p), NumberTitle="off");
    hold on;
    for i = 1:nBasis
        plot(xplot, B(i,:), LineWidth=1.2);
    end
    title(sprintf("Open uniform B-spline basis functions for p=%d",p));
    xlabel("x"); 
    ylabel(sprintf("N_{i,%d}(x)",p));
    grid on;

    % Check Partition of Unity
    % Verify partition of unity
    partitionUnity = sum(B, 1); % Sum of basis functions at each x
    fprintf("max sum is %.3e\n", max(partitionUnity)); % should be 1
    plot(xplot,partitionUnity,'r--',LineWidth=2);
end
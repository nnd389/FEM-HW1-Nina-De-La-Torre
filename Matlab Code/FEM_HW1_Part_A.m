a = 0; b = 10; N = 10;
h = (b-a)/N;

for p = 1:3
    if p==1
        global_nodes = linspace(a,b,N+1); % [0,1]
        numBasis = length(global_nodes) % p=1 has 11 bases
    elseif p==2
        endpoints = linspace(a,b,N+1);
        midpoints = (endpoints(1:end-1)+endpoints(2:end))/2;
        global_nodes = sort([endpoints midpoints]); % [0,.5,1]
        numBasis = length(global_nodes);
    elseif p==3
        endpoints = linspace(a,b,N+1);
        internal_nodes = [];
        for e = 1:N
            xL = endpoints(e);
            xR = endpoints(e+1);
            internal_nodes = [internal_nodes, xL + (xR-xL)/3, xL + 2*(xR-xL)/3];
        end
        global_nodes = sort([endpoints internal_nodes]);
        numBasis = length(global_nodes);
    end

    xx = linspace(a,b,1000);
    phi = zeros(numBasis, length(xx));

    for i = 1:numBasis
        for j = 1:length(xx)
            x = xx(j);
            val = 0;
            for e = 1:N
                xL = a + (e-1)*h;
                xR = a + e*h;
                if x >= xL && x <= xR
                    xi = (x-xL)/h; % map to [0,1]
                    if p==1
                        % Linear shape functions
                        L1 = 1-xi;
                        L2 = xi;
                        local_vals = [L1 L2];
                        elem_nodes = [e,e+1];
                    elseif p == 2
                        % quadratic shape functions
                        L1 = 2*(xi-0.5).*(xi-1);
                        L2 = 4*xi.*(1-xi);
                        L3 = 2*xi.*(xi-0.5);
                        local_vals = [L1 L2 L3];
                        elem_nodes = [2*(e-1)+1, 2*(e-1)+2, 2*(e-1)+3];
                    elseif p == 3
                        % cubic shape functions
                        L1 = -4.5*(xi-1/3).*(xi-2/3).*(xi-1);
                        L2 = 13.5*xi.*(xi-2/3).*(xi-1);
                        L3 = -13.5*xi.*(xi-1/3).*(xi-1);
                        L4 = 4.5*xi.*(xi-1/3).*(xi-2/3);
                        local_vals = [L1 L2 L3 L4];
                        elem_nodes = [3*(e-1)+1, 3*(e-1)+2, 3*(e-1)+3, 3*(e-1)+4];
                    end
                    % pick out the correct basis value
                    for k = 1:length(elem_nodes)
                        if elem_nodes(k) == i
                            val = local_vals(k);
                        end
                    end
                end
            end
            phi(i,j) = val;
        end
    end

    % Plot
    figure; hold on;
    for i = 1:numBasis
        plot(xx,phi(i,:));
    end
    title(sprintf("Global Basis Functions for p = %d", p));
    grid on;

    % Part E: check the partition of Unity
    sum_phi = sum(phi,1);
    plot(xx,sum_phi,'r--',LineWidth=2);


    % Part F
    % Confirm Local Support
    for i = 1:numBasis
        support = xx(phi(i,:) > 1e-12); % small tolerance
        fprintf("Basis %d is nonzero on [%.2f, %.2f]\n",i,min(support),max(support));
    end

    % Confirm Non-Negativity
    min_val = min(phi(:));
    if min_val >=0
        disp("ALL basis functions are non-negative!");
    else
        disp("WARNING: some basis functions are negative! AAH!")% are quadratic and cubis basis functions allowed to be negative? the internet says yes?
    end
end



% Plot Lebesgue functions for equally spaced points
% and for Chebyshev points.
for part=1:2 % Equally spaced nodes for part=1, Cheby nodes for part=2.
    for size=1:3
        if size==1, n = 5; end
        if size==2, n = 10; end
        if size==3, n = 20; end
        x = zeros(n+1,1);
        if part==1
            for i=0:n, x(i+1) = -1 + (2*i)/n; end % Equally spaced nodes
        elseif part==2
            for i=0:n, x(i+1) = cos((2*i+1)*pi/(2*n+2)); end % Chebyshev nodes
        end
        pts = x(1); % These are the pts where Lebesgue function will be evaluated.
        for i=1:n
            pts = [pts; [x(i) + ((x(i+1)-x(i))/20)*[1:20]']];
        end

        Lebesguefun = zeros(length(pts),1);
        for i=1:n+1
            Lebesguefun = Lebesguefun + abs(lsubi(i,x,pts));
        end
        figure(part)
        if part==1
            plot(pts,log10(Lebesguefun),'-','LineWidth',2); hold on; shg, pause
        else
            plot(pts,Lebesguefun,'-','LineWidth',2); hold on; shg, pause
        end
    end
    if part==1
        title('Log10 of Lebesgue functions for equally spaced nodes')
    else
        title('Lebesgue functions for Chebyshev nodes')
    end
end

function lsubivals = lsubi(i,nodes,pts)
    % Computes l_i for the given nodes and evaluates it at pts.
    nnodes = length(nodes);
    npts = length(pts);
    lsubivals = ones(npts,1);
    for j=1:nnodes
        if j~=i
            lsubivals = lsubivals.*((pts - nodes(j))./(nodes(i)-nodes(j)));
        end
    end
end
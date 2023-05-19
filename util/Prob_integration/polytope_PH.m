function [PA, Pb, C] = polytope_PH(Z)
% polytope - Converts a zonotope to a H-representation polytope
% (only for 2D or 3D)
% this function is based off of the CORA_2018 function which is at
% .../CORA_2018/contSet/@zonotope/polytope.m
%
%obtain number of generators, dimensions
c = Z(:, 1);
G = Z(:, 2:end);

[dim,nrOfGenerators]=size(G);

if dim > 1
    comb=combinator(nrOfGenerators,dim-1,'c');

    %build C matrices
    C=[];
    if dim == 2
        C = G;
        %... get perpendicular vector
        C = [-C(2, :); C(1, :)];
    elseif dim == 3
        % basically a cross product in matrix form
        Q = [G(:, comb(:, 1)); G(:, comb(:, 2))];
        C = [Q(2, :).*Q(6, :) - Q(3, :).*Q(5, :); -(Q(1, :).*Q(6, :) - Q(3, :).*Q(4, :)); Q(1, :).*Q(5, :) - Q(2, :).*Q(4, :)];
    else
        error('Dimension not supported.');
    end
    %remove NaN rows due to rank deficiency
    index = find(sum(isnan(C),2));
    if ~isempty(index)
        C(index,:) = [];
    end
else
    C = G;
end

% normalize normal vectors
C = (C./sqrt(sum(C.^2)))';

%build d vector
deltaD = sum(abs((C*G)'))';

%compute dPos, dNeg
d = C*c;

PA = [C; -C];
Pb = [d+deltaD; -d+deltaD];

%------------- END OF CODE --------------



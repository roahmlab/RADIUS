function [Zred]=deleteAligned_noslice(Z,slice_dim)
% deleteAligned - combines aligned generators
%
% Syntax:  
%    [Zred]=deleteAligned(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Zred - reduced zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Matthias Althoff
% Written: 15-January-2009
% Last update: ---
% Last revision: ---

% modified: JL
% update: do not delete sliceable generator

%------------- BEGIN CODE --------------


%get Z-matrix from zonotope Z
Zmatrix=get(Z,'Z');

%extract generator matrix
c=Zmatrix(:,1);
G=Zmatrix(:,2:end);

%Delete zero-generators
G=nonzeroFilter(G);


% JL added
if nargin < 2
    slice_dim = [];
    sliceg_idx = [];
    g_hold = [];
else
    % find sliceg_idx
    gtemp = sum(G(slice_dim,:));
    sliceg_idx = find(gtemp~=0);
    g_hold = G(:,sliceg_idx);
    G(:,sliceg_idx) = [];
end


%normalize generators
for i=1:length(G(1,:))
    G_norm(:,i) = G(:,i)/norm(G(:,i));
end

%find equal generators
i = 1;
while i < length(G(1,:))
    
    G_act = G_norm(:,i);
    ind = find(abs(G_act'*G_norm(:,(i+1):end)) > 1-1e-3);
    if ~isempty(ind)
        ind = ind+i;
        for iAdd = 1:length(ind)  
            %add generators
            G(:,i) = G(:,i) + sign(G_act'*G_norm(:,ind(iAdd)))*G(:,ind(iAdd));
        end
        for iAdd = 1:length(ind)
            %remove generators
            G(:,ind(iAdd)) = [];
            G_norm(:,ind(iAdd)) = [];
            %change ind to correct for cancellation
            ind = ind - 1;
        end
    end  
    %increase i
    i = i + 1;
end

Zred=zonotope([c,G,g_hold]);


function [G]=nonzeroFilter(G)
%Delete zero-generators
i=1;
while i<=length(G(1,:))
    if G(:,i)==0*G(:,i)
        G(:,i)=[];
    else
        i=i+1;
    end
end

%------------- END OF CODE --------------
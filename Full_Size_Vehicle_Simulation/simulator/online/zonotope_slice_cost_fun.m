function [newzono,slice_lambda,slice_G,slice_G_full] = zonotope_slice_cost_fun(zono,slice_dim, slice_pt)
% function for taking a slice of a zonotope along the corresponding
% dimensions...
% only works if there's one generator in each slice_dim

if size(slice_pt, 2) ~= 1
    error('Slice point should be a column vector');
end

Z = get(zono, 'Z');
c = Z(:, 1);
G = Z(:, 2:end);

slice_idx = [];
for i = 1:length(slice_dim)
   myidxs = find(G(slice_dim(i), :) ~= 0);
   if length(myidxs) ~= 1
       if isempty(myidxs)
           error('No generator for slice index');
       else
           error('More than one generator for slice index');
       end
   end
   slice_idx(i, 1) = myidxs;
end

slice_c = c(slice_dim, 1);
slice_G = G(slice_dim, slice_idx);
slice_G_full = G(:,slice_idx);
slice_lambda = slice_G\(slice_pt - slice_c);
if size(slice_lambda, 2) > 1
    error('slice_lambda is not 1D');
end
if any(abs(slice_lambda) > 1)
    
    if any(abs(slice_lambda) > 1.01)
        slice_lambda
        error('Slice point is outside bounds of reach set, and therefore is not verified');
    else
        replace_idx = abs(slice_lambda) > 0.999;
        slice_lambda(replace_idx) = 0.999*sign(slice_lambda(replace_idx));
    end
end

newG = G;
newG(:, slice_idx) = [];
newc = c + G(:, slice_idx)*slice_lambda;

newzono = zonotope([newc, newG]);

end


function [newzono] = zonotope_slice(zono,slice_dim, slice_pt)
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
       if length(myidxs) == 0
           slice_dim(i)
           error('No generator for slice index');
       else
           error('More than one generator for slice index');
       end
   end
   slice_idx(i, 1) = myidxs;
end

slice_c = c(slice_dim, 1);
slice_G = G(slice_dim, slice_idx);
slice_lambda = slice_G\(slice_pt - slice_c);
%basically to slice a range, you'll get a "slice_lambda_1" for the lower bound of the range, and "slice_lambda_2" for the upper bound
if size(slice_lambda, 2) > 1
    error('slice_lambda is not 1D');
end
if any(abs(slice_lambda) > 1)
     if any(abs(slice_lambda) > 1.02) % avoid numerical issue 
         error(num2str(slice_lambda(:)')+"Slice point is outside bounds of reach set, and therefore is not verified");
     end
    
    slice_lambda(slice_lambda > 1) = 1;
    slice_lambda(slice_lambda < -1) = -1;
end

newG = G;
newG(:, slice_idx) = [];
newc = c + G(:, slice_idx)*slice_lambda;

newzono = zonotope([newc, newG]);

end


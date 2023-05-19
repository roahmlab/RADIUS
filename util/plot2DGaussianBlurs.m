%Script to plot gaussian as faded blurs
%sets a range and plots the maximum weight for multiple gaussians

function h = plot2DGaussianBlurs(mu, Sigma,range,ptol,N, plot_flag)

if nargin<6
    plot_flag = 1;
end
if nargin<5
    N = 100;
if nargin<4
    ptol = 0.99;
end
end

x = linspace(range(1,1),range(1,2),N);
y = linspace(range(2,1),range(2,2),N);


[XX,YY] = meshgrid(x,y);
XXcol = XX(:);
YYcol = YY(:);
pts = [XXcol YYcol];

s = -2 * log(1 - ptol);

w_function = zeros(N);
L_alpha = false(1,N^2);

for i = 1:size(mu,2)
    w_function = max(w_function,2*pi*sqrt(det(Sigma{i}(1:2,1:2)))*reshape(mvnpdf([XXcol YYcol],mu(1:2,i)',Sigma{i}(1:2,1:2)),[N N]));
    xr = pts'-mu(1:2,i);
    Ltmp = arrayfun(@(j) xr(:,j)'*(Sigma{i}(1:2,1:2)\xr(:,j)) <= s ,1:N^2);
    L_alpha = L_alpha | Ltmp;
end

L_alpha = reshape(L_alpha,[N N]);

alpha = double(L_alpha);

if plot_flag
    h = imagesc(x,y,w_function,'AlphaData',alpha*0.45,'AlphaDataMapping','none'); % image with alpha
    colormap cool
else
    h = [];
end

end
% function D = distmat(xs,ys)
%
%  computes distance matrix D(i,j) = |x_i - y_j|
% 
%  Inputs
%   xs,ys Nx by 2 and Ny by 2 coordinate matrices
%
%  Output
%   D Nx by Ny distance matrix

% Fernando Guevara Vasquez 2009
function D = distmat(xs,ys)
% get sizes
Ny=size(ys,1); Nx=size(xs,1);
% compute norms squared of xs and ys
nxs = xs(:,1).^2 + xs(:,2).^2;
nys = ys(:,1).^2 + ys(:,2).^2;
% compute the norm of ys minus xr
D = sqrt(nxs*ones(1,Ny) -2*xs*ys' + ones(Nx,1)*nys');

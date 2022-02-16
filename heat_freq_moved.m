%Solution of the heat equation in the frequeuncy domain with sources moved
%to exterior points

%Inputs:
%   s: frequency
%   k0: thermal diffusivity
%   X1: x-coordinates grid points
%   X2: y-coordinates grid points
%   Y1: x-coordinate source
%   Y2: y-coordinate source
%   a: distance to move new sources
%   n_trunc: number to terms to truncate Graff's addition formula
%   n_move: number of new sources
%Outputs:
%   U: Array of fields at given freqs

function[U_f] = heat_freq_moved(s,k0,X1,X2,Y1,Y2,geo,a,n_trunc,n_move)
w = 1i*s;
ks = 1i*w/k0; %for ifft
k = sqrt(ks);
p = [X1(:),X2(:)];
%[pmc1, pmc2] = point_diff(p,geo.ctr);
[new_srcs, Y1m, Y2m] = move_src(a, geo, n_move); 
for j = 1:length(s)   
    [U, UDX1, UDX2] = field_point_src(geo.ctr(:,1),geo.ctr(:,2),Y1,Y2,k(j));
    F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2);
    %Apply SLP & DLP
    U1_moved = reshape(slp_move(F,p,geo, k(j),n_trunc,new_srcs,Y1m,Y2m),size(X1));
    U2_moved = reshape(dlp_move(U,p,geo, k(j),n_trunc,new_srcs,Y1m,Y2m),size(X1));
    U_rec_ext = U1_moved - U2_moved;
    U_f(:,j,:) = U_rec_ext/k0;
end
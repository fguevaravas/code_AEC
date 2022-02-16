%Solution of the helmholtz equation with sources moved

%Inputs:
%   k:wave number
%   X1: x-coordinates grid points
%   X2: y-coordinates grid points
%   Y1: x-coordinate source
%   Y2: y-coordinate source
%   a: distance to move new sources
%   n_trunc: number to terms to truncate Graff's addition formula
%   n_move: number of new sources
%Outputs:
%   U: Array of fields at given freqs

function[U_f] = helm_moved(k,X1,X2,Y1,Y2,geo,a,n_trunc,n_move)
%[pmc1, pmc2] = point_diff(p,geo.ctr);
p = [X1(:),X2(:)];
[new_srcs, Y1m, Y2m] = move_src(a, geo, n_move); 
[U, UDX1, UDX2] = field_point_src(geo.ctr(:,1),geo.ctr(:,2),Y1,Y2,k);
F = UDX1.*geo.nv(:,1) +UDX2.*geo.nv(:,2);
%Apply SLP & DLP
U1_moved = reshape(slp_move(F,p,geo, k,n_trunc,new_srcs,Y1m,Y2m),size(X1));
U2_moved = reshape(dlp_move(U,p,geo, k,n_trunc,new_srcs,Y1m,Y2m),size(X1));
U_f = U1_moved - U2_moved;
end
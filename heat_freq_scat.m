%Solution of the scattering problem in the frequeuncy domain

%Inputs:
%   s: frequency
%   k0: thermal diffusivity
%   X1: x-coordinates grid points
%   X2: y-coordinates grid points
%   geo: geometry of inclusion
%   dir_data: function to give dirichlet data at a given frequeuncy
%Outputs:
%   U: Array of fields at given freqs
function[U_f] = heat_freq_scat(s,k0,X1,X2,geo, dir_data)
w = 1i*s;
ks = 1i*w/k0; %for ifft
k = sqrt(ks);
p = [X1(:),X2(:)];
n = length(geo.ctr(:,1));
%Point differences
[pmc1, pmc2] = point_diff(p,geo.ctr); % grid points minus inclusion centers
dir = dir_data(s);
for j = 1:length(s)   
    if imag(k(j)) > 1/2
        %eta = 1/2;
        eta = imag(k(j));
    else
        eta = 1/2;
    end
    phi = Kapur_scat(geo,k(j),dir(:,j),eta);
    U1 = slp(phi,pmc1,pmc2,geo.h,k(j));
    U2 = dlp(phi,geo.nv, pmc1,pmc2,geo.h,k(j));
    U1 = reshape(U1,size(X1));
    U2 = reshape(U2,size(X1));
    U_rec_ext = U2-1i*eta*U1;
    U_f(:,j,:) = U_rec_ext/k0;
end
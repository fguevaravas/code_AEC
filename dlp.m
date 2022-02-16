%Double Layer Potential for the 2d Helmholtz equation 
%
%Input:
%   f: Dirichlet data
%   pmc1: x-coordinates of difference between points on surface and grid
%   points
%   pmc2: y-coordinates of difference between points on surface and grid
%   points
%   h: difference between points on surface
%   nv: normal vectors to surface
%   k: 1i*freq/thermal diffusivity
%Output:
%   U: single layer component of recreated field

function[U] = dlp(f,nv, pmc1,pmc2,h,k)
    R = sqrt(pmc1.^2+pmc2.^2);
    DG1 = k*(1i/4)*besselh(1,k*R).*pmc1./R; 
    DG2 = k*(1i/4)*besselh(1,k*R).*pmc2./R;
    %approximate integral with trapezoidal rule
    U = (DG1.*nv(:,1)'+DG2.*nv(:,2)')*(f.*h);
end%function
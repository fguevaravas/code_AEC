%Single Layer Potential for the 2d Helmholtz equation w/ imaginary
%frequency
%
%Input:
%   f: Neumann data
%   pmc1: x-coordinates of difference between points on surface and grid
%   points
%   pmc2: y-coordinates of difference between points on surface and grid
%   points
%   h: difference between points on surface
%   k: 1i*freq/thermal diffusivity
%Output:
%   U: single layer component of recreated field

function[U] = slp(f,pmc1,pmc2,h,k)
    R = sqrt(pmc1.^2+pmc2.^2);
    G = (1i/4)*besselh(0,k*R); %Green's function
    %approximate integral with trapezoidal rule
    U =G*(f.*h);
end%function
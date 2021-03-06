%Field generated by point source w/ imaginary frequency 
%
%
%Inputs:
%   X1: x-coordinate of points to evaluate
%   X2: y-coordinate of points to evaluate
%   Y1: x-coordinate of point source
%   Y2: y-coordinate of point source
%   k: 1i*freq/thermal diffusivity
%Ouputs:
%   U: field at given points
%   UDX1: x-component of gradient of field
%   UDX2: y-component of gradient of field
function[U, UDX1, UDX2] = field_point_src(X1,X2,Y1,Y2,k)
    R = sqrt((X1-Y1).^2+(X2-Y2).^2);
    U = (1i/4)*besselh(0,k.*R);
    UDX1 = -k.*(1i/4).*besselh(1,k.*R).*(X1-Y1)./R; 
    UDX2 = -k.*(1i/4).*besselh(1,k.*R).*(X2-Y2)./R;
end%function
%Solution of the heat equation in the Laplace domain 
%
%Input: 
%   s: frequency
%   k0: thermal diffusivity
%   X1: X-coordinates of points to eval
%   X2: Y-coordinates of points to eval
%   Y1: X-coordinate of point source
%   Y2: Y-coordinate of point source
%   res: reshape option (=0 matrix)
%   d: d = 0 temp; d=1 x-partial; d=2 y-partial;
%Ouputs:
%   U: temperature field or derivative
function[U] = heat_freq(s,k0,X1,X2,Y1,Y2,res,d)
w = 1i*s;
ks = 1i*w/k0; %for ifft
k = sqrt(ks);
if res == 0
    for j = 1:length(s)
        if d ==0
            [U_t, ~, ~] = field_point_src(X1,X2,Y1,Y2,k(j));
        elseif d==1
            [~, U_t, ~] = field_point_src(X1,X2,Y1,Y2,k(j));
        elseif d==2
            [~, ~, U_t] = field_point_src(X1,X2,Y1,Y2,k(j));
        end
        U(:,j,:) = U_t/k0;
    end
else
    [U, ~, ~] = field_point_src(X1,X2,Y1,Y2,k);
    U = U/k0;
end
end %function

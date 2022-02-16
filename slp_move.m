%Moves SLP to exterior points
%
%Note: it is possible to get rid of the for loops it would require defining
%bessel functions again
%
%Input:
%   f: Neumann data
%   p: coordinates of grid   
%   geo: geometry of cloaking region (including og src location)
%   k: 1i*freq/thermal diffusivity
%   N: truncation for Graffs addition formulae
%   new_src: new source location
%   Y1m: x-cord old sources col corresponds to new src loc
%   Y2m: y-cord old sources col corresponds to new src loc
%Output:
%   U: single layer component of recreated field

function[U] = slp_move(f,p,geo, k,N,new_src,Y1m,Y2m)
    X1 = p(:,1);
    X2 = p(:,2);
    i_tot = 1; %manual iteration
    U=zeros(size(X1)); %initialize solution
    for i = 1:length(new_src)
        Y1move = new_src(i,1);
        Y2move = new_src(i,2);
        Y1 = Y1m(:,i);
        Y2 = Y2m(:,i);
        %Move each source on th e segment to the new
        for j = 1:length(Y1)
            v = sqrt((X1-Y1move).^2+(X2-Y2move).^2);
            u = sqrt((Y1(j)-Y1move)^2+(Y2(j)-Y2move).^2);
            theta1 = atan2(X2-Y2move, X1-Y1move);
            theta2 = atan2(Y2(j)-Y2move, Y1(j)-Y1move);
            moved_field = zeros(size(X1));
            for ii = -N:N
                moved_field = moved_field+besselh(ii,1,k.*v)*besselj(ii,k.*u,1).*exp(1i*ii*(theta1- theta2));
            end
            G_n = 1i/4*moved_field.*exp(abs(imag(k.*u)));
            U = U+G_n*f(i_tot).*geo.h(i_tot);
            i_tot = i_tot+1;
        end
    end
end%function
%Moves DLP to exterior points
%
%Input:
%   f: Dirichlet data
%   p: coordinates of grid   
%   geo: geometry of cloaking region (including og src location)
%   k: 1i*freq/thermal diffusivity
%   N: truncation for Graffs addition formulae
%   new_src: new source location
%   Y1m: x-cord old sources col corresponds to new src loc
%   Y2m: y-cord old sources col corresponds to new src loc
%Output:
%   U: single layer component of recreated field

function[U] = dlp_move(f,p,geo, k,N,new_src,Y1m,Y2m)
    X1 = p(:,1);
    X2 = p(:,2);
    i_tot = 1; %manual iteration
    U=zeros(size(X1)); %initialize solution
    for i = 1:length(new_src)
        Y1move = new_src(i,1);
        Y2move = new_src(i,2);
        Y1 = Y1m(:,i);
        Y2 = Y2m(:,i);
        %Move each source on the segment to the new 
        for j = 1:length(Y1)
            %Graffs addition formula
            v = sqrt((X1-Y1move).^2+(X2-Y2move).^2);
            u = sqrt((Y1(j)-Y1move)^2+(Y2(j)-Y2move).^2);
            w = sqrt((X1-Y1(j)).^2+(X2-Y2(j)).^2);
            theta1 = atan2(X2-Y2move, X1-Y1move);
            theta2 = atan2(Y2(j)-Y2move, Y1(j)-Y1move);
            theta3 = atan2(Y2(j) - X2, Y1(j) - X1);
            theta4 = atan2(Y2move-X2, Y1move - X1);
            DG_move = zeros(size(X1));
            for ii = -N:N
                DG_move = DG_move+besselh(1+ii,1,k*v)*besselj(ii,k*u,1).*exp(1i*ii*(theta1- theta2)).*exp(-1i*(theta3 - theta4));
            end
            DG_move = DG_move.*exp(abs(imag(k.*u)));
            DG1_move = k*(1i/4)*DG_move.*(X1-Y1(j))./w;
            DG2_move = k*(1i/4)*DG_move.*(X2-Y2(j))./w;
            G_n = DG1_move*geo.nv(i_tot,1)+DG2_move*geo.nv(i_tot,2);
            %Integrate
            U = U+G_n*f(i_tot).*geo.h(i_tot);
            i_tot = i_tot+1;
        end    
    end
end%function
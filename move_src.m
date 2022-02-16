%Map sources on a circle to new source locations

%Inputs:
%   a: distance to move
%   geo: geometry of cloaking region
%   n_move: number of new points to move to
%Ouputs:
%   new_srcs: location of new sources
%   Y1m: x-coordinates of old sources, col corresponds to move location
%   Y2m: y-coordinates of old sources, col corresponds to move location
function[new_srcs, Y1m, Y2m] = move_src(a, geo, n_move)
    if rem(length(geo.ctr),n_move) ==0
        seg =  length(geo.ctr)/n_move; %number of points between segments
        Y1m = zeros(seg,n_move);
        Y2m = zeros(seg, n_move);
        new_srcs = zeros(n_move,2); %initalize vector for new source locations
        for i = 1:n_move
            %get exterior points by moving away from the circle in the normal
            %direction at equispaced points
            if rem(seg,2) == 0 %even case
                Y1loc = geo.ctr((i-1)*seg+seg/2,1)/2+geo.ctr((i-1)*seg+seg/2+1,1)/2;
                Y2loc = geo.ctr((i-1)*seg+seg/2,2)/2+geo.ctr((i-1)*seg+seg/2+1,2)/2;
                Ypert = [geo.Rcenter(1)-Y1loc, geo.Rcenter(2)-Y2loc];
                Ypert = Ypert/norm(Ypert);
                Y1move = Y1loc-a*Ypert(1);
                Y2move = Y2loc-a*Ypert(2);
            else %odd case
               seg_mp = (i-1)*seg+(seg+1)/2;
               Y1move = geo.ctr(seg_mp,1)+a*geo.nv(seg_mp,1);
               Y2move = geo.ctr(seg_mp,2)+a*geo.nv(seg_mp,2); 
            end
            new_srcs(i,:) = [Y1move,Y2move];
            Y1m(:,i) = geo.ctr((i-1)*seg+1:i*seg,1);
            Y2m(:,i) = geo.ctr((i-1)*seg+1:i*seg,2);         
        end
%             clf;
%             hold on;
%             plot(new_srcs(1,1),new_srcs(1,2),'rx')
%             plot(Y1m(:,1),Y2m(:,1), 'rx')
%              plot(new_srcs(2,1),new_srcs(2,2),'kx')
%             plot(Y1m(:,2),Y2m(:,2), 'kx')
%              plot(new_srcs(3,1),new_srcs(3,2),'rx')
%             plot(Y1m(:,3),Y2m(:,3), 'rx')
%              plot(new_srcs(4,1),new_srcs(4,2),'kx')
%             plot(Y1m(:,4),Y2m(:,4), 'kx')
%             hold off
%             axis square
            
    else
                disp('cp must be divisible by n_move');
    end
end%function
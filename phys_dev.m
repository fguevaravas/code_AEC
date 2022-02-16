%Calculate the physical size of the disc devices 
%Inputs:
%   U_rec: recreated solution from the frequeuncy domain
%   new_src_locs: location of new sources
%   maxt: maximum realizable temperature 
%   X1: y-cord grid points
%   X2: x-cord grid points
%   N1: discretization for x
%Ouput:
%   U_phys: physically realizable U_rec
%   all_rads: radius of disc devices
function[U_phys, all_rads] = phys_dev(U_rec,new_src_locs, maxt, X1, X2,N1)
    all_rads = zeros(1,length(new_src_locs)); %initialize vector to store radii
    r_max = norm(new_src_locs(1,:)-new_src_locs(2,:))/2; %maximum radius for any circle
    rad_step = ((max(max(X1))-min(min(X1)))+1)/N1/10; %grid point size for uniform grid
    U_phys = U_rec;
    for i = 1:length(new_src_locs)
        mask = (X1-new_src_locs(i,1)).^2 + (X2 - new_src_locs(i,2)).^2 < r_max.^2;
        current_max = max(max(abs(U_rec).*mask));
        rad = rad_step;
        while current_max>maxt
            inmask = (X1-new_src_locs(i,1)).^2 + (X2 - new_src_locs(i,2)).^2 > rad.^2;
            full_mask = inmask.*mask;
            current_max = max(max(abs(U_rec).*full_mask));
            rad = rad+rad_step;
        end
        all_rads(i) = rad;
        U_phys = U_phys.*inmask;
    end
   
end%function
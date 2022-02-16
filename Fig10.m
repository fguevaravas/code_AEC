%Exterior cloaking of a kite object 

%%%%Parameters
%Physical
k0 = 1.5; %thermal diffusivity
Y1 = 8; Y2 =5; %location of point source
max_temp_scale = 100; %factor to allow maximum temp of pelt devices

%NILT 
t_max = 2;
alfa =0;
M = 2^10; %number of time steps
dt = t_max/M;

%Reproduction
rep = 0; %0 for full cloaking

%Temperature function
temp = @(X1,X2,t) temp_source(X1,X2,t,0, Y1, Y2, 1, k0);


%%%%Geometry
%Grid 
L1 = 10; L2 = 10; % dimensions of domain
N1 = 200; N2 = 200; % number of points on grid
[X1,X2] = ndgrid(linspace(0,L1,N1),linspace(0,L2,N2));
p = [X1(:),X2(:)];

%Cloaking region
cp =256; % number of points on circle
Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/6; % center and radius
a = (sqrt(2)-1)*radius; %distance to new point sources 
geo = circ(cp, Rcenter, radius); 

%Inclusion
cp_k  = 512; %must be even
j = (0:(cp_k-1))';
k_p = 2*j*pi/cp_k;
geo_k = kite(k_p);

%Addition formula 
n_move = 4; %number of exterior points to use for cloaking 
n_trunc = 22; %ceil(15*(a+radius-radius/2)); %truncation for Graff's addition formula

%%%Reproduction
%Full field
fun = @(s) heat_freq_moved(s,k0,X1,X2,Y1,Y2,geo,a,n_trunc,n_move);
[ft, t] = fftilt(fun,t_max,M,alfa); %numerical inverse laplace transform
U_rec = reshape(ft(:,end,:),size(X1)); %final time soln



if rep ==0
%On cloaking region boundary
fun_bd = @(s) heat_freq_moved(s,k0,geo_k.ctr(:,1),geo_k.ctr(:,2),Y1,Y2,geo,a,n_trunc,n_move);
scat_eval_cl  = @(s)heat_freq_scat(s,k0,X1,X2,geo_k,fun_bd);
[fscat_cl, ~] = fftilt(scat_eval_cl,t_max,M,alfa);
U_scat_cl = reshape(fscat_cl(:,end,:),size(X1));

fun_bd_actual = @(s) heat_freq(s,k0,geo_k.ctr(:,1),geo_k.ctr(:,2),Y1,Y2,0,0);
scat_eval  = @(s)heat_freq_scat(s,k0,X1,X2,geo_k,fun_bd_actual);
[fscat, ~] = fftilt(scat_eval,t_max,M,alfa);
U_scat = reshape(fscat(:,end,:),size(X1));

end

%%%%Errors
%Reproduction
field = temp(X1,X2,t_max);
%field = reshape(field(:,end,:), size(X1));
eps = 0;
ext_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius+eps)^2;
int_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius+eps)^2;
int_err = int_mask.*log10(abs(field-U_rec));
ext_err = ext_mask.*log10(abs(U_rec));
tot_err = int_err+ext_err;

%calculate maximum allowable temperature
max_temp  = 0; %initiate

%take max over x in omega from [0,t_max]
for i = 1:M
    if max(max(temp(X1,X2,dt*i).*int_mask.*max_temp_scale))>max_temp
        max_temp = max(max(temp(X1,X2,dt*i).*int_mask)).*max_temp_scale;
    end
end
save('ext_cloaking_data.mat')

%calculate size of "disc devices"
[new_src_locs, Y1m, Y2m] = move_src(a, geo, n_move); 
[U_rec_phys,rads] = phys_dev(U_rec,new_src_locs, max_temp, X1, X2,N1);
%%%%Plotting
%time 1
thickLines;
cmax = max(max(field));
cmin = -cmax;

figure(1); clf; 
hold on;
[~,h]=contourf(X1,X2,field+U_scat,[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo_k.ctr(:,1),geo_k.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis square;
axis tight;
axis xy;
axis off;
colorbar;
colormap(redblue)
caxis([cmin,cmax]);
hold off;

figure(2); clf; 
hold on;
[~,h] = contourf(X1,X2,(field-U_rec+U_scat-U_scat_cl),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
h= patch(geo_k.ctr(:,1),geo_k.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
axis square;
axis tight;
axis xy;
axis off;
colorbar;
caxis([cmin,cmax]);
colormap(redblue)
viscircles(new_src_locs,rads, 'Color', 'k')
hold off;

figure(3); clf; 
hold on;
[~,h]=contourf(X1,X2,(field),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
axis square;
axis tight;
axis xy;
axis off;
colorbar;
colormap(redblue)
caxis([cmin,cmax]);
hold off;
% 
% %time 2
% new_m= 2^9/2; %number of time steps
% t_new = dt*new_m;
% U_rec = reshape(ft(:,new_m,:),size(X1));
% U_scat_cl = reshape(fscat_cl(:,new_m,:),size(X1));
% U_scat = reshape(fscat(:,new_m,:),size(X1));
% field = temp(X1,X2,t_new);
% figure(4); clf; 
% hold on;
% [~,h]=contourf(X1,X2,field+U_scat,[-2e300,linspace(cmin,cmax,22)]);
% set(h,'edgecolor','none');
% h= patch(geo_k.ctr(:,1),geo_k.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
% axis square;
% axis tight;
% axis xy;
% axis off;
% colorbar;
% colormap(redblue)
% caxis([cmin,cmax]);
% hold off;
% 
% figure(5); clf; 
% hold on;
% [~,h]=contourf(X1,X2,(field-U_rec+U_scat-U_scat_cl),[-2e300,linspace(cmin,cmax,22)]);
% set(h,'edgecolor','none');
% h= patch(geo_k.ctr(:,1),geo_k.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
% axis square;
% axis tight;
% axis xy;
% axis off;
% colorbar;
% caxis([cmin,cmax]);
% colormap(redblue)
% viscircles(new_src_locs,rads, 'Color', 'k')
% hold off;
% 
% figure(6); clf; 
% hold on;
% [~,h]=contourf(X1,X2,(field),[-2e300,linspace(cmin,cmax,22)]);
% set(h,'edgecolor','none');
% axis square;
% axis tight;
% axis xy;
% axis off;
% colorbar;
% colormap(redblue)
% caxis([cmin,cmax]);
% hold off;
% 
% 
% 
% %time 3
% new_m= 2^9/2^3; %number of time steps
% t_new = dt*new_m;
% U_rec = reshape(ft(:,new_m,:),size(X1));
% U_scat_cl = reshape(fscat_cl(:,new_m,:),size(X1));
% U_scat = reshape(fscat(:,new_m,:),size(X1));
% field = temp(X1,X2,t_new);
% figure(7); clf; 
% hold on;
% [~,h]=contourf(X1,X2,field+U_scat,[-2e300,linspace(cmin,cmax,22)]);
% set(h,'edgecolor','none');
% h= patch(geo_k.ctr(:,1),geo_k.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
% axis square;
% axis tight;
% axis xy;
% axis off;
% colorbar;
% colormap(redblue)
% caxis([cmin,cmax]);
% hold off;
% 
% figure(8); clf; 
% hold on;
% [~,h]=contourf(X1,X2,(field-U_rec+U_scat-U_scat_cl),[-2e300,linspace(cmin,cmax,22)]);
% set(h,'edgecolor','none');
% h= patch(geo_k.ctr(:,1),geo_k.ctr(:,2),'k'); set(h, 'edgecolor', 'none')
% axis square;
% axis tight;
% axis xy;
% axis off;
% colorbar;
% caxis([cmin,cmax]);
% colormap(redblue)
% viscircles(new_src_locs,rads, 'Color', 'k')
% hold off;
% 
% figure(9); clf; 
% hold on;
% [~,h]=contourf(X1,X2,(field),[-2e300,linspace(cmin,cmax,22)]);
% set(h,'edgecolor','none');
% axis square;
% axis tight;
% axis xy;
% axis off;
% colorbar;
% colormap(redblue)
% caxis([cmin,cmax]);
% hold off;
% 
% 
% 






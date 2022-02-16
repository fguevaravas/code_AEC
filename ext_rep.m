%Exterior reproduction in the frequeuncy domain 
%Inputs:
%   k: wave number
%Ouputs:
%   figures: exterior field reproduction and baseline 
function[] = ext_rep(k)
%%%%Parameters
%Physical
Y1 = 5; Y2 =2; %location of point source

%%%%Geometry
%Grid 
L1 = 10; L2 = 10; % dimensions of domain
N1 = 200; N2 = 200; % number of points on grid
[X1,X2] = ndgrid(linspace(0,L1,N1),linspace(0,L2,N2));

%Cloaking region
cp =128; % number of points on circle
Rcenter = [L1/2,L2/2]; radius = min (L1,L2)/6; % center and radius
a = (sqrt(2)-1)*radius; %distance to new point sources 
geo = circ(cp, Rcenter, radius); 

%Addition formula 
n_move = 4; %number of exterior points to use for cloaking 
n_trunc = 22; %truncation for Graff's addition formula

%%%Reproduction
%Full field
[field,~,~] = field_point_src(X1,X2, Y1,Y2,k);
U_rec = helm_moved(k,X1,X2,Y1,Y2,geo,a,n_trunc,n_move);

thickLines;
cmax = max(max(abs(field)))/10;
cmin = -cmax;
figure; clf;
hold on;
[~,h]=contourf(X1,X2,(real(field.')),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
axis xy;
axis square;
axis tight;
axis off;
caxis([cmin,cmax]);
colorbar
hold off;

figure; clf;
hold on;
[~,h]=contourf(X1,X2,real(U_rec.'),[-2e300,linspace(cmin,cmax,22)]);
set(h,'edgecolor','none');
axis xy;
axis square;
axis tight;
axis off;
caxis([cmin,cmax]);
colorbar
hold off;

%Masks for errors
[new_src_locs, ~, ~] = move_src(a, geo, n_move); 
urad= norm(new_src_locs(1,:)-new_src_locs(2,:))/2-.1;

%urchin mask
mask1 = (new_src_locs(1,1)-X1).^2+(new_src_locs(1,2)-X2).^2>urad.^2;
mask2 = (new_src_locs(2,1)-X1).^2+(new_src_locs(2,2)-X2).^2>urad.^2;
mask3 = (new_src_locs(3,1)-X1).^2+(new_src_locs(3,2)-X2).^2>urad.^2;
mask4 = (new_src_locs(4,1)-X1).^2+(new_src_locs(4,2)-X2).^2>urad.^2;
urchmask = mask1.*mask2.*mask3.*mask4;

%interior error mask
mask1 = X2>new_src_locs(1,2);
mask2 = X1>new_src_locs(1,1);
mask3 = X1<new_src_locs(2,1);
mask4 = X2<new_src_locs(3,2);
intmask = mask1+mask2+mask3+mask4==0;
extmask = not(intmask);

tot_err = log10(abs(U_rec.'-field.')).*intmask + log10(abs(U_rec.')).*extmask;
tot_err = tot_err.*urchmask;

%C_i radius
cirad= a+radius - norm(new_src_locs(1,:)-new_src_locs(2,:))/2-.05;
figure; clf;
hold on;
imagesc(linspace(0,L1,N1),linspace(0,L2,N2),tot_err);
axis xy;
axis square;
axis tight;
axis off;
%caxis([cmin,cmax]);
cmap = colormap;
[N,~] = size(cmap);
cmap(N,:) = 1;
colormap(cmap);
colorbar
viscircles(Rcenter,cirad, 'Color', 'k', 'LineStyle', '--','LineWidth', 1)
hold off;
end%function








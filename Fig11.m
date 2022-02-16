%%%%%Device radius vs. radius of green circ

%%%%Parameters
%Physical
k0 = 1.3; %thermal diffusivity
Y1 = 10; Y2 =1; %location of point source
max_temp_scale = 100; %factor to allow maximum temp of pelt devices

%NILT 
t_max = 1;
alfa =0;
M = 2^6; %number of time steps
dt = t_max/M;

%Reproduction
rep = 0; %0 for full cloaking

%Temperature function
temp = @(X1,X2,t) temp_source(X1,X2,t,0, Y1, Y2, 1, k0);


%%%%Geometry
%Grid 
L1 = 20; L2 = 20; % dimensions of domain
N1 = 100; N2 = 100; % number of points on grid
[X1,X2] = ndgrid(linspace(0,L1,N1),linspace(0,L2,N2));
p = [X1(:),X2(:)];

%Cloaking region
cp =128 ; % number of points on circle
Rcenter = [L1/2,L2/2]; all_radius =linspace(1,8,40)/sqrt(2);

deltas = zeros(length(all_radius),1);
devr = zeros(length(all_radius),1);

for j= 1:length(all_radius)
radius = all_radius(j);
eps = 0;
ext_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 > (radius+eps)^2;
int_mask = (X1-Rcenter(1)).^2 + (X2 - Rcenter(2)).^2 < (radius+eps)^2;

geo = circ(cp, Rcenter, radius); 

%Addition formula 
a = (sqrt(2)-1)*radius; %distance to new point sources 
n_move = 4; %number of exterior points to use for cloaking 
n_trunc = 22; %truncation for Graff's addition formula

%%%Reproduction
%Full field
  fun = @(s) heat_freq_moved(s,k0,X1,X2,Y1,Y2,geo,a,n_trunc,n_move);
[ft, t] = fftilt(fun,t_max,M,alfa); %numerical inverse laplace transform
U_rec = reshape(ft(:,end,:),size(X1)); %final time soln

%calculate maximum allowable temperature
max_temp  = 0; %initiate

%take max over x in omega from [0,t_max]
for i = 1:M
    if max(max(temp(X1,X2,dt*i).*int_mask.*max_temp_scale))>max_temp
        max_temp = max(max(temp(X1,X2,dt*i).*int_mask)).*max_temp_scale;
    end
end

[new_src_locs, Y1m, Y2m] = move_src(a, geo, n_move); 
[U_rec_phys,rads] = phys_dev(U_rec,new_src_locs, max_temp, X1, X2,N1);
deltas(j) = a+radius;
devr(j) = max(rads)/deltas(j);
end

max_rad = zeros(size(t));
for i=1:M
    current = reshape(ft(:,i,:),size(X1));
    [~,rads] = phys_dev(current,new_src_locs, max_temp, X1, X2,N1);
    max_rad(i) = max(rads)/(radius+a);
end

save('size_study.mat')
%plot

%value that circs would touch 
touch = 0.7059029628545967494090973559650592505931854248046875;

thickLines;
figure(1); clf;
hold on
plot(deltas, devr)
plot(deltas, touch*ones(length(deltas)), 'LineStyle', '--', 'Color', 'k')
xlim([min(deltas),max(deltas)])
xlabel('\delta_C')
ylabel('device radius/\delta_C')
hold off;

figure(2); clf;
hold on
plot(t, max_rad)
plot(t, touch*ones(length(t)), 'LineStyle', '--', 'Color', 'k')
ylim([0,.8])
hold off
xlabel('time')
ylabel('device radius/\delta_C')

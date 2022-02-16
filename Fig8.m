%generate the grafs truncation vs frequency 
thickLines;
%original src location 
Y1 = 8; Y2 =5;

%frequencies
theta = linspace(0.05,20,100);
freqs = (1+sqrt(2)*1i)*1/sqrt(3)*theta;

%initailize vectors to store errors
act_er = zeros(length(freqs),1);
est_er = act_er;

M = 22;
n_move = 4; %number of sources for exterior cloak 
emp_n = 3; %truncation for empirical estimation 

%geometry of interior cloak
cp =128; % number of points on circle
Rcenter = [5,5]; radius = min (10,10)/6; % center and radius
a = (sqrt(2)-1)*radius; %distance to new point sources 
geo = circ(cp, Rcenter, radius); 

%find X1s & X2s to evaluate at 
eps = radius*.1; %buffer for evaluating field;
[new_srcs, ~, ~] = move_src(a, geo, n_move); 
ext_rad = norm(new_srcs(1,:)-new_srcs(2,:))/2;
Xs = zeros(2*n_move,2);
for i =1:n_move
    v_dif = new_srcs(i,:)-Rcenter;
    v_dif = v_dif/norm(v_dif);
    Xs(i,:) = Rcenter + v_dif*(a+radius-ext_rad-eps);
    Xs(4+i,:) = new_srcs(i,:) + v_dif*(ext_rad+eps);
end


%find truncations
for i = 1: length(freqs)
    k = freqs(i);
    [act_er(i),est_er(i)] = graf_est(k,Y1,Y2,Xs(:,1),Xs(:,2),geo,a,n_move, M, emp_n,eps,ext_rad);
end

%plot 
thickLines;
figure(1); clf; 
semilogy(theta,act_er)
hold on
plot(theta,est_er,'--')
xlabel('\theta')
ylabel('absolute error')
legend('error','bound')
hold off


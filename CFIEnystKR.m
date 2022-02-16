%Originally from: Alex Barnett ``Boundary integral equations for BVPs, and their
%high-order Nystr¨om quadratures: a tutorial" https://math.dartmouth.edu/~fastdirect/notes/quadr.pdf
%Modifed by TD for geometry we use.

%%% Kapur-Rokhlin scheme for solving exterior dirichlet problem 
%Inputs:
%   geo: geometry of the scatterer
%   i/j: position in the paramterization of the scatterer 
%   k: wave number
%   eta: coupling parameter
%Ouputs:
%   a: ijth component of matrix of (D-ietaS) (D/S are double/single
%   boundary layer operators)
function a = CFIEnystKR(geo,i,j,k,eta)
g6 = [4.967362978287758 -16.20501504859126 25.85153761832639 ...
-22.22599466791883 9.930104998037539 -1.817995878141594]; % 6th order
sw = 2*pi/length(geo.ctr(:,1))*sqrt(geo.dx(j)^2+geo.dy(j)^2);
N = numel(geo.ctr(:,1));
l = mod(j-i,N); 
if l>N/2
    l=N-l; 
end % index distance i to j

if l>0 && l<=6
    sw = sw * (1 + g6(l));
end % apply correction to sw
% sp1 = geo.dy(j)*(geo.ctr(j,1)-geo.ctr(i,1));
% sp2 = geo.dx(j)*(geo.ctr(j,2)-geo.ctr(i,2));
% R = norm(geo.ctr(i,:)-geo.ctr(j,:));
% K2 = 1i/2*k*(sp1-sp2)*besselh(1,k*R)/R*sw;
% K1 = 1i/2*k*besselh(0,k*R)*norm([geo.dx(j),geo.dy(j)])*sw;
% a = K2+1i*eta*K1;
pmc1 = geo.ctr(i,1)-geo.ctr(j,1);
pmc2 = geo.ctr(i,2)-geo.ctr(j,2);
R = sqrt(pmc1.^2+pmc2.^2);
DG1 = k*(1i/2)*besselh(1,k*R)*pmc1/R;
DG2 = k*(1i/2)*besselh(1,k*R)*pmc2/R;
K2 = DG1*geo.nv(j,1)+DG2.*geo.nv(j,2);
K1= (1i/2)*besselh(0,k*R);
a = (K2-1i*eta*K1)*sw;
if i==j
    a = 0;
end % kill diagonal




% g6 = [4.967362978287758 -16.20501504859126 25.85153761832639 ...
% -22.22599466791883 9.930104998037539 -1.817995878141594]; % 6th order
% sw = G.sp(j)*G.w; % speed weight
% N = numel(G.x);
% l = mod(i-j,N); 
% if l>N/2
%     l=N-l; 
% end % index distance i to j
% 
% if l>0 && l<=6
%     sw = sw * (1 + g6(l));
% end % apply correction to sw
% 
% d = G.x(i)-G.x(j); kr = k*abs(d); % last 3 lines do CFIE kernel:
% costheta = real(conj(G.nx(j)).*d)./abs(d); % theta angle between x-y & ny
% a = (1i/4) * (k*costheta*besselh(1,kr) - 1i*eta*besselh(0,kr)) * sw;
% 
% if i==j
%     a = 0;
% end % kill diagonal
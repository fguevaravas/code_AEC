%Originally from: Alex Barnett ``Boundary integral equations for BVPs, and their
%high-order Nyström quadratures: a tutorial" https://math.dartmouth.edu/~fastdirect/notes/quadr.pdf
%Modified by TD for the geometry we use

function[phi] = Kapur_scat(geo,k,dir,eta)
rhs = -2*dir; % bdry data at nodes
N = length(dir);
A = zeros(N);
for i = 1:N
    for j=1:N
        A(i,j) = CFIEnystKR(geo,i,j,k,eta);
    end
end
phi = (eye(N)+A) \ rhs; % dense solve
end%function
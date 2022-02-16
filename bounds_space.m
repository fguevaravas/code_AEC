%%% Finds the monopole/dipole errors & bounds for moving one source as a
%%% function of space
%Inputs:
%   k: wave number
%   Y1: og position of point source (x)
%   Y2: og position of point source (y)
%   Y1move: position to move point source to (x)
%   Y2move: position to move point source to (y)
%   X1: position to evaluate field at (x)
%   X2s: positions to evaluate field at (y)
%Ouputs:
%   mono_errors: actual error on monopoles
%   mono_bound: bound for monopole error
%   di_errors: actual error for dipoles
%   di_bound: bound for dipole error
%   dists: |y-xj|/|x-xj|

function[mono_errors, mono_bound, di_errors, di_bound, dists] = bounds_space(k, Y1, Y2, Y1move, Y2move, X1, X2s, N)
%Compute dist that we moved the source
u = sqrt((Y1-Y1move)^2+(Y2-Y2move)^2);

%Initalize vectors
soln = zeros(size(X2s));
truncsoln = zeros(size(X2s));
soln_di = zeros(size(X2s));
truncsoln_di = zeros(size(X2s));
errbd1 = zeros(size(X2s));
errbd1_di = zeros(size(X2s));

%%%Empirically approximate error bound constants
C1 = 0;
C2 = 0;
for i = 1:length(X2s)
    %components for Graf
    X2  = X2s(i);
    w = sqrt((X1-Y1).^2+(X2-Y2).^2);
    v = sqrt((X1-Y1move).^2+(X2-Y2move).^2);
    theta1 = atan2(X2-Y2move, X1-Y1move);
    theta2 = atan2(Y2-Y2move, Y1-Y1move);
    theta3 = atan2(Y2- X2, Y1- X1);
    theta4 = atan2(Y2move-X2, Y1move - X1);
    N_est =4;
    
    %Aproximate monopole and dipole field
    moved_field = 0;
    moved_field_di = 0;
    for j = -N_est:N_est
        moved_field = moved_field+besselh(j,k*v)*besselj(j,k*u).*exp(1i*j*(theta1- theta2));
        moved_field_di = moved_field_di + besselh(1+j,k*v)*besselj(j,k*u).*exp(1i*j*(theta1- theta2)).*exp(-1i*(theta3 - theta4));
    end
    
    %Compute actual fields and remainders
    original_field = 1i/4*besselh(0,k*w);
    original_field_di  = 1i/4*besselh(1,k*w);
    remainder = abs(original_field- 1i/4*moved_field);
    remainder_di = abs(original_field_di-1i/4*moved_field_di);
    
    %Compute error bound series
    z1 = abs(u/v);
    S1 = (-log(1-z1));
    for j = 1:N_est
        S1 = S1 - z1^j/j;
    end
    
    %Approximate C
    C = remainder/abs(S1);
    if C>C1
        C1 = C;
    end
    
    C = remainder_di/abs((z1)^(N_est+1)/(1-z1));
    if C>C2
        C2 = C;
    end
end


for j = 1:length(X2s)
    %Compute components for Graf
    X2 = X2s(j);
    w = sqrt((X1-Y1).^2+(X2-Y2).^2);
    v = sqrt((X1-Y1move).^2+(X2-Y2move).^2);
    theta1 = atan2(X2-Y2move, X1-Y1move);
    theta2 = atan2(Y2-Y2move, Y1-Y1move);
    theta3 = atan2(Y2- X2, Y1- X1);
    theta4 = atan2(Y2move-X2, Y1move - X1);
    
    %Approximate fields
    moved_field = 0;
    moved_field_di = 0;
    for i = -N:N
        moved_field = moved_field+besselh(i,k*v)*besselj(i,k*u).*exp(1i*i*(theta1- theta2));
        moved_field_di = moved_field_di + besselh(1+i,k*v)*besselj(i,k*u).*exp(1i*i*(theta1- theta2)).*exp(-1i*(theta3 - theta4));
    end
    
    %Calculate actual fields and remainders 
    original_field = 1i/4*besselh(0,k*w);
    original_field_di  = 1i/4*besselh(1,k*w);
    
    truncsoln(j) = original_field;
    truncsoln_di(j) = original_field_di;
    
    soln_di(j) = 1i/4*moved_field_di;
    soln(j) = 1i/4*moved_field;
    
     %Monopole error bounds 
    z1 = abs(u/v);
    S1 = (-log(1-z1));
    for i = 1:N
       S1 = S1 - z1^i/i;
    end
    errbd1(j) = C1*S1;

    
    %Dipole Error bounds ;
    errbd1_di(j) = C2*(z1)^(N+1)/(1-z1);
end
%Return errors, error bounds, and |y-x_j|/|x-x_j|
dists =u./sqrt((X1-Y1move).^2 + (X2s-Y2move).^2);
di_bound = errbd1_di;
mono_bound = errbd1;
mono_errors = abs(soln-truncsoln);
di_errors = abs(soln_di-truncsoln_di);
end%function

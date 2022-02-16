%%% Finds the monopole/dipole errors & bounds for moving one source as a
%%% function of frequency 
%Inputs:
%   k: wave number
%   Y1: og position of point source (x)
%   Y2: og position of point source (y)
%   Y1move: position to move point source to (x)
%   Y2move: position to move point source to (y)
%   X1: position to evaluate field at (x)
%   X2: position to evaluate field at (y)
%   trunc: vector of truncations 
%Ouputs:
%   mono_errors: actual error on monopoles
%   mono_bound: bound for monopole error
%   di_errors: actual error for dipoles
%   di_bound: bound for dipole error

function[mono_errors, mono_bound, di_errors, di_bound] = bounds_freq(ks, Y1, Y2, Y1move, Y2move, X1, X2, N)
%All of the components necessary for Grafs addition formula
w = sqrt((X1-Y1).^2+(X2-Y2).^2);
v = sqrt((X1-Y1move).^2+(X2-Y2move).^2);
u = sqrt((Y1-Y1move)^2+(Y2-Y2move)^2);

theta1 = atan2(X2-Y2move, X1-Y1move);
theta2 = atan2(Y2-Y2move, Y1-Y1move);
theta3 = atan2(Y2- X2, Y1- X1);
theta4 = atan2(Y2move-X2, Y1move - X1);

%Initalize error vectors
soln = zeros(size(ks));
truncsoln = zeros(size(ks));
soln_di = zeros(size(ks));
truncsoln_di = zeros(size(ks));
errbd1 = zeros(size(ks));
errbd1_di = zeros(size(ks));


%%%%Estimate constant for error bounds
C1 = 0;
C2 = 0;
%Take max over all frequencies
for i = 1:length(ks)
    k = ks(i);
    N_est =4;
    
    %Estimate field using Grafs
    moved_field = 0;
    moved_field_di = 0;
    for j = -N_est:N_est
        moved_field = moved_field+besselh(j,k*v)*besselj(j,k*u).*exp(1i*j*(theta1- theta2));
        moved_field_di = moved_field_di + besselh(1+j,k*v)*besselj(j,k*u).*exp(1i*j*(theta1- theta2)).*exp(-1i*(theta3 - theta4));
    end
    
    %Find actual field & remainders
    original_field = 1i/4*besselh(0,k*w);
    original_field_di  = 1i/4*besselh(1,k*w);
    remainder = abs(original_field- 1i/4*moved_field);
    remainder_di = abs(original_field_di-1i/4*moved_field_di);
    
    %Compute bounds
    z1 = abs(u/v);
    S1 = (-log(1-z1));
    for j = 1:N_est
        S1 = S1 - z1^j/j;
    end
    
    %Estimate C for monopole and dipole
    C= remainder/abs(S1);
    if C>C1
        C1 = C;
    end
    C = remainder_di/abs((z1)^(N_est+1)/(1-z1));
    if C>C2
        C2 = C;
    end
end

%%%Find errors and error bounds
for j = 1:length(ks)
    k = ks(j);
    moved_field = 0;
    moved_field_di =0;
    for i = -N:N %grafs addition formula
        moved_field = moved_field+besselh(i,k*v)*besselj(i,k*u).*exp(1i*i*(theta1- theta2));
        moved_field_di = moved_field_di + besselh(1+i,k*v)*besselj(i,k*u).*exp(1i*i*(theta1- theta2)).*exp(-1i*(theta3 - theta4));
    end
    
    %Find actual fields
    original_field = 1i/4*besselh(0,k*w);
    original_field_di  = 1i/4*besselh(1,k*w);
    truncsoln(j) = original_field;
    truncsoln_di(j) = original_field_di;
    soln(j) = 1i/4*moved_field;
    soln_di(j) = 1i/4*moved_field_di;
    %Monopole error bounds
    S1 = (-log(1-z1));
    for i = 1:N
        S1 = S1 - z1^i/i;
    end
    errbd1(j) = C1*S1;
    %Dipole Error bounds ;
    errbd1_di(j) = C2*(z1)^(N+1)/(1-z1);
end
%Return vectors of errors and error bounds
di_bound = errbd1_di;
mono_bound = errbd1;
mono_errors = abs(soln-truncsoln);
di_errors = abs(soln_di-truncsoln_di);
end%function





%Actual error vs error bound at a fixed M
%Inputs:
%   k: wave number
%   Y1: position of point source (x)
%   Y2: og position of point source (y)
%   X1: positions to evaluate field at (x)
%   X2: positions to evaluate field at (y)
%   geo: geometry of the interior cloak
%   a: distance for new point sources
%   n_move: number of new point sources
%   M: actual truncation
%   emp_n: empirical truncation to estimate error constant
%   eps: buffer for x
%   ext_rad: radius of non-convergence region 
%Ouputs:
%   
function[acter,erest] = graf_est(k,Y1,Y2,X1,X2,geo,a,n_move, M, emp_n,eps,ext_rad)
z = ext_rad/(ext_rad+eps); %ratio for geo series 
%empirically guess C 
soln = field_point_src(X1(1:4),X2(1:4),Y1,Y2,k);
field_est = helm_moved(k,X1,X2,Y1,Y2,geo,a,emp_n,n_move);

%Find errors at closest z
Rm(1:4) = abs(soln-field_est(1:4)); 
Rm(5:8) =  abs(field_est(5:8));
%Estimate C
C = max(Rm/(z)^(emp_n+1)/(1-z));
%Estimate error
erest=C*z^(M+1)/(1-z);

%Find actual error
est = helm_moved(k,X1,X2,Y1,Y2,geo,a,M,n_move);
er_in = max(abs(est(1:4)-soln));
er_out = max(abs(est(5:8)));
acter = max(er_in,er_out);
end%function
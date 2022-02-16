% function [x] = kite(t)
%
% Kite scatterer in parametric form, as defined in Fig 3.1 of Colton and Kress,
% Inverse acoustic and electromagnetic scattering.
% 
% Input:
%   t: parameter
% Ouput:
%   geo.h: length of segments
%   geo.ctr: center points of segments
%   geo.nv: normal vectors to center points of segments
function [geo] = kite(t)
  %curve
  x=  .25*(cos(t)  +   0.65*cos(2*t)) + 4.95;
  y =  .5*sin(t)+5;
  p = [x,y];
  s = [diff(p,1,1); p(1,:) - p(end,:)];
  nv = [ s(:,2) , -s(:,1) ];
  for i=1:length(nv)
    nv(i,:) = nv (i,:)/norm(nv(i,:));
  end
  geo.h = sqrt(s (:,1).^2 + s (:,2).^2);
  geo.ctr = [x,y];
  geo.nv = nv;
  geo.dx = -.25*(sin(t)  +   0.65*2*sin(2*t));
  geo.dy = .5*cos(t);
end

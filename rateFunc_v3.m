function [ fx,gx ] = rateFunc_v3( xi,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x
%   and parameter values. It interpolates the rate functions from libraries
% computed with function rateFunc_v2 using a nearest neighbours approach
% inspired by the qinterp1 function
xi = xi - parms.xLib(1);      % subtract minimum of xLib

rxi = round(xi/parms.dxLib)+1;   % find indices of nearest neighbours

gx = parms.gxLib(rxi);
fx = parms.fxLib(rxi);
end

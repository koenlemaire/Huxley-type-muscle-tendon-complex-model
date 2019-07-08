function [ fx,gx ] = rateFunc_v4( xi,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x
%   and parameter values. It interpolates the rate functions from libraries
% computed with function rateFunc_v2 using a linear interpolation approach
% inspired by the qinterp1 function written by N. Brahms (2006)

% Koen Lemaire 11/2014

ndxLib=parms.ndxLib;
fxLib=parms.fxLib;
gxLib=parms.gxLib;
% 
xi = xi(:) - parms.xLib(1);      % subtract minimum of xLib

fxi=floor(xi*ndxLib)+1; % get index of lowest nearest neighbours
x1=fxi-xi*ndxLib;       % help vars
x2=1-fxi+xi*ndxLib;

fx = x1.*fxLib(fxi) + x2.*fxLib(fxi+1);
gx = x1.*gxLib(fxi) + x2.*gxLib(fxi+1); 
end


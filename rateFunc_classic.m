function [ fx,gx ] = rateFunc_classic( x,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x
%   and parameter values, for the absolute classic Huxley model

% Koen Lemaire 2014
f1=parms.f1;
g1=parms.g1;
g2=parms.g2;

fx= (x>=0 & x<=1)*f1.*x;
gx= (x>=0)*g1.*x + (x<0)*g2;
end

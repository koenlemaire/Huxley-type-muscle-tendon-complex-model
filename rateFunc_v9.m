function [ fx,gx ] = rateFunc_v9( x,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x 
%   and parameter values. 

% Koen Lemaire
f1=parms.f1; % asymptote of fx
g1=parms.g1; % lowest value of gx
g2=parms.g2; % left horizontal asymptote of gx
g3=parms.g3; % right horizontal asymptote of gx

%% indices
idx_f1=x>0&x<=1;
idx_g2=x<0;
idx_g3=x>1;

%% gx 
gx=g1*x.*(idx_f1 | idx_g3) + g2*idx_g2 + (g3-g1).*idx_g3;

%% fx
fx=x.*f1.*idx_f1;
end

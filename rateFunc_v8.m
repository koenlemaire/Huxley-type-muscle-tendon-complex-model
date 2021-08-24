function [ fx,gx ] = rateFunc_v8( x,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x 
%   and parameter values. 

% Koen Lemaire
f1=parms.f1; % asymptote of fx
g1=parms.g1; % lowest value of gx
g2=parms.g2; % left horizontal asymptote of gx
g3=parms.g3; % right horizontal asymptote of gx

%% indices
idx_f1=x>0&x<2;
idx_g1=x>0;

%% hvs
% calculate only once; expensive!
hvs_steepness=50;
hvs_f1=1./(1+exp(hvs_steepness*(x-1))); % centered around x=1, left asymptote is 1
hvs_g2=1./(1+exp(hvs_steepness*x)); % centered around x=0, left asymptote is 1
hvs_g3=1./(1+exp(-hvs_steepness*(x-1.2))); % centered around x=1.2, right asymptote is 1

%% gx 
% halfway point, as measured from minimum value: x_halfway = a^(1/n)
gx=g1*x.*idx_g1 + g2*hvs_g2 + (g3+g1)*hvs_g3;

%% fx
fx=x.*f1.*hvs_f1.*idx_f1;
end

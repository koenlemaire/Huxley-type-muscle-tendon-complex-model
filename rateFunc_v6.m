function [ fx,gx ] = rateFunc_v6( x,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x 
%   and parameter values. it deploys a piecewise linear ratefunction
%   connected by sinosoids and exponential in the outer regions. The
%   functions have continuous first derivatives ...

% Koen Lemaire 2016
f1=parms.f1; % asymptote of fx
g1=parms.g1; % lowest value of gx
g2=parms.g2; % left horizontal asymptote of gx
g3=parms.g3; % right horizontal asymptote of gx


%% gx 
% halfway point, as measured from minimum value: x_halfway = a^(1/n)
x_halfway=0.5;
n=5; % related to the slope at the halfway point ... 
a=x_halfway^n;

% first make GX. gx is constructed as two sigmoids attached to eachother
% which originate in (0.5,g1) and which have horizontal asymptotes g2 and g3.
x1=x(x<=0.5); % left s curve
x2=x(x>0.5); % right s curve 

% first calculate powers (expensive), help variables
hgx1=(-x1+0.5).^n; 
hgx2=(x2-0.5).^n; 

% division, scaling and moving in y direction
gx = [(g2-g1).*(hgx1 ./ (hgx1 + a)); (g3-g1).*(hgx2 ./ (hgx2 + a))] + g1;

%% fx
% fx is a sigmoid which is truncated by a smooth hvs function
fidx=x>0&x<=2;
x1=x(fidx);

hfx=x1.^n; % help variable x^n

% smooth heaviside step function:
hvs=1./(1+exp(-40*(1-x1)));

fx=zeros(size(gx));
fx(fidx)= f1*hvs.*hfx./(hfx+a);
%fx(fidx)= f1.*hfx./(hfx+a);
end

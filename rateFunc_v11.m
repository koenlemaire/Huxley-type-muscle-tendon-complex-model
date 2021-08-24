function [ fx,gx ] = rateFunc_v11( x,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x
%   and parameter values. it deploys a piecewise linear ratefunction
%   connected by sinosoids

% Koen Lemaire 2014
f1=parms.f1;
g1=parms.g1;
g2=parms.g2;
g3=parms.g3;

x1=((x>-.1)&(x<=.1));
x2=((x>.1)&(x<=.9));
x3=((x>.9)&(x<=1.1));

fx= x1*f1.*.5.*(1+cos(pi/.2*(x-.1))) + ...
    x2*f1 + ...
    x3*f1.*.5.*(1+cos(pi/.2*(x-.9)));

gx= (x<-.1).*(g2*(x+.1).^2+g2) + ...
    x1.*.5.*(g2+g2.*cos(pi/.2*(x+.1))) + ...
    x3.*.5.*(g3-g3.*cos(pi/.2*(x-.9))) + ...
    (x>1.1).*(g3*(x-1.1).^2+g3);
end

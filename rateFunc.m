function [ fx,gx ] = rateFunc( x,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x
%   and parameter values. it deploys a piecewise linear ratefunction
%   connected by sinosoids
f1=parms.f1;
g1=parms.g1;
g2=parms.g2;
g3=parms.g3;

fx= ((x>-.1)&(x<=.1))*f1.*.5.*(1+cos(pi/.2*(x-.1))) + ...
    ((x>.1)&(x<=.9))*f1 + ...
    ((x>.9)&(x<=1.1))*f1.*.5.*(1+cos(pi/.2*(x-.9)));

gx= (x<-.1)*g2 + ...
    ((x>-.1)&(x<=.1)).*.5.*((g2+g1)+(g2-g1).*cos(pi/.2*(x+.1))) + ...
    ((x>.1)&(x<=.9))*g1 + ...
    ((x>.9)&(x<=1.1)).*.5.*((g1+g3)+(g1-g3).*cos(pi/.2*(x-.9))) + ...
    (x>1.1)*g3;

end


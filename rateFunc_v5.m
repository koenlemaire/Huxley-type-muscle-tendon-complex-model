function [ fx,gx ] = rateFunc_v5( x,parms )
%[ fx,gx ] = rateFunc( x,parms )
%   this function computes the values of the rate functions given current x 
%   and parameter values. it deploys a piecewise linear ratefunction
%   connected by sinosoids and polynomial in the outer regions. The
%   functions have continuous first derivatives ...

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

f1=parms.f1;
g2=parms.g2;
g3=parms.g3;

if isfield(parms,'g1')
    g1=parms.g1;
else
    g1=0;
end

int=.05; % range of values amongst transition region

x1=(x<-int); % lower than 0
x2=((x>-int)&(x<=int)); % transition region
x3=((x>int)&(x<=1-int)); % between 0 and h
x4=((x>1-int)&(x<=1+int)); % transition region
x5=(x>1+int); % higher than h

fx= x2.*f1.*.5.*(1+cos(pi/(2*int).*(x-int))) + ...
    x3.*f1 + ...
    x4.*f1.*.5.*(1+cos(pi/(2*int).*(x-1+int)));

gx= x1.*(g2*(x+int).^2+g2) + ...
    x2.*.5.*((g2+g1)+(g2-g1).*cos(pi/(2*int).*(x+int))) + ...
    x3.*g1 + ...
    x4.*.5.*((g3+g1)-(g3-g1).*cos(pi/(2*int).*(x-1+int))) + ...
    x5.*(g3*(x-(1+int)).^2+g3);
end

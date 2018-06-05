function [ q ] = activeState( gamma,com )
%function [ q ] = activeState( gamma,parms )
%   input: free calcium concentration gamma and parms
%   output: active state q
% formula according to Curtin (1998)

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

% bottom parameter chosen to more or less match: 
% kernell (1983)
% mannard (1973)
% Rack & westbury (1969)
% stephenson (1984)
n=com.n; %=2 CAREFULL!!
k=com.k; %=.35
q = (1+k.^n).*(gamma.^n)./(gamma.^n + k.^n);
qmin=com.qmin;
q(q<qmin)=qmin;
return


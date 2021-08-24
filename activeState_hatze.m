function [ q ] = activeState_hatze( gamma,lcerel,parms )
%function [ q ] = activeState_hatze( gamma,lcerel,parms )
%   input: free calcium concentration gamma, relative CE length lcerel, parms
%   output: fraction of troponin with calcium bound (potential isometric force); active state q
% formula from Hatze 1977/1978, as described in Rockenfeller and Guenther
% (2017), eqns B4 and figure B3

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

% parameter values
v=parms.v; % 3 according to Hatze/Rockenfeller
lp=parms.lp; % 2.9 according to rockenfeller
% rho_c=[Ca]max * rho0 (ie molar volume); 1.373e-4 * 55700 (for v=3)
% according to Hatze. But rho0 ranges from 10e-6 to 16e-6 according to
% Rockenfeller (and some extent Kistemaker 2007). rho_c = 7 according to
% rockenfeller figure B4 but seems more likely to be 0.7 given the above
% and the discussion in Rockenfeller appendix C
rho_c=parms.rho_c; % 7?? 0.7?? 

% equations as eqns B4 in rockenfeller 2017
rho=rho_c.*(lp-1)./(lp./lcerel-1);
q=(gamma.^v)./(1/(rho.^v) + gamma.^v); 
q0=1e-10; % arbitrary number
q(q<q0)=q0;

return


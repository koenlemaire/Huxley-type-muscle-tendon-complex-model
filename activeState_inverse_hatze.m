function [ Err ] = activeState_inverse_hatze( gamma,lcerel0,parms )
%function [ q ] = activeState_hatze( gamma,lcerel,parms )
%   input: free calcium concentration gamma, relative CE length lcerel, parms
%   output: fraction of troponin with calcium bound (potential isometric force); active state q
% formula from Hatze 1977/1978, as described in Rockenfeller and Guenther
% (2017), eqns B4 and figure B3

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)
[ q ] = activeState_hatze( gamma,lcerel0,parms );
q0=parms.q0; % desired initial q
Err=q0-q;
return

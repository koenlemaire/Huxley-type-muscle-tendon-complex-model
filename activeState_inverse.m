function [ Err ] = activeState_inverse( gamma,parms )
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
[ q ] = activeState( gamma,parms );
q0=parms.q0;
Err=q0-q;
return

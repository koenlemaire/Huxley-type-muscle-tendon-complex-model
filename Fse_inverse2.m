
function [ Err ] = Fse_inverse2( lse,parms )
% function used to iteratively find lse corresponding to fse0

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)
fse0=parms.fse0;

lpe=parms.lceopt; % not important anyway ... 

[fse, ~, ~, ~] = CEEC_simple2(lse,lpe,parms);

Err=fse0-fse;
return

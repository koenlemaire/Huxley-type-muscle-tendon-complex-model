function [ fisomrel,forcedotpar ] = ce_fl_simple( lcerel,parms )
%function [ fisomrel,forcedotpar ] = forceLengthRelation( lcerel,parms )
%
% forcelength relation of CE described by a power amplified gaussian curve
%
%   INPUT:
% lcerel:       relative CE length
% parms:        parameter struct
%   OUTPUT:
% fisomrel:     isometric force normalised to fmax with q = 1
% forcedotpar:  paramter as such that dforcelength/dt =
%                forcedotpar*lcereldot

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

C = parms.C;

fisomrel=C.*(lcerel-1).^4 + 1;
forcedotpar = 4.*C.*(lcerel-1).^3;

fisomrel(fisomrel<1e-6)=1e-6;
forcedotpar(forcedotpar<0)=0;
return


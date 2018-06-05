function [ gammad ] = gammadot( gamma,stim,parms )
%function [ gammad ] = gammadot( gamma,stim,parms )
%   input: freecalcium concentration gamma, muscle stimulation stim and parms
%   output: time derivative of gamma

% KKL 2/2013
%gamma = gamma(:);stim = stim(:);
% read out parameters
% tau_act=parms.tau_act;
% tau_deact=parms.tau_deact;

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

% simplified form, use only for scalar input
if stim>=gamma
    gammad = (stim-gamma)/parms.tau_act;
else
    gammad = (stim-gamma)/parms.tau_deact;
end
% use bottom code for vectorized inputs!!!
% a=stim>=gamma; % ascending part
% d=~a; % descending part
% gammad(a) = (stim(a)-gamma(a))/tau_act;
% gammad(d) = (stim(d)-gamma(d))/tau_deact;
return


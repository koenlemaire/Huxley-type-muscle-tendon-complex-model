function [fse, fpe, kse, kpe] = CEEC_simple(lse,lpe,parms)
% [fse, fpe, kse, kpe] = CEEC_simple(lse,lpe,parms)
%
% Calculate Elastic Element Characteristic (CEEC) calculates either SE or
% PE force at current CE  and MTC length, or SE and PE instantanious
% stiffness at current CE and MTC length. Output can be selected via
% 'option'.
%
% INPUT
% lse:      absolute SE length [m]
% lpe:      absolute CE length [m]
% parms:    paramater struct

% PROCESS
% CEEC calculates Stiffness as the length derivative of the
% force-length relation specified by parms. CEEC
% calculates current Force via the same force length relation: use

% OUTPUT
% fse/kse: Series elastic Element Force or Stiffness (wrt to CE length!!) at lse [N]
% fpe/kpe: Series elastic Element Force or Stiffness (wrt to CE length!!) at lpe [N]
%
% CURRENT SETTING: Force length relation: second order polynomial function.

% Koen Lemaire, februari 2011
%% read out parameters
se_shape = parms.se_shape;
lse_slack = parms.lse_slack;
pe_shape = parms.pe_shape;
lpe_slack = parms.lpe_slack;
%% calculate SE/PE Force/Stiffness
% Force
fse = se_shape*(lse-lse_slack).^2;
fpe = pe_shape*(lpe-lpe_slack).^2;

% stiffness, derivative to lse (kse) and lpe (kpe)
kse = 2*se_shape*(lse-lse_slack); 
kpe = 2*pe_shape*(lpe-lpe_slack);

% correct for slack tendon:
fse(lse<=lse_slack)=0;
fpe(lpe<=lpe_slack)=0;
kse(lse<=lse_slack)=0;
kpe(lpe<=lpe_slack)=0;
return


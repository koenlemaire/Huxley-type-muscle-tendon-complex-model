function [fse, fpe, kse, kpe] = CEEC_simple2(lse,lpe,parms)
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
%
% PROCESS
% CEEC calculates Stiffness as the length derivative of the
% force-length relation specified by parms. CEEC
% calculates current Force via the same force length relation: use
%
% OUTPUT
% fse/kse: Series elastic Element Force or Stiffness (wrt to CE length!!) at lse [N]
% fpe/kpe: Series elastic Element Force or Stiffness (wrt to CE length!!) at lpe [N]

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)
%% read out parameters
Fmax = parms.Fmax;
se_strain = parms.se_strain; % [] relative strain @ Fmax
lse_slack = parms.lse_slack; % [m] slack length
pe_strain = parms.pe_strain; % [] relative strain @ Fmax
lpe_slack = parms.lpe_slack; % [m]
lceopt = parms.lceopt; % [m]

%% calculate force and stiffness
% the characteristics are modelled as the sum of a linear and a quadratic
% spring ...

k1_se=.001*Fmax./lse_slack; % between 0 and l_slack
k1_pe=.001*Fmax./lceopt; % stiffness such that force is 0.005 at lceopt...

k2_se=(Fmax-(1+se_strain).*k1_se)./((se_strain.*lse_slack).^2);
k2_pe=(Fmax-(1+pe_strain).*k1_pe)./((pe_strain.*lpe_slack).^2);

tmp=ones(size(lse));
% se characteristics
idx=lse>lse_slack;
fse = k1_se.*lse;
kse = k1_se.*tmp;
fse(idx)=fse(idx) + k2_se.*(lse(idx)-lse_slack).^2; 
kse(idx)=kse(idx) + 2*k2_se.*(lse(idx)-lse_slack); 

% pe characteristics
idx=lpe>lpe_slack;
fpe = k1_pe.*lpe;
kpe = k1_pe.*tmp;
fpe(idx)=fpe(idx) + k2_pe.*(lpe(idx)-lpe_slack).^2;
kpe(idx)=kpe(idx) + 2*k2_pe.*(lpe(idx)-lpe_slack);
return


function [stim,lmtc,lmtcd]=kinematic_model_input(t,parms)
% stim
if t<=.1
    stim=parms.gamma0; % we start in steady state
elseif t>.1 && t<2 % full activation
    stim=1; 
else
    stim=0.2; % relaxation to low value
end
lmtc0=parms.lmtc0;
lceopt=parms.lceopt;
% mtc length
% isometric at lmtc0
lmtc=lmtc0+0.2*lceopt*sin(2*pi*t);
lmtcd=0.2*2*pi*lceopt*cos(2*pi*t);
end
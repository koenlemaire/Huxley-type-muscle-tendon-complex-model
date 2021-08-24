function [ dstatedt,varargout ] = hux_tutorial_dynamic_quasi( t,state,parms )
%function [ dstatedt,out,check,x,n,dndt ] = hux_tutorial( t,state,parms )
%   INPUT: state [n gamma lce]
%   OUTPUT: dstatedt: time derivitave of state
%   (output args 2-6 are optional)
%   out: [stim gamma fce fpe fse lmtc lmtcdot Edot];
%   check: [dfcedt dfpedt dfsedt regime clrX(1) clrX(end)]
%   x, n and dndt: only the relevant entries of the respective vectors
%
% Huxley cross bridge dynamics modified from Zahalak (1981) eq. 3:
% (dndt)x - v*(dndx)t = q*fisom*f(x) - [f(x) + g(x)]*n                      (1)
% parametrisation by methods of characteristics with definitions:
% a(x,t) = v(t);
% b(x,t) = 1;
% c(x,t) = -[f(x,t) + g(x)];
% d(x,t) = q*fisom*f(x)
% the system is transformed into 3 ODE's:
% dxds = a(x,t) = v(t), x(0) = x0;
% dtds = b(x,t) = 1, t(0) = t0;
% dnds = c(x,t)*n + d(x,t) = -[f(x,t) + g(x)]*n + q*fisom*f(x), n(0) = u0.
% the second equation is dumped because t = s. x0 is the initial domain of
% n, n0 is the inital state of n. For details not listed here, please see
% the functions called within this function.
%
% The complete model is described in:
% Lemaire, K. K., Baan, G. C., Jaspers, R. T., & van Soest, A. J. (2016).
% Comparison of the validity of Hill and Huxley muscle tendon complex
% models using experimental data obtained from rat m. Soleus in situ.
% Journal of Experimental Biology, 219, 977-987. DOI: 10.1242/jeb.128280   

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

%% unravel state vector
state=state(:);
n = state(1:end-4); % [] fraction of bound cross bridges
gamma = state(end-3); % [] relative free Ca2+ concentration
lce = state(end-2); % [m] contractile element length
lmtc=state(end-1); % [m] muscle tendon complex length
lmtcd=state(end); % [m/s] time derivative 
%% read out parameters
scale_factor = parms.scale_factor; % [] scale factor between lcerel and x
x0 = parms.x0; % initial vector x0 at t=0
lce0 = parms.lce0; % [m] ce length at t=0
lceopt = parms.lceopt; % [m] optimum ce length
dndt = parms.dndt; % =zeros(size(x0))
rateFun=parms.rateFun; % fx/gx rate function
c_act=parms.c_act;
c_cb=parms.c_cb;
%% model input
if t<=.1
    stim=parms.gamma0; % we start in steady state
elseif t>.1 && t<2 % full activation
    stim=1; 
else
    stim=0.2; % relaxation to low value
end
%% calculate muscle components lengths
lcerel=lce./lceopt; % [] relative CE length
lse = lmtc - lce;  % [m] SE length
%% calculate gammad and q
gammad = gammadot(gamma,stim,parms); 
q = activeState(gamma,parms); % [] relative Ca2+ bound to troponin
%% create current bond lengths (x) vector
dlcerel = (lce - lce0)/lceopt; % [] difference of current lce to lce0 scaled to lcerel
x = x0 + dlcerel*scale_factor; % [] update x0 to current x 
%% select relevant part of x/n vector
iRel = (x<2 & x>-1) | abs(n)>0; % these are the values where dndt~=0 AND/OR n~=0
xRel = x(iRel); % define xRel and nRel
nRel = n(iRel);
%% check sparsity assumption
% we need to make sure the distribution does not 'run off', ie the lowest
% and highest values of the relevant part of x may not exceed the highest
% and lowest values of the total x vector
clrX = [xRel(1)-x(1) x(end)-xRel(end)]; % clearance between edges of xRel and x
if min(clrX) < .1 % now we are too close to the edge!
    err=true; 
else
    err=false;
end
%% calculate fisomrel
[fisomrel] = ce_fl_simple(lcerel,parms); % [] relative isometric CE force
%% calculate SE and PE force and instantanious stiffness
[fse, fpe, kse, kpe] = CEEC_simple(lse,lce,parms); % [N N N/m N/m]
%% calculate f(x), g(x) and dndt
[fx,gx]=rateFun(xRel(:)); % see function for details
dndtRel=fisomrel*q*fx-(fx+gx).*nRel; % see Lemaire et al. 2016 for details, now incorporates both q and fisom
dndt(iRel)=dndtRel; % update dndt with new (nonzero) values
%% calculate lced
% help variables:
kf=parms.k_f; % [N/h] scaling between distribution and force 
int_nx  = sum(xRel.*nRel)*kf; % CE force [N]
int_n   = sum(nRel)*kf; % CE stifness [N/h]
int_dnx = sum(xRel.*dndtRel)*kf; % [N/s]

% now calculate lced (see Lemaire et al. 2016 for details)
lced = (-int_dnx + lmtcd*kse) ./ (int_n*scale_factor/lceopt + kse + kpe); % [m/s]

%% calculate lmtcdd
Fz=-9.81*parms.mass; % [N] gravitational force
lmtcdd=-(fse+Fz)/parms.mass; % [m/s^2] Newtons law, minus sign because moving up means muscle shortening

%% complete stated
dstatedt = [dndt; gammad; lced; lmtcd; lmtcdd];
%% calculate optional output parameters and error handling
if nargout > 1 || err == true
    % NOTE: for all output parameters same notes as in "calucalate lced"
    % section holds!!
    fce = int_nx; % [N]
    p_act=c_act*gamma; % [W] metabolic power for calcium pumping
    p_cb=c_cb*sum(gx.*nRel); % [W] metabolic power for CB cycling
    
    varargout{1}=[stim q fce fpe fse fisomrel p_cb p_act];
    if nargout>2
        dfcedt = int_dnx + scale_factor*int_n*lced/lceopt; % [N/s]
        dfsedt = kse*(lmtcd-lced); % [N/s]
        dfpedt = kpe*lced; % [N/s]
        varargout{2}=[dfcedt dfpedt dfsedt];% clrX(1) clrX(end)];
        if nargout>3
            varargout{3}=xRel;
            if nargout>4
                varargout{4}=nRel;
                if nargout>5
                    varargout{5}=dndtRel;
                end
            end
        end
    end
end
if err == true && nargout > 1 % we are not evaluating jacobian!!
    warning('sparsity assumption possibly voilated, check results')
    keyboard
end
return

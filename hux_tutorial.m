function [ dstatedt,varargout ] = hux_tutorial( t,state,parms )
%function [ dstatedt,out,check,x,n,dndt ] = hux_tutorial( t,state,parms )
%   INPUT: state [n q lce]
%   OUTPUT: time derivitave of state
%   (output args 2-6 are optional)
%   out: [stim gamma fce fpe fse lmtc lmtcdot Edot];
%   check: [dfcedt dfpedt dfsedt regime clrX(1) clrX(end)]
%   x, n and dndt: only the relevant entries of the respective vectors
%
% Huxley model modified from Zahalak (1981) eq. 3:
% (dndt)x - v*(dndx)t = q*fisom*f(x) - [f(x) + g(x)]*n                      (1)
% parametrisation by methods of characteristics results in
% a(x,t) = v(t);
% b(x,t) = 1;
% c(x,t) = -[f(x,t) + g(x)];
% d(x,t) = q*fisom*f(x)
% the system is then transformed into 3 ODE's:
% dxds = a(x,t) = v(t), x(0) = x0;
% dtds = b(x,t) = 1, t(0) = t0;
% dnds = c(x,t)*n + d(x,t) = -[f(x,t) + g(x)]*n + q*fisom*f(x), n(0) = u0.
% the second one is dumped because t = s. x0 is the initial domain of n, n0
% is the inital state of n. For details not listed here, please see the
% functions called within this function.
%
% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

% The complete model is described in:
% Lemaire, K. K., Baan, G. C., Jaspers, R. T., & van Soest, A. J. (2016).
% Comparison of the validity of Hill and Huxley muscle tendon complex
% models using experimental data obtained from rat m. Soleus in situ.
% Journal of Experimental Biology, 219, 977-987. DOI: 10.1242/jeb.128280   


%% unravel state
state=state(:);
n = state(1:end-4);
gamma = state(end-3);
lce = state(end-2);
lmtc=state(end-1);
lmtcd=state(end);
%% read out parameters
scale_factor = parms.scale_factor; % [] scale factor between lcerel and x
x0 = parms.x0; % initial vector x0 at t=0
lce0 = parms.lce0; % [m] ce length at t=0
lceopt = parms.lceopt; % [m] optimum ce length
dndt = parms.dndt; % =zeros(size(x0))

%% calculate model inputs
if t<=.2
    stim=parms.gamma0; % we start in steady state
elseif t>.2 && t<1.5 % full activation
    stim=1; 
else
    stim=0.2; % relaxation to low value
end
%% calculate muscle components lengths
%lse_slack = parms.lse_slack; % [m] SEE slack length
lcerel=lce./lceopt; % [] relative CE length
lse = lmtc - lce;  % [m] SE length
%% calculate gammad and q
gammad = gammadot(gamma,stim,parms);
q = activeState(gamma,parms);
%% create current x vector
dlcerel = (lce - lce0)/lceopt; % [] difference of current lce to lce0 scaled to lcerel
x = x0 + dlcerel*scale_factor; % [] update x0 to current x
%% select relevant part of x/n vector
iRel = x<1.2 & x>-.2 | abs(n)>1e-16; % these are the values where dndt~=0 AND/OR n~=0
xRel = x(iRel); % define xRel and nRel
nRel = n(iRel);
%% check sparsity assumption
clrX = [xRel(1)-x(1) x(end)-xRel(end)]; % clearance between edges of xRel and x
if min(clrX) < 1.5 % now we are too close to the edge!
    err=true; 
else
    err=false;
end
%% calculate fisomrel
[fisomrel] = ce_fl_simple(lcerel,parms);
%% calculate SE and PE instantanious stiffness
% see function file for details
[fse, fpe, kse, kpe] = CEEC_simple2(lse,lce,parms); % [N] and [N/m]
%% calculate f(x), g(x) and dndt
%parms.hux.f1=parms.hux.f1*q; % enable this line and remove factor q from line dndtRel = ... to set this function to bogus form
[fx,gx]=rateFunc_v5(xRel(:),parms); % see function for details
dndtRel=fisomrel*q*fx-(fx+gx).*nRel; % see thesis Koen Lemaire for details, now incorporates both q and fisom
dndt(iRel)=dndtRel; % update dndt with new (nonzero) values
%% calculate lced
kf=parms.k_f;
int_nx  = sum(xRel.*nRel)*kf; % CE force [N]
int_n   = sum(nRel)*kf; % CE stifness [N/h]
int_dnx = sum(xRel.*dndtRel)*kf; % [N/s]

% now calculate lced (see thesis Koen Lemaire for math details)
lced = (-int_dnx + lmtcd*kse) ./ (int_n*scale_factor/lceopt + kse + kpe); % [m/s]

%% calculate lmtcdd
Fz=-9.81*parms.mass; % [N]
lmtcdd=-(fse+Fz)/parms.mass; % [m/s^2] Newton, minus sign as moving up means muscle shortening

%% complete stated
dstatedt = [dndt; gammad; lced; lmtcd; lmtcdd];
%% calculate optional output parameters and error handling
if nargout > 1 || err == true
    % NOTE: for all output parameters same notes as in "calucalate lced"
    % section holds!!
    fce = int_nx; % [N]
    varargout{1}=[stim q fce fpe fse];
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
    figure
    plot(xRel,nRel,xRel,xRel.*nRel,'--')
    legend('n over x','n*x over x')
    title('n over x at current integration step')
    xlabel('x')
    ylabel('n')
    warning('sparsity assumption possibly voilated, see figure')
    keyboard
end
return

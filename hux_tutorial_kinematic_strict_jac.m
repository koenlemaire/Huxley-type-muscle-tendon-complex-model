function [ jac ] = hux_tutorial_kinematic_strict_jac( t,state,parms )
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
n = state(1:end-1); % [] fraction of bound cross bridges
gamma = state(end); % [] relative free Ca2+ concentration

%% read out parameters
scale_factor = parms.scale_factor; % [] scale factor between lcerel and x
x0 = parms.x0; % initial vector x0 at t=0
lce0 = parms.lce0; % [m] ce length at t=0
lceopt = parms.lceopt; % [m] optimum ce length
dndt = parms.dndt; % =zeros(size(x0))
rateFun=parms.rateFun; % fx/gx rate function
c_act=parms.c_act;
c_cb=parms.c_cb;
kf=parms.k_f; % [N/h] scaling between distribution and force 

%% model input
[stim,lmtc,lmtcd]=kinematic_model_input(t,parms);

%% calculate lce given current state
% already done in other function, in which the value of lce is overridden,
% we don't want to do that here because of potential unexpected outcome ...
global CURRENT_LCE
lce=CURRENT_LCE;
lcerel=lce./lceopt; % [] relative CE length
%% calculate gammad and q
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
%% calculate f(x), g(x) and dndt
[fx,gx]=rateFun(xRel(:)); % see function for details
dndtRel=fisomrel*q*fx-(fx+gx).*nRel; % see Lemaire et al. 2016 for details, now incorporates both q and fisom
dndt(iRel)=dndtRel; % update dndt with new (nonzero) values

%% calculate Jacobian, using sparse matrix
% initialize sparse jacobian; generate NxN sparse matrix with at most 2N
% entries. We might need to compute the complete jacobian (without skipping
% places where ndot is zero, i.e. using the iRel trick), but I'm not sure
% about that, my gut tells me this is fine. Test later. 
N=length(state);
jac=spalloc(N,N,2*sum(iRel)+1);
idx=1:N; % help variable for indexing

% dndt wrt n
dndt_wrt_n=-(fx+gx);
row_idx_dndt_wrt_n=idx(iRel);
col_idx_dndt_wrt_n=idx(iRel);

% dndt wrt gamma
n=parms.n; 
k=parms.k; 
dq_dgamma=((1+k.^n).*n.*gamma.^(n-1).*k.^n) ./ (gamma.^n + k.^n).^2; % simplified version KKL
dndt_wrt_gamma=fisomrel.*fx.*dq_dgamma;
row_idx_dndt_wrt_gamma=idx(iRel);
col_idx_dndt_wrt_gamma=ones(1,sum(iRel))*N(end);

% gammad wrt gamma
% simplified form, use only for scalar input
if stim>=gamma
    gammad_wrt_gamma = -1/parms.tau_act;
else
    gammad_wrt_gamma = -1/parms.tau_deact;
end
row_idx_gammad_wrt_gamma = N;
col_idx_gammad_wrt_gamma = N;
% gammad wrt n
% this is zero, and thus does not need to be replaced in the sparse
% jacobian ... 

% combine all parts of jacobian
total_row_idx=[row_idx_dndt_wrt_n row_idx_dndt_wrt_gamma row_idx_gammad_wrt_gamma];
total_col_idx=[col_idx_dndt_wrt_n col_idx_dndt_wrt_gamma col_idx_gammad_wrt_gamma];
total_jac_vec=[dndt_wrt_n; dndt_wrt_gamma; gammad_wrt_gamma];

jac=sparse(total_row_idx,total_col_idx,total_jac_vec,N,N);
return

% we are going to simulate a huxley MTC suspended from the ceiling with a
% mass attached to it. Positive direction is upward (1 mechanical DOF).
% this version resembles a standard neuromusculoskeletal simulation in
% terms of its initial condition settings. 

% This file released under the terms of the GNU General Public License,
% version 3. See http://www.gnu.org/licenses/gpl.html 
% Author: KK Lemaire (kklemaire_edu@posteo.nl)

% The complete model is described in:
% Lemaire, K. K., Baan, G. C., Jaspers, R. T., & van Soest, A. J. (2016).
% Comparison of the validity of Hill and Huxley muscle tendon complex
% models using experimental data obtained from rat m. Soleus in situ.
% Journal of Experimental Biology, 219, 977-987. DOI: 10.1242/jeb.128280   



clear
close all
clc

addpath(genpath(cd)) % make sure all subfunctions are included
%% set parameters
% simulation processing:
diagnostics=true; % figures and data for sanity checking results
normFig=true; % figures for viewing data
Animate=true; % do animation

% general:
parms.Fmax=1000; % [N]
parms.lceopt=0.05; % [m] CE optimum length
parms.lpe_slack=1.1*parms.lceopt; % [m] PE slack length
parms.lse_slack=0.25; % [m] SE slack length
parms.se_strain=.05; % [N/m^2] SE shape, Fse=Fmax at 4% strain
parms.pe_strain=.5; % [N/m^2] PE shape, Fpe=0.5*Fmax at 4% strain????
width=0.56; % [] width of force length relation
parms.width=width;
parms.C=-(1/width)^4;

% activation dynamics:     %% DONT MATTER FOR INITIAL SS CONDITION
parms.k=.35; % 50% point in q-gamma relation
parms.n=2; % curvature in q-gamma relation
parms.tau_act=0.080; %[s] rising time constant of gamma(stim) dynamics
parms.tau_deact=0.100; %[s] falling time constant of gamma(stim) dynamics
parms.qmin=1e-10; % minimum possible value for q
parms.gamma_min=parms.qmin;

% huxley model:
h = 1e-8;           % attachment 'range' for myosin head [m]
s = 2.6e-6;         % sarcomere length [m]
parms.scale_factor = s/(2*h); % [] scaling between x and lcerel
parms.dx=.05; % [h] stepsize in x
parms.g1=200; % [Hz] detachment rate parameter
parms.f1=800; % [Hz] attachment rate parameter
parms.g2=3000; % [Hz] detachment rate parameter
parms.g3=1400; % [Hz] detachment rate parameter

% rate function form; see different supplied functions
parms.rateFun=@(x)rateFunc_v8(x,parms);
%% Initial conditions

% -calculate (or set, in this case) lmtc0, fse0 and lse0
% -calculate lce0 and fisom0
% -calculate q0 (note that F_SS = fisom*q*Fmax, just as in Hill mod ...)
% -invert q(gamma) to find gamma0
% -set n0(x0) to SS distribution, corresponding to q0 and fisom0, thus:
% n0 = nSS*fisom0*q0

% set F0 and mass:
parms.fse0=0.5*parms.Fmax; % [N]
parms.mass=parms.fse0/9.81; % [kg]

% set lmtc0, activation follows from this
lmtc0=(1+parms.se_strain)*parms.lse_slack+1.1*parms.lceopt; % [m]
lmtcd0=0; % [m/s]

[lse0]=fzero(@Fse_inverse,parms.lse_slack*[1 (1+parms.se_strain)],[],parms); % [m]

% compute lce0
lce0 = lmtc0-lse0; % [m]
parms.lce0=lce0;
lcerel0=lce0/parms.lceopt; % []

% calculate fisom:
[ fisomrel0,~] = ce_fl_simple( lcerel0,parms ); % []

% bottom code if we wish to start in steady state, otherwise choose gamma0
% freely
% calculate q0:
%parms.q0=activeState(gamma0,parms);
[fse0, fpe0, ~, ~] = CEEC_simple2(lse0,lce0,parms);
parms.q0=(fse0-fpe0)/(parms.Fmax*fisomrel0);

% invert q(gamma) to find gamma0 (fzero has tolX = 1e-16 ...):
[gamma0,~]=fzero(@activeState_inverse,[parms.qmin 1],[],parms);
parms.gamma0=gamma0;

% now set the domain for x0

% first determine boundaries in terms of lce lengths, because If
% u=scale_factor*lcereldot than xi=scale_factor*lcerel. We are going to use
% this to calculate maximum values for x1 and x2, so that during
% optimization the distribution can not 'run off'.

% NOTE: positive direction of lce is to the right!
% NOTE: The smallest value lcerel can attain is the left bound of the lce
% domain. The largest is the right bound. The LEFT bound of the xi domain
% is given by the difference between lcerel0 and the RIGHT bound of the
% lcerel domain, scaled to the xi domain. Vica versa for the right bound of
% the xi domain. See bottom script!!

lce_L = 0.2*parms.lceopt; % [m] smallest value lce is expected to attain
lce_R = 1.8*parms.lceopt; % [m] largest value lce is expected to attain

x1 = round((lce0-lce_R)*parms.scale_factor/parms.lceopt); % [h] left bound x
x2 = round((lce0-lce_L)*parms.scale_factor/parms.lceopt); % [h] right bound x

% initial state vector for x, with stepsize dx
x0 = (x1:parms.dx:x2)';
parms.x0 = x0; % [h]

% NOTE: n0 is assumed Steady State
% first compute f0 and g0
[fx0,gx0]=parms.rateFun(x0(:)); % [Hz]

% ss condition is given by f(f+g)|f(x)>g(x)
nSS = zeros(length(x0),1); % initiate n0 for indexing by find
parms.dndt=nSS; % set base for dndt vector!

i = find(fx0>0);
nSS(i) = fx0(i)./(fx0(i)+gx0(i));  % inital state vector for n0, isometric ss condition is given by f/(f+g) | f(x)>0!!
n0 = nSS*parms.q0*fisomrel0; % n0 such that fce0 + fpe0 = fse0

figure;plot(x0,fx0,x0,gx0)
title('rate parameter functions')
xlabel('x [h]')
ylabel('rate [Hz]')
xlim([-2 2])
legend('attachment','detachment')
figure;plot(x0,nSS); axis([-2 2 0 1])
title('steady state isometric n(x)')
xlabel('x [h]')
ylabel('n []')

% Scaling factor for the moments of n(x,t) to force
k_f=parms.Fmax/(sum(x0.*nSS));
parms.k_f=k_f; % [N/h]
%% run simulation
state0 = [n0' gamma0 lce0 lmtc0 lmtcd0];

[stated0,y0,check0,xRel0,nRel0,dndtRel0] = hux_tutorial(0,state0,parms);
t_end=4; % [s] simulation time, starting at t=0 ...

tSpan=[0:.01:t_end]; % chop up time to save memory (reduce state vector)
ode_fun=@(t,state)hux_tutorial(t,state,parms);
odeparms=odeset('abstol',1e-8,'reltol',1e-8,'maxstep',.02);
[t,state] = ode45(ode_fun,tSpan,state0,odeparms);

% unravel state
n = state(:,1:end-4);
gamma = state(:,end-3);
lce= state(:,end-2); % [m]
lmtc= state(:,end-1); % [m]
lmtcd= state(:,end); % [m]
lse=lmtc-lce; % [m]

% initialise stated and optional output
check = zeros(length(check0),length(t));
stated=zeros(size(state'));
y = zeros(length(y0),length(t));
% output block
if diagnostics == true
    figure
end
for i=1:length(t)
    [stated(:,i),y(:,i),check(:,i),x,n,dndt] = hux_tutorial(t(i),state(i,:)',parms);
    if diagnostics==true
        plot3(x,ones(size(x))*t(i),n);hold on
        xlabel x; ylabel t; zlabel n
    end
end
% handle output
stated=stated'; y=y';
gammad = stated(:,end-1);
lced = stated(:,end); %[m/s]


% unravel y
Fmax=parms.Fmax;
stim = y(:,1);
q = y(:,2);
fce=y(:,3); %[N]
fpe=y(:,4); %[N]
fse=y(:,5); %[N]
fcerel = fce/parms.Fmax;
fperel = fpe/parms.Fmax;
fserel = fse/parms.Fmax;
clear y
lcerel = lce/parms.lceopt;
lserel = lse/parms.lse_slack;
lmtcrel=lmtc/(parms.lceopt+(1+.05)*parms.lse_slack);

% calculate energy/work terms:
ce_work=cumtrapz(-lce,fce); %[J]
se_work=cumtrapz(-lse,fse); %[J]
pe_work=cumtrapz(-lce,fpe); %[J]
lmtc_work=cumtrapz(-lmtc,fse); %[J]
Ekin = .5*parms.mass*lmtcd.^2; %[J]
Ekin = Ekin-Ekin(1);
W_grav = -parms.mass*(lmtc-lmtc(1))*-9.81; %[J]



%% diagnostics
if diagnostics == true
    % error and other checks
    check=check';
    dforceErrordt = -check(:,1)-check(:,2)+check(:,3);
    forceError = fse-fpe-fce;
    lengthError = lmtc - lse - lce;
    clear check
    figure
    subplot(221)
    plot (t,forceError)
    title(['force error'])
    xlabel ('time [s]')
    ylabel ('F_c_e + F_p_e - F_s_e [N]')
    subplot(222)
    plot (t,dforceErrordt)
    title('derivative of force error')
    xlabel ('time [s]')
    ylabel ('dF_c_e/dt + dF_p_e/dt - dF_s_e/dt [N/s]')
    subplot(223)
    plot (t,[ce_work+se_work+pe_work-lmtc_work Ekin-W_grav-lmtc_work])
    title('energy errors')
    xlabel ('time [s]')
    ylabel ('energy [J]')
    legend('w_c_e + w_p_e + w_s_e - w_m_t_c','Ekin - w_m_t_c - w_F_z')
    subplot(224)
    plot (t,lengthError)
    title('length error')
    xlabel ('time [s]')
    ylabel ('l_mt_c - l_c_e - l_s_e [m]')
end

%% standard output figs
if normFig==true
    figure
    subplot(221)
    plot(t,[stim gamma q])
    axis([0 t(end) 0 1])
    legend('STIM','gamma','q')
    title('Stimulation') %'stimulation as set in this protocol')
    xlabel ('Time [s]')
    ylabel ('Normalized STIM/Active state')
    %         grid
    subplot(222)
    plot(t,[lce lse lmtc])
    title('Muscle components lengths')
    legend('ce','se','mtc')
    xlabel ('Time [s]')
    ylabel ('length [m]')
    grid
    subplot(223)
    plot(t,[fse fpe fce])
    title('forces')
    legend('se','pe','ce')
    xlabel('Time [s]')
    ylabel('Force [N]')
    grid
    subplot(224)
    plot(t,[ce_work pe_work se_work lmtc_work Ekin W_grav])
    title('mechanical energy work')
    ylabel [J]
    xlabel('Time [s]')
    legend('w_c_e','w_p_e','w_s_e','w_m_t_c','E_k_i_n','w_F_z')
end

%% animation
if Animate
    % for video ...
    % writerObj = VideoWriter('hux_tutorial.avi');
    % writerObj.FrameRate = 100;
    % open(writerObj)
    
    figure
    for iSample=1:length(t)
        plot([-.5 .5]*lmtc0,[0 0],'k','linewidth',1);hold on
        nline=12;
        tmp=linspace(-.5*lmtc0,.5*lmtc0,nline);
        for iLine=1:nline
            plot([tmp(iLine) tmp(iLine)+.04*lmtc0],[0 .04*lmtc0],'k','linewidth',1)
        end
        set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        box off
        ylim([-.4 0.1])
        axis equal
        text(.01,-.5*lce(iSample),'CE')
        text(.01,-.5*(lmtc(iSample)+lce(iSample)),'SEE')
        
        text(-.13,-.1,['t=',num2str(t(iSample),3)])
        text(-.17,-.125,['STIM=',num2str(stim(iSample),3)])
        
        plot([0 0],[0 -lce(iSample)],'r','linewidth',2);
        plot([0 0],[-lce(iSample) -lmtc(iSample)],'b','linewidth',2);
        plot(0,-lmtc(iSample),'ko','linewidth',8); hold off
        drawnow
        
        % for video
        % currFrame=getframe(gcf);
        % writeVideo(writerObj,currFrame)
    end
    % close(writerObj)
end
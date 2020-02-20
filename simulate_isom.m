%% Huxley MTC model tutorial
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

% what have we learned today:
% -Why is activation cost constant??
% since there is an algebraic relation between gamma and q, the dynamics
% between stim and gamma don't matter the only thing that matters is the
% dynamics between q and force. This is because there is a one to one
% algebraic relation between stim and gamma, and gamma is the only thing
% that determines calcium pumping cost, ie gamma=c*calcium cost. stim does
% not really matter as long as you can find stim that satisfies q(t) given
% gamma ... above is true for hill and huxley, where in Huxley model there
% is additional dynamics (apart from CE-SE interacion) in the relationship 
% between q and force. Apparently, the CE-SE interaction does not add much
% dynamics (in fast Huxley dynamics simulation the q-force relation is
% straight line ...) in the Hill case. 
% -why is CB cost constant??
% 

clear
close all
clc

% sort out pathing and stuff
tmp = mfilename('fullpath');
tmp = tmp(1:length(tmp)-length(mfilename));
cd(tmp)
addpath(genpath(cd)) % make sure all subfunctions are included
%% set parameters
% simulation processing:
parms.diagnostics=true; % figures and data for sanity checking results
parms.normFig=true; % figures for viewing data
makeFig=true; % false runs optimization; true shows figures
parms.opt=true;

% general:
parms.Fmax=1; %[N]
parms.lceopt=0.081; % [m] CE optimum length
parms.lpe_slack=1.1*parms.lceopt; % [m] PE slack length
parms.lse_slack=0.34; % [m] SE slack length
parms.se_strain=.05; % [N/m^2] SE shape, Fse=Fmax at 4% strain
parms.pe_strain=.2; % [N/m^2] PE shape, Fpe=0.5*Fmax at 4% strain????
width=0.56; % [] width of force length relation
parms.width=width;
parms.C=-(1/width)^4;

% activation dynamics curtin:     %% DONT MATTER FOR INITIAL SS CONDITION
% parms.k=.35; % 50% point in q-gamma relation
% parms.n=2; % curvature in q-gamma relation
% parms.tau_act=0.080; %[s] rising time constant of gamma(stim) dynamics
% parms.tau_deact=0.100; %[s] falling time constant of gamma(stim) dynamics
% parms.qmin=1e-10; % minimum possible value for q
% parms.gamma_min=parms.qmin;

% activation dynamics hatze
parms.v=3; % 3 according to Hatze/Rockenfeller
parms.lp=2.9; % 2.9 according to rockenfeller
% rho_c=[Ca]max * rho0 (ie molar volume); 1.373e-4 * 55700 (for v=3)
% according to Hatze. But rho0 ranges from 10e-6 to 16e-6 according to
% Rockenfeller (and some extent Kistemaker 2007). rho_c = 7 according to
% rockenfeller figure B4 but seems more likely to be 0.7 given the above
% and the discussion in Rockenfeller appendix C
parms.rho_c=1.373e-4 * 55700; % ~7?? ??
parms.m=11.25; % [fast twitch, 1/3.67 for slow twitch]
parms.qmin=1e-10; % minimum possible value for q
%parms.gamma_min=parms.qmin;

% huxley model:
h = 1e-8;           % attachment 'range' for myosin head [m]
s = 2.6e-6;         % sarcomere length [m]
tmpscale=1; % for slow simulations
parms.scale_factor = s/(2*h)/tmpscale; % [] scaling between x and lcerel
parms.dx=.05; % [h] stepsize in x
tmp=1.0e+03 * [0.8890    0.4275    3.1703    0.7796];
parms.g1=tmp(1)/tmpscale; % [Hz] detachment rate parameter
parms.f1=tmp(2)/tmpscale; % [Hz] attachment rate parameter
parms.g2=tmp(3)/tmpscale; % [Hz] detachment rate parameter
parms.g3=tmp(4)/tmpscale; % [Hz] detachment rate parameter

% energetics parms scaling: need to fit/optimize these

% rate function form; see different supplied functions
parms.rateFun=@(x)rateFunc_v8(x,parms);
%% Initial conditions

% -calculate (or set, in this case) lmtc0, fse0 and lse0
% -calculate lce0 and fisom0
% -calculate q0 (note that F_SS = fisom*q*Fmax, just as in Hill mod ...)
% -invert q(gamma) to find gamma0
% -set n0(x0) to SS distribution, corresponding to q0 and fisom0, thus:
% n0 = nSS*fisom0*q0

% set Fse0 :
parms.fse0=0.2*parms.Fmax; % [N]

% set lmtc0, activation follows from this
lmtc0=1.02*parms.lse_slack+parms.lceopt; % [m]
parms.lmtc0=lmtc0;
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
zeroFun=@(x)activeState_inverse_hatze(x,lcerel0,parms);
[gamma0,~]=fzero(zeroFun,[0 1]);
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

lcerel_L = 0.7; % [m] smallest value lce is expected to attain
lcerel_R = 1.4; % [m] largest value lce is expected to attain

x1 = (lcerel0-lcerel_R)*parms.scale_factor-1; % [h] left bound x
x2 = (lcerel0-lcerel_L)*parms.scale_factor+2; % [h] right bound x

% initial state vector for x, with stepsize dx
x0 = (x1:parms.dx:x2)';
parms.x0 = x0; % [h]

% NOTE: n0 is assumed Steady State
% first compute f0 and g0
[fx0,gx0]=parms.rateFun(x0(:)); % [Hz]

% ss condition is given by f/(f+g)|f(x)>g(x)
nSS = zeros(length(x0),1); % initiate n0 for indexing by find
parms.dndt=nSS; % set base for dndt vector!

i = find(fx0>0);
nSS(i) = fx0(i)./(fx0(i)+gx0(i));  % inital state vector for n0, isometric ss condition is given by f/(f+g) | f(x)>0!!
n0 = nSS*parms.q0*fisomrel0; % n0 such that fce0 + fpe0 = fse0

% figure;plot(x0,fx0,x0,gx0)
% title('rate parameter functions')
% xlabel('x [h]')
% ylabel('rate [Hz]')
% xlim([-2 2])
% legend('attachment','detachment')
% figure;plot(x0,nSS); axis([-2 2 0 1])
% title('steady state isometric n(x)')
% xlabel('x [h]')
% ylabel('n []')

% Scaling factor for the moments of n(x,t) to force
k_f=parms.Fmax/(sum(x0.*nSS));
parms.k_f=k_f; % [N/h]
parms.c_act=1;
parms.c_cb=(sum(x0.*nSS.*gx0));

state0 = [n0' gamma0 lce0];
parms.state0=state0;

% scaling factor for energetics
% loop over frequencies

freqs=[0.5 1 2];
for iFreq=1:length(freqs)    
    parms.freq=freqs(iFreq);  % [stim frequency Hz]
    parms.iFreq=iFreq;    
    parms.n_cycli=2;
    tspan=linspace(0,parms.n_cycli/parms.freq,100*parms.n_cycli);
    parms.tspan=tspan;
    n_control=11;
    parms.desiredForce=parms.Fmax*0.2*(sin(2*pi*parms.freq*tspan)+1);
    
    % bounds and initial conditions    
    
    %% figures
    if ~makeFig % if makeFig==false
        lb=zeros(1,n_control);
        ub=ones(1,n_control)*.4;
        %initial_guess=(lb+ub)/2;
        
        tmp=load('control');
        initial_guess=tmp.control(iFreq,:);
                
        optimFun=@(x)simulate_isometric_contraction(x,parms);
                
        fmins_opt=optimset('fminsearch');
        n_hours=.5;
        fmins_opt.TolFun=1e-10;
        fmins_opt.TolX=1e-10;
        fmins_opt.MaxTime=n_hours*3600; % [s]
        fmins_opt.MaxFunEvals=500;
        fmins_opt.Display='iter';
        
        [control(iFreq,:),fval(iFreq,1),flag,out] = fminsearch(optimFun,initial_guess,fmins_opt);
        save('tmpLastControl','control')
    else
        load control
        parms.opt=false;
        x=control(iFreq,:);
        parms.nFreqs=length(freqs);
        
        err = simulate_isometric_contraction(x,parms);
    end
    
end

function err = simulate_isometric_contraction(x,parms)
%% run simulation

control=abs(x);
ode_fun=@(t,state)hux_tutorial_kinematic(t,state,control,parms);
state0=parms.state0;
tspan=parms.tspan;
% test odefun
[stated0,y0,check0,xRel0,nRel0,dndtRel0] = ode_fun(0,state0);

odeparms=odeset('abstol',1e-8,'reltol',1e-8,'maxstep',.02);
%tic
[t,state] = ode45(ode_fun,tspan,state0,odeparms);
%toc

% unravel state
n = state(:,1:end-2);
gamma = state(:,end-1);
lce= state(:,end); % [m]
%lmtc= state(:,end-1); % [m]
%lmtcd= state(:,end); % [m]


% initialise stated and optional output
check = zeros(length(check0),length(t));
stated=zeros(size(state'));
y = zeros(length(y0),length(t));
% output block

if parms.opt==true % optimizing, minimal output
    for i=1:length(t)
        [~,y(:,i)] = ode_fun(t(i),state(i,:)');
    end
    fse=y(5,:); %[N]
else % not optimizing, full scale output
    if parms.diagnostics == true
        figure
    end    
    for i=1:length(t)
        [stated(:,i),y(:,i),check(:,i),x,n,dndt] = ode_fun(t(i),state(i,:)');
        if parms.diagnostics==true && mod(i-1,5)==0 % every 5 samples
            plot3(x,ones(size(x))*t(i),n);hold on
            xlabel x; ylabel t; zlabel n
        end
    end
    xlim([-3 3])
    ylim([0 t(end)])
    zlim([0 1])
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
    fisomrel=y(:,6); %[]
    p_cb=y(:,7); %[W]
    p_act=y(:,8); %[W]
    lmtc=y(:,9); %[m]
    lmtcd=y(:,10); %[m/s]
    lse=lmtc-lce; % [m]
    
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
    %Ekin = .5*parms.mass*lmtcd.^2; %[J]
    %Ekin = Ekin-Ekin(1);
    %W_grav = -parms.mass*(lmtc-lmtc(1))*-9.81; %[J]
    
    %% diagnostics
    if parms.diagnostics == true
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
        plot (t,[ce_work+se_work+pe_work-lmtc_work])
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
    if parms.normFig==true
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
        %yyaxis left
        plot(t,[ce_work pe_work se_work lmtc_work])
        title('mechanical energy/work')
        ylabel ('mechanical work [J]')
        xlabel('Time [s]')
        legend('w_c_e','w_p_e','w_s_e','w_m_t_c','E_k_i_n','w_F_z')
        
        figure
        plot(t,[p_cb p_act])
        title('metabolic power')
        ylabel ('metabolic power [normalized]')
        xlabel('Time [s]')
        legend('P_c_b','P_a_c_t')
        
        
        figure
        subplot(221)
        plot(lced/parms.lceopt,fce/parms.Fmax,'.')
        title('normalized CE force vs normalized CE velocity, from simulation result')
        ylabel ('CE force [Fmax]')
        xlabel('CE velocity [lceopt/s]')
        
        subplot(222)
        plot(lced/parms.lceopt,fce./q./fisomrel/parms.Fmax,'.')
        title('normalized CE force / (q*Fisom) vs normalized CE velocity, from simulation result')
        ylabel ('norm CE force [Fce/(q*Fisom)]')
        xlabel('CE velocity [lceopt/s]')
        % ! this would be a single curve in a Hill model !
        
        subplot(223)
        plot(gamma,q)
        title('q vs gamma, from simulation result')
        ylabel ('q []')
        xlabel('gamma []')
        
        subplot(224)
        plot(lce/parms.lceopt,fisomrel)
        title('normalized CE isometric force vs normalized CE length, from simulation result')
        ylabel ('norm isometric CE force []')
        xlabel('normalized CE length []')
        % ! this would be a single curve in a Hill model !
        
        %% isometric-frequency figs
        figure(50)
        iFreq=parms.iFreq;
        nFreqs=parms.nFreqs;
        subplot(4,nFreqs,iFreq); % stim/gamma
        plot(t,[stim gamma q])
        axis([0 t(end) 0 1])
        legend('STIM','calcium','active state')
        title(['Hatze act. dyn., freq=',num2str(parms.freq)]) %'stimulation as set in this protocol')
        xlabel ('Time [s]')
        ylabel ('normalized units')
        
        subplot(4,nFreqs,iFreq+nFreqs); % forces
        plot(tspan,parms.desiredForce,tspan,fse);
        axis([0 t(end) 0 .6])
        legend('desired','simulation')
        title(['Forces, freq=',num2str(parms.freq)]) %'stimulation as set in this protocol')
        xlabel ('Time [s]')
        ylabel ('Force [Fmax]')        
        
        subplot(4,nFreqs,iFreq+2*nFreqs); % power / work
        idx=tspan>tspan(end)-1/parms.freq; %only consider last (ss) cycle for evaluation
        p_avg_cb=trapz(t(idx),p_cb(idx))*parms.freq;
        p_avg_act=trapz(t(idx),stim(idx))*parms.freq;        
        plot(t,[p_cb p_act])
        axis([0 t(end) 0 1])
        title(['Average power: cb=',num2str(p_avg_cb,2),' act=',num2str(p_avg_act,2)])
        ylabel ('metabolic power [au]')
        xlabel('Time [s]')
        legend('P_c_b','P_a_c_t')
        
        subplot(4,nFreqs,iFreq+3*nFreqs); % power / work        
        plot(q(idx),fse(idx))        
        title('force vs q in ss')
        ylabel ('force [Fmax]')
        xlabel('active state [au]')        
    end
end
%figure(99)
%plot(tspan,parms.desiredForce,tspan,fse); drawnow

idx=tspan>tspan(end)-1/parms.freq; %only consider last (ss) cycle for evaluation
err = sqrt(mean((fse(idx)-parms.desiredForce(idx)).^2));
end
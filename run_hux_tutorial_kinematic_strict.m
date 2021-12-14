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

clear
close all
clc

% sort out pathing and stuff
tmp = mfilename('fullpath');
tmp = tmp(1:length(tmp)-length(mfilename));
cd(tmp)
addpath(genpath(cd)) % make sure all subfunctions are included

%% set parameters
tutorial_parms
%% Initial conditions

% -calculate (or set, in this case) lmtc0, fse0 and lse0
% -calculate lce0 and fisom0
% -calculate q0 (note that F_SS = fisom*q*Fmax, just as in Hill mod ...)
% -invert q(gamma) to find gamma0
% -set n0(x0) to SS distribution, corresponding to q0 and fisom0, thus:
% n0 = nSS*fisom0*q0

% set Fse0 :
parms.fse0=0.1*parms.Fmax; % [N]

% set lmtc0, activation follows from this
lmtc0=(1+se_strain)*parms.lse_slack+.95*parms.lceopt; % [m]
parms.lmtc0=lmtc0;
lmtcd0=0; % [m/s]

[lse0]=fzero(@Fse_inverse,parms.lse_slack*[1 (1+se_strain)],[],parms); % [m]

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
[fse0, fpe0, ~, ~] = CEEC_simple(lse0,lce0,parms);
parms.q0=(fse0-fpe0)/(parms.Fmax*fisomrel0);

% invert q(gamma) to find gamma0 (fzero has tolX = 1e-16 ...):
zeroFun=@(x)activeState_inverse(x,parms);
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

lce_L = 0.4*parms.lceopt; % [m] smallest value lce is expected to attain
lce_R = 1.6*parms.lceopt; % [m] largest value lce is expected to attain

x1 = round((lce0-lce_R)*parms.scale_factor/parms.lceopt); % [h] left bound x
x2 = round((lce0-lce_L)*parms.scale_factor/parms.lceopt); % [h] right bound x

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

%% test calc_lce function
global CURRENT_LCE
CURRENT_LCE = lce0;

[lcetmp]=calc_lce_strict(n0,lmtc0,parms);

%% run simulation
state0 = [n0(:); gamma0];

sim_opt='normal'; % stiff or normal

[stated0,y0,xRel0,nRel0,dndtRel0] = hux_tutorial_kinematic_strict(0,state0,parms);

tSpan=[0 t_end]; % chop up time to save memory (reduce state vector)
ode_fun=@(t,state)hux_tutorial_kinematic_strict(t,state,parms);
odeparms=odeset('abstol',1e-5,'reltol',1e-5,'maxstep',.02,'Stats','on');
%odeparms.JPattern=sparsity_pattern;
tic
switch sim_opt
    case 'normal'
        [t,state] = ode45(ode_fun,tSpan,state0,odeparms);
    case 'stiff'
        jacFunc=@(t,state)hux_tutorial_kinematic_strict_jac2(t,state,parms);
        odeparms.Jacobian=jacFunc;
        %odeparms.Jpattern=eye(length(state0)); % test with numerical case
        %odeparms.Jpattern(:,end)=1;
        [t,state] = ode23s(ode_fun,tSpan,state0,odeparms);
end
toc

% unravel state
%n = state(:,1:end-4);
gamma = state(:,end);

% initialise stated and optional output
stated=zeros(size(state'));
y = zeros(length(y0),length(t));
% output block
if diagnostics == true
    figure
end

for i=1:length(t)
    [stated(:,i),y(:,i),x,n,dndt] = ode_fun(t(i),state(i,:)');
    if diagnostics==true %&& mod(i-1,5)==0 % every 5 samples
        subplot(121)
        plot3(x,ones(size(x))*t(i),n);hold on
        xlabel x; ylabel t; zlabel n
        subplot(122)
        plot3(x,ones(size(x))*t(i),dndt./n);hold on
        xlabel x; ylabel t; zlabel('local eigenvalue')
    end
end
subplot(121)
xlim([-5 5])
ylim([0 t(end)])
zlim([0 1])
subplot(122)
xlim([-5 5])
ylim([0 t(end)])

%handle output
stated=stated'; y=y';
gammad = stated(:,end-2);

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
lce=y(:,9); %[W]
lmtc=y(:,10); %[W]
lmtcd=y(:,11); %[W]

lced=central_diff(lce,1/diff(t(1:2)));
lse=lmtc-lce; % [m]

fcerel = fce/parms.Fmax;
fperel = fpe/parms.Fmax;
fserel = fse/parms.Fmax;
%clear y
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
    figure
    plot (t,[ce_work+se_work+pe_work-lmtc_work])
    title('energy errors')
    xlabel ('time [s]')
    ylabel ('energy [J]')
    legend('w_c_e + w_p_e + w_s_e - w_m_t_c','Ekin - w_m_t_c - w_F_z')
end

%% standard output figs
if normFig==true
    figure
    subplot(221)
    plot(t,[stim gamma q]); hold on
    plot(t,gamma,'.')
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
    plot(t,[p_cb/trapz(t,p_cb) p_act/trapz(t,p_act)])
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
end

%% animation
if Animate
    % for video ...
    % writerObj = VideoWriter('hux_tutorial.avi');
    % writerObj.FrameRate = 100;
    % open(writerObj)
    
    figure
    for iSample=1:10:length(t)
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
        ylim([-.27 0.02])
        axis equal
        text(.01,-.5*lce(iSample),['CE (STIM=',num2str(stim(iSample),3),')'])
        text(.01,-.5*(lmtc(iSample)+lce(iSample)),'SEE')
        
        text(-.13,-.07,['t=',num2str(t(iSample),3)])
        %text(-.15,-.115,['STIM=',num2str(stim(iSample),3)])
        
        plot([0 0],[0 -lce(iSample)],'r','linewidth',2);
        plot([0 0],[-lce(iSample) -lmtc(iSample)],'b','linewidth',2);
        plot([-.015 .015],[-lmtc(iSample) -lmtc(iSample)],'k','linewidth',1.5); hold off
        drawnow
        
        % for video
        % currFrame=getframe(gcf);
        % writeVideo(writerObj,currFrame)
    end
    % close(writerObj)
end

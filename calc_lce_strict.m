function [lce, varargout]=calc_lce_strict(n,lmtc,parms)
% make function, something like: 
% [lce, fse, fpe, fce]=calc_lce_strict(n,lmtc,parms)
% we'll need:
% lce0, lceopt, kf, pe_shape, se_shape, lse_slack, lpe_slack, scale_factor, x0
% then compute a, b, c as in text
% somehow deal with slack cases?
% easy to calculate jacobian in case of kinematically driven model, only
% states n and gamma, so use of stiff integrator well possible ... 
do_fig=false;
% read out parameters
scale_factor = parms.scale_factor; % [] scale factor between lcerel and x
x0 = parms.x0; % initial vector x0 at t=0
lce0 = parms.lce0; % [m] ce length at t=0
lceopt = parms.lceopt; % [m] optimum ce length
se_shape= parms.se_shape; % [] relative strain @ Fmax
lse_slack = parms.lse_slack; % [m] slack length
pe_shape= parms.pe_shape; % [] relative strain @ Fmax
lpe_slack = parms.lpe_slack; % [m]
kf=parms.k_f; % [N/h] scaling between distribution and force 

% read out current lce ....
% this is a suuuper ugly solution but it is what it is at this point
global CURRENT_LCE

% compute integrals of x and n
int_n=trapz(n);
int_nx0=trapz(x0.*n);

%% compute possible lce 
% we are going to systematically evaluate all possible outcomes, which are:
% 1) both pe and se slack -> lce at Fce(lce)=0 satisfies this condition -> DONE 
% 2) only pe slack
% 3) neither pe nor se slack
% 4) only se slack, will never occur??

% initiate control variable
possible_lce=[]; % initiate

% condition 1; both pe and se slack
lce_tmp=(kf*scale_factor/lceopt*lce0*int_n - kf*int_nx0)/(kf*scale_factor/lceopt*int_n);
if lce_tmp<lpe_slack && (lmtc-lce_tmp)<lse_slack % this condition is feasible
    possible_lce=[possible_lce lce_tmp];
end

% condition 2; only pe is slack
% problem in terms of ABC formula
if isempty(possible_lce)
A=se_shape;
B=2*se_shape*(lse_slack-lmtc) - kf*scale_factor/lceopt*int_n;
C=se_shape*(lmtc.^2 + lse_slack.^2 -2*lmtc*lse_slack) + kf*(scale_factor*lce0/lceopt*int_n - int_nx0);
 
D=B.^2-4*A*C; % determinant
if D>0 % otherwise we don't even consider this   
    lce_tmp=(-B+[sqrt(D) -sqrt(D)])./(2*A); % two solutions to square equation
    lce_tmp=lce_tmp(lce_tmp<lpe_slack); % enforce PE condition
    lce_tmp=lce_tmp((lmtc-lce_tmp)>lse_slack); % enforce SE condition
    % update possible lce
    possible_lce=[possible_lce lce_tmp];
end

end
% condition 3; neither lpe nor lse slack
if isempty(possible_lce)
    A=se_shape-pe_shape;
    B=2*se_shape*(lse_slack-lmtc) - kf*scale_factor/lceopt*int_n + 2*pe_shape*lpe_slack;
    C=se_shape*(lmtc.^2 + lse_slack.^2 -2*lmtc*lse_slack) - pe_shape*lpe_slack.^2 + kf*(scale_factor*lce0/lceopt*int_n - int_nx0);
    D=B.^2-4*A*C; % determinant
    if D>0
        lce_tmp=(-B+[sqrt(D) -sqrt(D)])./(2*A); % two solutions to square equation
        lce_tmp=lce_tmp(lce_tmp>lpe_slack); % enforce PE condition
        lce_tmp=lce_tmp((lmtc-lce_tmp)>lse_slack); % enforce SE condition
        possible_lce=[possible_lce lce_tmp];
    end
end

%% select lce closest to current value
if ~isempty(possible_lce) % this value is filled, we are in business
    [~,iPick]=min(abs(CURRENT_LCE-possible_lce));
    lce=possible_lce(iPick);
else % check of the situation, something has gone wrong
 % visualize check of problem formulation
    lcetmp=(.5:.001:1.5)*parms.lceopt;
    Fse=se_shape*(lmtc-lcetmp-lse_slack).^2;
    Fse(lmtc-lcetmp<lse_slack)=0;
    Fpe=pe_shape*(lcetmp-lpe_slack).^2;
    Fpe(lcetmp<lpe_slack)=0;
    
    Fce=kf*int_nx0 - kf*scale_factor/lceopt*lce0*int_n + kf*scale_factor/lceopt*int_n*lcetmp;
    figure;
    subplot(311)
    plot(lcetmp,Fse); title('Fse vs lce')
    xlabel lce(m)
    ylabel Force(N)    
    subplot(312)
    plot(lcetmp,Fpe); title('Fpe vs lce')
    xlabel lce(m)
    ylabel Force(N)    
    subplot(313)
    plot(lcetmp,Fce); title('Fce vs lce')
    xlabel lce(m)
    ylabel Force(N)
    
    
    figure;
    subplot(211)
    plot(lcetmp,Fse-Fce); title('Fse-Fce-Fpe vs lce. Black 0 lpe_slack, red 0 lse_slack'); hold on
    plot(lcetmp,zeros(size(lcetmp)),'k--')
    plot([lpe_slack lpe_slack],[0 0],'k>')
    plot([lmtc-lse_slack lmtc-lse_slack],[0 0],'r<')
    plot(possible_lce,[0 0],'k+')
    plot(lce,0,'ro')
    xlabel lce(m)
    ylabel Force(N)
    
    subplot(212)
    plot(lcetmp,A*lcetmp.^2+B*lcetmp+C); title('ABC formula vs lce'); hold on
    plot(lcetmp,zeros(size(lcetmp)),'k--')
    plot([lpe_slack lpe_slack],[0 0],'k>')
    plot([lmtc-lse_slack lmtc-lse_slack],[0 0],'r<')
    plot(possible_lce,[0 0],'k+')
    plot(lce,0,'ro')
    xlabel lce(m)
    ylabel Force(N)    
    keyboard
end

% parse solution; feasibility and choose between lce1 and lce2
% pick solution that is closest to lceopt? this might fail at some point
% idea: the parabola really isn't a parabola, but half a parabola; look for
% the feasible solution. If lmtc > lce_slack + lse_slack, choose the point
% for which lse>lse_slack, otherwise choose the other point. We now found
% that it is possible for both solutions to be at lse>lse_slack and can
% then also be either slack or not. This means there is no principled way
% to distinguish such solutions. For now try with global variable (very
% ugly), other option is to do linear elastic elements instead of quadratic
% (will give only 1 solution), but that might also not be ideal... 

% Jettisen PE?? If lce<lpe_slack we don't want PE force to contribute to
% the solution; perhaps make two cases or so. 


%% update current LCE
CURRENT_LCE=lce;

% optional output; fce, fse and fpe
if nargout>1    
    % compute forces
    Fse=se_shape*(lmtc-lce-lse_slack).^2;
    Fse(lmtc-lce<lse_slack)=0;
    Fpe=pe_shape*(lce-lpe_slack).^2;
    Fpe(lce<lpe_slack)=0;    
    Fce=kf*int_nx0 - kf*scale_factor/lceopt*lce0*int_n + kf*scale_factor/lceopt*int_n*lce;
    varargout{1}=Fse;
    varargout{2}=Fce;
    varargout{3}=Fpe;
end
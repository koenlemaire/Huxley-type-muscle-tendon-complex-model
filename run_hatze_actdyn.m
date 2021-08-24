% hatze activation model according to the description in Rockenfeller and
% Gunther 2017
clear
close all
clc 

parms.v=3; % 3 according to Hatze/Rockenfeller
parms.lp=2.9; % 2.9 according to rockenfeller
% rho_c=[Ca]max * rho0 (ie molar volume); 1.373e-4 * 55700 (for v=3)
% according to Hatze. But rho0 ranges from 10e-6 to 16e-6 according to
% Rockenfeller (and some extent Kistemaker 2007). rho_c = 7 according to
% rockenfeller figure B4 but seems more likely to be 0.7 given the above
% and the discussion in Rockenfeller appendix C
parms.rho_c=1.373e-4 * 55700; % 7?? ?? 
parms.m=11.25; % [fast twitch, 1/3.67 for slow twitch]

gamma=linspace(0,1);
n=10;
lcerel_vec=linspace(.5,1.5,n);
figure
for i=1:n
    lcerel=lcerel_vec(i);
    q=activeState_hatze(gamma,lcerel,parms);
    plot(gamma,q,'k'); hold on
    title('q(gamma) for lcerel=linspace(0.5,1.5,10)')
    xlabel('relative free [Ca2+]')
    ylabel('active state, ie potential isometric force')
end


figure;
for i=1:5
    gamma0=0;
    parms.stim=1/i;  
    odeopt=odeset('abstol',1e-9,'reltol',1e-9);
    [t,gamma]=ode45(@simulate_hatze_ode,[0 .6],gamma0,odeopt,parms);
    lcerel=1;
    q=activeState_hatze(gamma,lcerel,parms);
    subplot(211);plot(t,gamma,'k'); hold on
    title('gamma(t) for stim=1/(1:5)')
    xlabel('time [s]')
    ylabel('relative free [Ca2+]')
    subplot(212);plot(t,q,'k'); hold on
    title('q(t) for stim=1/(1:5)')    
    xlabel('time [s]')
    ylabel('active state, ie potential isometric force')
   
end
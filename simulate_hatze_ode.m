function [ gammad ] = simulate_hatze_ode( t,gamma,parms )

stim=0;
if t<.4
    stim=parms.stim;
end
    
m=parms.m; % 1/time constant of hatze activation dynamics
gammad=m*(stim-gamma);
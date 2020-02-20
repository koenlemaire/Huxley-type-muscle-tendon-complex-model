function [ gammad ] = calciumdyn_hatze( gamma,stim,parms )

m=parms.m; % 1/time constant of hatze activation dynamics
gammad=m*(stim-gamma);
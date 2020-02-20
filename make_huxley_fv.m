clear all
%close all
%clc

% CONCLUSIONS BASED ON THIS CODE (03-2020 KKL):
% -effect of scaling on ss fv curve: 
% if vce=c*dx/dt and f1*=f1/c, g1*=g1/c, etc. the ss fv relationship is
% invariant for all c, for a given rate function and rate parameters. This
% is because the ode's are all linear (I think). However the speed at which
% the ss is reached in a dynamical contraction is of course dependent on c;
% the larger c, the quicker you reach equilibrium. In terms of
% computational cost, the cost for a dynamic contraction is lower for lower
% values of c, right untill the time scale of cb attachment is in the same
% ballpark as CE SEE interaction. 
% -effect of q and fisomrel:
% the effects of q and fisomrel on ss fv curve are the same; they scale the
% ss fv curve in the vertical direction. 

h = 1e-8;           % attachment 'range' for myosin head [m]
s = 2.6e-6;         % sarcomere length [m]
scaletmp=1; % drastically reduces computation time while hardly affecting results, very worthwhile implementing it seems ...
parms.scale_factor = s/(2*h)/scaletmp; % [] scaling between x and lcerel
parms.dx=.01; % [h] stepsize in x
tmp=1.0e+03 * [0.8890    0.4275    3.1703    0.7796];
parms.f1=tmp(1)/scaletmp; % [Hz] detachment rate parameter
parms.g1=tmp(2)/scaletmp; % [Hz] attachment rate parameter
parms.g2=tmp(3)/scaletmp; % [Hz] detachment rate parameter
parms.g3=tmp(4)/scaletmp; % [Hz] detachment rate parameter
% parms.g1=500/scaletmp; % [Hz] detachment rate parameter
% parms.f1=800/scaletmp; % [Hz] attachment rate parameter
% parms.g2=3000/scaletmp; % [Hz] detachment rate parameter
% parms.g3=1400/scaletmp; % [Hz] detachment rate parameter
parms.q=1;
parms.fisomrel=1;

rateFun=@(x)rateFunc_v8(x,parms);
parms.rateFun=rateFun;

% % domain for x of the initial curve
x1 =  -5; %[bond length au]
x2 =  5; %[bond length au]

% stepsize for initial x, default is 1000 steps
dx = 1/1000;

% create x for library
x = x1:dx:x2;
[fx,gx]=rateFun(x(:));

figure;plot(x,fx,x,gx)
title('fx and gx')
n_ss = fx./(fx+gx);  % inital state vector for n0, isometric ss condition is given by f(f+g)|f(x)>0!!
figure;plot(x,n_ss)
title('steady state nx for v=0')

hux_vce=-10:.1:10;
Fhux=zeros(size(hux_vce));
tic
for i = 1: length(hux_vce)
    hux_vce(i)
    [Fhux(i)]=Fv_huxley_simple(hux_vce(i),parms);
end
toc
figure;
plot(hux_vce,Fhux)


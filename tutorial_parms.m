%% set parameters
% simulation processing:
diagnostics=false; % figures and data for sanity checking results
normFig=true; % figures for viewing data
Animate=false; % do simple animation

% general:
parms.mass=1; % [kg]
parms.Fmax=2*parms.mass*9.81; %[N]
parms.lceopt=0.07; % [m] CE optimum length
parms.lpe_slack=1.1*parms.lceopt; % [m] PE slack length
parms.lse_slack=0.13; % [m] SE slack length
se_strain=.05; % [N/m^2] SE shape, Fse=Fmax at 4% strain
pe_strain=.2; % [N/m^2] PE shape, Fpe=0.5*Fmax at 4% strain????

% compute kse and kpe
parms.se_shape=parms.Fmax./((se_strain.*parms.lse_slack).^2);
parms.pe_shape=parms.Fmax./((pe_strain.*parms.lpe_slack).^2);

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
scaletmp=1; % drastically reduces computation time while hardly affecting results, very worthwhile implementing it seems ...
parms.scale_factor = s/(2*h)/scaletmp; % [] scaling between x and lcerel
parms.dx=.01; % [h] stepsize in x
parms.g1=200/scaletmp; % [Hz] detachment rate parameter
parms.f1=800/scaletmp; % [Hz] attachment rate parameter
parms.g2=3000/scaletmp; % [Hz] detachment rate parameter
parms.g3=1400/scaletmp; % [Hz] detachment rate parameter

% energetics parms scaling: need to fit/optimize these
parms.c_cb=1; % scale factor between cross bridge uncoupling and metabolic power
parms.c_act=1; % scale factor between free calcium and metabolic power 

% rate function form; see different supplied functions
parms.rateFun=@(x)rateFunc_v8(x,parms);

% simulation time
t_end=.15; % [s] simulation time, starting at t=0 ...


% input: parameters defining hill force velocity relation at q=1 and
% lce=lceopt, parameter struct parms, rateFunc (function handle to Huxley
% ratefunction)
% output: huxleyParms = [f1 g1 g2 g3]
% make sure parms.scale_factor exists for huxley model and parameters
% for isometric force length relation exist for Hill model.

% !! CRITICAL NOTE:
% the success of the optimization depends wholy on the chosen bounds (lb
% and ub) or region where the optmizer is allowed to look. If you make
% changes in the Hill F-v curve, it is likely you have to make some changes
% in the bounds to accomodate the new parameters; it is generally desired
% to choose these bounds as tightly as possible around the optimum because
% then the algorithm can and will converge quicker and more precisely. This
% is mostly a matter of trial and error, ie if you notice the optimal
% solution does not look very optimal and one the parameters runs into its
% bound you should increase the bound. Conversely if the algorithm is slow
% to converge or doesn't converge well, try tightening the bounds. In
% general, quite good correspondance (cost ~ 0.01) should be possible for
% a range of Hill F-v curves, note that this also depends on the rate
% functions and on the hill F-v curve. In general the huxley curve will not
% deal well large (>1.5) slopefac values, note though that the
% high frequency behaviour (stiffness) in the Huxley model is different
% from the F-V relation calculated here ... 
%% parameters and initial conditions
clear; clc; 
parms.Fmax=1;

% huxley model:
h = 1e-8;           % attachment 'range' for myosin head [m]
s = 2.6e-6;         % sarcomere length [m]
parms.scale_factor = s/(2*h); % [] scaling between xrel and lcerel!!

parms.Fmax=1;
parms.width=.56;
parms.C=-1/parms.width^4;
parms.C=-parms.width.^(1/4);
parms.qmin=1e-10;
Arel=0.41;
Brel=5.2;
Fasymp=1.5;
slopfac=1.2;

parms.rateFunName='rateFunc_v8';

% for the above settings parameters [f1 g1 g2 g3] are:
% [11    837    1236    368]

parms.Fasymp=Fasymp;
parms.slopfac=slopfac;
parms.Arel=Arel;
parms.Brel=Brel;

fce=linspace(0,0.95*Fasymp,2000)';
%fce=[linspace(.95,.995,10) linspace(1.005,1.05,10)]';
lcerel=1;
fisomrel=1;
q=1;
parms.fisomrel=1;
parms.q=1;
vcerel=zeros(size(fce));
for i=1:length(fce)
    [vcerel(i,1),~] = lcereldotnew_v2(lcerel,fce(i),fisomrel,q,parms);
end

parms.hillData = [vcerel fce];
parms.hux_vce=linspace(vcerel(2),vcerel(end-1),20)'; % 20 querie points
parms.hux_vce=parms.hux_vce(parms.hux_vce~=0); % make sure no zeros in here ...
lb=[5e2 1e2 1e3 1e2]; % 
ub=[1.5e3 7e2 5e3 1.5e3]; % f1 g1 g2 g3 

x_init=(lb+ub)/2; % just in the middle of our space ...
parms.x_init=x_init;
parms.spread_x=ub-lb; % complete space ...

%% Optimization
% Here I found that simulated annealing works better than
% fmincon/fminsearch, who tend to get trapped in local minima a lot. Still,
% also with sa you'll have to tweak the upper and lower bounds to guarantee
% convergence; if your Hill parms are far away from what the current
% settings are it may not converge successfully.  
sa_opt=optimoptions('simulannealbnd');
n_hours=.5;
sa_opt.MaxTime=n_hours*3600; % [s]
sa_opt.Display='diagnose';
sa_opt.InitialTemperature=sqrt(length(x_init)); % if space is normalized such that axes are within [-.5,.5], the maximum step length is sqrt(length(x_init)) ...
sa_opt.ReannealInterval=250;
sa_opt.AnnealingFcn='annealingboltz';
sa_opt.TemperatureFcn=@temperatureexp_kkl; % T=T_init*0.97^k
sa_opt.MaxFunEvals=inf;
sa_opt.MaxStallIterations=1000*length(x_init);

parms.dispFig=false;

fcost=@(x)calculateHuxleyError(x,parms); % DOESNT WORK WITH LSQNONLIN!!!

% scale design parameter input:
lbrel=scale_x(lb,'abs_to_rel',parms);
ubrel=scale_x(ub,'abs_to_rel',parms);

for iOpt=1%:10
    x0rel(iOpt,:)=rand(size(x_init))-.5+scale_x(x_init,'abs_to_rel',parms);
    %[xrel(iOpt,:),fval(iOpt),flag] = fminsearch(fcost,x0rel(iOpt,:),fminsearch_opt);
    %[xrel(iOpt,:),fval(iOpt),flag,out] = fmincon(fcost,x0rel,[],[],[],[],lbrel,ubrel,[],fmincon_opt);
end

%[xrel,fval,flag,out] = simulannealbnd(fcost,x0rel,lbrel,ubrel,sa_opt);
    
%[~,iMin]=min(fval);
%huxleyParms=scale_x(xrel(iMin,:),'rel_to_abs',parms);
%huxleyParms=scale_x(xrel,'rel_to_abs',parms);

parms.dispFig = true;
tmp=1.0e+03 * [0.8890    0.4275    3.1703    0.7796];
tmp=scale_x(1.0e+03 * [0.8890    0.4275    3.1703    0.7796],'abs_to_rel',parms);
[cost]=calculateHuxleyError(tmp,parms);
keyboard    % check stuff when finished ...


function [lcereld,regime] = lcereldotnew_v2(lcerel,fce,fisomrel,q,parms)
%function [lcereld,regime] = lcereldotnew(lcerel,fce,fisomrel,q,muspar)
%
% inputs relative ce length, absolute ce force, relative isometric force, active state and parms
% output time derivative of lcerel and regime, a test variable that is normally not used
% this function models the force velocity relation, given all other variables
% KvS 01/2003

% wijziging 1: fact vervangen door parabool met fact=0.1 bij q=qmin en fact=1 bij q=0.3
%              en dfact/dq=0 bij q=0.3 zodat continue 1e afg aldaar

% wijziging 2: concentrische fv relatie vereenvoudigd

% wijziging 3: nieuwe eccentrische fv relatie gemaakt, hyperbool met schuine asymptoot

% wijziging 4: sloplin een functie gemaakt van overige parameters, en wel zodanig dat
%              bij q=qmin en lcerel>1+0.95*width er geen eccentrische hyperbool meer kan zijn

% wijziging 5: in uitzonderingsgeval van 4 de ecc helling gelijk gemaakt aan concentrische
%              vermenigvuldigd met slopfac (ipv gelijk aan sloplin)

% wijziging 6: vfact in brel gestopt (dwz brel=brel.*vfact)

% 30/10/2002: lengte die sloplin definieert gezet op 1+0.8*width ipv
% 1+0.95*width
% en verder: defective cases anders:
% als helling in isom punt verkeerd dan lineaire F-v relatie zowel
% concentrisch als eccentrisch, met hellingverhouding slopfact
%
% this model is now used by a program also running a huxley model; folowing
% changes have been made:
% 'protocol' is used to elicit various contraction protocols
% Koen Lemaire 4/2011
%
% v2: BUG DISCOVERED, KKL 4/2013
% during regime 1: fcerel<0 can result in a positive lcereld, which can
% lead to instability if fcerel+q>0, is now mitigated by demanding fcerel>0
% in regime 1 and adding a (concentric) regime for fcerel<0.

% lcerel=lcerel(:)';
% fce=fce(:)';
% fisomrel=fisomrel(:)';
% q=q(:)';

% Bottom parameter read script slightly adapted for my own purposes, KKL
% 4/2011

% read out muscle parameters
lcerel=lcerel(:);
fce=fce(:);
fisomrel=fisomrel(:);
q=q(:);

fasymp    = parms.Fasymp;
fmax      = parms.Fmax(:);
arel      = parms.Arel(:);
brel      = parms.Brel(:);
%lceopt    = parms.lceopt;
%sloplin   = parms.hill.slope_linear;
slopfac   = parms.slopfac;
width     = parms.width(:);
hatze_q0  = parms.qmin(:);


qcrit=0.3;vfactmin=0.1; % see below for explanation; these might be made input parameters??

% new sloplin depends on minimal value of brel which equals vfactmin*brel (see also below)
% sloplin is calculated so that it equals the eccentric isometric dvcerel/dfcerel
% at q=qmin and lcerel=1+0.95*width
% Note: sloplin may be removed as an input parameter??

sloplin=vfactmin.*brel./(slopfac.*hatze_q0.*ce_fl_simple( 1+0.95*width,parms ).*(1+arel));

% calculate fcerel
fcerel = fce./fmax;	    					%

% if lcerel>1, arel scales with fisomrel so that lcereldmax is unaffected by length
i       = find(lcerel>1);
arel(i) = arel(i).*fisomrel(i);

% if q<qcrit, brel depends on q so that at low q lcereldmax is affected by q
% NOTE THIS REPLACES VFACT IN EARLIER VERSIONS
brel=brel.*(1-(q<qcrit).*(1-vfactmin).*((q-qcrit)./(hatze_q0-qcrit)).^2);

dvdf_isom_con=brel./(q.*(fisomrel+arel)); % desired eccentric isometric slope
dvdf_isom_ecc=dvdf_isom_con./slopfac;

% concentric part; simplified form 8/2002
c       = find( ((fcerel./q) <= fisomrel) & (dvdf_isom_con<=sloplin) & (fcerel>=0));
lcereld(c) = brel(c).*(fcerel(c)-q(c).*fisomrel(c))...
    ./(fcerel(c)+q(c).*arel(c));
regime(c)=1;
% partial derivative of concentric lcereld wrt fcerel; for checking purposes
%dvreldfrel(c)=brel(c).*q(c).*(arel(c)+fisomrel(c))./((fcerel(c)+q(c).*arel(c)).^2);

%NEW ECCENTRIC PART
% the hyperbola:
% (fcerel-fasymp*fisom*q-vcerel/sloplin)*(vcerel+p1)=p2
e     = find( ((fcerel./q) > fisomrel) & (dvdf_isom_ecc<=(sloplin./slopfac)) );
k1=q.*fisomrel.*(1-fasymp); % fasymp>1!!
k2=slopfac.*q.*(fisomrel+arel)./brel;
p1=k1.*sloplin./(1-k2.*sloplin); % kan hier ooit iets misgaan???
p2=k1.*p1;
r1=sloplin.*(fcerel-fasymp.*fisomrel.*q); % kan hier ooit iets misgaan???

lcereld(e)=(r1(e)-p1(e)+sqrt((p1(e)-r1(e)).^2-4*(p2(e).*sloplin(e)-r1(e).*p1(e))))/2; % altijd de hoogste van de twee wortels!!
regime(e)=2;
% partial derivative of lcereld wrt fcerel; for checking purposes
%dvreldfrel(e)=1./((1./sloplin(e))-p2(e)./((lcereld(e)+p1(e)).^2));

% now for the concentric defective case
% ie dvdf in isometric point higher than allowed
c2     = find( ((fcerel./q) <= fisomrel) & (dvdf_isom_con>sloplin) & fcerel>=0);
lcereld(c2)=sloplin(c2).*(fcerel(c2)-q(c2).*fisomrel(c2));
regime(c2)=3;
% and now for the eccentric defective case
e2     = find( ((fcerel./q) > fisomrel) & (dvdf_isom_ecc>(sloplin./slopfac)) & fcerel>=0 );

lcereld(e2)=(sloplin(e2)./slopfac(e2)).*(fcerel(e2)-q(e2).*fisomrel(e2));
regime(e2)=4;

% See note in file history 4/2013, new concentric part when fcerel<0:
% lcereld = -brel/arel. Note this is not an ideal solution! Idea for now is
% to compute
c3 = find(fcerel<0);
lcereld(c3) = -brel(c3)./arel(c3);
regime(c3)=5;

% if ~isreal(lcereld)==1
%     save('errorLog')
%     error('lcereld ~isreal')
% end
% partial derivative of lcereld wrt fcerel; for checking purposes
%dvreldfrel(e2)=dvdf_isom_ecc(e2);

end
function [ error ] = calculateHuxleyError( huxleyParms,parms )
%[ ERROR ] = calculateHuxleyError( huxleyParms0,dataPoints )
%   calculateHuxleyError computes the residual of points predicted on a
%   huxley force velocity curve given by huxleyParms and dataPoints. This
%   function can be used to fit huxleyParms to datapoints by means of
%   optimization
%
% INPUT
%   huxleyParms: [f1 (g1) g2 g3] Huxley's parameters g1 is optional
%   parms:
%
% OUTPUT
%   error:      scalar value
%
% Koen Lemaire 01/2013

%read out parameters
huxleyParms=abs(scale_x(huxleyParms,'rel_to_abs',parms));
parms.f1 = huxleyParms(1);
parms.g1 = huxleyParms(2);
parms.g2 = huxleyParms(3);
parms.g3 = huxleyParms(4);
hillData = parms.hillData;
hux_vce=parms.hux_vce;
eval(['rateFun=@(x)',parms.rateFunName,'(x,parms);'])
parms.rateFun=rateFun;

% domain for x of the initial curve
x1 =  -2; %[bond length au]
x2 =  3; %[bond length au]

% stepsize for initial x, default is 1000 steps
dx = 1/1000;

% create x for library
x = x1:dx:x2;
[fx,gx]=rateFun(x(:));



Fhux=zeros(size(hux_vce));
for i = 1: length(hux_vce)
    [Fhux(i)]=Fv_huxley_simple(hux_vce(i),parms);
end

huxData=[hux_vce Fhux];

[idx,D]=knnsearch(hillData,huxData); % gives distance to nearest neighbours..

error = sum(D.^2);

if parms.dispFig==true
    figure(98)
    plot(x,[fx gx])
    figure(99)
    plot(hillData(:,1),hillData(:,2),'k',huxData(:,1),huxData(:,2),'ko',hillData(idx,1),hillData(idx,2),'ro')
    drawnow
    keyboard
end

%error=norm(huxleyParms);
end



function [ x_out ] = scale_x( x_in,direction,parms )
% this function deals with scaling our design vector. direction specifies
% the change. direction<0: from x_abs to x_rel. Direction>0: from x_rel to
% x_abs.
%[m,~]=size(x_in);

% x=[f1 g1 g2 g3 tau_act tau_deact]
mean_x=parms.x_init;
spread_x=parms.spread_x;
%spread_x=[1500 500 2000 2000 .02 .04];

%mean_x=repmat(mean_x,m,1);
%spread_x=repmat(spread_x,m,1);

switch direction
    case 'abs_to_rel' % we go from absolute to relative x
        x_out=(x_in-mean_x)./spread_x;
    case 'rel_to_abs' % we go from relative to absolute
        x_out=x_in.*spread_x + mean_x;
end
end

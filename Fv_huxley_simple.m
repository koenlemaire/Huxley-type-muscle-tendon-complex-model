function [F,varargout] = Fv_huxley_simple (vce,parms)
% program for calculating F as a function of vcerel, using the huxley model

% extract parameters:

scale_factor=parms.scale_factor; % scaling between vce and u
rateFun=parms.rateFun; % rate function to be used here ...
% define u, shortening velocity in terms of xb bond lengths
parms.u = vce*scale_factor; % [h/s]

% domain for x of the initial curve
x1 =  -5; %bond length [h]
x2 =  5; %bond length [h]

x_ss=linspace(x1,x2,5000)';
[fx0,gx0]=rateFun(x_ss);

n_ss = fx0./(fx0+gx0);  % inital state vector for n0, isometric ss condition is given by f(f+g)|f(x)>0!!
k_f=trapz(x_ss,x_ss.*n_ss); % [m] % note trapz obligatory!! because dx~=constant

i=find(n_ss); % find indices where n>0

if vce < 0
    state0 = [x_ss(i(end)+1) n_ss(i(end)+1)]; %[x0 n0]
elseif vce == 0
    F=parms.q*parms.fisomrel;
    if nargout>1
        varargout{1}=x_ss;
        varargout{2}=n_ss;
        if nargout>3
            varargout{3}=trapz(x_ss,gx0.*n_ss); % metabolic energy [au/s]
        end        
    end
    return
else
    state0 = [x_ss(i(1)-1) n_ss(i(1)-1)]; %[x0 n0]
end

% simulate the solution
tspan = [0 inf]; % inf: we should wait for termination criterium
%odeparms=odeset('abstol',1e-6,'reltol',1e-6,'maxstep',.5,'events',@huxleyEvent);
odeparms=odeset('abstol',1e-6,'reltol',1e-6,'maxstep',.1,'events',@huxleyEvent);
[~,state] = ode113(@huxley_SS,tspan,state0,odeparms,parms);
if vce < 0
    x=flipud(state(:,1));
    n=flipud(state(:,2));
else
    x=state(:,1);
    n=state(:,2);
end
F=trapz(x,x.*n)/k_f; %[Fisomrel] % trapz is obligatory!!

if nargout>1
    varargout{1}=x;
    varargout{2}=n;
    if nargout>3
        [~,gx]=rateFun(x(:));
        varargout{3}=trapz(x,gx.*n); % metabolic energy [au/s]
    end
end


% figure(8)
% subplot(211)
% plot(xSS,nSS,xSS,xSS.*nSS)
% title('SS n over x at current v')
% xlabel('xSS')
% ylabel('nSS')
% subplot(212)
% plot(xSS,f,xSS,g)
% title('f & g')
% xlabel('x')
% ylabel('rate')
% figure(99)
% plot(x,n)
% title(['n over x at v=',num2str(vce),' F=',num2str(F)])
% \xlabel('x')
% ylabel('n')
% xlim([-5 5])
% ylim([0 1])
% drawnow

end
function [ dstatedt ] = huxley_SS( t,state,parms )
%function [ dstatedt ] = huxley( t,state )
%   input: state (x n)
%   output: time derivitave of state
% Huxleymodel as described in Zahalac (1981) eq. 3:
% (dndt)x - v*(dndx)t = f(x) - [f(x) + g(x)]*n                      (1)
% parametrisation by methods of characteristics results in
% a(x,t) = v(t);
% b(x,t) = 1;
% c(x,t) = -[f(x) + g(x)];
% d(x,t) = f(x).
% the system is then transformed into 3 ODE's:
% dxds = a(x,t) = v(t), x(0) = x0;
% dtds = b(x,t) = 1, t(0) = t0;
% dnds = c(x,t)*n + d(x,t) = -[f(x) + g(x)]*n + f(x), n(0) = u0.
% the second one is dumped because t = s. x0 is the initial domain of n, n0
% is the inital state of n

u = parms.u; % [lceopt/s * scale_factor] current velocity in xb bond length!!
q=parms.q; % active state
fisomrel=parms.fisomrel; % norm
rateFun=parms.rateFun;

x = state(1);
n = state(2);
% calculate f(x) and g(x) and determine dndt. NOTE: ratefunc_v4 requires
% previously calculated input parms...
[fx,gx]=rateFun(x(:));

dndt=q*fisomrel*fx-(fx+gx).*n;
dxdt=u;
dstatedt = [dxdt;dndt];
end
function [ VALUE,ISTERMINAL,direction ] = huxleyEvent( t,state,parms )
%[ VALUE,ISTERMINAL,DIRECTION ] = huxleyEvent( t,state,parms )
%   this event function calculates at which point n is zero while x<0, at
%   this point we assume that the integration can be terminated, since the
%   complete curve has been granted
ISTERMINAL = 1;
direction = 0;
x=state(1);
n=state(2);

if parms.u<0    
    if x<0
        VALUE=n;
    else
        VALUE=1;
    end
else
    if x>1
        VALUE=n;
    else
        VALUE=1;
    end    
end

end

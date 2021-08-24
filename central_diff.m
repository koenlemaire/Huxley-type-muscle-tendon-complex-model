function [ dxdt ] = central_diff( x,fs)
%[ dxdt ] = diff1d( x,dt)
%   differentiates N-Dimensional signal using central difference. Only
%   signals with equidistant time sampling are suited for this function.
%   Note: the first and last sample are unreliable estimates of the
%   derivative. function such that size(dxdt)=size(x) for all x
% input: x signal, fs sample frequency.
% output: dxdt 
%
% example:
% fs=1000;
% t=0:(1/fs):1;
% x=[sin(2*pi*t); t.^2+t-1];
% dxdt=diff1d(x,fs);
% figure;subplot(211); plot(t,x)
% subplot(212);plot(t,dxdt)

% koen lemaire 02/2021

% give same output dimensions as input dimensions, assuming largest
% dimension is nr of samples:

[d1,d2]=size(x);

if d1>d2 % nr of samples in first dimension
    x=x';
else % do nothing, algorithm already accomodates this
end
% central difference algorithm:
dxdt=(x(:,3:end)-x(:,1:end-2))*.5*fs;% derivative
dxdt=[dxdt(:,1) dxdt dxdt(:,end)]; % correct dimension
if d1>d2
    dxdt=dxdt'; % transpose back
end
end


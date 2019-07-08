% script file for testing various implementations of the ratefunc function
% clear;
% close all
% clc
% first the classical approach:

opt=1.0e+04 * [0.031274633000000   2.662285600000000   0.067457802000000   0.000003625249400   0.000008013419200   0.000493615020000];

parms.g2=opt(2);
parms.f1=opt(1);
parms.g3=opt(3);
parms.g1=100;
h=1;
%h=1e-8;
parms.h=h;
x=linspace(-3*h,3*h,600);

checkFunc=@(x)rateFunc_v5(x,parms);

[fx,gx]=checkFunc(x(:));
figure;plot(x,[fx gx], x(1:end-1),[diff(fx)./diff(x) diff(gx)./diff(x)])
figure;plot(x,[fx gx])
%legend('fx','dfxdx','gx','dgxdx')
 
% tstart=cputime;
% for i=1:1000
%     [ fx,gx ] = rateFunc_v5( x,parms );
% end
% tspend1=cputime-tstart
% now the new appraoch: first make libraries large enough for
% interpolation:
% dx=.001;
% parms.ndxLib=1/dx;
% xLib=-16:dx:16;
% 
% [ fxLib,gxLib ] = rateFunc_v5( xLib,parms );
% figure;plot(xLib,fxLib,xLib,gxLib)
% 
% parms.fxLib=fxLib(:);
% parms.gxLib=gxLib(:);
% parms.xLib=xLib(:);

% tstart=cputime;
% for i=1:1000
%     [ fx,gx ] = rateFunc_v6( x,parms );
% end
% tspend12=cputime-tstart

% 



% % rateFunc figure!
% load par
% x=par(3).hux.xLib;
% fx=par(3).hux.fxLib;
% gx=par(3).hux.gxLib;
% figure;plot(x,fx,'k','LineWidth',2); hold on
% plot(x,gx,'k--','LineWidth',2)
%         set(gca, 'XLim', [-1, 2]);        % Sets start and end of displayed part of axis
%         set(gca, 'XTick', [-1, 0,1, 2]); % Sets the numbers at which small indicator lines and numbers appear along axis
%         set(gca, 'YLim', [0, 900]);
%         set(gca, 'YTick', [0 200 400 600 800]);
%         set(gca, 'FontSize',12, 'FontWeight','demi');
%         
%         xlabel('Bond length [h]')
%         ylabel('Rate [$s^{-1}$]')
%         legend('Attachment','Detachment')
%         legend('Location','northeast')    % Or 'NorthWest', 'SouthEast' etc.
%         box off 
%         legend boxoff
%         
%         figureHandle = gca;
%         set(findall(figureHandle,'type','text'),'fontSize',12, 'FontWeight','demi')    % Sets rest of text also at fontsize 14
%         set(figureHandle, 'Units', 'centimeters', 'OuterPosition', [0 0 16 7], 'LineWidth', 1)       % Sets the dimensions of the plot (i
% % store result in parms:
% break
% parms.xLib=xLib(:);
% parms.gxLib=gxLib(:);
% parms.fxLib=fxLib(:);
% parms.rateLib=[fxLib(:) gxLib(:)];
% 
% tstart=cputime;
% for i=1:1000
%     [ fx2,gx2 ] = rateFunc_v4( x,parms );
% end
% tspend=cputime-tstart
% 
% % now run new approach:
% x=x(:);
% tstart=cputime;
% for i=1:1000
%     [ fx3,gx3 ] = rateFunc_v3( x,parms );
% end
% tspend=cputime-tstart
% 
% 
% 
% figure;plot(x,fx,x,fx2);title('fx and fx2')
% figure;plot(x,gx,x,gx2);title('gx and gx2')
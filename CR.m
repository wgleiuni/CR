F=1e-2;
proc=1; % 1: semi-Dirac; 2: massive; 3: massless
px0=0;
py0=1e-4;
f=@(x,y)CR_diff(x,y,F,proc);
[T,Y]=ode45(f,[0 5000*2*pi],[px0 py0]);
figure
set(gcf,'color','w','position',[2000 280 560 630])
subplot(3,1,1)
plot(T/(2*pi),Y)
if proc==1
    syst='semi-Dirac';
elseif proc==2
    syst='massive';
elseif proc==3
    syst='massless';
end
title([syst,' F=',num2str(F),' px0=',num2str(px0),' py0=',num2str(py0)])
subplot(3,1,2)
plot(T/(2*pi),sqrt(Y(:,1).^2+Y(:,2).^2))
title('total energy')
subplot(3,1,3)
plot(T/(2*pi),atan2(Y(:,2),Y(:,1))/pi)
title('phase')
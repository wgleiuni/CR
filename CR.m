F=1e-0;
proc=3; % 1: semi-Dirac; 2: massive; 3: massless
px0=-1;
py0=0;
f=@(x,y)CR_diff(x,y,F,proc);

options=['RelTol',1e-16,'AbsTol',[1e-16 1e-16]];
[T,Y]=ode113(f,[0 300*2*pi],[px0 py0],options);
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
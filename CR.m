F=1e-2;
proc=1; % 1: semi-Dirac; 2: massive; 3: massless; 4: semi-Dirac with delta; 5: semi-Dirac in p-phi space
px0=-1;
py0=0;
delta=1e3*0; % delta=Delta/(m*v^2)
f=@(x,y)CR_diff(x,y,F,proc,delta);

options=['RelTol',1e-16,'AbsTol',[1e-16 1e-16]];
[T,Y]=ode113(f,[0 300*2*pi],[px0 py0],options);
figure
set(gcf,'color','w','position',[2000 280 560 630])
if proc==1
    syst='semi-Dirac';
elseif proc==2
    syst='massive';
elseif proc==3
    syst='massless';
elseif proc==4
    syst='D semi-Dirac';
elseif proc==5
    syst='semi-Dirac in phase space';
end
if proc~=5
    subplot(3,1,1)
    plot(T/(2*pi),Y)
    title([syst,' F=',num2str(F),' \Delta=',num2str(delta),' px0=',num2str(px0),' py0=',num2str(py0)])
    subplot(3,1,2)
    plot(T/(2*pi),sqrt(Y(:,1).^2+Y(:,2).^2))
    title('total energy')
    subplot(3,1,3)
    plot(T/(2*pi),atan2(Y(:,2),Y(:,1))/pi)
    title('phase')
elseif proc==5
    subplot(3,1,1)
    plot(T/(2*pi),Y(:,1))
    title([syst,' F=',num2str(F),' \Delta=',num2str(delta),' px0=',num2str(px0),' py0=',num2str(py0)])
    subplot(3,1,2)
    plot(T/(2*pi),Y(:,2))
    title('phase')
    subplot(3,1,3)
    plot(Y(end-1000:end,1),Y(end-1000:end,2))
    title('phase')
end
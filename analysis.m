step=1000;
t=zeros(step,50000);
Kx=t;
Ky=t;
E=t;
for i=1:1:step
    eval(['load data/data4/D',num2str(i),'.mat']);
    temp=eval(['D',num2str(i)]);
    if ~isempty(temp)
        t(i,:)=temp(:,1);
        Kx(i,:)=temp(:,2);
        Ky(i,:)=temp(:,3);
        E(i,:)=temp(:,4);
        eval(['clear D',num2str(i)])
    end
end

%%
step=1000;
lya=zeros(step,2);
for i=1:1:step
    eval(['load data/data6/L',num2str(i),'.mat']);
    ltemp=eval(['L',num2str(i)]);
    if ~isempty(ltemp)
        lya(i,:)=ltemp;
        eval(['clear L',num2str(i)])
    end
end
%%
figure
set(gcf,'color','w')
tI=-4+9/1000*(1:1000);
para=10.^tI;
plot(para,lya)
%%
cut=2000;

%%
figure
set(gcf,'color','w')
tI=-6+7/10000*(1:10000);
para=10.^tI;
para=2/1000*(1:1000);
hold on
for i=1:1:step
    y=findpeaks(Kx(i,end-cut:end));
    x=ones(length(y),1)*para(i);
    scatter(x,y,2,'b')
end

for i=1:1:step
    y=findpeaks(Ky(i,end-cut:end));
    x=ones(length(y),1)*para(i);
    scatter(x,y,2,'r')
end

%%
total_Kx=mean(Kx,1);
total_Ky=mean(Ky,1);


%%
D=D1;
figure
set(gcf,'color','w','position',[2000 280 560 420])
subplot(2,1,1)
plot(D(:,1),D(:,2),'b')
hold on
plot(D(:,1),D(:,3),'r')
subplot(2,1,2)
plot(D(:,1),D(:,4))

%% FFT
ftemp=Ky(900,end-cut:end);
dt=0.001*100;
Fs=1/dt;
L=length(ftemp);
t=(0:L-1)*dt;
NFFT=2^nextpow2(L);
Y = fft(ftemp,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure
semilogy(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
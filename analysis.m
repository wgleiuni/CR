step=1000;
t=zeros(step,50000);
Kx=t;
Ky=t;
E=t;
for i=1:1:step
    eval(['load data/data1/D',num2str(i),'.mat']);
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
cut=2000;

%%
figure
set(gcf,'color','w')
tI=-6+7/1000*(1:1000);
para=10.^tI;
%para=2/1000*(1:1000);
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
D=D1;
figure
set(gcf,'color','w','position',[2000 280 560 420])
subplot(2,1,1)
plot(D(:,1),D(:,2),'b')
hold on
plot(D(:,1),D(:,3),'r')
subplot(2,1,2)
plot(D(:,1),D(:,4))
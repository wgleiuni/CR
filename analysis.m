load('D1.mat');
D=D1;
figure
set(gcf,'color','w','position',[2000 280 560 420])
subplot(2,1,1)
plot(D(:,1),D(:,2),'b')
hold on
plot(D(:,1),D(:,3),'r')
subplot(2,1,2)
plot(D(:,1),D(:,4))
function dy=CR_diff(t,y,F,i)

dy=zeros(2,1);
w=1;
p=sqrt(y(1)^2+y(2)^2);
if i==1 % semi-Dirac
    temp=sqrt(y(1)^4+y(2)^2);
    dy(1)=-y(2)/(2*temp)-1*F*cos(w*t);
    dy(2)=y(1)^3/temp-1*F*sin(w*t);
elseif i==2 % massive
    dy(1)=-y(2)-F*cos(w*t);
    dy(2)=y(1)-F*sin(w*t);
elseif i==3 % massless
    dy(1)=-y(2)/p-F*cos(w*t);
    dy(2)=y(1)/p-F*sin(w*t);
end
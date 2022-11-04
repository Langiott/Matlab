 clear all;
close all;
clc;
global omega;
n=10;
ta=0;
tb=pi/2;

err=inf;
toll=1e-6;
temp=0;
omega=1;
u(1)=0.; v(1)=1.0; 
K(1)=fK(u(1),v(1));
while err>toll
    t=linspace(ta,tb,n);
    h=t(2)-t(1);
    h2=0.5*h;
% fai runge kutta
    for i=1:n-1
        k1=f(t(i),u(i),v(i));
        l1=g(t(i),u(i),v(i));
        k2=f(t(i)+h2,u(i)+h2*k1,v(i)+h2*l1);
        l2=g(t(i)+h2,u(i)+h2*k1,v(i)+h2*l1);
        k3=f(t(i)+h2,u(i)+h2*k2,v(i)+h2*l2);
        l3=g(t(i)+h2,u(i)+h2*k2,v(i)+h2*l2);
        k4=f(t(i)+h,u(i)+h*k3,v(i)+h*l3);
        l4=g(t(i)+h,u(i)+h*k3,v(i)+h*l3);
        u(i+1)=u(i)+h*(k1+2*k2+2*k3+k4)/6;
        v(i+1)=v(i)+h*(l1+2*l2+2*l3+l4)/6;
        K(i+1)=fK(u(i+1),v(i+1));
        %fprintf ("%e\n", K(i+1)-K(1));
    end
    %posPI2=round((pi/2)/h);
    err=abs(v(end)-temp);
    temp=v(end);
    n=n*2;
end
disp(err);
figure
plot(t,u,"r",t,v,"b")
figure
plot(u,v,"k")
figure
plot(t,K-K(1),"k")
axis([0 tb -0.001 0.001]);
function z=f(t,u,v)
    z=v;
end
function z=g(t,u,v)
global omega;
    z=-omega^2*sin(u);
end
function z=fK(u,v)
global omega;
    z=v^2/2-omega^2*cos(u);
end
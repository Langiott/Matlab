clear all;
close all;
clc;
global omega eps;
n=10001;
ta=0;
tb=2*pi;
h=(tb-ta)/(n-1);
h2=0.0049;
t=linspace(ta,tb,n);
omega=1;
eps=0.;
u(1)=0.0; v(1)=1% aumentando la velocità si ottengono le oscillazione intorno allo zero 
K(1)=fK(u(1),v(1));
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
ut=abs(u(n)-u(1))
vt=abs(v(n)-v(1))
figure
plot(t,u,"r",t,v,"b")
hold on
figure
plot(u,v,"k")
hold on
figure
plot(t,K,"k")
axis([0 tb 0 1]);
hold on
function z=f(t,u,v)
    z=v;
end
function z=g(t,u,v)
global omega eps;
    z=-omega^2*u-eps*u^3;
end
function z=fK(u,v)
global omega eps;
    z=v^2/2-omega^2*u^2/2+eps*u^4/4;
end

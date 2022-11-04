close all;
clear all;
clc;
% analizziamo il caso in cui v e D sono costanti 
global v D ya yb;
n=11;
nm2=n-2;
v=-0.2;
D=0.1;
xa=0;
xb=1;
ya=0;
yb=1;
h=(xb-xa)/(n-1);
hsq=h*h;
x=linspace(xa,xb,n);
% il resto rimane uguale
for i=1:n
    p(i)=fp(x(i));
    q(i)=fq(x(i));
    r(i)=fr(x(i));
end
a=zeros(n,n);
a(1,1)=1;
a(n,n)=1;
for i=2:n-1
    a(i,i)=-2*(2+hsq*q(i));
    a(i,i-1)=2+h*p(i);
    a(i,i+1)=2-h*p(i);
end
b(1,1)=ya;
b(n,1)=yb;
for i=2:n-1
    b(i,1)=2*hsq*r(i);
end
for i=1:n
    z(i)=zsol(x(i));
end
l(1,1)=a(1,1);
u(1,2)=a(1,2)/l(1,1);
y(1)=b(1)/l(1,1);
for i=2:n
    l(i,i-1)=a(i,i-1);
end
for i=2:n-1
    l(i,i-1)=a(i,i-1);
    l(i,i)=a(i,i)-l(i,i-1)*u(i-1,i);
    u(i,i+1)=a(i,i+1)/l(i,i);
    y(i)=(b(i)-l(i,i-1)*y(i-1))/l(i,i);
end
l(n,n)=a(n,n)-l(n,n-1)*u(n-1,n);
y(n)=(b(n)-l(n,n-1)*y(n-1))/l(n,n);
ww(n)=y(n);
ww(n-1)=y(n-1)-u(n-1,n)*ww(n);
for i=n-1:-1:2
    ww(i-1)=y(i-1)-u(i-1,i)*ww(i);
    ww(1)=y(1)-u(1,2)*ww(2);
end
w=linsolve(a,b);
plot(x,w,"xr",x, ww,"ob",x,z)
function z=fp(x)% costruisco la variabile P
global v D; 
    z=v/D;
end
function z=fq(x) % z è uguale a zero
    z=0;
end
function z=fr(x)% r è uguale a zero
    z=0;
end
function z=zsol(x) % in questo caso è disponibile a soluzione esatta 
global v D ya yb;
    b=(ya-yb)/(1-exp(v/D));
    a=ya-b;
    z=a+b*exp(v/D*x);
end

close all;
clear all;
clc;


%% dati
a=0;
b=4*pi;
zr=11;
toll=1e-5;

%% esercizio 1 

[x,runge]=RK4(a,b,toll);
nr=length(runge);
errore=abs(esatta(4*pi)-runge(end));
hold on;
plot(x,runge,"o");
yline(0)


%% interpolazione runge 

% mi costruisco l'asse x per i nodi
z1=linspace(a,b,zr);

% costruisco i vettre dei punti
z2=linspace(a,b,101);

% costruiscoo i nodi della funzione

for s=1:zr
    fattore = round((nr-1)/(zr-1));
    xs(s)=x((s-1)*fattore+1);
    zs(s)=runge((s-1)*fattore+1);   
end 
%plot(z1,zs,"or");
for i=1:101
    inter(i) =interpol(xs,zs,z2(i),zr);
end
%plot(z2,inter)

%% esercizio 2 
%interpolazione lineare
zrunge=runge-0.5; % scalo runge di 0,5
plot(x,zrunge);% lo plotto
count =0;
xint=[];
for j=1:70
    if zrunge(j)*zrunge(j+1)<0
        count=j;
        xint(end+1)=interx(zrunge(j),zrunge(j+1),x(j),x(j+1));%mi salvo i punti
    end
end
xint

%% funzioni

function w=interpol(x,y,z,n)
    sum=0;
    for k=1:n
        num=1;
        den=1;
        for j=1:k-1
            num=num*(z-x(j));
            den=den*(x(k)-x(j));
        end
        for j=k+1:n
            num=num*(z-x(j));
            den=den*(x(k)-x(j));
        end
        lnk=num/den;
        sum=sum+y(k)*lnk;
    end
    w=sum;
end


function [xt,rk4]=RK4(ta,tb,toll)
err=inf;
appo=2;
n=11;
    while err>toll
    t=linspace(ta,tb,n);
    h=(tb-ta)/(n);
    h2=0.5*h;
    u=t;
    u(1)=2;
    for i=1:n-1
        k1=f(t(i),u(i));
        k2=f(t(i)+h2,u(i)+h2*k1);
        k3=f(t(i)+h2,u(i)+h2*k2);
        k4=f(t(i)+h,u(i)+h*k3);
        u(i+1)=u(i)+h*(k1+2*k2+2*k3+k4)/6;
    end
    err=abs(u(end)-appo);
    appo=u(end);
    n=n+10;
    end
    rk4=u;
    xt=t;
    function  h=f(x,y)
        h=-y*sin(x);
    end
end
function g=esatta(x)
g=2*exp(cos(x)-1);
end
function retta=interx(y1,y2,x1,x2)
    retta=(x2-x1)*(0-y1)/(y2-y1)+x1;
end
    
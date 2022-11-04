clear all;
close all;
clc;
a=0;
b=4*pi;
toll=1e-5;
[runge,t]=RK4(a,b,toll,2);% chiamata della funzione 
lung=length(runge);% ora posso creare la variabile x per plottare 
x=linspace(a,b,lung);
plot(x,runge);


function [rk4,tn]=RK4(ta,tb,toll,s)
%input:ta,tb=estremi,toll=tolleranza,s=condizione iniziale della y
err=inf;
appo=0;
n=2;
    while err>toll
    t=linspace(ta,tb,n);
    h=(tb-ta)/(n);
    h2=0.5*h;
    u=t;
    u(1)=s;% indicare la condizione iniziale
    for i=1:n-1
        k1=f(t(i),u(i));
        k2=f(t(i)+h2,u(i)+h2*k1);
        k3=f(t(i)+h2,u(i)+h2*k2);
        k4=f(t(i)+h,u(i)+h*k3);
        u(i+1)=u(i)+h*(k1+2*k2+2*k3+k4)/6;
    end
    err=abs(u(end)-appo);
    appo=u(end);
    n=n+1;
    end
    rk4=u;
    tn=t;

    function  h=f(x,y)
        h=-y*sin(x);
    end
end
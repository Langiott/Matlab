clear all;
close all;
clc;
a=0;
b=6;
toll=1e-5;
[runge,t]=RK4(a,b,toll,5,0);% chiamata della funzione 
lung=length(runge);% ora posso creare la variabile x per plottare 
x=linspace(a,b,lung);
plot(x,runge);


function [rk4,tn]=RK4(ta,tb,toll,s,r)
%input:ta,tb=estremi,toll=tolleranza,s=condizione iniziale della y
err=inf;
appo=0;
n=2;
    while err>toll
    t=linspace(ta,tb,n);
    h=(tb-ta)/(n);
    h2=0.5*h;
    u=zeros(1,n); % questo rappresenta la x/y
    v=zeros(1,n);% questo rappresenta la x'/y'
    u(1)=s; v(1)=r; 

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
    
    end
    err=abs(u(end)-appo);
    appo=u(end);
    n=n+1;
    end
    rk4=u;
    tn=t;

    function z=f(t,u,v)
        z=v;
    end
    function z=g(t,u,v)
        z=-1-v;
    end
end
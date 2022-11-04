clear all;
close all; 
clc;

%% esercizio 1)

% introduzione delle variabili 
n1=11;
ta=0;
tb=2*pi;
tb2=3;
toll=1.e-5;
s=0;
r=1;
[rungevet,t]=RK4(ta,tb,toll,s,r);
[rungevet2,t2]=RK4(ta,tb2,toll,s,r);
hold on;
plot(t,rungevet,"or");

[yy,dyy,ddyy,xx1,a1,b1,c1,d1]=splin(t,rungevet,n1);
[yy2,dyy2,ddyy2,xx12,a2,b2,c2,d2]=splin(t2,rungevet2,n1);
plot(xx1,yy);
format short;
dyy2;


%% bisezione 

tolass=toll;% definisco la tolleranza 
va=ta;
vb=tb2;
na(1)=va; % inizializzo i vettori per il calcolo della funzione  
nb(1)=vb;
nmax=25;
nstim=log2((vb-va)/tolass);% definisco la stima del numero di iterazioni
fprintf ("n stimato = %d\n", nstim);% lo stampo in output

for n=1:nmax
    ic=0;
    cc(n)=0.5*(na(n)+nb(n))% valore intermedio 
    %disp(n);
    %disp(c(n));
    % if (abs(f(c(n)))<=tolf)
    if (abs(nb(n)-na(n))<=tolass)% verifica della tolleranza 
        % if (abs((sb(n)-sa(n))/c(n))<=tolrel)% tolleranza relativa
        z=cc(n);% mi salvo l'ultimo valore di c preima di fermare tutto
        break
    end

    fa=ds(8,na(n),a2,b2,c2,d2,t2); % mi calcolo il valore di f in a 
    fc=ds(8,cc(n),a2,b2,c2,d2,t2);% mi calcolo il valore di f in b
    if (fa*fc < 0.0) % procedo con l'analisi degli intervalli
        na(n+1)=na(n);% preparo il valore successivo di sa
        nb(n+1)=cc(n);% preparo il valore successivo di sb
    else
        na(n+1)=cc(n);
        nb(n+1)=nb(n);

    end
end
% zzero della derivata 
z

% stampo lo zero della eq.differenziale
y= ds(8,z,a2,b2,c2,d2,t2);
plot(xx1,dyy);
plot(z,y,"Og");


nrunge=abs(z-t);
[app,rku]=min(nrunge)
plot(z,rungevet(rku),"Og");



%% funzioni

%% RK4 vettoriale

function [rk4,tn]=RK4(ta,tb,toll,s,r)
%input:ta,tb=estremi,toll=tolleranza,s=condizione iniziale della y
err=inf;
appo=0;
n=11;

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
    
    rk4=u;
    tn=t;
    function z=f(~,~,v)
        z=v;
    end
    function z=g(t,~,~)
        z=-t*sin(t);
    end

end

%% metodo interpolazione

function[yy,dyy,ddyy,xx1,a,b,c,d]=splin(x,vel,n1)

    n=n1-1; % lo faccio perchè l'indice nesimo me lo darà il termine i+1
    
    for i=1:n1
        y(i)=vel(i);% mi costruisco la funzione originale
    end
 
    n=n1-1; % lo faccio perchè l'indice nesimo me lo darà il termine i+1
   
    for i=1:n
        h(i)=x(i+1)-x(i);% costruisco il vettore degli h
    end
    
    bb(1)=0;% setto a zero l'elemento1
    bb(n1)=0;% setto a zero l'elemento n
    
    for i=2:n
        bb(i)=3*(y(i+1)-y(i))/h(i)-3*(y(i)-y(i-1))/h(i-1); %calcolo b dell' equazione risolutiva
    end
    
    a=y; % mi sono determinato i coefficienti di s(k)=f(k)
    
    for i=1:n1
        aa(i,i)=0;% predispongo la diagonale principale 
    end
    
    aa(1,1)=1;% il primo elemento della diagonale lo setto a 1
    aa(n1,n1)=1;
    
    for i=2:n
        % mi costruisco la matrice tri-diagonale
        aa(i,i)=2*(h(i-1)+h(i));% mi costruisco la matrice tridiagonale
        aa(i,i-1)=h(i-1);
        aa(i,i+1)=h(i);
    end
    
    c=linsolve(aa,bb'); % souzione del sistema, mi serve per determinare le c
    
    for i=1:n
        d(i)=(c(i+1)-c(i))/(3*h(i));% mi ricavo il coefficiente d
        b(i)=(a(i+1)-a(i))/h(i)-h(i)*(2*c(i)+c(i+1))/3;% mi ricavo il coefficiente b
    end
    
    npts=1000; %numero punti in ciascun intervallo
    ki=1;
   
    for k=1:n
        xx=linspace(x(k),x(k+1),npts);% costruisco la funzione polinomiale a tratti
        for i=1:npts % interpolo con 51 punti
            xx1(ki)=xx(i); % determino il primo coefficiente
            yy(ki)=s1(k,xx(i),a,b,c,d,x);% interpolazione della funzione semplice
            dyy(ki)=ds(k,xx(i),a,b,c,d,x);% interpolazione della derivata della funzione 
            ddyy(ki)=dds(k,xx(i),a,b,c,d,x);% interpolazione della derivata seconda
            ytrue(ki)=vel(k);% 
            ki=ki+1;
        end
    end
    % mi sono cstruito la funzione interpolante della funzione semplice  
end
%% funzioni di supporto alla interpolazione

function w=s1(k,z,a,b,c,d,x) % costruisco s(k)
    w=a(k)+b(k)*(z-x(k))+c(k)*(z-x(k))^2+d(k)*(z-x(k))^3; % mi calcolo s(k)
end

function w=ds(k,z,a,b,c,d,x)% faccio la derivata di s(k)
    w=b(k)+2*c(k)*(z-x(k))+3*d(k)*(z-x(k))^2;
end

function w=dds(k,z,a,b,c,d,x)% faccio la derivata seconda si s(k)
    w=2*c(k)+6*d(k)*(z-x(k));
end

function w=ints(k,z,a,b,c,d,x) % integrale con il metodo di simpson     
    w=a(k)*z+0.5*b(k)*(z-x(k))^2+(c(k)/3)*(z-x(k))^3+(d(k)/4)*(z-x(k))^4;
end







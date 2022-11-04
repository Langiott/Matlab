close all;
clear all;
clc;

%definisco gli estremi
a=0;
b=2;
toll=1.e-6;
n=11;
np=10;
t=linspace(0,2,n);

%passo per runge kutta
h1=(b-a)/(n-1);
h2=(0.5*h1);

%passo integrale
h3=(t(2)-t(1));

%% esercizio 1)

% soluzione esatta eq.differenziale
yse=g(t);

%plot(x,runge);% abbiamo graficato la soluzione esatta

u(1)=1;
for i=1:n-1
k1=f(t(i),u(i));
k2=f(t(i)+h2,u(i)+h2*k1);
k3=f(t(i)+h2,u(i)+h2*k2);
k4=f(t(i)+h1,u(i)+h1*k3);
u(i+1)=u(i)+h1*(k1+2*k2+2*k3+k4)/6;
end

% individuazine dell'indice con i valore più vicinno all'1

tn1 = abs(1-t);
[tmin,kl]=min(tn1);

errn=abs(u(kl)-g(1));

%% esercizio 2)

% interpolo
z=linspace(a,b,100); % vettore per graficare 

inter(1)=1;
for j=1:100
    inter(j)=interpol(t,u,z(j));
end

hold on
plot(t,yse);% abbiamo graficato la soluzione esatta 
plot(t,u,"or");% nodi per rungekutta
plot(z,inter);% plot interpolazione
legend("soluzione esatta","valore ai nodi della funzione","interpolazione")


%% esercizio 3)

% ho calcolato semplicemente l'integrale di u
% in questo caso il passo è dispari
 ws(1)=1.0;% inizializzo a 1 il primo termine dei coefficienti
 ws(np+1)=1.0;% inizializzo a 1 gli elementi successivi

 for i=2:2:np-2
    ws(i)=4.0% inserisco 4 nei coefficienti parti
    ws(i+1)=2.0% inserisco i 2 nei coefficienti pari
 end
 ws(np)=4% questo è il penultimo elemento
 ints=h3*sum(ws.*u)/3.0% formula per simpson



%% funzioni 

function w=interpol(x,y,z)
%x,y vettori dei nodi
%z punto in cui interpolare
    sum=0;
    n=length(x);
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


function  h=f(x,y)
    h=(1-x)/y;
end

function seqd= g(x)
    seqd=sqrt(-x.^2+2*x+1);
end
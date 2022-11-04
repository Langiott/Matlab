clear all;
close all;
clc;
n=10;% numero di nodi
a=-1;
b=1; 
q=0;%inizializzo la variabile
e=linspace(a,b,n);% creo un vettore di dieci elementi che mi rappresenteranno i 
% nodi della funzione 
for i=1:n
     q(i)=f(e(i));% costruisco i punti della funzione
end
plot(e,q,"or");
hold on;grid 
xx=linspace(a,b,10000);% costruisco i vettore delle x della funzione interpolante
for i=1:10000
    %y1(i)=f(xx(i));
    yy(i)=interpol(e,q,xx(i));% salvo la funzione interpolata 
end
plot(xx,yy);
%axis([a b -0.5 1.5]);
function w=f(z)
   %w=tanh(z);
   %w=sin(z);
   w=2*exp(cos(z)-1);
   %w=alpha*z;
end

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


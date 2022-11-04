clear all;
close all;
clc;
% costruzione della funzione
a=0; b=1; %estremi dell' integrale
toll=1.e-5;
h=0.1
x=a:h:b
for i=1:length(x)
    y(i)=sims(a,x(i),toll)
end

hold on
plot(x,f(x))
plot(x,y,'o')

n=101;
z=linspace(a,b,n);
for i=1:n
    inter(i)=interpol(x,y,z(i));
end
plot(z,inter)
legend("e^{-x^2}","\Phi(x)","\Phi(x) interpol")

function Simpson=sims(a,b,toll)
ap=0;
np=14;
er=1;
    while er>toll
        h=(b-a)/np;
        x=linspace(a,b,np+1);
        for i=1:np+1
            y(i)=f(x(i));% creo i pu ti della funzione 
        end
        lung=length(y);
        % creo il vettore dei coefficienti per il metodo di simpson per i termini 
        % pari e dispari
        ws(1)=1.0;% inizializzo a 1 il primo termine dei coefficienti
        ws(np+1)=1.0;% inizializzo a 1 gli elementi successivi

        for i=2:2:np-2
            ws(i)=4.0;% inserisco 4 nei coefficienti parti
            ws(i+1)=2.0;% inserisco i 2 nei coefficienti pari
        end
        ws(np)=4;% questo è il penultimo elemento
        ints=h*sum(ws.*y)/3.0;% formula per simpson
        er=abs(ints-ap);
        ap=ints;
        np=np+2;
    end
    Simpson=ints;
   
end

function  g=f(x)
    g=exp(-x.^2);
end


function w=interpol(x,y,z)
%x,y vettori dei nodi
%z punto in cui interpolare
%restituisce f(z)
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




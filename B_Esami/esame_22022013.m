clear all;
close all;
clc;
%% esercizio 2)
% calcolia due intervalli ifferenti e poi ne facciamo la somma sfruttando
% la propietà di additività

% DATI
a0=0; b0=1/3;
a1=b0; b1=1;
tol=1.e-6;
n=6 ;

% costruisco il primo integrale 

% il linspace lo devo costruire in modo furbo sapendo
% che il numero degli elementi della funzione devono essere dispari 
t1=linspace(a0,b0,n+2);
h1=t1(2)-t1(1);
y=zeros(1,n+1);

% ciclo per il primo integrale 
for i=1:n+2
    if t1(i)< b0
        y(i)=w(t1(i));
    end
end

% calcolo secondo integrale

t2=linspace(a1,b1,n+1);
h2=t2(2)-t2(1);
y1=zeros(1,n+1);

% ciclo per il secondo integrale 
for i=1:n+1
    if t2(i)< b1
        y1(i)=s(t2(i));
    end
end
y1(end)=[];

% ora calcolo l'integrale
ws(1)=1.0;% inizializzo a 1 il primo termine dei coefficienti  
ws(n+1)=1;
for i=2:2:n-2
   ws(i)=4.0% inserisco 4 nei coefficienti parti
   ws(i+1)=2.0% inserisco i 2 nei coefficienti pari
end 
ws(n)=4
% questo è il penultimo elemento
ints1=h1*sum(ws.*y)/3.0;% formula per simpson

% calcolo i secondo integrale 
ws1(1)=1.0;% inizializzo a 1 il primo termine dei coefficienti
ws1(n)=1.0;    
for i=2:2:n-2
   ws1(i)=4.0;% inserisco 4 nei coefficienti parti
   ws1(i+1)=2.0;% inserisco i 2 nei coefficienti pari
end 
ws1(n-1)=4;% questo è il penultimo elemento
ints2=h2*sum(ws1.*y1)/3.0;% formula per simpson

fprintf("l'integrale tra 0 e 1/3 vale %f\n mentre il secondo compreso tra 1/3 e 1 vale %f\t\n ",ints1,ints2);
ints=ints1+ints2;
fprintf("l'integrale complessivamente vale %f\n",ints)


%% funzioni

function g=w(x)
    g=exp(x);
end
function h=s(x)
    h=exp(-x)+2;
end
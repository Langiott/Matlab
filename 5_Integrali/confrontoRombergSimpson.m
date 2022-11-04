clear all;
close all;
clc;

x=linspace(1,5,5);
y=[2.4142,2.6734,2.8974,3.0976,3.2804];

[dum,n]=size(x);
a=x(1);
b=x(n);
h(1)=b-a;
r(1,1)=0.5*h(1)*(y(1)+y(n));
k=0;
for k=2:n-2
    m=2^(k-1);
    h(k)=0.5*h(k-1);
    sumo=0.0;
    for l=1:2:m
        sumo=sumo+y(1+l);
    end
    r(k,1)=0.5*r(k-1,1)+h(k)*sumo;
    for j=2:n
        r(k,j)=r(k,j-1)+(r(k,j-1)-r(k-1,j-1))/(4^(j-1)-1);
    end
    %fprintf ("k, %d\t%e\n",k, abs(r(k,k)-r(k-1,k-1)));
   
end
 fprintf ("k, %f\n",k, r(k,k));
 
 %% applichimao il metodo di simpson
 ns = 5;
 ws(1)=1.0;
 ws(ns)=1.0;
 us=[2.4142,2.6734,2.8974,3.0976,3.2804];

for s=2:2:ns-2
    ws(s)=4.0;
    ws(s+1)=2.0;
end

%il passo lho determinato facendo la differenza tra due 
%t successive t(1),t(2)
ws(ns-1)=4.0;
h6=1;
int(1)=h6*sum(us.*ws)/3


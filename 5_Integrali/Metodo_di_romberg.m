clear all
clc

f=@(x) (1/5)*x.^4-(5/4)*x.^3+(3/2)*x.^2;
a=0; 
b=4.5488; %punti a e b dall'esercizio 2
tol=1.e-8;
itmax=50;
If=romberg(f,a,b,tol,itmax) %calcolo gli integrali con romberg

function[I]=romberg(func,a,b,tol,kmax) %function che implementa il metodo di romberg

kmax=abs(kmax);
R=zeros(1,kmax+1);
err=1;
Ip=	0;
R(1)=((b-a)/2)*(func(a)+func(b));
k=1;
while(err>tol*abs(Ip))
	R(k+1)=trapezi(func,a,b,k+1,R(k));
    for j=k:-1:1
        p=4^(k-j+1);
        R(j)=(p*R(j+1)-R(j))/(p-1);
    end
  	err=abs(R(1)-Ip);
    Ip=R(1);
    k=k+1;
    if k==kmax
        warning('non converge, massime iterazioni raggiunte');
        break;
    end
end
I=R(1);

function[I]=trapezi(func,a,b,k,Ip) %function che implementa il metodo dei trapezi
I               =	0;
H               =	b-a;
for i=1:2^(k-2)
	I           =	I+func(a+(2*i-1)*H/(2^(k-1)));
end
I               =   0.5*Ip+(H/(2^(k-1)))*I;
end
end
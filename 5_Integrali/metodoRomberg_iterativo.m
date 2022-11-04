clear all;
clc;
tol=10^(-6);
a=0.0;
b=1.0;
n=10;
h(1)=b-a;
true=fp(b)-fp(a);
r(1,1)=0.5*h(1)*(f(a)+f(b));
k=2;
test=1.0;
while test>tol & k<=n
    m=2^(k-1);
    h(k)=0.5*h(k-1);
    sum=0.0;
    for l=1:2:m
        sum=sum+f(a+l*h(k));
    end
    r(k,1)=0.5*r(k-1,1)+h(k)*sum;
    for j=2:n
        r(k,j)=r(k,j-1)+(r(k,j-1)-r(k-1,j-1))/(4^(j-1)-1);
    end
    test=abs(r(k,k)-r(k-1,k-1));
    fprintf ("k, %d\t%e\t%e\n",k, r(k,k),true-r(k,k));
    k=k+1;
end
fprintf ("true = %e\n error = %e\n", true, r(k-1,k-1)-true);
function y=f(x)
y=x^2*exp(-x);
end
function y=fp(x)
y=(-x^2-2*x-2)*exp(-x);
end

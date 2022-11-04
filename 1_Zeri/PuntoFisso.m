clear all;
close all;
clc;
tolass=1e-6;
x= linspace(1,2,100);
nmax=100;
x0=0.8;
x(1)=x0;
for n=1:nmax
    x(n+1)=g1(x(n));
        %if (abs(f(x(n+1)))<=tolf)
        if (abs(x(n+1)-x(n))<=tolass)
            break;
        end
        
end
x(n+1)
g1(x(n+1))
%printf ( z, g1(z));

function y=g1(x)
y=sqrt((x+3-x.^4)./2);
end


% Cuadratura Gaussiana
syms x
f=(pi/24*x+pi/8)^2*sin(pi/24*x+pi/8)*pi/24
p=x^6-(15/11)*x^4+5/11*x^2-5/231
n=6
X=double(solve(p))
Y=double(subs(f,X))
h=diff(p)
for i=1:n
d=int(p/(x-X(i)))
m=subs(d,1)-subs(d,-1)
W(i)=double(norm(1/subs(h,X(i))*m))
end
double(W)
integral=0
for i=1:n
integral=integral+W(i)*Y(i)
end
vpa(integral)
clear
% Fin cuadratura gaussiana

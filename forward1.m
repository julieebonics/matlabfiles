Matlab Function in Euler's Forward Method for Solving Ordinary Differential Equations

clear;

L = 1;
T = .5;
n = 10;
maxk = 50;
K = T/maxk;
h = L/16;
alpha = 1;
r = (alpha*K)/(h^2);



for i = 1:n+1
    x(i) = (i-1)*h;
    u(i,1) = sin(pi*x(i));
end
for k = 1:maxk+1
    u(1,k) = 0;
    u(n+1, k) = 0;
    time(k) = (k-1)*K;
end

for k=1:maxk
    for i=2:n;
        u(i, k+1) = u(i,k) + .5*r*(u(i-1,k) + u(i+1,k) -2.*u(i,k));
    end
end

[x,t] = meshgrid(x,time);

utrue = exp(-pi.^2*t).*sin(pi*x);



nrm = norm(u,inf)
err = norm((utrue-u'),inf)


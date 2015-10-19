Matlab function for Backward Step Method for solving Ordinary Differential Equations
clear;

L = 1;
T = .5;
n = 10;
maxk = 50;
K = T/maxk;
h = L/64;
alpha = 1;
r = 2*(alpha*K)/(h^2);


for i = 1:n+1
    x(i) = (i-1)*h;
    u(i,1) = sin(pi*x(i));
end
for k = 1:maxk+1
    u(1,k) = 0;
    u(n+1, k) = 0;
    time(k) = (k-1)*K;
end

aa(1:n-2) = -r;
bb(1:n-1) = 1+2*r;
cc(1:n-2) = -r;
MM = inv(diag(bb,0)+diag(aa,-1)+diag(cc,1));

for k=2:maxk
    uu = u(2:n,k-1);
    u(2:n,k) = MM*uu;

end

[x,t] = meshgrid(x,time);

utrue = exp(-pi.^2*t).*sin(pi*x);

%figure(1)
%surf(x,time,u')
%title('Approximation at h =1/10 and k=1/100')

%figure(2)
%surf(x,time, utrue)
%title('True Solution')

%figure(3)
%mesh(x,time, u'-utrue)

nrm = norm(u',inf)
err = norm((utrue-u'),inf)
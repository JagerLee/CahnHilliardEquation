function E = energy(u,n,t)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
global  k X derta
N_x=X*k/t;
h=X/N_x;
u_1=(u(2:end,:)-u(1:end-1,:))/h;
if n==1
    I=F(u);
else
    I=(u.^2-1).^2/(4*derta);
end
I=I(2:end,:)+0.5*u_1.^2;
E=h*sum(I,1);
end


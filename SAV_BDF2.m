function u=SAV_BDF2(t)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
global k X beta derta C_0 T
count=round(T/t);
N_x=round(X*k/t);
h=X/N_x;
u0=initial(h);
[Lap_1,Lap_2]=Laplace(h);
tmp_1=Lap_1(3:N_x-1,3:N_x-1);
tmp_1(1,1)=-0.5/h^2;tmp_1(1,2)=0.5/h^2;tmp_1(end,end-1)=0.5/h^2;tmp_1(end,end)=-0.5/h^2;
u=zeros(N_x+1,count+1);

A=eye(N_x-3)+t*Lap_2(3:N_x-1,3:N_x-1)-(beta*t/derta)*tmp_1;
B=eye(N_x-3)+t*Lap_2(3:N_x-1,3:N_x-1)-(t*beta/derta)*tmp_1;
u1=u0;
    r_n=sqrt(E_1(u1,h)+C_0);
    u_n_half=ones(N_x+1,1);
    tmp=Lap_1*dF(u1);
    u_n_half(3:N_x-1)=B\(u1(3:N_x-1)+t*tmp(3:N_x-1));
    u_n_half(1)=3/2*u_n_half(3)-1/2*u_n_half(4);u_n_half(2)=u_n_half(1);
    u_n_half(N_x)=3/2*u_n_half(N_x-1)-1/2*u_n_half(N_x-2);u_n_half(N_x+1)=u_n_half(N_x);
    b_n=dF(u_n_half)/sqrt(E_1(u_n_half,h)+C_0);
    L_b=Lap_1*b_n;
    tmp_A=A;
    tmp_A(:,1)=A(:,1)-1.5*0.5*t*h*(b_n(1)+b_n(2))*L_b(3:N_x-1);
    tmp_A(:,2)=A(:,2)+0.5*0.5*t*h*(b_n(1)+b_n(2))*L_b(3:N_x-1);
    tmp_A(:,end)=A(:,end)-1.5*0.5*t*h*(b_n(end)+b_n(end-1))*L_b(3:N_x-1);
    tmp_A(:,end-1)=A(:,end-1)+0.5*0.5*t*h*(b_n(end)+b_n(end-1))*L_b(3:N_x-1);
    g_n=u1+t*(r_n-1/2*inner_product(b_n,u1,h))*L_b;
    tmp_b=tmp_A\L_b(3:N_x-1);
    gamma_n=-inner_product(b_n(3:N_x-1),tmp_b,h);
    tmp_g=tmp_A\g_n(3:N_x-1);
    inner_p=inner_product(b_n(3:N_x-1),tmp_g,h)/(1+t*gamma_n/2);
    u1(3:N_x-1)=0.5*t*inner_p*tmp_b+tmp_g;
    u1(1)=1.5*u1(3)-0.5*u1(4);u1(2)=u1(1);
    u1(end)=1.5*u1(N_x-1)-0.5*u1(N_x-2);u1(end-1)=u1(end);

u(:,1:2)=[u0 u1];
A=eye(N_x-3)+(2*t/3)*Lap_2(3:N_x-1,3:N_x-1)-(2*t*beta/(3*derta))*tmp_1;
B=eye(N_x-3)+2*t/3*Lap_2(3:N_x-1,3:N_x-1)-(2*t*beta/derta/3)*tmp_1;
disp(['------------t=',num2str(t),'-----------']);
tic
for n=2:count
    r_n=sqrt(E_1(u(:,n),h)+C_0);
    r_n_1=sqrt(E_1(u(:,n-1),h)+C_0);
    u_n_half=ones(N_x+1,1);
    tmp=Lap_1*dF(u(:,n));
    u_n_half(3:N_x-1)=B\(4/3*u(3:N_x-1,n)-1/3*u(3:N_x-1,n-1)+2*t/3*tmp(3:N_x-1));
    u_n_half(1)=3/2*u_n_half(3)-1/2*u_n_half(4);u_n_half(2)=u_n_half(1);
    u_n_half(N_x)=3/2*u_n_half(N_x-1)-1/2*u_n_half(N_x-2);u_n_half(N_x+1)=u_n_half(N_x);
    b_n=dF(u_n_half)/sqrt(E_1(u_n_half,h)+C_0);
    L_b=Lap_1*b_n;L_b(1)=L_b(2);L_b(end)=L_b(end-1);
    tmp_A=A;
    tmp_A(:,1)=A(:,1)-1.5/3*t*h*(b_n(1)+b_n(2))*L_b(3:N_x-1);
    tmp_A(:,2)=A(:,2)+0.5/3*t*h*(b_n(1)+b_n(2))*L_b(3:N_x-1);
    tmp_A(:,end)=A(:,end)-1.5/3*t*h*(b_n(end)+b_n(end-1))*L_b(3:N_x-1);
    tmp_A(:,end-1)=A(:,end-1)+0.5/3*t*h*(b_n(end)+b_n(end-1))*L_b(3:N_x-1);
    g_n=4/3*u(:,n)-1/3*u(:,n-1)+2*t/9*(4*r_n-r_n_1-1/2*inner_product(b_n,4*u(:,n)-u(:,n-1),h))*L_b;
    tmp_b=tmp_A\L_b(3:N_x-1);
    gamma_n=-inner_product(b_n(3:N_x-1),tmp_b,h);
    tmp_g=tmp_A\g_n(3:N_x-1);
    inner_p=inner_product(b_n(3:N_x-1),tmp_g,h)/(1+t*gamma_n/3);
    u(3:N_x-1,n+1)=(t/3*inner_p)*tmp_b+tmp_g;
    u(1,n+1)=3/2*u(3,n+1)-1/2*u(4,n+1);
    u(2,n+1)=u(1,n+1);
    u(N_x,n+1)=3/2*u(N_x-1,n+1)-1/2*u(N_x-2,n+1);
    u(N_x+1,n+1)=u(N_x,n+1);
    if mod(n,2000)==0
        disp(['progress ',num2str(n*100/count),'%']);
        toc
        disp('///');
    end
end
disp('progress 100%');
toc
end


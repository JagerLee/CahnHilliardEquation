function  u= SAV_BDF1(t)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
global k X beta derta C_0 T
count=round(T/t);
N_x=round(X*k/t);
h=X/N_x;
u0=initial(h);
[Lap_1,Lap_2]=Laplace(h);
tmp_1=Lap_1(3:N_x-1,3:N_x-1);
tmp_1(1,1)=-0.5/h^2;tmp_1(1,2)=0.5/h^2;tmp_1(end,end-1)=0.5/h^2;tmp_1(end,end)=-0.5/h^2;
A=eye(N_x-3)+t*Lap_2(3:N_x-1,3:N_x-1)-(beta*t/derta)*tmp_1;
B=eye(N_x-3)+t*Lap_2(3:N_x-1,3:N_x-1)-(t*beta/derta)*tmp_1;
u1=u0;
u=ones(N_x+1,count+1);
% w=u;
u(:,1)=u0; 
disp(['------------t=',num2str(t),'-----------']);
tic
% tmp_A=inv(A);
% inv_B=inv(B);
r_n=sqrt(E_1(u1,h)+C_0);
for i=1:count

    u_n_half=ones(N_x+1,1);
    tmp=Lap_1*dF(u(:,i));
    u_n_half(3:N_x-1)=B\(u(3:N_x-1,i)+t*tmp(3:N_x-1));
    u_n_half(1)=3/2*u_n_half(3)-1/2*u_n_half(4);u_n_half(2)=u_n_half(1);
    u_n_half(N_x)=3/2*u_n_half(N_x-1)-1/2*u_n_half(N_x-2);u_n_half(N_x+1)=u_n_half(N_x);
    b_n=dF(u_n_half)/sqrt(E_1(u_n_half,h)+C_0);
    L_b=Lap_1*b_n;
    tmp_A=A;
    tmp_A(:,1)=A(:,1)-1.5*0.5*t*h*(b_n(1)+b_n(2))*L_b(3:N_x-1);
    tmp_A(:,2)=A(:,2)+0.5*0.5*t*h*(b_n(1)+b_n(2))*L_b(3:N_x-1);
    tmp_A(:,end)=A(:,end)-1.5*0.5*t*h*(b_n(end)+b_n(end-1))*L_b(3:N_x-1);
    tmp_A(:,end-1)=A(:,end-1)+0.5*0.5*t*h*(b_n(end)+b_n(end-1))*L_b(3:N_x-1);
    g_n=u(:,i)+t*(r_n-1/2*inner_product(b_n,u(:,i),h))*L_b;
    tmp_b=tmp_A\L_b(3:N_x-1);
    gamma_n=-inner_product(b_n(3:N_x-1),tmp_b,h);
    tmp_g=tmp_A\g_n(3:N_x-1);
    inner_p=inner_product(b_n(3:N_x-1),tmp_g,h)/(1+t*gamma_n/2);
    u(3:N_x-1,i+1)=0.5*t*inner_p*tmp_b+tmp_g;
    u(1,i+1)=1.5*u(3,i+1)-0.5*u(4,i+1);u(2,i+1)=u(1,i+1);
    u(end,i+1)=1.5*u(N_x-1,i+1)-0.5*u(N_x-2,i+1);u(end-1,i+1)=u(end,i+1);
    r_n=r_n+0.5*inner_product(b_n,u(:,i+1)-u(:,i),h);
%     w(:,i+1)=-Lap_1*u1+(r_n )
    if mod(i,200)==0
        disp(['progress ',num2str(i*100/count),'%']);
        toc
        disp('///');
    end
end
disp('progress 100%');
toc
end

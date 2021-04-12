clear
global derta beta C_0 k X T
k=0.05/pi;
derta=0.01;
beta=4;
C_0=0;
T=0.1;
X=2*pi;

t=1.5*1e-4;
N_t=round(T/t);
N_x=round(X*k/t);
h=X/N_x;
xx=linspace(0,X,N_x+1);
tt=linspace(0,T,N_t+1);

u=SAV_BDF1(t);
% [xxx,ttt]=meshgrid(xx,tt);
% mesh(xxx',ttt', u); 
% xlabel('x')
% ylabel('t')
% zlabel('¦Õ(x,t)')
% title('Solution in [0,2¦Ð]X[0,T]')

E0=energy(u,0,t);
plot(tt,E0)
xlabel('t')
ylabel('E(¦Õ)')
title('Energy dissipation of SAV')
M=h*sum(u,1); 
plot(tt,M)
xlabel('t')
ylabel('mass')
title('Mass conservation of SAV')

tn=5;
derta_t=0.1*2^-6*(2/3).^(1:tn);

u_1=SAV_BDF1(derta_t(1));
u_2=SAV_BDF1(derta_t(2));
u_3=SAV_BDF1(derta_t(3));
u_4=SAV_BDF1(derta_t(4));
u_5=SAV_BDF1(derta_t(5));
% u_6=SAV_BDF2(derta_t(6));

u_r=SAV_BDF1(0.1*2^-6*(2/3)^8);
error=ones(1,tn);
error_ref=norm(u_r(:,end));
nn=length(u_r(:,end));
error(1)=norm(funfit(u_1(:,end),nn)-u_r(:,end))/error_ref;
error(2)=norm(funfit(u_2(:,end),nn)-u_r(:,end))/error_ref;
error(3)=norm(funfit(u_3(:,end),nn)-u_r(:,end))/error_ref;
error(4)=norm(funfit(u_4(:,end),nn)-u_r(:,end))/error_ref;
error(5)=norm(funfit(u_5(:,end),nn)-u_r(:,end))/error_ref;
% error(6)=norm(funfit(u_6(:,end),nn)-u_r(:,end))/error_ref;

p=polyfit(log(derta_t),log(error),1);
x1=linspace(1.2*derta_t(1),0.8*derta_t(end));
y1=exp(polyval(p,log(x1)));
loglog(derta_t,error,'*',x1,y1);
xlabel('log(¡÷t)');
ylabel('log(err)');
title('SAV/BDF1');
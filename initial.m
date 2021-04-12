function u0 = initial(h)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
global X 
N_x=round(X/h);
xx=linspace(0,X,N_x+1);
u0=0.6*sin(xx)+0.05*2*(rand(1,N_x+1)-0.5);
u0=u0';
end


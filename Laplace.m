function [Lap_1,Lap_2] = Laplace(h)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
global X
N_x=round(X/h);
v=ones(1,N_x);
Lap_1=-2*eye(N_x+1) ;
B=diag(v,-1);
Lap_1=Lap_1+B+B';
Lap_2=Lap_1^2;Lap_2(1,1)=6;
Lap_2(3,3)=6-4.5;Lap_2(3,4)=-4+1.5;
Lap_2(4,3)=-4+1.5;Lap_2(4,4)=6-0.5;
Lap_2(N_x-1,N_x-1)=6-4.5;Lap_2(N_x-1,N_x-2)=-4+1.5;
Lap_2(N_x-2,N_x-1)=-4+1.5;Lap_2(N_x-2,N_x-2)=6-0.5;
Lap_1=Lap_1/h^2;
Lap_2=Lap_2/h^4;

end


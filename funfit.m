function u = funfit(u0,n)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
n0=length(u0);
u=interp1(linspace(0,1,n0),u0,linspace(0,1,n));
u=u';
end


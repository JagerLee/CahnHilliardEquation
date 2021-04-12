function f = F(x)
%UNTITLED2 此处显示有关此函数的摘要
%   F(x)=1/4*(x^2-1)^2
global derta beta
f=(x.^2-1-beta).^2/(4*derta);
end


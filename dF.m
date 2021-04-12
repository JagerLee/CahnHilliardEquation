function df = dF(x)
%UNTITLED4 此处显示有关此函数的摘要
global derta beta
df=(x.^3-(1+beta)*x)/derta;

end


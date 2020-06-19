function error=fun(x,s1,lenchrom)
%与标准s1之间的MSE均方误差作为判断依据
N=length(lenchrom)
mse=0
for i=1:N
    mse=mse+(x(i)-s1(i))^2;
end
error=mse/N;
end


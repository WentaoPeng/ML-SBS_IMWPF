function error=fun(x,s1,lenchrom)
%���׼s1֮���MSE���������Ϊ�ж�����
N=length(lenchrom)
mse=0
for i=1:N
    mse=mse+(x(i)-s1(i))^2;
end
error=mse/N;
end


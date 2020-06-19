function ret=Code(lenchrom,bound)
%CODE 生成每个个体N个数据，模仿系统传播，加入噪声，产生初代群体
%全部返回到频域比较
global s1

S=s1;

flag=0
while flag==0
    N=length(lenchrom);
    for i=1:N
        S(i)=S(i)*exp(1i*rand()*2*pi);
    end
    iftx=ifftshift(S);
    iftxs=ifft(iftx);
    T=awgn(iftxs,30,'measured');
    R=abs(fftshift(fft(T)));
    for i=1:N
        ret(i)=R(i);
    end
    flag=test(lenchrom,bound,ret);
end
function ret=Code(lenchrom,bound)
%CODE ����ÿ������N�����ݣ�ģ��ϵͳ������������������������Ⱥ��
%ȫ�����ص�Ƶ��Ƚ�
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
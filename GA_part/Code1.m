function ret=Code(lenchrom,bound)
%添加优良基因仿真
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
    R=abs(fftshift(fft(iftxs)));
    pick=randperm(N,100);
    for i=1:N
%         pick=round(rand(1,100)*N);      
        for j=1:100 
             f=0
%             while pick(j)==0
%             pick(j)=round(rand(1,1)*N);
%             end
            if i==pick(j)
                ret(i)=R(i)+rand(1,1)*0.3
                f=1
                break;
            end
            
%             if f==1
%                 break;
%             end
        end
        
        if f==1
            i=i+1;
        else
            ret(i)=R(i);
        end
    end
    flag=test(lenchrom,bound,ret);
end


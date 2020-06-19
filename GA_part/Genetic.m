% ��ջ�������
clc
clear
% 
global s1 N V

%% ���ɱ�׼��������%%��������Ƶ�����岨��
V=2.^13;     %���ݷ�Χ
N=2.^14;%��������,��������
df=1.526e5;     %�ֱ���0.1MHz
f=(-N/2:1:N/2-1)*df;%Ƶ�ʷ�Χ
Fs=mean(diff(f));       %Sampling Frequency
Fn=Fs/2;        %Nyuquist Frequency
t=(0:0.00004/Fs:(0.00004/Fs)*(N-1));
w=1e9;%����1GHz

s1=rectpuls(f,w);
%s1=tripuls(f,w);

%���峤���뷶Χ
lenchrom=ones(1,N);         %����1*N��ȫ1����
bound=[-2*ones(N,1) 3*ones(N,1)];        %���ݷ�Χ
%% �Ŵ��㷨������ʼ��
maxgen=600;             %��������������������
sizepop=500;             %��Ⱥ��ģ

%����Ӧ���������ʶ��峣��,��ʼ������ֵ
pcross=[0.85];           %������ʣ�0~1֮��
pmutation=[0.05];        %������ʣ�0~1֮��

%% ��Ⱥ��ʼ�������ȵ����ɸ��岢�Ҽ����������ģ�»���
individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]);  %����Ⱥ��Ϣ����Ϊһ���ṹ��
avgfitness=[];                      %ÿһ����Ⱥ��ƽ����Ӧ��
bestfitness=[];                     %ÿһ����Ⱥ�������Ӧ��
bestchrom=[];                       %��Ӧ����õ�Ⱦɫ��
%��ʼ����Ⱥ
for i=1:sizepop
    %�����������������ʼ��Ⱥ
    if i>10
        individuals.chrom(i,:)=Code(lenchrom,bound);
        x=individuals.chrom(i,:);
    %������Ӧ��
        individuals.fitness(i)=fun(x,s1,lenchrom);        %������Ӧ��
    else
        individuals.chrom(i,:)=Code1(lenchrom,bound);
        x=individuals.chrom(i,:);
        individuals.fitness(i)=fun(x,s1,lenchrom);
    end
end

figure(4)
subplot(3,1,1)
plot(f,individuals.chrom);
title(['��ʼ��Ⱥ  ' '����=' num2str(sizepop)]);
xlabel('Ƶ��/Hz');ylabel('dBm');
hold on

[bestfitness bestindex]=min(individuals.fitness);
[worestfitness,worestindex]=max(individuals.fitness);
bestchrom=individuals.chrom(bestindex,:);  %��õ�Ⱦɫ��
avgfitness=sum(individuals.fitness)/sizepop; %Ⱦɫ���ƽ����Ӧ��
% ��¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
trace=[avgfitness bestfitness]; 


%% ��������ѳ�ʼ��ֵ��Ȩֵ
%��ʼ����

for j=1:maxgen
    j
    %��̬�������������ʣ����ÿһ������
    Mcmax=0.85;
    Mcmin=0.35;
    Mmmax=0.05;
    Mmmin=0.01;
    if avgfitness==worestfitness
        pcross=Mcmax-Mcmin;
        pmutation=Mmmax-Mmmin;
    else
        pcross=Mcmax-Mcmin*((worestfitness-avgfitness)/(worestfitness-bestfitness));
        pmutation=Mmmax-Mmmin*((worestfitness-avgfitness)/(worestfitness-bestfitness));
    end
    
    %ѡ�񣬶�����
    individuals=Select(individuals,sizepop);
    avgfitness=sum(individuals.fitness)/sizepop;
    
    %����
    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound,bestfitness,avgfitness,individuals.fitness);
   
    %����
    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,j,maxgen,bound,bestfitness,avgfitness,individuals.fitness)
    
    %����õ���������Ӧ��
    for k=1:sizepop
        x=individuals.chrom(k,:);
        individuals.fitness(k)=fun(x,s1,lenchrom);
    end
    
      %�ҵ���С�������Ӧ�ȵ�Ⱦɫ�弰��������Ⱥ�е�λ��
    [newbestfitness,newbestindex]=min(individuals.fitness);
    [worestfitness,worestindex]=max(individuals.fitness);
    % ������һ�ν�������õ�Ⱦɫ��
    if bestfitness>newbestfitness
        bestfitness=newbestfitness;
        bestchrom=individuals.chrom(newbestindex,:);
    end
    
    individuals.chrom(worestindex,:)=bestchrom;
    individuals.fitness(worestindex)=bestfitness;
    
    avgfitness=sum(individuals.fitness)/sizepop;
    
%     figure(1);
%     hold on;
%     plot(j,avgfitness);
%     plot(j,bestfitness);
    
    
    trace=[trace;avgfitness bestfitness]; %��¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
    
    
    
    
    
end


%% �Ŵ��㷨������� 
 figure(3)
[r c]=size(trace);
subplot(2,1,1)
plot([1:r]',trace(:,1),'b--');
title(['��Ӧ������  ' '��ֹ������' num2str(maxgen)]);
xlabel('��������');ylabel('��Ӧ��');
legend('ƽ����Ӧ��');
disp('��Ӧ��                   ����');
subplot(2,1,2)
plot([1:r]',trace(:,2),'b--');
title(['��Ӧ������  ' '��ֹ������' num2str(maxgen)]);
xlabel('��������');ylabel('��Ӧ��');
legend('�����Ӧ��');
disp('��Ӧ��                   ����');

x=bestchrom;
figure(4)
subplot(3,1,3)
plot(f,x);
title(['���Ÿ��� ']);
xlabel('Ƶ��/Hz');ylabel('dBm');
subplot(3,1,2)
plot(f,individuals.chrom);
title(['������Ⱥ  ' '����=' num2str(sizepop)]);
xlabel('Ƶ��/Hz');ylabel('dBm');




%% �ֲ��Ż��㷨��������������Ҫ֪����Ƶ��ĵ�ѹֵ�����Ҫ֪��ǰ�������⹦��ֵ�����ݱ�������������
%for k=1:sizepop
    
% 清空环境变量
clc
clear
% 
global s1 N V

%% 生成标准波形数据%%生产矩形频域脉冲波形
V=2.^13;     %数据范围
N=2.^14;%采样点数,总数剧量
df=1.526e5;     %分辨率0.1MHz
f=(-N/2:1:N/2-1)*df;%频率范围
Fs=mean(diff(f));       %Sampling Frequency
Fn=Fs/2;        %Nyuquist Frequency
t=(0:0.00004/Fs:(0.00004/Fs)*(N-1));
w=1e9;%带宽1GHz

s1=rectpuls(f,w);
%s1=tripuls(f,w);

%定义长度与范围
lenchrom=ones(1,N);         %生成1*N的全1矩阵
bound=[-2*ones(N,1) 3*ones(N,1)];        %数据范围
%% 遗传算法参数初始化
maxgen=600;             %进化代数，即迭代次数
sizepop=500;             %种群规模

%自适应交叉变异概率定义常数,初始化概率值
pcross=[0.85];           %交叉概率，0~1之间
pmutation=[0.05];        %变异概率，0~1之间

%% 种群初始化，首先得生成个体并且加上随机噪声模仿环境
individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]);  %将种群信息定义为一个结构体
avgfitness=[];                      %每一代种群的平均适应度
bestfitness=[];                     %每一代种群的最佳适应度
bestchrom=[];                       %适应度最好的染色体
%初始化种群
for i=1:sizepop
    %利用随机噪声产生初始种群
    if i>10
        individuals.chrom(i,:)=Code(lenchrom,bound);
        x=individuals.chrom(i,:);
    %计算适应度
        individuals.fitness(i)=fun(x,s1,lenchrom);        %计算适应度
    else
        individuals.chrom(i,:)=Code1(lenchrom,bound);
        x=individuals.chrom(i,:);
        individuals.fitness(i)=fun(x,s1,lenchrom);
    end
end

figure(4)
subplot(3,1,1)
plot(f,individuals.chrom);
title(['初始种群  ' '数量=' num2str(sizepop)]);
xlabel('频率/Hz');ylabel('dBm');
hold on

[bestfitness bestindex]=min(individuals.fitness);
[worestfitness,worestindex]=max(individuals.fitness);
bestchrom=individuals.chrom(bestindex,:);  %最好的染色体
avgfitness=sum(individuals.fitness)/sizepop; %染色体的平均适应度
% 记录每一代进化中最好的适应度和平均适应度
trace=[avgfitness bestfitness]; 


%% 求解最初最佳初始阈值和权值
%开始进化

for j=1:maxgen
    j
    %动态调整交叉变异概率，针对每一代调整
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
    
    %选择，赌轮盘
    individuals=Select(individuals,sizepop);
    avgfitness=sum(individuals.fitness)/sizepop;
    
    %交叉
    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound,bestfitness,avgfitness,individuals.fitness);
   
    %变异
    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,j,maxgen,bound,bestfitness,avgfitness,individuals.fitness)
    
    %计算该迭代代数适应度
    for k=1:sizepop
        x=individuals.chrom(k,:);
        individuals.fitness(k)=fun(x,s1,lenchrom);
    end
    
      %找到最小和最大适应度的染色体及它们在种群中的位置
    [newbestfitness,newbestindex]=min(individuals.fitness);
    [worestfitness,worestindex]=max(individuals.fitness);
    % 代替上一次进化中最好的染色体
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
    
    
    trace=[trace;avgfitness bestfitness]; %记录每一代进化中最好的适应度和平均适应度
    
    
    
    
    
end


%% 遗传算法结果分析 
 figure(3)
[r c]=size(trace);
subplot(2,1,1)
plot([1:r]',trace(:,1),'b--');
title(['适应度曲线  ' '终止代数＝' num2str(maxgen)]);
xlabel('进化代数');ylabel('适应度');
legend('平均适应度');
disp('适应度                   变量');
subplot(2,1,2)
plot([1:r]',trace(:,2),'b--');
title(['适应度曲线  ' '终止代数＝' num2str(maxgen)]);
xlabel('进化代数');ylabel('适应度');
legend('最佳适应度');
disp('适应度                   变量');

x=bestchrom;
figure(4)
subplot(3,1,3)
plot(f,x);
title(['最优个体 ']);
xlabel('频率/Hz');ylabel('dBm');
subplot(3,1,2)
plot(f,individuals.chrom);
title(['最终种群  ' '数量=' num2str(sizepop)]);
xlabel('频率/Hz');ylabel('dBm');




%% 局部优化算法————首先需要知道电频梳的电压值，其次要知道前后的增益光功率值，根据比例，反馈设置
%for k=1:sizepop
    
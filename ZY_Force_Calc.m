


%% 各尺寸参数初始化

Js=889;     %永磁体矫顽力
R1=1;       %动子内环内周半径
R2=3.5;     %定子内周半径
R3=6.5;    %动子外环内周半径
ws=1.5;     %定子宽度
w1=1.5;     %动子内环宽度
w2=2;     %动子外环宽度
h1=7.5;     %1、2号动子高度
h2=5;       %3、4号动子高度
h3=7.5;     %5、6号动子高度

hm=6.5;    %7号定子上表面距离x轴；8号定子下表面距离x轴
n=20;       %zss的步数
zss=linspace(-1.5,1.5,n);      %动子行程
hs=1.4;       %定子高度



%% 为计算结果Force预先分配空间
Force=zeros(1,n);

%% 计算7号定子和8号定子的受力 

%function Br=ZY_Br_Calc(Js,R1,R3,w1,w2,h1,h2,h3,hm,zss,hs)
%Br=[B1256_7u_8l,B1256_7l_8u,Bradial_7u_8l,Bradial_7l_8u]

syms r 
for j=1:length(zss)
    Br=ZY_Br_Calc(Js,R1,R3,w1,w2,h1,h2,h3,hm,zss(j),hs); %hm，hs为正时计算7号

    F7_upper_1256=str2func(['@(thiv0,r)',vectorize(Br(1)*(-Js)*r)]);
    F7_lower_1256=str2func(['@(thiv0,r)',vectorize(Br(2)*Js*r)]);

    F7_upper_radial=str2func(['@(thiv0,R,r)',vectorize(Br(3)*(-Js)*r)]);
    F7_lower_radial=str2func(['@(thiv0,R,r)',vectorize(Br(4)*Js*r)]); %符号函数转化函数句柄

    F7=2*pi*integral2(F7_upper_1256,0,2*pi,R2,R2+ws)...
        +2*pi*integral2(F7_lower_1256,0,2*pi,R2,R2+ws)...%这两行计算轴向充磁的磁环产生的力
        +2*pi*integral3(F7_upper_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...
        +2*pi*integral3(F7_lower_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...%这两行计算3号(+r)磁环产生的力
        +2*pi*integral3(F7_upper_radial,0,2*pi,R3,R3+w2,R2,R2+ws)...
        +2*pi*integral3(F7_lower_radial,0,2*pi,R3,R3+w2,R2,R2+ws); %这两行计算4号(+r)磁环产生的力

    
    Br=ZY_Br_Calc(Js,R1,R3,w1,w2,h1,h2,h3,-hm,zss(j),-hs);%hm，hs为负时计算8号

    F8_lower_1256=str2func(['@(thiv0,r)',vectorize(Br(1)*(-Js)*r)]);
    F8_upper_1256=str2func(['@(thiv0,r)',vectorize(Br(2)*Js*r)]);

    F8_lower_radial=str2func(['@(thiv0,R,r)',vectorize(Br(3)*(-Js)*r)]);
    F8_upper_radial=str2func(['@(thiv0,R,r)',vectorize(Br(4)*Js*r)]); %符号函数转化函数句柄

    F8=2*pi*integral2(F8_lower_1256,0,2*pi,R2,R2+ws)...
        +2*pi*integral2(F8_upper_1256,0,2*pi,R2,R2+ws)...%这两行计算轴向充磁的磁环产生的力
        +2*pi*integral3(F8_lower_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...
        +2*pi*integral3(F8_upper_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...%这两行计算3号(+r)磁环产生的力
        +2*pi*integral3(F8_lower_radial,0,2*pi,R3,R3+w2,R2,R2+ws)...
        +2*pi*integral3(F8_upper_radial,0,2*pi,R3,R3+w2,R2,R2+ws); %这两行计算4号(+r)磁环产生的力
    
    
    Force(j)=F7+F8;
end

%% 绘制图像
plot(zss,abs(Force),'-*');

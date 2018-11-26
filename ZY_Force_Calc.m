


%% ���ߴ������ʼ��

Js=889;     %�����������
R1=1;       %�����ڻ����ܰ뾶
R2=3.5;     %�������ܰ뾶
R3=6.5;    %�����⻷���ܰ뾶
ws=1.5;     %���ӿ��
w1=1.5;     %�����ڻ����
w2=2;     %�����⻷���
h1=7.5;     %1��2�Ŷ��Ӹ߶�
h2=5;       %3��4�Ŷ��Ӹ߶�
h3=7.5;     %5��6�Ŷ��Ӹ߶�

hm=6.5;    %7�Ŷ����ϱ������x�᣻8�Ŷ����±������x��
n=20;       %zss�Ĳ���
zss=linspace(-1.5,1.5,n);      %�����г�
hs=1.4;       %���Ӹ߶�



%% Ϊ������ForceԤ�ȷ���ռ�
Force=zeros(1,n);

%% ����7�Ŷ��Ӻ�8�Ŷ��ӵ����� 

%function Br=ZY_Br_Calc(Js,R1,R3,w1,w2,h1,h2,h3,hm,zss,hs)
%Br=[B1256_7u_8l,B1256_7l_8u,Bradial_7u_8l,Bradial_7l_8u]

syms r 
for j=1:length(zss)
    Br=ZY_Br_Calc(Js,R1,R3,w1,w2,h1,h2,h3,hm,zss(j),hs); %hm��hsΪ��ʱ����7��

    F7_upper_1256=str2func(['@(thiv0,r)',vectorize(Br(1)*(-Js)*r)]);
    F7_lower_1256=str2func(['@(thiv0,r)',vectorize(Br(2)*Js*r)]);

    F7_upper_radial=str2func(['@(thiv0,R,r)',vectorize(Br(3)*(-Js)*r)]);
    F7_lower_radial=str2func(['@(thiv0,R,r)',vectorize(Br(4)*Js*r)]); %���ź���ת���������

    F7=2*pi*integral2(F7_upper_1256,0,2*pi,R2,R2+ws)...
        +2*pi*integral2(F7_lower_1256,0,2*pi,R2,R2+ws)...%�����м��������ŵĴŻ���������
        +2*pi*integral3(F7_upper_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...
        +2*pi*integral3(F7_lower_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...%�����м���3��(+r)�Ż���������
        +2*pi*integral3(F7_upper_radial,0,2*pi,R3,R3+w2,R2,R2+ws)...
        +2*pi*integral3(F7_lower_radial,0,2*pi,R3,R3+w2,R2,R2+ws); %�����м���4��(+r)�Ż���������

    
    Br=ZY_Br_Calc(Js,R1,R3,w1,w2,h1,h2,h3,-hm,zss(j),-hs);%hm��hsΪ��ʱ����8��

    F8_lower_1256=str2func(['@(thiv0,r)',vectorize(Br(1)*(-Js)*r)]);
    F8_upper_1256=str2func(['@(thiv0,r)',vectorize(Br(2)*Js*r)]);

    F8_lower_radial=str2func(['@(thiv0,R,r)',vectorize(Br(3)*(-Js)*r)]);
    F8_upper_radial=str2func(['@(thiv0,R,r)',vectorize(Br(4)*Js*r)]); %���ź���ת���������

    F8=2*pi*integral2(F8_lower_1256,0,2*pi,R2,R2+ws)...
        +2*pi*integral2(F8_upper_1256,0,2*pi,R2,R2+ws)...%�����м��������ŵĴŻ���������
        +2*pi*integral3(F8_lower_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...
        +2*pi*integral3(F8_upper_radial,0,2*pi,R1,R1+w1,R2,R2+ws)...%�����м���3��(+r)�Ż���������
        +2*pi*integral3(F8_lower_radial,0,2*pi,R3,R3+w2,R2,R2+ws)...
        +2*pi*integral3(F8_upper_radial,0,2*pi,R3,R3+w2,R2,R2+ws); %�����м���4��(+r)�Ż���������
    
    
    Force(j)=F7+F8;
end

%% ����ͼ��
plot(zss,abs(Force),'-*');

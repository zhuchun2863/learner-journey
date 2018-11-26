function Br=ZY_Br_Calc(Js,R1,R3,w1,w2,h1,h2,h3,hm,zss,hs)  
%% 文件说明
%这个函数文件只是用来运算，参与ZY_Force_Calc.m文件中受力积分式，的积分函数，并不是计算磁感应强度
%因为暂时没有发现直接数值计算Br的方法，只能直接数值计算受力。但是本文件仍以Br命名

%% 真空磁导率
u0=4*pi*10^-7; 

%% 辅助函数构造
syms r z R thiv0 z0
gz=cos(thiv0)*R/sqrt(r^2+R^2-2*r*R*cos(thiv0)+(z-z0)^2);            %轴向充磁计算的辅助函数
gr=(z-z0)*cos(thiv0)*R/((r^2+R^2-2*r*R*cos(thiv0)+(z-z0)^2)^(3/2)); %径向充磁计算的辅助函数


%% -----------------轴向充磁磁环在7号定子上表面（或者8号定子下表面）--------------------------
B1256_7u_8l=u0*Js/4/pi*(-subs(gz,[R,z,z0],[R1,hm+zss,h1+h2/2])...               
                   +subs(gz,[R,z,z0],[R1,hm+zss,h2/2])...       %1号动子(-z)内周半径的面电流  
                   +subs(gz,[R,z,z0],[R1+w1,hm+zss,h1+h2/2])...                             
                   -subs(gz,[R,z,z0],[R1+w1,hm+zss,h2/2])...    %1号动子(-z)外周半径的面电流  
                   +subs(gz,[R,z,z0],[R3,hm+zss,h1+h2/2])...
                   -subs(gz,[R,z,z0],[R3,hm+zss,h2/2])...       %2号动子(+z)内周半径的面电流
                   -subs(gz,[R,z,z0],[R3+w2,hm+zss,h1+h2/2])...
                   +subs(gz,[R,z,z0],[R3+w2,hm+zss,h2/2])...    %2号动子(+z)外周半径的面电流
               +subs(gz,[R,z,z0],[R1,hm+zss,-h2/2])...               
                   -subs(gz,[R,z,z0],[R1,hm+zss,-h2/2-h3])...   %5号动子(+z)内周半径的面电流  
                   -subs(gz,[R,z,z0],[R1+w1,hm+zss,-h2/2])...                             
                   +subs(gz,[R,z,z0],[R1+w1,hm+zss,-h2/2-h3])...%5号动子(+z)外周半径的面电流  
                   -subs(gz,[R,z,z0],[R3,hm+zss,-h2/2])...
                   +subs(gz,[R,z,z0],[R3,hm+zss,-h2/2-h3])...   %6号动子(-z)内周半径的面电流
                   +subs(gz,[R,z,z0],[R3+w2,hm+zss,-h2/2])...
                   -subs(gz,[R,z,z0],[R3+w2,hm+zss,-h2/2-h3])); %6号动子(-z)外周半径的面电流
               
               
%% -----------------轴向充磁磁环在7号定子下表面（或者8号定子上表面）-------------------------               
B1256_7l_8u=u0*Js/4/pi*(-subs(gz,[R,z,z0],[R1,hm-hs+zss,h1+h2/2])...              
                   +subs(gz,[R,z,z0],[R1,hm-hs+zss,h2/2])...    %1号动子(-z)内周半径的面电流  
                   +subs(gz,[R,z,z0],[R1+w1,hm-hs+zss,h1+h2/2])...                             
                   -subs(gz,[R,z,z0],[R1+w1,hm-hs+zss,h2/2])... %1号动子(-z)外周半径的面电流  
                   +subs(gz,[R,z,z0],[R3,hm-hs+zss,h1+h2/2])...
                   -subs(gz,[R,z,z0],[R3,hm-hs+zss,h2/2])...    %2号动子(+z)内周半径的面电流
                   -subs(gz,[R,z,z0],[R3+w2,hm-hs+zss,h1+h2/2])...
                   +subs(gz,[R,z,z0],[R3+w2,hm-hs+zss,h2/2])... %2号动子(+z)外周半径的面电流
               +subs(gz,[R,z,z0],[R1,hm-hs+zss,-h2/2])...               
                   -subs(gz,[R,z,z0],[R1,hm-hs+zss,-h2/2-h3])...%5号动子(+z)内周半径的面电流  
                   -subs(gz,[R,z,z0],[R1+w1,hm-hs+zss,-h2/2])...                             
                   +subs(gz,[R,z,z0],[R1+w1,hm-hs+zss,-h2/2-h3])... %5号动子(+z)外周半径的面电流  
                   -subs(gz,[R,z,z0],[R3,hm-hs+zss,-h2/2])...
                   +subs(gz,[R,z,z0],[R3,hm-hs+zss,-h2/2-h3])...    %6号动子(-z)内周半径的面电流
                   +subs(gz,[R,z,z0],[R3+w2,hm-hs+zss,-h2/2])...
                   -subs(gz,[R,z,z0],[R3+w2,hm-hs+zss,-h2/2-h3]));  %6号动子(-z)外周半径的面电流             

%% -----------------径向充磁磁环在7号定子上表面（或者8号定子下表面）-----------
%gr=(z-z0)*cos(thiv0)/((r^2+R^2-2*r*R*cos(thiv0)+(z-z0)^2)^(3/2));%径向充磁计算的辅助函数

Bradial_7u_8l=u0*Js/4/pi*(+subs(gr,[z,z0],[hm+zss,-h2/2])...               
                   -subs(gr,[z,z0],[hm+zss,h2/2]));   %动子(+r)的面电流    


%% -----------------径向充磁磁环在7号定子下表面（或者8号定子上表面）-----------
%gr=(z-z0)*cos(thiv0)/((r^2+R^2-2*r*R*cos(thiv0)+(z-z0)^2)^(3/2));%径向充磁计算的辅助函数

Bradial_7l_8u=u0*Js/4/pi*(+subs(gr,[z,z0],[hm-hs+zss,-h2/2])...              
                   -subs(gr,[z,z0],[hm-hs+zss,h2/2]));   %动子(+r)的面电流    

%% 将符号函数建立符号矩阵输出
Br=[B1256_7u_8l,B1256_7l_8u,Bradial_7u_8l,Bradial_7l_8u];

end
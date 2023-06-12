% Standard script file for input variables of CT intepretation system.
% Please change the values of these parameters freely for specific problem.
% meshInterval:The interval of ONE cell (= 0.5 or 1.0,in general);it will
%              be constant for the whole intepretation process.
% nodeNum:The total number of nodes on ONE cell (multiplier of 4,in general)
% Clb:Left-bottom (x,y)
% Clt:Left-top (x,y)
% Crb:Right-bottom (x,y)
% Crt:Right-top (x,y)
% draw:Control the outputs of figures(后期用数组控制的更细一点)
%      =1 
%      =2
%      =3
% s_x,s_y:only provide s_y. s_x can be generated automatically.
% t_x,t_y:only provide t_y. t_x can be generated automatically.
%   Sources and receivers arrangement rule:
%CANNOT be placed at angle points of cells simultaneously!

% -----By Chenglong Duan,Nanjing University,2015.-----

clear;
meshInterval=0.5; 
nodeNum=16;     % nodeNum和源点/接收点有关，若接收点在边界中点，则不能采用12节点，这会导致搜索不到接收点
Clb=[0 -12]; 
Clt=[0 0]; 
Crb=[7 -12]; 
Crt=[7 0];  
% 震源点坐标定义
s_y=-92.5:0.5:-80.5;  %对跨孔来说，y比x重要
s_x=ones(1,length(s_y))*Clb(1);
% 接收点坐标定义
t_y=-92.25:0.5:-80.75;  %对跨孔来说，y比x重要
t_x=ones(1,length(t_y))*Crb(1);
T_file='D:\GW39-37(time).xlsx'; % T文件的路径(注意：Excel 2010版本的后缀是xlsx)
                      % T文件编制的规定：第1炮点为第1行，依次放置1#...12#检波器接收的时间；
                      %                 第2炮点为第2行，依次放置1#...12#检波器接收的时间；依次类推。

file_XYZ = 'GW39-37.xlsx';  % file_XYZ为最终保存的Excel文件名，用法如'GB3-2.xls'(注意：Excel 2010版本的后缀是xlsx)                      
% Controlling handle
initial_fig=0; % Plot initial velocity figure
draw=0; % 控制射线追踪过程中的一系列图，一般都设为0不画
damp_fig=1; % 画damp_find里的图

% 合并fbk文件






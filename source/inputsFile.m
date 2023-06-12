% Standard script file for input variables of CT intepretation system.
% Please change the values of these parameters freely for specific problem.
% meshInterval:The interval of ONE cell (= 0.5 or 1.0,in general);it will
%              be constant for the whole intepretation process.
% nodeNum:The total number of nodes on ONE cell (multiplier of 4,in general)
% Clb:Left-bottom (x,y)
% Clt:Left-top (x,y)
% Crb:Right-bottom (x,y)
% Crt:Right-top (x,y)
% draw:Control the outputs of figures(������������Ƶĸ�ϸһ��)
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
nodeNum=16;     % nodeNum��Դ��/���յ��йأ������յ��ڱ߽��е㣬���ܲ���12�ڵ㣬��ᵼ�������������յ�
Clb=[0 -12]; 
Clt=[0 0]; 
Crb=[7 -12]; 
Crt=[7 0];  
% ��Դ�����궨��
s_y=-92.5:0.5:-80.5;  %�Կ����˵��y��x��Ҫ
s_x=ones(1,length(s_y))*Clb(1);
% ���յ����궨��
t_y=-92.25:0.5:-80.75;  %�Կ����˵��y��x��Ҫ
t_x=ones(1,length(t_y))*Crb(1);
T_file='D:\GW39-37(time).xlsx'; % T�ļ���·��(ע�⣺Excel 2010�汾�ĺ�׺��xlsx)
                      % T�ļ����ƵĹ涨����1�ڵ�Ϊ��1�У����η���1#...12#�첨�����յ�ʱ�䣻
                      %                 ��2�ڵ�Ϊ��2�У����η���1#...12#�첨�����յ�ʱ�䣻�������ơ�

file_XYZ = 'GW39-37.xlsx';  % file_XYZΪ���ձ����Excel�ļ������÷���'GB3-2.xls'(ע�⣺Excel 2010�汾�ĺ�׺��xlsx)                      
% Controlling handle
initial_fig=0; % Plot initial velocity figure
draw=0; % ��������׷�ٹ����е�һϵ��ͼ��һ�㶼��Ϊ0����
damp_fig=1; % ��damp_find���ͼ

% �ϲ�fbk�ļ�






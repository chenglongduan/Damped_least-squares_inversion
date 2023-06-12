function  [t_print, L]=rayTracing( meshInterval, nodeNum, Clb, Clt, Crb, Crt, draw, s_x, s_y, t_x, t_y, velocity, m, n )
%UNTITLED2 Summary of this function goes here
%   速度放在网格
% source的x坐标放在向量s_x中，y坐标放在s_y中
% terminus的x坐标放在向量t_x中，y坐标放在t_y中
% Define cell-velocity
%
% m,n来自initialField.m，因为假定网格划分自始至终都用一套。
% 其中：m--L的行数，等于射线总数；n---L的列数，等于未知数总数。

%-----By Chenglong Duan,Nanjing University,2015.-----

xLength=abs(Clb(1)-Crb(1));
yLength=abs(max(Clt(2),Crt(2))-min(Clb(2),Crb(2)));
if mod(xLength,meshInterval)~=0 || mod(yLength,meshInterval)~=0
    error('mesh error!');
end

%-------------预留一个输入初始速度文件的模块---------------------------------
%velocity=input('Please input velocity file:');
%--------------------------------------------------------------------------

%% Velocity map (尽量不画)
% 当速度值放在中心时，四周就出现了meshInterval/2的空缺，需要填补（左右没问题，但是上下在不规则情况下要重新考虑！！！）
%if Clb(2)~=Crb(2) || Clt(2)~=Crt(2)
if draw==1
    figure;
    imagesc(velocity);
    axis equal;
    set(gca,'YDir','normal');
    %-------------------------------虚线框中的部分第一次画速度图时先不运行，等画完后根据自然坐标，填入xlim,ylim,xtick,ytick再运行
    xlim([0.5 3.5]);
    ylim([0.5 2.5]);
    set(gca,'xtick',0.5:3.5,'xticklabel',Clb(1):meshInterval:Crb(1)); 
    set(gca,'ytick',0.5:2.5,'yticklabel',min(Clb(2),Crb(2)):meshInterval:max(Clt(2),Crt(2)));
    %-----------------------------------------------------------------------------------------------------------------------
    colorbar;
    xlabel('Horizontal distance(m)');ylabel('Elevation(m)');
end

%% 计算关联矩阵
if draw==1
    figure;
end
[s_row,s_col,s_val,node,coordcell,N1,N2,N3] = relateMatrix( meshInterval, nodeNum, Clb, Clt, Crb, Crt, draw, velocity );
fprintf('Iter %d:The related matrix has created.\n',1); % 进度语句1

%% 计算SRP
%------先利用输入的源、检位置找出对应的节点------
source=zeros(1,length(s_x));terminus=zeros(1,length(t_x));
if any(s_x~=Clb(1)) || any(t_x~=Crb(1))
    error('Please Check the y-boundary!');
end
count_source=1;count_terminus=1;
for i=1:(N1*N2+N2+1)
    if any(coordcell{i}(2)==s_y)
        source(count_source)=i;
        count_source=count_source+1;
    end
end
for i=(N1*N2+N2+1)*(N3+1)+(N2+1)*N1*N3-N2*(N1+1):(N1*N2+N2+1)*(N3+1)+(N2+1)*N1*N3
    if any(coordcell{i}(2)==t_y)
        terminus(count_terminus)=i;
        count_terminus=count_terminus+1;
    end
end
clear count_source count_terminus; % count_source,count_terminus是计数标记
fprintf('Iter %d:The source and terminus points have been found.\n',1); % 进度语句2
%把coordcell保存到默认文件夹下，需要时再调出，以减小空间
% save coordcell;
% clear coordcell;
%---------------------------------------------
L_source=length(source);L_terminus=length(terminus);
t_print=zeros(L_source,L_terminus);
pos=1;  % 构造稀疏矩阵元素的位置索引
figure;
for i=1:L_source
    tmin=zeros(1,L_terminus);
    for j=1:L_terminus
        xPath=[];yPath=[];
        [t,path]=dijkstra(s_row,s_col,s_val,node,source(i),terminus(j));
        fprintf('Finish source %d,terminus %d. Drawing now...\n',i,j); % 进度语句3
        tmin(j)=t;
        %------------------------------------------------------
        for k=1:length(path)
            xPath=[xPath coordcell{path(k)}(1)]; 
            yPath=[yPath coordcell{path(k)}(2)];
        end
        for k=1:length(xPath)-1
            judge1=yPath(k)==yPath(k+1) && yPath(k)==min(Clb(2),Crb(2));
            judge2=xPath(k)==xPath(k+1) && xPath(k)==Crb(1);
            judge3=xPath(k)==xPath(k+1) && xPath(k)~=Clb(1) && xPath(k)~=Crb(1) && abs(yPath(k)-yPath(k+1))~=meshInterval;
            judge4=yPath(k)==yPath(k+1) && yPath(k)~=min(Clb(2),Crb(2)) && yPath(k)~=max(Clt(2),Crt(2)) && abs(xPath(k)-xPath(k+1))~=meshInterval;
            if judge1==1  %两节点连线是下边界
                condition=1;
            end
            if judge2==1  %两节点连线是右边界
                condition=2;
            end
            if judge3==1  %竖直连线
                condition=3;
            end
            if judge4==1  %水平连线
                condition=4;
            end
            if judge1==0 && judge2==0 && judge3==0 && judge4==0
                condition=0;
            end
            switch condition
                case 1
                    meshNO_c=floor((xPath(k)-Clb(1))/meshInterval)*N2;
                    meshNO_f=1;
                    L_col(pos)=meshNO_c+meshNO_f;
                case 2
                    meshNO_c=((xPath(k)-Clb(1))/meshInterval-1)*N2;
                    meshNO_f1=ceil((yPath(k)-min(Clb(2),Crb(2)))/meshInterval);
                    meshNO_f2=ceil((yPath(k+1)-min(Clb(2),Crb(2)))/meshInterval);
                    meshNO_f=max(meshNO_f1,meshNO_f2);
                    L_col(pos)=meshNO_c+meshNO_f;
                case 3
                    meshNO_c=floor((xPath(k)-Clb(1))/meshInterval-1)*N2;
                    meshNO_f1=ceil((yPath(k)-min(Clb(2),Crb(2)))/meshInterval);
                    meshNO_f2=ceil((yPath(k+1)-min(Clb(2),Crb(2)))/meshInterval);
                    meshNO_f=max(meshNO_f1,meshNO_f2);
                    pool=[meshNO_c+meshNO_f,meshNO_c+meshNO_f+N2];
                    [~,nn]=max([velocity(pool(1)),velocity(pool(2))]);
                    L_col(pos)=pool(nn);
                case 4
                    meshNO_c=floor((xPath(k)-Clb(1))/meshInterval)*N2;
                    meshNO_f=ceil((yPath(k)-min(Clb(2),Crb(2)))/meshInterval);
                    pool=[meshNO_c+meshNO_f,meshNO_c+meshNO_f+1];
                    [~,nn]=max([velocity(pool(1)),velocity(pool(2))]);
                    L_col(pos)=pool(nn);
                otherwise
                    meshNO_c=floor((xPath(k)-Clb(1))/meshInterval)*N2;
                    meshNO_f1=ceil((yPath(k)-min(Clb(2),Crb(2)))/meshInterval);
                    meshNO_f2=ceil((yPath(k+1)-min(Clb(2),Crb(2)))/meshInterval);
                    meshNO_f=max(meshNO_f1,meshNO_f2);
                    L_col(pos)=meshNO_c+meshNO_f;
            end
            L_row(pos)=(i-1)*length(t_x)+j;
            L_val(pos)=euclid_dist(xPath(k:k+1),yPath(k:k+1),2);
            pos=pos+1;
        end
        %--------------------------------------------------------
        plot(xPath,yPath,'b-');
        axis equal;
        xlim([Clb(1) Crb(1)]);ylim([min(Clb(2),Crb(2)) max(Clt(2),Crt(2))]);
        hold on;
    end
    t_print(i,:)=tmin;
end
hold off;

L=sparse(L_row,L_col,L_val,m,n); % Generate final sparse matrix: L


end


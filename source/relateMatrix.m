function [s_row, s_col, s_val, node, coordcell, N1, N2, N3] = relateMatrix( meshInterval, nodeNum, Clb, Clt, Crb, Crt, draw, velocity )
%*******************************EXPLANATION********************************
% 求解关联矩阵，并返回每个节点的坐标
% meshInterval 网格间距(m). 
% Clb 左下角坐标. Clt 左上角坐标. Crb 右下角坐标. Crt 右上角坐标.
% nodeNum 指的是一个网格的节点数,4的倍数
% flag控制是否画网格、节点图
% 还缺一个“打死”模块，用于应对不规则剖面（@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@）
%****************************BY CHENGLONG DUAN*****************************

% MinInterval = 1/1000; % 设置一个最小误差范围，用于浮点数比较（何琨修改）
%% Initialize the geometry of the profile
xLength=abs(Clb(1)-Crb(1));
yLength=abs(max(Clt(2),Crt(2))-min(Clb(2),Crb(2)));
if mod(xLength,meshInterval)~=0 || mod(yLength,meshInterval)~=0
    error('mesh error!');
end

%% Establish the relation among Node, Mesh, Coordinate
meshNum=(xLength/meshInterval)*(yLength/meshInterval);
[xmesh,ymesh]=meshgrid(Clb(1):meshInterval:Crb(1),min(Clb(2),Crb(2)):meshInterval:max(Clt(2),Crt(2)));

if draw==1  % 画网格图
    plot(xmesh,ymesh,'k',xmesh',ymesh','k'); % plot all the meshes in the plane
%     axis equal;
%     xlim([Clb(1) Crb(1)]);ylim([min(Clb(2),Crb(2)) max(Clt(2),Crt(2))]);
    hold on;
end

N1=(nodeNum-4)/4; % Number of nodes in each mesh edge
N2=yLength/meshInterval; % Number of edges in y-direction
N3=xLength/meshInterval; % Number of edges in x-direction
row=N1*N2+(N2+1);
column=N1*N3+(N3+1);
nodemat=zeros(row,column);coordcell=cell(1,row*column-N1*N1*meshNum);
node=1;i=1;
for x=Clb(1):xLength/(column-1):Crb(1)
    for y=min(Clb(2),Crb(2)):yLength/(row-1):max(Clt(2),Crt(2))
        if all(x~=xmesh(1,:))==1 && all(y~=ymesh(:,1))==1 % 何琨修改后的语句(适用于网格特别小的时候)：min(abs(x - xmesh(1,:)))>MinInterval && min(abs(y - ymesh(:,1)))>MinInterval ///段成龙的语句：all(x~=xmesh(1,:))==1 && all(y~=ymesh(:,1))==1
            nodemat(i)=0;
        else
            nodemat(i)=node;
            coordcell{node}=[x y];
            node=node+1;
        end
        i=i+1;
    end
end
%nodemat=flipud(nodemat);
memory=cell(meshNum,2);
for mesh=1:meshNum
    memory{mesh,1}=mesh;
    minC=N1*(ceil(mesh/N2)-1)+ceil(mesh/N2);
    maxC=minC+(N1+1);
    minR=N1*((mesh-N2*(ceil(mesh/N2)-1))-1)+(mesh-N2*(ceil(mesh/N2)-1));
    maxR=minR+(N1+1);
    %---把目标节点列向量化并剔除0----
    n=nodemat(minR:maxR,minC:maxC);
    n=n(:);
    n=n(n~=0);
    %-------------------------------
    memory{mesh,2}=n;
end
%--------------------------------------------------------------------------
% memory存储网格号和对应的节点号，节点调用其坐标可用coordcell{node}来实现
%--------------------------------------------------------------------------
if draw==1 % %画节点图
    for i=1:node-1
        plot(coordcell{i}(1),coordcell{i}(2),'o','color','k','Markerfacecolor','k','markeredgecolor','k');
        hold on;
    end
    hold off;
end

%% Establish the relate matrix
s_row=[];s_col=[];s_val=[];
for mesh=1:meshNum
    node_array=memory{mesh,2}; % 列向量
    for node1=node_array'  %遍历该网格下属的每一个节点
        waitingNode=node_array(node_array~=node1);
        C1=coordcell{node1}; %node1的坐标
        for node2=waitingNode'
            flag=0; %控制“不能越点传”
            flag2=0; %当横着或竖着的两点是一个大X或者大Y时，该标记就亮了
            C2=coordcell{node2}; %node2的坐标
            %if node2>node1  排除回传情况     && C2(1)~=C1(1)竖直传
            %----------------------------------------------------------------------------------
            if (C2(2)==C1(2) && C2(1)>(C1(1)+xLength/(N3*(1+N1)))) || (C2(2)==C1(2) && C2(1)<(C1(1)-xLength/(N3*(1+N1)))) %检验两点距离是一个小x的情况(+/-)
                for veri_node=waitingNode(waitingNode~=node2)'
                    if (coordcell{veri_node}(2)==C2(2))  %中间有点就跳出
                        flag=1;
                        break;
                    else
                        continue;
                    end
                end
                if flag==1
                    continue;
                else
                    flag2=1; %横着或竖着的两点是一个大X
                end
            end

            if (C2(1)==C1(1) && C2(2)>(C1(2)+yLength/(N2*(1+N1)))) || (C2(1)==C1(1) && C2(2)<(C1(2)-yLength/(N2*(1+N1)))) %检验两点距离是一个小y的情况(+/-)
                for veri_node=waitingNode(waitingNode~=node2)'
                    if (coordcell{veri_node}(1)==C2(1))
                        flag=1;
                        break;
                    else
                        continue;
                    end
                end
                if flag==1
                    continue;
                else
                    flag2=1; %横着或竖着的两点是一个大Y
                end
            end
            %------------------------------------------------------------------------------------            
                %------------------搜索速度对应的网格号---------------------
                [Row1,Col1]=find(nodemat==node1);[Row2,Col2]=find(nodemat==node2);
                inf=(N2*(Col1-2-N1)+Row1-1)/(1+N1);sup=(N2*(Col1-1)+Row1+N1)/(1+N1);
                inf=max(1,inf);
                set1=ceil(inf):1:sup;
                inf=(N2*(Col2-2-N1)+Row2-1)/(1+N1);sup=(N2*(Col2-1)+Row2+N1)/(1+N1);
                inf=max(1,inf);
                set2=ceil(inf):1:sup;
                m=intersect(set1,set2);
                for i=1:length(m)
                    index=ismember([node1 node2],memory{m(i),2});
                    if all(index)
                        meshNO=m(i);
                        break;
                    end
                end
                %---------------------------------------------------------
                s_row=[s_row node1];s_col=[s_col node2];
                %--------------------对边界速度进行处理---------------------
                if C1(2)==C2(2)&&C1(2)~=min(Clb(2),Crb(2))&&C1(2)~=max(Clt(2),Crt(2))&&flag2~=1 %横着连且不是上下边界且中间没有点
                    speed=max(velocity(meshNO),velocity(meshNO+1)); %Alternative:speed=(velocity(meshNO)+velocity(meshNO+1))/2;speed=max(velocity(meshNO),velocity(meshNO+1));
                    t=euclid_dist(C1,C2,1)/speed;
                else if C1(1)==C2(1)&&C1(1)~=Clb(1)&&C1(1)~=Crb(1)&&flag2~=1 %竖着连且不是左右边界且中间没有点
                        speed=max(velocity(meshNO),velocity(meshNO+N2)); %Alternative:speed=(velocity(meshNO)+velocity(meshNO+N2))/2;speed=max(velocity(meshNO),velocity(meshNO+N2));
                        t=euclid_dist(C1,C2,1)/speed;
                %---------------------------------------------------------
                    else
                        t=euclid_dist(C1,C2,1)/velocity(meshNO);
                    end
                end
                s_val=[s_val t];%还原
                %W(node1,node2)=t;
            %else
                %continue;
            %end
        end
    end
end
%W=sparse(s_row,s_col,s_val,node-1,node-1);
    
end
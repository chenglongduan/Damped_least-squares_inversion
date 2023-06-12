function [s_row, s_col, s_val, node, coordcell, N1, N2, N3] = relateMatrix( meshInterval, nodeNum, Clb, Clt, Crb, Crt, draw, velocity )
%*******************************EXPLANATION********************************
% ���������󣬲�����ÿ���ڵ������
% meshInterval ������(m). 
% Clb ���½�����. Clt ���Ͻ�����. Crb ���½�����. Crt ���Ͻ�����.
% nodeNum ָ����һ������Ľڵ���,4�ı���
% flag�����Ƿ����񡢽ڵ�ͼ
% ��ȱһ����������ģ�飬����Ӧ�Բ��������棨@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@��
%****************************BY CHENGLONG DUAN*****************************

% MinInterval = 1/1000; % ����һ����С��Χ�����ڸ������Ƚϣ������޸ģ�
%% Initialize the geometry of the profile
xLength=abs(Clb(1)-Crb(1));
yLength=abs(max(Clt(2),Crt(2))-min(Clb(2),Crb(2)));
if mod(xLength,meshInterval)~=0 || mod(yLength,meshInterval)~=0
    error('mesh error!');
end

%% Establish the relation among Node, Mesh, Coordinate
meshNum=(xLength/meshInterval)*(yLength/meshInterval);
[xmesh,ymesh]=meshgrid(Clb(1):meshInterval:Crb(1),min(Clb(2),Crb(2)):meshInterval:max(Clt(2),Crt(2)));

if draw==1  % ������ͼ
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
        if all(x~=xmesh(1,:))==1 && all(y~=ymesh(:,1))==1 % �����޸ĺ�����(�����������ر�С��ʱ��)��min(abs(x - xmesh(1,:)))>MinInterval && min(abs(y - ymesh(:,1)))>MinInterval ///�γ�������䣺all(x~=xmesh(1,:))==1 && all(y~=ymesh(:,1))==1
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
    %---��Ŀ��ڵ������������޳�0----
    n=nodemat(minR:maxR,minC:maxC);
    n=n(:);
    n=n(n~=0);
    %-------------------------------
    memory{mesh,2}=n;
end
%--------------------------------------------------------------------------
% memory�洢����źͶ�Ӧ�Ľڵ�ţ��ڵ�������������coordcell{node}��ʵ��
%--------------------------------------------------------------------------
if draw==1 % %���ڵ�ͼ
    for i=1:node-1
        plot(coordcell{i}(1),coordcell{i}(2),'o','color','k','Markerfacecolor','k','markeredgecolor','k');
        hold on;
    end
    hold off;
end

%% Establish the relate matrix
s_row=[];s_col=[];s_val=[];
for mesh=1:meshNum
    node_array=memory{mesh,2}; % ������
    for node1=node_array'  %����������������ÿһ���ڵ�
        waitingNode=node_array(node_array~=node1);
        C1=coordcell{node1}; %node1������
        for node2=waitingNode'
            flag=0; %���ơ�����Խ�㴫��
            flag2=0; %�����Ż����ŵ�������һ����X���ߴ�Yʱ���ñ�Ǿ�����
            C2=coordcell{node2}; %node2������
            %if node2>node1  �ų��ش����     && C2(1)~=C1(1)��ֱ��
            %----------------------------------------------------------------------------------
            if (C2(2)==C1(2) && C2(1)>(C1(1)+xLength/(N3*(1+N1)))) || (C2(2)==C1(2) && C2(1)<(C1(1)-xLength/(N3*(1+N1)))) %�������������һ��Сx�����(+/-)
                for veri_node=waitingNode(waitingNode~=node2)'
                    if (coordcell{veri_node}(2)==C2(2))  %�м��е������
                        flag=1;
                        break;
                    else
                        continue;
                    end
                end
                if flag==1
                    continue;
                else
                    flag2=1; %���Ż����ŵ�������һ����X
                end
            end

            if (C2(1)==C1(1) && C2(2)>(C1(2)+yLength/(N2*(1+N1)))) || (C2(1)==C1(1) && C2(2)<(C1(2)-yLength/(N2*(1+N1)))) %�������������һ��Сy�����(+/-)
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
                    flag2=1; %���Ż����ŵ�������һ����Y
                end
            end
            %------------------------------------------------------------------------------------            
                %------------------�����ٶȶ�Ӧ�������---------------------
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
                %--------------------�Ա߽��ٶȽ��д���---------------------
                if C1(2)==C2(2)&&C1(2)~=min(Clb(2),Crb(2))&&C1(2)~=max(Clt(2),Crt(2))&&flag2~=1 %�������Ҳ������±߽����м�û�е�
                    speed=max(velocity(meshNO),velocity(meshNO+1)); %Alternative:speed=(velocity(meshNO)+velocity(meshNO+1))/2;speed=max(velocity(meshNO),velocity(meshNO+1));
                    t=euclid_dist(C1,C2,1)/speed;
                else if C1(1)==C2(1)&&C1(1)~=Clb(1)&&C1(1)~=Crb(1)&&flag2~=1 %�������Ҳ������ұ߽����м�û�е�
                        speed=max(velocity(meshNO),velocity(meshNO+N2)); %Alternative:speed=(velocity(meshNO)+velocity(meshNO+N2))/2;speed=max(velocity(meshNO),velocity(meshNO+N2));
                        t=euclid_dist(C1,C2,1)/speed;
                %---------------------------------------------------------
                    else
                        t=euclid_dist(C1,C2,1)/velocity(meshNO);
                    end
                end
                s_val=[s_val t];%��ԭ
                %W(node1,node2)=t;
            %else
                %continue;
            %end
        end
    end
end
%W=sparse(s_row,s_col,s_val,node-1,node-1);
    
end
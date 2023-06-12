function  [t_print, L]=rayTracing( meshInterval, nodeNum, Clb, Clt, Crb, Crt, draw, s_x, s_y, t_x, t_y, velocity, m, n )
%UNTITLED2 Summary of this function goes here
%   �ٶȷ�������
% source��x�����������s_x�У�y�������s_y��
% terminus��x�����������t_x�У�y�������t_y��
% Define cell-velocity
%
% m,n����initialField.m����Ϊ�ٶ����񻮷���ʼ���ն���һ�ס�
% ���У�m--L����������������������n---L������������δ֪��������

%-----By Chenglong Duan,Nanjing University,2015.-----

xLength=abs(Clb(1)-Crb(1));
yLength=abs(max(Clt(2),Crt(2))-min(Clb(2),Crb(2)));
if mod(xLength,meshInterval)~=0 || mod(yLength,meshInterval)~=0
    error('mesh error!');
end

%-------------Ԥ��һ�������ʼ�ٶ��ļ���ģ��---------------------------------
%velocity=input('Please input velocity file:');
%--------------------------------------------------------------------------

%% Velocity map (��������)
% ���ٶ�ֵ��������ʱ�����ܾͳ�����meshInterval/2�Ŀ�ȱ����Ҫ�������û���⣬���������ڲ����������Ҫ���¿��ǣ�������
%if Clb(2)~=Crb(2) || Clt(2)~=Crt(2)
if draw==1
    figure;
    imagesc(velocity);
    axis equal;
    set(gca,'YDir','normal');
    %-------------------------------���߿��еĲ��ֵ�һ�λ��ٶ�ͼʱ�Ȳ����У��Ȼ���������Ȼ���꣬����xlim,ylim,xtick,ytick������
    xlim([0.5 3.5]);
    ylim([0.5 2.5]);
    set(gca,'xtick',0.5:3.5,'xticklabel',Clb(1):meshInterval:Crb(1)); 
    set(gca,'ytick',0.5:2.5,'yticklabel',min(Clb(2),Crb(2)):meshInterval:max(Clt(2),Crt(2)));
    %-----------------------------------------------------------------------------------------------------------------------
    colorbar;
    xlabel('Horizontal distance(m)');ylabel('Elevation(m)');
end

%% �����������
if draw==1
    figure;
end
[s_row,s_col,s_val,node,coordcell,N1,N2,N3] = relateMatrix( meshInterval, nodeNum, Clb, Clt, Crb, Crt, draw, velocity );
fprintf('Iter %d:The related matrix has created.\n',1); % �������1

%% ����SRP
%------�����������Դ����λ���ҳ���Ӧ�Ľڵ�------
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
clear count_source count_terminus; % count_source,count_terminus�Ǽ������
fprintf('Iter %d:The source and terminus points have been found.\n',1); % �������2
%��coordcell���浽Ĭ���ļ����£���Ҫʱ�ٵ������Լ�С�ռ�
% save coordcell;
% clear coordcell;
%---------------------------------------------
L_source=length(source);L_terminus=length(terminus);
t_print=zeros(L_source,L_terminus);
pos=1;  % ����ϡ�����Ԫ�ص�λ������
figure;
for i=1:L_source
    tmin=zeros(1,L_terminus);
    for j=1:L_terminus
        xPath=[];yPath=[];
        [t,path]=dijkstra(s_row,s_col,s_val,node,source(i),terminus(j));
        fprintf('Finish source %d,terminus %d. Drawing now...\n',i,j); % �������3
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
            if judge1==1  %���ڵ��������±߽�
                condition=1;
            end
            if judge2==1  %���ڵ��������ұ߽�
                condition=2;
            end
            if judge3==1  %��ֱ����
                condition=3;
            end
            if judge4==1  %ˮƽ����
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


function velocity = neg_vel_rev( neg_pos, xNum, yNum, velocity )
% �Ը��ٶȽ���ƽ����������Ϊ��Dijkstra�㷨���ԣ������������������ζ�����·�������ڡ�
% ���⣬���ڴ���negative circuit�����磬Bellman-Ford�㷨Ҳ��ʧЧ��
% neg_pos---���ٶ����ٶ������е�λ�ã�velocity---ԭʼ�ٶ�
% xNum---������yNum---����
% Output:velocity ��������ٶ�

% -----By Chenglong Duan,Nanjing University,2015.-----


neg_col = ceil(neg_pos/yNum);         % ��neg_posת���ɸ������ڵڼ���
neg_row = neg_pos-(neg_col-1)*yNum;   % ��neg_posת���ɸ������ڵڼ���

for i=1:length(neg_pos)
    flag = 0;  % ����λ�õ�
    flag1 = 0;  % ������������һ����
    
    if neg_row(i) == 1 && neg_col(i) == 1  % ���ٶ�λ�����½�
        if velocity(neg_pos(i)+1)>0 && velocity(neg_pos(i)+yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)+1)+velocity(neg_pos(i)+yNum))/2;
            flag1 = 1;
        end
        if velocity(neg_pos(i)+1)*velocity(neg_pos(i)+yNum)<0
            velocity(neg_pos(i)) = max(velocity(neg_pos(i)+1),velocity(neg_pos(i)+yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in left bottom:cell %d\n',neg_pos(i));
        end
        flag = 1; 
    end
    
    if neg_row(i) == yNum && neg_col(i) == 1  % ���ٶ�λ�����Ͻ�
        if velocity(neg_pos(i)-1)>0 && velocity(neg_pos(i)+yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)-1)+velocity(neg_pos(i)+yNum))/2;
            flag1 = 1;
        end
        if velocity(neg_pos(i)-1)*velocity(neg_pos(i)+yNum)<0
            velocity(neg_pos(i)) = max(velocity(neg_pos(i)-1),velocity(neg_pos(i)+yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in left top:cell %d\n',neg_pos(i));
        end
        flag = 1;
    end
    
    if neg_row(i) == yNum && neg_col(i) == xNum  % ���ٶ�λ�����Ͻ�
        if velocity(neg_pos(i)-1)>0 && velocity(neg_pos(i)-yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)-1)+velocity(neg_pos(i)-yNum))/2;
            flag1 = 1;
        end
        if velocity(neg_pos(i)-1)*velocity(neg_pos(i)-yNum)<0
            velocity(neg_pos(i)) = max(velocity(neg_pos(i)-1),velocity(neg_pos(i)-yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in right top:cell %d\n',neg_pos(i));
        end
        flag = 1;
    end
    
    if neg_row(i) == 1 && neg_col(i) == xNum  % ���ٶ�λ�����½�
        if velocity(neg_pos(i)+1)>0 && velocity(neg_pos(i)-yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)+1)+velocity(neg_pos(i)-yNum))/2;
            flag1 = 1;
        end
        if velocity(neg_pos(i)+1)*velocity(neg_pos(i)-yNum)<0
            velocity(neg_pos(i)) = max(velocity(neg_pos(i)+1),velocity(neg_pos(i)-yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in right bottom:cell %d\n',neg_pos(i));
        end
        flag = 1;
    end
    
    if neg_col(i) == 1 && flag == 0   % ���ٶ�λ����߽�
        if velocity(neg_pos(i)-1)>0 && velocity(neg_pos(i)+1)>0 && velocity(neg_pos(i)+yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)-1)+velocity(neg_pos(i)+1)+velocity(neg_pos(i)+yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(neg_pos(i)-1),velocity(neg_pos(i)+1),velocity(neg_pos(i)+yNum)];
            pool_positive = pool(pool>0);
            if length(pool_positive)==2
                velocity(neg_pos(i)) = sum(pool_positive)/2;
                flag1 = 1;
            end
            if length(pool_positive)==1
                velocity(neg_pos(i)) = pool_positive;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in left boundary:cell %d\n',neg_pos(i));
        end
        flag = 1;
    end
    
    if neg_col(i) == xNum && flag == 0   % ���ٶ�λ���ұ߽�
        if velocity(neg_pos(i)-1)>0 && velocity(neg_pos(i)+1)>0 && velocity(neg_pos(i)-yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)-1)+velocity(neg_pos(i)+1)+velocity(neg_pos(i)-yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(neg_pos(i)-1),velocity(neg_pos(i)+1),velocity(neg_pos(i)-yNum)];
            pool_positive = pool(pool>0);
            if length(pool_positive)==2
                velocity(neg_pos(i)) = sum(pool_positive)/2;
                flag1 = 1;
            end
            if length(pool_positive)==1
                velocity(neg_pos(i)) = pool_positive;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in right boundary:cell %d\n',neg_pos(i));
        end
        flag = 1;
    end
    
    if neg_row(i) == yNum && flag == 0   % ���ٶ�λ���ϱ߽�
        if velocity(neg_pos(i)-1)>0 && velocity(neg_pos(i)-yNum)>0 && velocity(neg_pos(i)+yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)-1)+velocity(neg_pos(i)-yNum)+velocity(neg_pos(i)+yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(neg_pos(i)-1),velocity(neg_pos(i)-yNum),velocity(neg_pos(i)+yNum)];
            pool_positive = pool(pool>0);
            if length(pool_positive)==2
                velocity(neg_pos(i)) = sum(pool_positive)/2;
                flag1 = 1;
            end
            if length(pool_positive)==1
                velocity(neg_pos(i)) = pool_positive;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in up boundary:cell %d\n',neg_pos(i));
        end
        flag = 1;
    end
    
    if neg_row(i) == 1 && flag == 0   % ���ٶ�λ���±߽�
        if velocity(neg_pos(i)+1)>0 && velocity(neg_pos(i)-yNum)>0 && velocity(neg_pos(i)+yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)+1)+velocity(neg_pos(i)-yNum)+velocity(neg_pos(i)+yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(neg_pos(i)+1),velocity(neg_pos(i)-yNum),velocity(neg_pos(i)+yNum)];
            pool_positive = pool(pool>0);
            if length(pool_positive)==2
                velocity(neg_pos(i)) = sum(pool_positive)/2;
                flag1 = 1;
            end
            if length(pool_positive)==1
                velocity(neg_pos(i)) = pool_positive;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in down boundary:cell %d\n',neg_pos(i));
        end
        flag = 1;
    end
    
    if flag == 0    % ��ͨ���
        if velocity(neg_pos(i)+1)>0 && velocity(neg_pos(i)-1)>0 && velocity(neg_pos(i)-yNum)>0 && velocity(neg_pos(i)+yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)+1)+velocity(neg_pos(i)-1)+velocity(neg_pos(i)-yNum)+velocity(neg_pos(i)+yNum))/4;
            flag1 = 1;
        end
        if flag1 ==0 && velocity(neg_pos(i)+1)>0 && velocity(neg_pos(i)-1)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)+1)+velocity(neg_pos(i)-1))/2;
            flag1 = 1;
        end
        if flag1 ==0 && velocity(neg_pos(i)-yNum)>0 && velocity(neg_pos(i)+yNum)>0
            velocity(neg_pos(i)) = (velocity(neg_pos(i)-yNum)+velocity(neg_pos(i)+yNum))/2;
            flag1 = 1;
        end
        if flag1 ==0
            pool = [velocity(neg_pos(i)+1),velocity(neg_pos(i)-1),velocity(neg_pos(i)-yNum),velocity(neg_pos(i)+yNum)];
            pool_positive = pool(pool>0);
            if length(pool_positive)>=2
                velocity(neg_pos(i)) = sum(pool_positive)/length(pool_positive);
                flag1 = 1;
            end
            if length(pool_positive)==1
                velocity(neg_pos(i)) = pool_positive;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in the middle plane:cell %d',neg_pos(i));
        end
    end
    
end

end


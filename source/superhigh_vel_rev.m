function velocity = superhigh_vel_rev( superhigh_pos, xNum, yNum, velocity )
% ���ٶ���������С�
% superhigh����Ϊ�ٶȴ���5 km/s��
% superhigh_pos---���ٶ����ٶ������е�λ�ã�velocity---ԭʼ�ٶ�
% xNum---������yNum---����
% Output:velocity ��������ٶ�

% -----By Chenglong Duan,Nanjing University,2015.-----


superhigh_col = ceil(superhigh_pos/yNum);         % ��superhigh_posת���ɸ������ڵڼ���
superhigh_row = superhigh_pos-(superhigh_col-1)*yNum;   % ��superhigh_posת���ɸ������ڵڼ���

for i=1:length(superhigh_pos)
    flag = 0;  % ����λ�õ�
    flag1 = 0;  % ������������һ����
    
    if superhigh_row(i) == 1 && superhigh_col(i) == 1  % ���ٶ�λ�����½�
        if velocity(superhigh_pos(i)+1)<=5 && velocity(superhigh_pos(i)+yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)+1)+velocity(superhigh_pos(i)+yNum))/2;
            flag1 = 1;
        end
        if any([velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)+yNum)]<=5)
            velocity(superhigh_pos(i)) = min(velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)+yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in left bottom:cell %d\n',superhigh_pos(i));
        end
        flag = 1; 
    end
    
    if superhigh_row(i) == yNum && superhigh_col(i) == 1  % ���ٶ�λ�����Ͻ�
        if velocity(superhigh_pos(i)-1)<=5 && velocity(superhigh_pos(i)+yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)-1)+velocity(superhigh_pos(i)+yNum))/2;
            flag1 = 1;
        end
        if any([velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)+yNum)]<=5)
            velocity(superhigh_pos(i)) = min(velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)+yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in left top:cell %d\n',superhigh_pos(i));
        end
        flag = 1;
    end
    
    if superhigh_row(i) == yNum && superhigh_col(i) == xNum  % ���ٶ�λ�����Ͻ�
        if velocity(superhigh_pos(i)-1)<=5 && velocity(superhigh_pos(i)-yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)-1)+velocity(superhigh_pos(i)-yNum))/2;
            flag1 = 1;
        end
        if any([velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)-yNum)]<=5)
            velocity(superhigh_pos(i)) = min(velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)-yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in right top:cell %d\n',superhigh_pos(i));
        end
        flag = 1;
    end
    
    if superhigh_row(i) == 1 && superhigh_col(i) == xNum  % ���ٶ�λ�����½�
        if velocity(superhigh_pos(i)+1)<=5 && velocity(superhigh_pos(i)-yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)+1)+velocity(superhigh_pos(i)-yNum))/2;
            flag1 = 1;
        end
        if any([velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)-yNum)]<=5)
            velocity(superhigh_pos(i)) = min(velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)-yNum));
            flag1 = 1;
        end
        if flag1 == 0
            fprintf('No choice in right bottom:cell %d\n',superhigh_pos(i));
        end
        flag = 1;
    end
    
    if superhigh_col(i) == 1 && flag == 0   % ���ٶ�λ����߽�
        if velocity(superhigh_pos(i)-1)<=5 && velocity(superhigh_pos(i)+1)<=5 && velocity(superhigh_pos(i)+yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)-1)+velocity(superhigh_pos(i)+1)+velocity(superhigh_pos(i)+yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)+yNum)];
            pool_normal = pool(pool<=5);
            if length(pool_normal)==2
                velocity(superhigh_pos(i)) = sum(pool_normal)/2;
                flag1 = 1;
            end
            if length(pool_normal)==1
                velocity(superhigh_pos(i)) = pool_normal;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in left boundary:cell %d\n',superhigh_pos(i));
        end
        flag = 1;
    end
    
    if superhigh_col(i) == xNum && flag == 0   % ���ٶ�λ���ұ߽�
        if velocity(superhigh_pos(i)-1)<=5 && velocity(superhigh_pos(i)+1)<=5 && velocity(superhigh_pos(i)-yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)-1)+velocity(superhigh_pos(i)+1)+velocity(superhigh_pos(i)-yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)-yNum)];
            pool_normal = pool(pool<=5);
            if length(pool_normal)==2
                velocity(superhigh_pos(i)) = sum(pool_normal)/2;
                flag1 = 1;
            end
            if length(pool_normal)==1
                velocity(superhigh_pos(i)) = pool_normal;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in right boundary:cell %d\n',superhigh_pos(i));
        end
        flag = 1;
    end
    
    if superhigh_row(i) == yNum && flag == 0   % ���ٶ�λ���ϱ߽�
        if velocity(superhigh_pos(i)-1)<=5 && velocity(superhigh_pos(i)-yNum)<=5 && velocity(superhigh_pos(i)+yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)-1)+velocity(superhigh_pos(i)-yNum)+velocity(superhigh_pos(i)+yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)-yNum),velocity(superhigh_pos(i)+yNum)];
            pool_normal = pool(pool<=5);
            if length(pool_normal)==2
                velocity(superhigh_pos(i)) = sum(pool_normal)/2;
                flag1 = 1;
            end
            if length(pool_normal)==1
                velocity(superhigh_pos(i)) = pool_normal;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in up boundary:cell %d\n',superhigh_pos(i));
        end
        flag = 1;
    end
    
    if superhigh_row(i) == 1 && flag == 0   % ���ٶ�λ���±߽�
        if velocity(superhigh_pos(i)+1)<=5 && velocity(superhigh_pos(i)-yNum)<=5 && velocity(superhigh_pos(i)+yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)+1)+velocity(superhigh_pos(i)-yNum)+velocity(superhigh_pos(i)+yNum))/3;
            flag1 = 1;
        end
        if flag1 == 0
            pool = [velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)-yNum),velocity(superhigh_pos(i)+yNum)];
            pool_normal = pool(pool<=5);
            if length(pool_normal)==2
                velocity(superhigh_pos(i)) = sum(pool_normal)/2;
                flag1 = 1;
            end
            if length(pool_normal)==1
                velocity(superhigh_pos(i)) = pool_normal;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in down boundary:cell %d\n',superhigh_pos(i));
        end
        flag = 1;
    end
    
    if flag == 0    % ��ͨ���
        if velocity(superhigh_pos(i)+1)<=5 && velocity(superhigh_pos(i)-1)<=5 && velocity(superhigh_pos(i)-yNum)<=5 && velocity(superhigh_pos(i)+yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)+1)+velocity(superhigh_pos(i)-1)+velocity(superhigh_pos(i)-yNum)+velocity(superhigh_pos(i)+yNum))/4;
            flag1 = 1;
        end
        if flag1 ==0 && velocity(superhigh_pos(i)+1)<=5 && velocity(superhigh_pos(i)-1)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)+1)+velocity(superhigh_pos(i)-1))/2;
            flag1 = 1;
        end
        if flag1 ==0 && velocity(superhigh_pos(i)-yNum)<=5 && velocity(superhigh_pos(i)+yNum)<=5
            velocity(superhigh_pos(i)) = (velocity(superhigh_pos(i)-yNum)+velocity(superhigh_pos(i)+yNum))/2;
            flag1 = 1;
        end
        if flag1 ==0
            pool = [velocity(superhigh_pos(i)+1),velocity(superhigh_pos(i)-1),velocity(superhigh_pos(i)-yNum),velocity(superhigh_pos(i)+yNum)];
            pool_normal = pool(pool<=5);
            if length(pool_normal)>=2
                velocity(superhigh_pos(i)) = sum(pool_normal)/length(pool_normal);
                flag1 = 1;
            end
            if length(pool_normal)==1
                velocity(superhigh_pos(i)) = pool_normal;
                flag1 = 1;
            end
        end
        if flag1 == 0
            fprintf('No choice in the middle plane:cell %d',superhigh_pos(i));
        end
    end
    
end

end

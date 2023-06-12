function velocity = neg_vel_rev( neg_pos, xNum, yNum, velocity )
% 对负速度进行平滑修正。因为对Dijkstra算法而言，如果包含负环，则意味着最短路径不存在。
% 此外，对于存在negative circuit的网络，Bellman-Ford算法也会失效。
% neg_pos---负速度在速度向量中的位置，velocity---原始速度
% xNum---列数，yNum---行数
% Output:velocity 修正后的速度

% -----By Chenglong Duan,Nanjing University,2015.-----


neg_col = ceil(neg_pos/yNum);         % 把neg_pos转化成该网格在第几列
neg_row = neg_pos-(neg_col-1)*yNum;   % 把neg_pos转化成该网格在第几行

for i=1:length(neg_pos)
    flag = 0;  % 控制位置的
    flag1 = 0;  % 控制优先做哪一步的
    
    if neg_row(i) == 1 && neg_col(i) == 1  % 负速度位于左下角
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
    
    if neg_row(i) == yNum && neg_col(i) == 1  % 负速度位于左上角
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
    
    if neg_row(i) == yNum && neg_col(i) == xNum  % 负速度位于右上角
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
    
    if neg_row(i) == 1 && neg_col(i) == xNum  % 负速度位于右下角
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
    
    if neg_col(i) == 1 && flag == 0   % 负速度位于左边界
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
    
    if neg_col(i) == xNum && flag == 0   % 负速度位于右边界
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
    
    if neg_row(i) == yNum && flag == 0   % 负速度位于上边界
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
    
    if neg_row(i) == 1 && flag == 0   % 负速度位于下边界
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
    
    if flag == 0    % 普通情况
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


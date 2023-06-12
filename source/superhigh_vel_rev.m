function velocity = superhigh_vel_rev( superhigh_pos, xNum, yNum, velocity )
% 负速度修正后进行。
% superhigh定义为速度大于5 km/s。
% superhigh_pos---高速度在速度向量中的位置，velocity---原始速度
% xNum---列数，yNum---行数
% Output:velocity 修正后的速度

% -----By Chenglong Duan,Nanjing University,2015.-----


superhigh_col = ceil(superhigh_pos/yNum);         % 把superhigh_pos转化成该网格在第几列
superhigh_row = superhigh_pos-(superhigh_col-1)*yNum;   % 把superhigh_pos转化成该网格在第几行

for i=1:length(superhigh_pos)
    flag = 0;  % 控制位置的
    flag1 = 0;  % 控制优先做哪一步的
    
    if superhigh_row(i) == 1 && superhigh_col(i) == 1  % 高速度位于左下角
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
    
    if superhigh_row(i) == yNum && superhigh_col(i) == 1  % 负速度位于左上角
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
    
    if superhigh_row(i) == yNum && superhigh_col(i) == xNum  % 负速度位于右上角
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
    
    if superhigh_row(i) == 1 && superhigh_col(i) == xNum  % 负速度位于右下角
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
    
    if superhigh_col(i) == 1 && flag == 0   % 负速度位于左边界
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
    
    if superhigh_col(i) == xNum && flag == 0   % 负速度位于右边界
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
    
    if superhigh_row(i) == yNum && flag == 0   % 负速度位于上边界
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
    
    if superhigh_row(i) == 1 && flag == 0   % 负速度位于下边界
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
    
    if flag == 0    % 普通情况
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

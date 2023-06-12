function [s0, m, n, xNum, yNum, rnorm, xnorm, R, var] = initialField( T, meshInterval, Clb, Clt, Crb, Crt, s_x, s_y, t_x, t_y )
%为计算初始场Ls(0)=t提供L & s0(v0).
%   L为straight-ray tracing情况下的系数矩阵.
%   网格数和SPM中的网格一致，以避免SPM计算中重新识别波速.
%   Sources and receivers arrangement rule:
%CANNOT be placed at angle points of cells simultaneously!
%   网格排布规则见论文.
%   Inputs:
%   T; Observed traveltime data.
%   meshInterval:网格间隔长度 (一般取0.5或1.0)
%   Clb:Left-bottom (x,y)
%   Clt:Left-top (x,y)
%   Crb:Right-bottom (x,y)
%   Crt:Right-top (x,y)
%
%   Outputs:
%   v0: A vector! (排布规则见论文)
%   m,n: size(L)

%-----By Chenglong Duan,Nanjing University,2015.-----


%% Extract geometry info. of the profile
xNum=(Crb(1)-Clb(1))/meshInterval;
yNum=abs(max(Clt(2),Crt(2))-min(Clb(2),Crb(2)))/meshInterval;
n=xNum*yNum;  % L的列数
m=length(s_x)*length(t_x); % L的行数

%% Check if sources and receivers are placed at angle points of cells simultaneously
if any(s_x~=Clb(1)) || any(t_x~=Crb(1))
    error('Sources and receivers should be side boundaries of the profile!');
end
left=s_y/meshInterval;right=t_y/meshInterval;
if any(left==fix(left)) && any(right==fix(right))
    error('Sources and receivers CANNOT be placed at angle points of cells simultaneously!');
end


%% Compute matrix L
% 构建射线长度和网格号之间的关系
L_row=zeros(1,m*n);L_col=zeros(1,m*n);L_val=zeros(1,m*n);
i=1;
for s=1:length(s_x)
    for t=1:length(t_x)
        slope=(s_y(s)-t_y(t))/(s_x(s)-t_x(t));
        % 给出射线与网格的交点坐标
        x1=s_x(s):meshInterval:t_x(t); % 与竖线一定有交点,且左右边界确定
        y1=s_y(s)+slope*(x1-s_x(s));
        
        if slope>=0  
            y2Inf=floor(s_y(s)/meshInterval)*meshInterval;
            y2Sup=ceil(t_y(t)/meshInterval)*meshInterval;
        end
        if slope<0 
            y2Inf=floor(t_y(t)/meshInterval)*meshInterval;
            y2Sup=ceil(s_y(s)/meshInterval)*meshInterval;
        end
        
        y2=y2Inf:meshInterval:y2Sup;
        x2=(y2-s_y(s))/slope+s_x(s);
        x2_logic=(x2>=s_x(s) & x2<=t_x(t));
        x2=x2(x2_logic);
        y2=y2(x2_logic);
        
        [~,~,index]=intersect(x1,x2);
        x2(index)=[];y2(index)=[];
        xComb=[x1,x2]; yComb=[y1,y2];
        [xComb,order]=sort(xComb,'ascend');yComb=yComb(order);
        
        % 对应网格号赋值
        for k=1:length(xComb)-1
            meshNO_c=floor((xComb(k)-Clb(1))/meshInterval)*yNum;
            meshNO_f1=ceil((yComb(k)-min(Clb(2),Crb(2)))/meshInterval); 
            meshNO_f2=ceil((yComb(k+1)-min(Clb(2),Crb(2)))/meshInterval); 
            meshNO_f=max(meshNO_f1,meshNO_f2);
            L_row(i)=(s-1)*length(t_x)+t;
            L_col(i)=meshNO_c+meshNO_f;
            L_val(i)=euclid_dist(xComb(k:k+1),yComb(k:k+1),2);
            i=i+1;
        end
    end
end
L_row(L_row==0)=[]; L_col(L_col==0)=[]; L_val(L_val==0)=[];
L=sparse(L_row,L_col,L_val,m,n); % Generate final sparse matrix: L


%% Compute v0
% References:
% W. Yang. 1997. Theory and Methods of Geophysical Inversion.
% Beijing:Geological publishing house.

msg=['The equation has an exact solution                                    '
     'We have a least square solution                                       '
     'Ill-posed equation(unstable). Compromised rule is employed            '
     'Under-determined equations. We are trying to get a min-length solution'];

L=full(L);
rule1=rank(L); % rank deficiency
rule2=rcond(L'*L); % ill-posed 接近0.0,也叫奇异
rule3= m==n;
rule4= m<n;
rule5= m>n;

if rule3 && abs(rule2-1)<1e-3 % 精确解
    disp(msg(1,:));
    s0=L\T;
    rnorm=0;
    xnorm=norm(s0);
    R=eye(n,n);  % R=eye(m,n);
    var='none';
end
if rule5 && rule1==n && rule2 >= 1e-3 % 超定，排除病态情况
    disp(msg(2,:));
    s0=pinv(L)*T;  % Moore-Penrose pseudoinverse method
    r=L*s0-T;
    rnorm=norm(r);
    xnorm=norm(s0);
    R=eye(n,n);  % High resolution  R=eye(m,n);
    var='none';
end
if rule4 && rule2 >= 1e-3 % 欠定，排除病态情况
    disp(msg(4,:));
    s0=pinv(L)*T;  % Moore-Penrose pseudoinverse method
    r=L*s0-T;
    rnorm=norm(r);
    xnorm=norm(s0);
    [~,~,V]=svd(L);
    Vp=V(:,1:rule1);
    R=Vp*Vp';  % Low resolution    R=L'*pinv(L'*L)*L;
    var='none';
end
if rule2 < 1e-3
    disp(msg(3,:));
    [U,S,V]=svd(L);
    % Yang's method(Wiggins method + Schochastic inverse)
    %-----------------I.Wiggins method-------------------
    % 单位方差
    s=diag(S);  s_length=length(s);
    varW=zeros(size(V,1),s_length);
    for k=1:size(V,1)
        for q=1:s_length
            ele=0;
            for i=1:q
                ele=ele+(V(k,i)/s(i))^2;
            end
            varW(k,q)=ele;
        end
    end
    % 分辨率
    spreadR=zeros(s_length,1);
    for q=1:s_length
        [~,~,V]=svd(L);
        if q < n
            V(:,q+1:end)=0;
        end
        RW=V*V';
        sr=RW-eye(size(RW));
        spreadR(q)=norm(sr,'fro')^2;
    end
    % 由于是单位方差，先选哪一行
    sampleNum=10; varw_length=size(varW,1); sum=0;
    for j=1:ceil(varw_length/sampleNum):varw_length  % 抽样
        var_vector=varW(j,:)';  % 1-->1 ~ n
        %figure;
        %plot(spreadR,var_vector,'*'); 
        %xlabel('Resolution');ylabel('Variation');hold on;
        delta=diff(var_vector);
        for i=1:length(delta)
            if abs(delta(i))>1e-2
                q_thr=i; % q_thr为保留奇异值的个数,即矩阵L的有效自由度
                break;
            end
        end
        %plot(spreadR(q_thr),var_vector(q_thr),'ro');
        sum=sum+q_thr;
    end
    q_thr=floor(sum/sampleNum);
    t=s(q_thr); % 保留的最小奇异值
    
    %-----------------II.Yang's method--------------------
    for k=1:s_length
        if k>=q_thr
            s(k)=t^2/s(k); % 对小奇异值做的修改
        end
        if k<q_thr
            s(k)=s(k)+(k/q_thr)*(t^2/s(q_thr+1)-s(q_thr)); % 对大奇异值做的修改
        end
    end
    
    S_size=size(S); lnComple=S_size(1)-s_length; colComple=S_size(2)-s_length;
    if lnComple
        S_mod=[diag(s);zeros(lnComple,s_length)];
    end 
    if colComple
        S_mod=[diag(s),zeros(s_length,colComple)];
    end
    if lnComple == 0 && colComple == 0
        S_mod=diag(s);
    end
    %-----------------------------------------------------
    s0=V*(S_mod\(U'*T));
    r=L*s0-T;
    rnorm=norm(r);
    xnorm=norm(s0);
    Vp=V(:,1:s_length);
    R=Vp*Vp';  % R=V*(S_mod\(S*V'));
    var=zeros(n,1);
    for k=1:n
        for j=1:s_length
            var(k)=var(k)+(V(k,j)/s(j))^2; % 单位方差(列向量：网格数*1),即每一未知数对应一个方差
        end
    end
end

% Last step
% v0=1./s0;

% 误差棒图
% n=1:256;
% errorbar(n,v0,var,'.');

end

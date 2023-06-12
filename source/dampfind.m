function damp = dampfind( L, T, n, num, p, delta, damp_fig )

%该函数用于找最优正则化因子(阻尼因子)
%   Inputs:
%   L,T,n---含义同之前。
%   num---离散点的个数
%   p---指数的初始值
%   delta---指数两点之间的间隔,e.g.lamda=[10^p,10^(p+delta)].
%   最终，用于试验的阻尼因子的范围在[10^p, 10^(p+(num-1)*delta)]
%   Outputs:
%   damp---用于迭代计算的最优阻尼因子

% -----By Chenglong Duan,Nanjing University,2015.-----

L=sparse(L);
I=speye(n);
z=zeros(n,1);

damp_pool=zeros(num,1);
xnorm=zeros(num,1);
rnorm=zeros(num,1);

for k=1:num
    lamda=10^p;
    B=[L;lamda*I];
    c=[T;z];
    x=B\c;
    r=T-L*x;
    damp_pool(k)=lamda;
    xnorm(k)=norm(x);
    rnorm(k)=norm(r);
    p=p+delta;
end

% para_de=polyfit(rnorm,xnorm,5); % 4阶多项式拟合
% poly_x=min(rnorm):0.0001:max(rnorm);
% poly_y=polyval(para_de,poly_x);
% [~, rnorm_op, xnorm_op] = curvature (para_de, min(rnorm):0.01:max(rnorm));

h=figure;
loglog(rnorm,xnorm,'o');
[rnorm_op, xnorm_op] = ginput(1);
close(h);

flag=0;
if ~isempty(find(rnorm_op==rnorm, 1)) % 验证拾取的点是否和原数据重合
    index= rnorm==rnorm_op;
    damp=damp_pool(index);
    flag=1;
end

% if flag==0
%     if rnorm(1)-rnorm_op > 0
%         damp=(damp_pool(i)-damp_pool(i+1))/(rnorm(i)-rnorm(i+1))*(rnorm_op-rnorm(i))+damp_pool(i); % 线性(向前)外插计算damp
%         flag=1;
%     end
%     if rnorm(end)-rnorm_op < 0
%         damp=; % 线性(向后)外插计算damp
%         flag=1;
%     end
% end

if flag==0
    for i=1:length(rnorm)-1
        in=(rnorm(i)-rnorm_op)*(rnorm(i+1)-rnorm_op);
        if in < 0
            damp=(damp_pool(i)-damp_pool(i+1))/(rnorm(i)-rnorm(i+1))*(rnorm_op-rnorm(i))+damp_pool(i);% 线性内插计算damp
            break;
        end
    end
end

if damp_fig==1
    figure(51);
    loglog(rnorm,xnorm,'o',rnorm_op,xnorm_op,'r*');
    hold on;
end

end


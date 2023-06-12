function [d,path]=dijkstra(s_row,s_col,s_val,node,s,t)

% W---关联矩阵；s---起点；t---终点
% d---最短距离；path---最短路径

%   给节点编号（可以画一张图）；
%   指定源、检节点号；
%   权重矩阵（关联矩阵）---Square matrix

%-----By Chenglong Duan,Nanjing University,2015.-----

%% relate matrix预处理
n=node-1;  % [n,m]=size(W);
spamat=sparse(s_row,s_col,s_val,n,n);
W=full(spamat); % full()  ans(ans==0)=inf
clear spamat;
W(W==0)=inf;
% ix=(W==0); % 返回一个逻辑判断矩阵(稀疏矩阵)-->存成sparse形式
% W(ix)=inf;
% if n~=m
%     error('Square W required');
% end

%% 初始化
visited(1:n)=0;
dist(1:n)=inf;
parent(1:n)=0;
dist(s)=0;
d=inf;

%% 每个节点与起始节点的关系
for i=1:n-1
    ix=(visited==0); % 未访问过的节点的位置
    vec(1:n)=inf;
    vec(ix)=dist(ix);
    [~,u]=min(vec);  % u是靶点
    visited(u)=1;
    for v=1:n   % 全部验证点
        if (W(u,v)+dist(u) < dist(v))  % dist(u)是上一轮的值,W(u,v)是从u出发的
                                       %所有可能值，二者之和为累加值。dist(v)是
                                       %所有已经计算过的值(对未计算过的都是inf，无意义)
            dist(v) = dist(u)+W(u,v);
            parent(v) = u;
        end
    end
end

%% 由终点回溯最短路径
if parent(t)~=0
    path=t;
    d=dist(t);
    while t~=s
        p=parent(t);
        path=[p path];
        t=p;
    end
end


end


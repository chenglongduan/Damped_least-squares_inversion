function [d,path]=dijkstra(s_row,s_col,s_val,node,s,t)

% W---��������s---��㣻t---�յ�
% d---��̾��룻path---���·��

%   ���ڵ��ţ����Ի�һ��ͼ����
%   ָ��Դ����ڵ�ţ�
%   Ȩ�ؾ��󣨹�������---Square matrix

%-----By Chenglong Duan,Nanjing University,2015.-----

%% relate matrixԤ����
n=node-1;  % [n,m]=size(W);
spamat=sparse(s_row,s_col,s_val,n,n);
W=full(spamat); % full()  ans(ans==0)=inf
clear spamat;
W(W==0)=inf;
% ix=(W==0); % ����һ���߼��жϾ���(ϡ�����)-->���sparse��ʽ
% W(ix)=inf;
% if n~=m
%     error('Square W required');
% end

%% ��ʼ��
visited(1:n)=0;
dist(1:n)=inf;
parent(1:n)=0;
dist(s)=0;
d=inf;

%% ÿ���ڵ�����ʼ�ڵ�Ĺ�ϵ
for i=1:n-1
    ix=(visited==0); % δ���ʹ��Ľڵ��λ��
    vec(1:n)=inf;
    vec(ix)=dist(ix);
    [~,u]=min(vec);  % u�ǰе�
    visited(u)=1;
    for v=1:n   % ȫ����֤��
        if (W(u,v)+dist(u) < dist(v))  % dist(u)����һ�ֵ�ֵ,W(u,v)�Ǵ�u������
                                       %���п���ֵ������֮��Ϊ�ۼ�ֵ��dist(v)��
                                       %�����Ѿ��������ֵ(��δ������Ķ���inf��������)
            dist(v) = dist(u)+W(u,v);
            parent(v) = u;
        end
    end
end

%% ���յ�������·��
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


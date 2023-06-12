function resultXY = xyGen (meshInterval, Clb, Clt, Crb, Crt, xNum, yNum)
% resultXY������Surfer�����ƽ������
% Clb�����ճ����������½ڵ�������ꣻ
% Crb�����ճ����������½ڵ�������ꣻ
% Clt�����ճ����������Ͻڵ�������ꣻ
% Crt�����ճ����������Ͻڵ�������ꡣ
% xNum������������
% yNum������������

% -----By Chenglong Duan,Nanjing University,2015.-----


distance = Crb(1)-Clb(1);
if (distance-floor(distance)~=0) && (distance-floor(distance)~=0.5)
    error('distance must be .0 or .5');
end
x = Clb(1)+meshInterval/2:meshInterval:(Crb(1)-meshInterval/2);
y = min(Clb(2),Crb(2))+meshInterval/2:meshInterval:(max(Clt(2),Crt(2))-meshInterval/2);
yLine1=(Crb(2)-Clb(2))*x/distance+Clb(2);
yLine2=(Crt(2)-Clt(2))*x/distance+Clt(2);
resultXY = zeros(xNum*yNum,2);
r2 = 0;
for i=1:length(x)
    count1 = round((yLine1(i)-min(Clb(2),Crb(2)))/meshInterval)+1;
    yStart = y(count1);  
    count2 = round((yLine2(i)-min(Clb(2),Crb(2)))/meshInterval)+1;
    yEnd = y(count2-1);
    
    yy=yStart:meshInterval:yEnd;
    r1 = length(yy);
    result1=zeros(r1,2);
    result1(:,1)=x(i);result1(:,2)=yy;
    resultXY(r2+1:r1+r2,:) = result1;
%     if i==1
%         result = result1;
%     else
%         result = [result;result1];
%     end
    r2 = r2+r1;
end

function  main( T_file, meshInterval, nodeNum, Clb, Clt, Crb, Crt, s_x, s_y, t_x, t_y, initial_fig, draw, damp_fig, file_XYZ )
% main����������(1)ʵ�ִ�ѭ����(2)���Ʋ�������������
%   ���д˳���ǰ���ȱ���inputsFile.m

%-----By Chenglong Duan,Nanjing University,2015.-----

% Read observed traveltime data.
T=xlsread(T_file,'sheet1');
disp('Observed traveltime data have been successfully loaded. Continue(1) or Cancel(2)?');
flag1=input('');
if flag1==2
    error('The program has been cancelled!');
end
clear flag1;

% Initialize.
T1 = T';
T_vector = T1(:); % Tת��������
[s0, m, n, xNum, yNum, rnorm_ini, xnorm_ini, R, var] ...
    = initialField( T_vector, meshInterval, Clb, Clt, Crb, Crt, s_x, s_y, t_x, t_y );
                      % Note:(1)�˴���m,n�Ǿ���L��size��������Ԫ��size;(2)s0��������
v0=1./s0;
fprintf(['Initial r_norm = %f\n' 'Initial x_norm = %f\n'],rnorm_ini,xnorm_ini);
if initial_fig == 1
    vfig1=reshape(v0,yNum,xNum); % Note:(1)xNum---ˮƽ����ĸ�����,yNum---��ֱ����ĸ�����.
                                 %      (2)xNum�������м���(����),yNum�������м���(����)
    vfig2=flipud(vfig1);
    figure('Name','Initial Pixel Map in Wave Velocity');
    imagesc(vfig2);
    axis equal;colorbar; % ����һ��axis equal
    
    % �����ͼ
    figure('Name','Initial Errorbar Map (slowness,Variance)');
    num=1:n;
    errorbar(num,s0,var,'.'); xlabel('k');ylabel('s_k');
    % �ֱ�����Ԫͼ
    figure('Name','Initial Pixel Map in Resolution');
    Rm1=reshape(diag(R),yNum,xNum);
    Rm2=flipud(Rm1);
    imagesc(Rm2);axis equal;colorbar;
end
% �ٶ�����
if any(v0<0)
    neg_pos = find(v0<0);
    v0 = neg_vel_rev( neg_pos, xNum, yNum, v0 );
end
if any(v0<0)
    neg_pos = v0<0;
    v0(neg_pos) = 0.1;
end
if any(v0>5)
    superhigh_pos = find(v0>5);
    v0 = superhigh_vel_rev( superhigh_pos, xNum, yNum, v0 );
end
if any(v0>5)
    superhigh_pos = v0>5;
    v0(superhigh_pos) = 5.0;
end
s0 = 1./v0; % s0��������
% Generate X-Y-Z
resultXY = xyGen (meshInterval, Clb, Clt, Crb, Crt, xNum, yNum);
surferXYZ = [resultXY,v0];
str_ini = strcat('Initial_',file_XYZ);
xlswrite(str_ini,surferXYZ);

% -------------------------Main iteration--------------------------------
bigIter = 1;
s=s0;  s1=s0;  % ����������ֵ(s1���м���)
big_itnlim=8;  % 500���������Ĵ�ѭ����
dampCom=zeros(big_itnlim,1); 
while (norm(s-s1)>=1e-3 || bigIter==1) && bigIter < big_itnlim  % ǰ�����μ�������
    s1=s;
    velocity=1./s1;
    
    [t_print, L]=rayTracing( meshInterval, nodeNum, Clb, Clt, Crb, Crt,...
        draw, s_x, s_y, t_x, t_y, velocity, m, n ); % ���ڲ���ͼ����draw
    
    dt1=abs(T-t_print)'; % T��t_print���Ǿ���
    dt=dt1(:); % ת��Ϊ������
    % Parameters in damp LSQR----------------------------------------------
    % (1)determine 'damp'
    number=17;
    p=-5;
    delta=0.5;
    damp = dampfind( L, dt, n, number, p, delta, damp_fig );
    dampCom(bigIter)=damp;
    % (2)stopping criteria for LSQR
    atol=0;  % A-tolϵ�������������(��׷�ٷ�ʽ�й�)
    btol=0;  % b-tol��ʱ���ݵ�������(��ʰȡ�����й�)
    conlim=0;% cond-limit�������������ļ���(ԽСԽ�ã�������С��0;����ܳ���1e+8)
    
    itnlim=400; % LSQR����������
    show=1; % gives an iteration log on the screen
    %----------------------------------------------------------------------
    [ ds, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, varLSQR ]...
        = lsqrSOL( m, n, L, dt, damp, atol, btol, conlim, itnlim, show );
    disp(varLSQR);
    
    fprintf('Big iteration %d has finished...\n',bigIter); % �����������
    
    bigIter = bigIter+1;
    s=s1+ds;
    v = 1./s;  % ���յ��ٶ�
    
    % �Ը��ٶȺͳ����ٶȽ�������-------------------
    if any(v<0)
        neg_pos = find(v<0);  % Find the negative position
        v = neg_vel_rev( neg_pos, xNum, yNum, v ); % Return the revised velocity
    end
    if any(v<0)    % Re-verify and re-revise if negative values still exist
        neg_pos = v<0;
        v(neg_pos) = 0.1;  % �����и��ٶȣ������ٶ�ͳһ����Ϊ0.1
    end
    
    if any(v>6)
        superhigh_pos = find(v>6);
        v = superhigh_vel_rev( superhigh_pos, xNum, yNum, v );
    end
    if any(v>6)
        superhigh_pos = v>6;
        v(superhigh_pos) = 6.0;  % �����г����ٶȣ������ٶ�ͳһ����Ϊ5.0
    end
    s = 1./v; % sҲ����һ������
    % ----------------------------------
    
    surferXYZ_final = [resultXY,v];
    ss=num2str(bigIter-1);
    sss=strcat(file_XYZ,ss);
    dlmwrite(sss,surferXYZ_final);
    if bigIter > 2
        figure('Name','Comparison between observed and estimated time');
        for i=1:length(s_x)
            plot(T(i,:),t_y,'b*-',t_print(i,:),t_y,'ro-');
            hold on;
        end
        xlabel('Arrival time');ylabel('Depth');
        legend;
    end
end
% v = 1./s;  % ���յ��ٶ�

% ����Ļ�����bigIter��dampCom
fprintf('The actual iteration number is: %d\n',bigIter-1);
dampCom=dampCom(1:bigIter-1);  % ��ÿһ�δ�����õ�damp������
fprintf('Damps for each main iteration:\n');
fprintf('       %e\n',dampCom);
% -------------------------iteration ends------------------------------


% ���ٿ��ӻ�����
% Pixel map
% v1=reshape(v,yNum,xNum);
% v2=flipud(v1); % v2����imagesc��ͼ
% figure('Name','Final Pixel Map in Wave Velocity');
% imagesc(v2); colorbar;

% % Surfer map
% surferXYZ_final = [resultXY,v];
% % Put X-Y-V into a xls.
% xlswrite(file_XYZ,surferXYZ_final,'sheet1');
% 
% �۲���ʱ��Ԥ����ʱ�ıȽ�(X;Arrival time;Y:Depth)
% figure('Name','Comparison between observed and estimated time');
% for i=1:length(s_x)
%     plot(T(i,:),t_y,'b*-',t_print(i,:),t_y,'ro-');
%     hold on;
% end
% xlabel('Arrival time');ylabel('Depth');
% legend;

% �����һ�����ļ�������һ�μ����������Ϣ
disp('END OF COMPUTATION. NOW GENERATING COMPUTATION LOG...');
dateinfo=floor(clock);
fileName=strcat('data','_',num2str(dateinfo(2)),'_',num2str(dateinfo(3)),...
    '_',num2str(dateinfo(1)),'_',num2str(dateinfo(4)),'-',...
    num2str(dateinfo(5)),'-',num2str(dateinfo(6)),'.txt');
fid=fopen(fileName,'wt');
head0='**************************Initial Estimation****************************'; % 80�ַ�
text1={'||r||_2 = ','||x||_2 = '};
fprintf(fid,[head0 '\n']);
fprintf(fid,'%s %f    %s %f',text1{1},rnorm_ini,text1{2},xnorm_ini);
head1='*************************Computed Traveltimes***************************'; % 80�ַ�
text2={'||r||_2 = ','sqrt(||r||^2+damp^2*||x||^2)','||A||_F','cond(A)',...
    '||Ar-damp^2*x||','||x||_2 = ','varLSQR','istop in LSQR (last)','itn in LSQR (last)'};
fprintf(fid,[head1 '\n']);
fprintf(fid,'%s %f\n    %s %f\n    %s %f\n    %s %f\n    %s %f\n    %s %f\n    %s %f\n',...
    text2{1},r1norm,text2{2},r2norm,text2{3},Anorm,text2{4},Acond,text2{5},...
    Arnorm,text2{6},xnorm,text2{7},varLSQR,text2{8},istop,text2{9},itn);
fclose(fid);
disp('DONE!');

end


function createfigure(X1, Y1, Y2)
%CREATEFIGURE1(X1, Y1, Y2)
%  X1:  x ���ݵ�ʸ��
%  Y1:  y ���ݵ�ʸ��
%  Y2:  y ���ݵ�ʸ��

%  �� MATLAB �� 15-Oct-2016 13:14:45 �Զ�����

% ���� figure
figure1 = figure;

% ���� subplot
subplot1 = subplot(2,1,1,'Parent',figure1);
%% ȡ�������е�ע���Ա���������� X ��Χ
% xlim(subplot1,[0.001 0.004]);
%% ȡ�������е�ע���Ա���������� Y ��Χ
% ylim(subplot1,[38 38.6]);
box(subplot1,'on');
hold(subplot1,'on');

% ���� ylabel
ylabel('PSNR (dB)');

% ���� xlabel
xlabel('\lambda');

% ���� plot
plot(X1,Y1,'Parent',subplot1,'MarkerFaceColor',[0 0 0],'MarkerSize',3,...
    'Marker','o',...
    'LineWidth',1,...
    'LineStyle','-.',...
    'Color',[0 0 1]);

% ���� subplot
subplot2 = subplot(2,1,2,'Parent',figure1);
%% ȡ�������е�ע���Ա���������� X ��Χ
% xlim(subplot2,[0.001 0.004]);
%% ȡ�������е�ע���Ա���������� Y ��Χ
% ylim(subplot2,[0.96 0.97]);
box(subplot2,'on');
hold(subplot2,'on');

% ���� ylabel
ylabel('SSIM');

% ���� xlabel
xlabel('\lambda');

% ���� plot
plot(X1,Y2,'Parent',subplot2,'MarkerFaceColor',[0 0 0],'MarkerSize',3,...
    'Marker','o',...
    'LineWidth',1,...
    'LineStyle','-.',...
    'Color',[1 0 0]);


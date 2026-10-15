%%%%%%%%%%%%% analysis of real data WPBC %%%%%%%%%%%
clear;
%% 加载并标准化数

yhats = importdata('results/yhats1.txt');
errors = importdata('results/errors1.txt');
N_plot = 2000; % 画图点数
N_plot2 = 1500;

methods = 6;
%compute plot messages
[fit_Xbeta, err, t2, err_reg2, err_con2] = compute_plot(yhats, errors, N_plot, methods);

%compute plot messages
%[fit_Xbeta, err, t2, err_reg2, err_con2] = compute_plot(yhat_pro, yhat_CC, yhat_Im_E, yhat_Im_W, ...
%    err_pro, err_CC, err_Im_E, err_Im_W, N_plot);
x_r = 7;
figure
i = 1;
subplot(2,3,i);
plot(fit_Xbeta(:,i),err(:,i),'o',t2(:,i),err_reg2(:,i),'-r',t2(:,i),err_con2(:,i),'--r',t2(:,i),err_con2(:,i + methods),'--r', [-0, x_r],[0,0],'-b')
xlabel({'${X}^\top\widehat{\bf\beta}_{n}(U)$'},'Interpreter','latex','FontSize',20)
ylabel({'$\widehat{\varepsilon}_{n}(Y, Z)$'},'Interpreter','latex','FontSize',20)
title('(IPW)','FontSize',14)
axis([-0, x_r, -2, 2])

i = 2;
subplot(2,3,i);
plot(fit_Xbeta(:,i),err(:,i),'o',t2(:,i),err_reg2(:,i),'-r',t2(:,i),err_con2(:,i),'--r',t2(:,i),err_con2(:,i + methods),'--r', [-0, x_r],[0,0],'-b')
xlabel({'${X}^\top\widehat{\bf\beta}_{n}(U)$'},'Interpreter','latex','FontSize',20)
ylabel({'$\widehat{\varepsilon}_{n}(Y, Z)$'},'Interpreter','latex','FontSize',20)
title('(CC)','FontSize',14)
axis([-0, x_r, -2, 2])

i = 3;
subplot(2,3,i);
plot(fit_Xbeta(:,i),err(:,i),'o',t2(:,i),err_reg2(:,i),'-r',t2(:,i),err_con2(:,i),'--r',t2(:,i),err_con2(:,i + methods),'--r', [-0, x_r],[0,0],'-b')
xlabel({'${X}^\top\widehat{\bf\beta}_{n}(U)$'},'Interpreter','latex','FontSize',20)
ylabel({'$\widehat{\varepsilon}_{n}(Y, Z)$'},'Interpreter','latex','FontSize',20)
title('(Naive)','FontSize',14)
axis([-0, x_r, -2, 2])


i = 4;
subplot(2,3,i);
plot(fit_Xbeta(:,i),err(:,i),'o',t2(:,i),err_reg2(:,i),'-r',t2(:,i),err_con2(:,i),'--r',t2(:,i),err_con2(:,i + methods),'--r', [-0, x_r],[0,0],'-b')
xlabel({'${X}^\top\widehat{\bf\beta}_{n}(U)$'},'Interpreter','latex','FontSize',20)
ylabel({'$\widehat{\varepsilon}_{n}(Y, Z)$'},'Interpreter','latex','FontSize',20)
title('(CMI\_E)','FontSize',14)
axis([-0, x_r, -2, 2])

i = 5;
subplot(2,3,i);
plot(fit_Xbeta(:,i),err(:,i),'o',t2(:,i),err_reg2(:,i),'-r',t2(:,i),err_con2(:,i),'--r',t2(:,i),err_con2(:,i + methods),'--r', [-0, x_r],[0,0],'-b')
xlabel({'${X}^\top\widehat{\bf\beta}_{n}(U)$'},'Interpreter','latex','FontSize',20)
ylabel({'$\widehat{\varepsilon}_{n}(Y, Z)$'},'Interpreter','latex','FontSize',20)
title('(CMI\_W)','FontSize',14)
axis([-0, x_r, -2, 2])

i = 6;
subplot(2,3,i);
plot(fit_Xbeta(:,i),err(:,i),'o',t2(:,i),err_reg2(:,i),'-r',t2(:,i),err_con2(:,i),'--r',t2(:,i),err_con2(:,i + methods),'--r', [-0, x_r],[0,0],'-b')
xlabel({'${X}^\top\widehat{\bf\beta}_{n}(U)$'},'Interpreter','latex','FontSize',20)
ylabel({'$\widehat{\varepsilon}_{n}(Y, Z)$'},'Interpreter','latex','FontSize',20)
title('(Mean)','FontSize',14)
axis([-0, x_r, -2, 2])


%saveas(gcf, '~Body_Fat.eps')
saveas(gcf, '~/Desktop/20251227模型验证/仿真与实例_Win/Section4/3ExampleAnalysis/PBC/PBC', 'epsc')
% saveas(figure, 'MEL.eps')



function [y_fit, err, t, err_reg, err_con] = compute_plot(yhats, errors, N_plot, methods)
%%
% err_reg：残差的常数核估计，err_con：残差常数核估计的置信区间，N_plot: 作图点数 
%%%%%%%%%%%%%%%%%% 残差图的散点图及其置信区间 %%%%%%%%%%%%%%%%%%
n = length(yhats(:,1));
n_con1 = 5;
n_con2 = 5;

err = errors % [err_pro err_CC err_Im_E err_Im_W]; % 残差

y_fit = yhats %[yhat_pro yhat_CC yhat_Im_E yhat_Im_W]; % 模型拟合值
t = zeros(N_plot, 4); % 作图网格
err_reg = zeros(N_plot,4);
err_con = zeros(N_plot,8);
for m = 1:methods
    t(:,m) = linspace(min(y_fit(:,m)),max(y_fit(:,m)),N_plot)';
    t_rep = repmat(t(:, m),1,n);
    h_fit = 2.34*std(y_fit(:,m))*n^(-1/n_con1); % 窗宽，n的次数可调
    y_fit_rep = repmat(y_fit(:,m),1,N_plot);
    Khy_fit = Epkf((t_rep'-y_fit_rep)/h_fit);
    err_rep = repmat(err(:,m),1,N_plot);
    err_reg(:,m) = sum(Khy_fit.*err_rep)./sum(Khy_fit);
    sig_err = sqrt(0.6*sum(err_rep.^2.*Khy_fit)./sum(Khy_fit).^2);
    err_con(:,m) = err_reg(:,m)-1.96*sig_err';
    err_con(:,methods+m) = err_reg(:,m)+1.96*sig_err';
end

%%%%%%%%%%%%%%%%%% 非线性部分 %%%%%%%%%%%%%%%%
%g_reg = zeros(N_plot1, 4);
%u = linspace(min(U),max(U),N_plot1)'; % 作图网格
%h_u = 2.34*std(U)*n^(-1/n_con2); % 窗宽，n的次数可调
%u_rep = repmat(u,1,n);
%U_rep = repmat(U,1,N_plot1);
%KhU = Epkf((u_rep' - U_rep)/h_u);
%for m = 1:4
%%   g_rep = repmat(gu(:,m),1,N_plot1);
%    g_reg(:,m) = sum(KhU.*g_rep)./sum(KhU);
%end

end


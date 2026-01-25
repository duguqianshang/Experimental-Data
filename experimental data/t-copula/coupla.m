%--------------------------------------------------------------------------
% 使用 Copula 分析风速与浪高的相关性
%--------------------------------------------------------------------------
clc;
clear;
close all;
warning off;
%% 导入数据
% 假设数据为两个 Excel 文件，每个文件 365 行 24 列，每行为一天，每列为每小时的数据
WindSpeedData = xlsread('wind_speed.xlsx'); % 风速数据文件
WaveHeightData = xlsread('wave_height.xlsx'); % 浪高数据文件
% 展平数据为一维向量
WindSpeed = reshape(WindSpeedData, [], 1); % 风速数据展平
WaveHeight = reshape(WaveHeightData, [], 1); % 浪高数据展平

%% 核密度估计法计算概率密度函数
% 计算带宽 h 使用 Silverman 经验法则
n = length(WindSpeed); % 样本数量
sigma_wind = std(WindSpeed);
sigma_wave = std(WaveHeight);
h_wind = 1.06 * sigma_wind * n^(-1/5); % 风速带宽
h_wave = 1.06 * sigma_wave * n^(-1/5); % 浪高带宽
% 输出带宽
fprintf('风速核密度估计带宽: %.4f\n', h_wind);
fprintf('浪高核密度估计带宽: %.4f\n', h_wave);
% 采用高斯核计算核密度估计值
[f_wind, xi_wind] = ksdensity(WindSpeed, 'Bandwidth', h_wind, 'Kernel', 'normal');
[f_wave, xi_wave] = ksdensity(WaveHeight, 'Bandwidth', h_wave, 'Kernel', 'normal');
% 风速概率密度函数图,线的
figure;
plot(xi_wind, f_wind, 'b-', 'LineWidth', 2);
title('概率密度', 'FontName', 'Times New Roman', 'FontSize', 7.5);
xlabel('风速', 'FontName', 'Times New Roman', 'FontSize', 7.5);
ylabel('密度', 'FontName', 'Times New Roman', 'FontSize', 7.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 7.5); % 坐标轴刻度字体
grid on;
% 浪高概率密度函数图
figure;
plot(xi_wave, f_wave, 'r-', 'LineWidth', 2);
title('概率密度', 'FontName', 'Times New Roman', 'FontSize', 7.5);
xlabel('浪高', 'FontName', 'Times New Roman', 'FontSize', 7.5);
ylabel('密度', 'FontName', 'Times New Roman', 'FontSize', 7.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 7.5); % 坐标轴刻度字体
grid on;

% 计算边缘分布函数
[F_wind, x_wind] = ecdf(WindSpeed);
[F_wave, x_wave] = ecdf(WaveHeight);

% 解决采样点重复的问题
[x_wind, unique_idx_wind] = unique(x_wind);
F_wind = F_wind(unique_idx_wind);
[x_wave, unique_idx_wave] = unique(x_wave);
F_wave = F_wave(unique_idx_wave);

% 确保 F_wind 和 F_wave 严格在 (0, 1) 之间
F_wind = max(min(F_wind, 1 - eps), eps);
F_wave = max(min(F_wave, 1 - eps), eps);

% 绘制边缘分布函数，曲线的用不上
figure;
plot(x_wind, F_wind, 'b-', 'LineWidth', 2);
title('风速边缘分布函数');
xlabel('风速'); ylabel('F(x)');
grid on;

figure;
plot(x_wave, F_wave, 'r-', 'LineWidth', 2);
title('浪高边缘分布函数');
xlabel('浪高'); ylabel('F(y)');
grid on;

%% 计算 Kendall 和 Spearman 秩相关系数
kendall_corr = corr(WindSpeed, WaveHeight, 'type', 'Kendall');
spearman_corr = corr(WindSpeed, WaveHeight, 'type', 'Spearman');

fprintf('Kendall 秩相关系数: %.4f\n', kendall_corr);
fprintf('Spearman 秩相关系数: %.4f\n', spearman_corr);

%% Copula 参数估计
% 将边缘分布函数转换为 [0,1] 的数据
U = interp1(x_wind, F_wind, WindSpeed, 'linear', 'extrap');
V = interp1(x_wave, F_wave, WaveHeight, 'linear', 'extrap');

% 确保 U 和 V 严格在 (0, 1) 之间
U = max(min(U, 1 - eps), eps);
V = max(min(V, 1 - eps), eps);
%% 经验 Copula 分布函数
% 定义经验 Copula 函数
C_empirical = @(u, v) mean((U <= u) .* (V <= v));

% 计算经验 Copula 值
[Udata, Vdata] = meshgrid(linspace(0, 1, 50));
C_empirical_values = arrayfun(C_empirical, Udata, Vdata);

% 绘制经验 Copula 分布函数图
figure;
surf(Udata, Vdata, C_empirical_values, 'EdgeColor', 'black'); % 添加网格线
title('经验 Copula 分布函数');
xlabel('U'); ylabel('V'); zlabel('C(u,v)');
zlim([0, 1]); % 确保 z 轴范围在 [0, 1]
colormap jet; % 设置颜色图
grid on; % 启用网格


% 估计不同 Copula 的参数
rho_gaussian = copulafit('Gaussian', [U, V]); % 正态 Copula
[rho_t, nu_t] = copulafit('t', [U, V]); % t-Copula
param_gumbel = copulafit('Gumbel', [U, V]); % Gumbel-Copula
param_clayton = copulafit('Clayton', [U, V]); % Clayton-Copula
param_frank = copulafit('Frank', [U, V]); % Frank-Copula

%% 计算样本与 Copula 的平方欧式距离和相关系数
% 定义经验 Copula 函数
C_empirical = @(u, v) mean((U <= u) .* (V <= v));

% 网格点计算经验 Copula 值
[Udata, Vdata] = meshgrid(linspace(0, 1, 50));
C_empirical_values = arrayfun(C_empirical, Udata, Vdata);

% 正态 Copula
C_gaussian = reshape(copulacdf('Gaussian', [Udata(:), Vdata(:)], rho_gaussian), size(Udata));
distance_gaussian = sum((C_empirical_values - C_gaussian).^2, 'all');
kendall_gaussian = copulastat('Gaussian', rho_gaussian);
spearman_gaussian = copulastat('Gaussian', rho_gaussian, 'type', 'Spearman');

% t-Copula
C_t = reshape(copulacdf('t', [Udata(:), Vdata(:)], rho_t, nu_t), size(Udata));
distance_t = sum((C_empirical_values - C_t).^2, 'all');
kendall_t = copulastat('t', rho_t);
spearman_t = copulastat('t', rho_t, nu_t, 'type', 'Spearman');

% Gumbel-Copula
C_gumbel = reshape(copulacdf('Gumbel', [Udata(:), Vdata(:)], param_gumbel), size(Udata));
distance_gumbel = sum((C_empirical_values - C_gumbel).^2, 'all');
kendall_gumbel = copulastat('Gumbel', param_gumbel);
spearman_gumbel = copulastat('Gumbel', param_gumbel, 'type', 'Spearman');

% Clayton-Copula
C_clayton = reshape(copulacdf('Clayton', [Udata(:), Vdata(:)], param_clayton), size(Udata));
distance_clayton = sum((C_empirical_values - C_clayton).^2, 'all');
kendall_clayton = copulastat('Clayton', param_clayton);
spearman_clayton = copulastat('Clayton', param_clayton, 'type', 'Spearman');

% Frank-Copula
C_frank = reshape(copulacdf('Frank', [Udata(:), Vdata(:)], param_frank), size(Udata));
distance_frank = sum((C_empirical_values - C_frank).^2, 'all');
kendall_frank = copulastat('Frank', param_frank);
spearman_frank = copulastat('Frank', param_frank, 'type', 'Spearman');

% 输出平方欧式距离和相关系数
fprintf('正态 Copula 的平方欧式距离: %.4f, Kendall: %.4f, Spearman: %.4f\n', distance_gaussian, kendall_gaussian, spearman_gaussian);
fprintf('t-Copula 的平方欧式距离: %.4f, Kendall: %.4f, Spearman: %.4f\n', distance_t, kendall_t, spearman_t);
fprintf('Gumbel-Copula 的平方欧式距离: %.4f, Kendall: %.4f, Spearman: %.4f\n', distance_gumbel, kendall_gumbel, spearman_gumbel);
fprintf('Clayton-Copula 的平方欧式距离: %.4f, Kendall: %.4f, Spearman: %.4f\n', distance_clayton, kendall_clayton, spearman_clayton);
fprintf('Frank-Copula 的平方欧式距离: %.4f, Kendall: %.4f, Spearman: %.4f\n', distance_frank, kendall_frank, spearman_frank);

%% 构建表格并显示结果
% 添加样本的实际值到表格中
CopulaResults = table(...
    {'Sample'; 'Gaussian'; 't-Copula'; 'Gumbel'; 'Clayton'; 'Frank'}, ...
    [NaN; distance_gaussian; distance_t; distance_gumbel; distance_clayton; distance_frank], ...
    [kendall_corr; kendall_gaussian(1); kendall_t(1); kendall_gumbel(1); kendall_clayton(1); kendall_frank(1)], ...
    [spearman_corr; spearman_gaussian(1); spearman_t(1); spearman_gumbel(1); spearman_clayton(1); spearman_frank(1)], ...
    'VariableNames', {'Copula', 'Euclidean_Distance', 'Kendall_Coefficient', 'Spearman_Coefficient'});

% 显示表格结果
disp('各 Copula 函数的评估结果:');
disp(CopulaResults);


%% 绘制 Copula 概率密度函数和分布函数
% 正态 Copula
Cpdf_gaussian = copulapdf('Gaussian', [Udata(:), Vdata(:)], rho_gaussian);
Cpdf_gaussian = reshape(Cpdf_gaussian, size(Udata));
Cpdf_gaussian_normalized = Cpdf_gaussian / max(Cpdf_gaussian(:)); % 归一化
figure;
surf(Udata, Vdata, Cpdf_gaussian_normalized);
title('Gaussian-Copula 密度函数');
xlabel('U'); ylabel('V'); zlabel('c(u,v)');

% t-Copula
Cpdf_t = copulapdf('t', [Udata(:), Vdata(:)], rho_t, nu_t);
Cpdf_t = reshape(Cpdf_t, size(Udata));
Cpdf_t_normalized = Cpdf_t / max(Cpdf_t(:)); % 归一化
figure;
surf(Udata, Vdata, Cpdf_t_normalized);
title('t-Copula Density Function')
title('t-Copula概率密度函数','FontName', 'Times New Roman')
xlabel('风速'); ylabel('浪高'); zlabel('c(u,v)');
figure;
% 绘制三维表面图
surf(Udata, Vdata, Cpdf_t_normalized);
% % 设置标题
% title('t-Copula Density Function', ...
%       'FontName', 'Times New Roman', 'FontSize', 12);
% 设置坐标轴标签
xlabel('U', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('V', 'FontName', 'Times New Roman', 'FontSize', 12);
zlabel('c(u,v)', 'FontName', 'Times New Roman', 'FontSize', 14);
% 设置坐标轴刻度字体
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
% Gumbel-Copula
 Cpdf_gumbel = copulapdf('Gumbel', [Udata(:), Vdata(:)], param_gumbel);
Cpdf_gumbel = reshape(Cpdf_gumbel, size(Udata));
Cpdf_gumbel_normalized = Cpdf_gumbel / max(Cpdf_gumbel(:)); % 归一化
figure;
surf(Udata, Vdata, Cpdf_gumbel_normalized);
title('Gumbel-Copula 密度函数');
xlabel('风速'); ylabel('浪高'); zlabel('联合密度函数');

% Clayton-Copula
Cpdf_clayton = copulapdf('Clayton', [Udata(:), Vdata(:)], param_clayton);
Cpdf_clayton = reshape(Cpdf_clayton, size(Udata));
Cpdf_clayton_normalized = Cpdf_clayton / max(Cpdf_clayton(:)); % 归一化
figure;
surf(Udata, Vdata, Cpdf_clayton_normalized);
title('Clayton-Copula 密度函数');
xlabel('U'); ylabel('V'); zlabel('c(u,v)');
% 
% Frank-Copula
Cpdf_frank = copulapdf('Frank', [Udata(:), Vdata(:)], param_frank);
Cpdf_frank = reshape(Cpdf_frank, size(Udata));
Cpdf_frank_normalized = Cpdf_frank / max(Cpdf_frank(:)); % 归一化
figure;
surf(Udata, Vdata, Cpdf_frank_normalized);
title('Frank-Copula 密度函数');
% 设置坐标轴标签
xlabel('U', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('V', 'FontName', 'Times New Roman', 'FontSize', 12);
zlabel('c(u,v)', 'FontName', 'Times New Roman', 'FontSize', 14);
% 设置坐标轴刻度字体
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
% xlabel('U'); ylabel('V'); zlabel('c(u,v)');

%% 选择最佳 Copula
% 根据平方欧式距离选择最佳 Copula
[~, best_copula_idx] = min([distance_gaussian, distance_t, distance_gumbel, distance_clayton, distance_frank]);
CopulaNames = {'正态 Copula', 't-Copula', 'Gumbel-Copula', 'Clayton-Copula', 'Frank-Copula'};
fprintf('最佳 Copula 是: %s\n', CopulaNames{best_copula_idx});

%% 选择最佳 Copula
% 修改最佳 Copula 为 t-Copula
BestCopula = 't';
BestParam = rho_t;
BestNu = nu_t;
% 绘制 t-Copula 拟合后 U-V 散点图，并标注 Kendall 和 Spearman 系数
% figure('Position', [100, 100, 800, 500])
% scatter(U, V, 20, 'filled', 'MarkerFaceAlpha', 0.5)
% grid on
% xlabel('Transformed Wind Speed (U)', 'FontSize', 12)
% ylabel('Transformed Wave Height (V)', 'FontSize', 12)
% title({'t-Copula Transformed Space', ...
%        sprintf('Kendall''s \\tau = %.3f, Spearman''s \\rho = %.3f', kendall_t, spearman_t)}, ...
%        'FontSize', 12)
% set(gca, 'FontSize', 12)
% box on
% === 绘图并统一字体 Times New Roman ===
figure('Position', [100, 100, 800, 500])

% 绘制散点图
scatter(U, V, 20, 'filled', 'MarkerFaceAlpha', 0.5)
grid on
box on

% 设置坐标轴标签
xlabel('Transformed Wind Speed (U)', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('Transformed Wave Height (V)', 'FontSize', 12, 'FontName', 'Times New Roman')

% 设置标题
title({'t-Copula Transformed Space', ...
       sprintf('Kendall''s \\tau = %.3f, Spearman''s \\rho = %.3f', kendall_t, spearman_t)}, ...
       'FontSize', 12, 'FontName', 'Times New Roman')

% 设置坐标轴
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'LineWidth', 1)

% 如果有 colorbar（色标条）
if exist('colorbar_on', 'var') && colorbar_on
    h_cb = colorbar;
    set(h_cb, 'FontSize', 12, 'FontName', 'Times New Roman')
end

% 如果有 legend（图例）
if exist('legend_on', 'var') && legend_on
    h_lgd = legend('show');
    set(h_lgd, 'FontSize', 12, 'FontName', 'Times New Roman')
end

% %% 使用最佳 Copula 生成场景
% % 从 t-Copula 生成 500 个场景
% BestNu = ceil(nu_t); % 确保自由度为整数
% SimData = copularnd(BestCopula, BestParam, 500, BestNu); % 生成 500 个场景
% 
% % 逆变换采样，生成风速和浪高场景
% SimWind = interp1(F_wind, x_wind, SimData(:, 1), 'linear', 'extrap');
% SimWave = interp1(F_wave, x_wave, SimData(:, 2), 'linear', 'extrap');
% 
% % 绘制生成的场景
% figure;
% scatter(SimWind, SimWave, 10, 'filled');
% title('t-Copula 生成的场景');
% xlabel('风速');
% ylabel('浪高');
% grid on;
% 
% 
% %% 使用后向缩减法生成 5 个场景
% % 使用 K-means 聚类将场景缩减为 5 类
% k = 5; % 最终场景数
% [ClusterIdx, ClusterCenters] = kmeans([SimWind, SimWave], k);
% 
% % 绘制聚类结果
% figure;
% gscatter(SimWind, SimWave, ClusterIdx);
% hold on;
% scatter(ClusterCenters(:, 1), ClusterCenters(:, 2), 100, 'kx', 'LineWidth', 2);
% title('聚类结果与质心 (t-Copula)');
% xlabel('风速');
% ylabel('浪高');
% grid on;
% 
% % 输出每类场景的概率
% ClusterProbs = histcounts(ClusterIdx, k) / length(ClusterIdx);
% for i = 1:k
%     fprintf('场景 %d 的概率为 %.4f\n', i, ClusterProbs(i));
% end
%% 绘制 Copula 联合分布函数
% 正态 Copula 联合分布函数
C_gaussian_cdf = copulacdf('Gaussian', [Udata(:), Vdata(:)], rho_gaussian);
C_gaussian_cdf = reshape(C_gaussian_cdf, size(Udata));
figure;
surf(Udata, Vdata, C_gaussian_cdf);
title('Gaussian-Copula联合分布函数');
xlabel('U'); ylabel('V'); zlabel('C(u,v)');
zlim([0, 1]); % 确保 z 轴范围在 [0, 1]

% % % t-Copula 联合分布函数
C_t_cdf = copulacdf('t', [Udata(:), Vdata(:)], rho_t, nu_t);
C_t_cdf = reshape(C_t_cdf, size(Udata));
% figure;
% surf(Udata, Vdata, C_t_cdf);
% title('t-Copula 联合分布函数');
% xlabel('风速'); ylabel('浪高'); zlabel('C(u,v)');
% zlim([0, 1]);
% t-Copula 联合分布函数
% C_t_cdf = copulacdf('t', [Udata(:), Vdata(:)], rho_t, nu_t);
% C_t_cdf = reshape(C_t_cdf, size(Udata));
figure;
surf(Udata, Vdata, C_t_cdf);
% title('t-Copula Distribution Function', 'FontName', 'Times New Roman','FontSize', 8);
xlabel('U', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('V', 'FontName', 'Times New Roman', 'FontSize', 12);
zlabel('C(u,v)', 'FontName', 'Times New Roman', 'FontSize', 14);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
zlim([0, 1]);
set(gca, 'FontName', 'Times New Roman');  % 设置坐标轴字体


% Gumbel-Copula 联合分布函数
C_gumbel_cdf = copulacdf('Gumbel', [Udata(:), Vdata(:)], param_gumbel);
% C_gumbel_cdf = reshape(C_gumbel_cdf, size(Udata));
% figure;
% surf(Udata, Vdata, C_gumbel_cdf);
% title('Gumbel-Copula 联合分布函数');
% xlabel('U'); ylabel('V'); zlabel('C(u,v)');
% zlim([0, 1]);
% 
% Clayton-Copula 联合分布函数
C_clayton_cdf = copulacdf('Clayton', [Udata(:), Vdata(:)], param_clayton);
C_clayton_cdf = reshape(C_clayton_cdf, size(Udata));
figure;
surf(Udata, Vdata, C_clayton_cdf);
title('Clayton-Copula 联合分布函数');
xlabel('U'); ylabel('V'); zlabel('C(u,v)');
zlim([0, 1]);
% 
% Frank-Copula 联合分布函数
C_frank_cdf = copulacdf('Frank', [Udata(:), Vdata(:)], param_frank);
C_frank_cdf = reshape(C_frank_cdf, size(Udata));
figure;
surf(Udata, Vdata, C_frank_cdf);
xlabel('U', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('V', 'FontName', 'Times New Roman', 'FontSize', 12);
zlabel('C(u,v)', 'FontName', 'Times New Roman', 'FontSize', 14);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
zlim([0, 1]);
% set(gca, 'FontName', 'Times New Roman');  % 设置坐标轴字体
% title('Frank-Copula 联合分布函数');
% xlabel('U'); ylabel('V'); zlabel('C(u,v)');
% zlim([0, 1]);
%拟合场景生成
scenarionum = 500;  % 初始场景数目，可修改
num_cluster = 5;     % 要削减到的场景数目，可修改
ntime = 24;  % 24小时
% X和Y分别存储风和光的24个时刻历史观测数据
X = []; Y = [];
for t = 1 : ntime
    X{t} = WindSpeedData(:, t);
    Y{t} = WaveHeightData(:, t);
end
%% t-Copula拟合
% t-Copula 函数用于捕捉变量之间的相关性，特别是尾部相关性
for i = 1 : ntime
    % 核密度估计：计算风速和浪高的经验分布函数（CDF）
    U = ksdensity(X{i}, X{i}, 'function', 'cdf');
    V = ksdensity(Y{i}, Y{i}, 'function', 'cdf');
    
    % 使用 t-Copula 拟合相关性
    [rho, nu] = copulafit('t', [U(:), V(:)]); % rho为相关矩阵，nu为自由度
    copulaparams.rho = rho; % 保存相关矩阵
    copulaparams.nu = nu;   % 保存自由度
    copModels(i) = copulaparams;       
end
%% 从t-Copula中采样
Data = cell(1, ntime);

for j = 1: ntime
    % 使用 t-Copula 的 rho 和 nu 生成场景
    Data{j} = copularnd('t', copModels(j).rho, copModels(j).nu, scenarionum); 
end
% 逆变换，转换为实际场景
w = cell(1, ntime);
for k = 1 : ntime
    [fwind, xwind] = ecdf(X{k});
    funwind = @(pdwind)spline(fwind, xwind, pdwind);  % 样条插值，用于逆变换
    
    [fsolar, xsolar] = ecdf(Y{k});
    funsolar = @(pdsolar)spline(fsolar, xsolar, pdsolar);
    
    % 对风和浪进行逆变换
    for i = 1 : 2
        for j = 1 : scenarionum
            if i == 1
                w{k}(j, i) = funwind(Data{k}(j, i));  % 风速逆变换
            else
                w{k}(j,i) = funsolar(Data{k}(j, i));  % 浪高逆变换
            end
        end
    end
end
% 将生成的风速和浪高场景数据分别保存为 Excel 文件
% 定义文件名
WindSpeedFileName = 'Generated_WindSpeed.xlsx';
WaveHeightFileName = 'Generated_WaveHeight.xlsx';

% 初始化存储矩阵
GeneratedWindSpeed = zeros(scenarionum, ntime);
GeneratedWaveHeight = zeros(scenarionum, ntime);

% 填充生成的数据
for t = 1 : ntime
    GeneratedWindSpeed(:, t) = w{t}(:, 1); % 提取风速场景
    GeneratedWaveHeight(:, t) = w{t}(:, 2); % 提取浪高场景
end

% 将数据写入 Excel 文件
xlswrite(WindSpeedFileName, GeneratedWindSpeed);
xlswrite(WaveHeightFileName, GeneratedWaveHeight);

fprintf('风速和浪高场景数据已保存为 Excel 文件:\n - %s\n - %s\n', WindSpeedFileName, WaveHeightFileName);

%% 场景削减
wind = [];
solar = [];

for i = 1 : ntime
    wind = [wind, w{i}(:, 1)];
    solar = [solar, w{i}(:, 2)];
end
% K-means削减
opts = statset('Display', 'final');
[idx, C] = kmeans([wind, solar], num_cluster, 'Distance', 'cityblock', 'Replicates', 20, 'Options', opts);

C_wind = C(:, 1:24);     % 最终的风电场景
C_solar = C(:, 25:end);  % 最终的光伏场景
%% 概率计算
for i = 1 : num_cluster
    p(i) =  length(find(idx == i)) / length(idx); % 计算每个场景的概率
end
disp(['各场景概率为: ', num2str(p)]);

P_wd2 = zeros(num_cluster, 24);
for j = 1:num_cluster
    for i = 1:24
         P_wd2(j,i) = p(j) * C_wind(j,i);
    end
end
P_WD = sum(P_wd2, 1); % 合并各场景的风速概率

P_pv2 = zeros(num_cluster, 24);
for j = 1:num_cluster
    for i = 1:24
         P_pv2(j,i) = p(j) * C_solar(j,i);
    end
end
P_PV = sum(P_pv2, 1); % 合并各场景的光伏概率

%% 绘图
% 原始场景画图
% figure;
% [ss, gg] = meshgrid(1:scenarionum, 1:24);
% plot3(ss, gg, wind, 'linewidth', 1);
% title(['基于t-Copula生成的风速场景（', num2str(scenarionum), '个）']);
% xlabel('场景编号'); ylabel('时刻（小时）'); zlabel('风速值');
% 
% figure;
% plot3(ss, gg, solar, 'linewidth', 1);
% title(['基于t-Copula生成的浪高场景（', num2str(scenarionum), '个）']);
% xlabel('场景编号'); ylabel('时刻（小时）'); zlabel('浪高值');
% % 绘制削减后的风速场景（三维图）
% 
% % 绘制削减后的风速场景（三维图）：以下内容未改变过
% figure;
% hold on;
% styles = {'--', '-', ':', '-.', '--'}; % 不同场景的线条样式
% colors = lines(num_cluster); % 不同场景的颜色
% [timeGrid, sceneGrid] = meshgrid(1:24, 1:num_cluster);
% 
% for i = 1:num_cluster
%     plot3(1:24, i * ones(1, 24), C_wind(i, :), 'LineStyle', styles{i}, ...
%           'Color', colors(i, :), 'LineWidth', 1.5);
% end
% legendLabels = cell(num_cluster, 1);
% for i = 1:num_cluster
%     legendLabels{i} = sprintf('场景%d (%.3f)', i, p(i)); % 图例显示概率值
% end
% legend(legendLabels, 'Location', 'best', 'FontSize', 10);
% xlabel('时间/h');
% ylabel('场景编号');
% zlabel('风电出力/p.u.');
% title('削减后的风速场景（三维）及其概率');
% grid on;
% view(45, 30); % 设置三维视角
% 
% % 绘制削减后的浪高场景（三维图）
% figure;
% hold on;
% for i = 1:num_cluster
%     plot3(1:24, i * ones(1, 24), C_solar(i, :), 'LineStyle', styles{i}, ...
%           'Color', colors(i, :), 'LineWidth', 1.5);
% end
% legendLabels = cell(num_cluster, 1);
% for i = 1:num_cluster
%     legendLabels{i} = sprintf('场景%d (%.3f)', i, p(i)); % 图例显示概率值
% end
% legend(legendLabels, 'Location', 'best', 'FontSize', 10);
% xlabel('时间/h');
% ylabel('场景编号');
% zlabel('浪高出力/p.u.');
% title('削减后的浪高场景（三维）及其概率');
% grid on;
% view(45, 30); % 设置三维视角
% 
% 
% % 场景概率分布
% figure;
% bar(p);
% xlabel('场景编号');
% ylabel('概率');
% title('削减后场景的概率分布');
% 
% % 风速不确定性
% figure;
% plot(P_WD, 'LineWidth', 1.5);
% xlabel('时间（小时）');
% ylabel('风速（m/s）');
% title('风速不确定性');
% 
% % 浪高不确定性
% figure;
% plot(P_PV, 'LineWidth', 1.5);
% xlabel('时间（小时）');
% ylabel('浪高（m）');
% title('浪高不确定性');
% 拟合部分
% %% 1. 风速频率分布直方图及拟合
% figure;
% hold on;
% % 绘制频率直方图（蓝色填充）
% histogram(WindSpeed, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.6, 1], 'EdgeColor', 'black');
% title('风速频率分布直方图及拟合效果');
% xlabel('风速 (m/s)');
% ylabel('概率密度');
% grid on;
% 
% % 核密度估计
% [f_kde, x_kde] = ksdensity(WindSpeed);
% plot(x_kde, f_kde, 'r-', 'LineWidth', 2, 'DisplayName', '核密度估计');
% 
% % 威布尔分布拟合
% paramWeibull = wblfit(WindSpeed); % 威布尔参数估计
% f_weibull = wblpdf(x_kde, paramWeibull(1), paramWeibull(2));
% plot(x_kde, f_weibull, 'y--', 'LineWidth', 2, 'DisplayName', '威布尔分布');
% 
% % 正态分布拟合
% mu = mean(WindSpeed);
% sigma = std(WindSpeed);
% f_normal = normpdf(x_kde, mu, sigma);
% plot(x_kde, f_normal, 'm-.', 'LineWidth', 2, 'DisplayName', '正态分布');
% 
% legend('频率直方图', '核密度估计', '威布尔分布', '正态分布', 'Location', 'NorthEast');
% hold off;
% 
% %% 2. 浪高频率分布直方图及拟合
% figure;
% hold on;
% % 绘制频率直方图（蓝色填充）
% histogram(WaveHeight, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.6, 1], 'EdgeColor', 'black');
% title('浪高频率分布直方图及拟合效果');
% xlabel('浪高 (m)');
% ylabel('概率密度');
% grid on;
% 
% % 核密度估计
% [f_kde_wave, x_kde_wave] = ksdensity(WaveHeight);
% plot(x_kde_wave, f_kde_wave, 'r-', 'LineWidth', 2, 'DisplayName', '核密度估计');
% 
% % 对数正态分布拟合
% paramLogNormal = lognfit(WaveHeight); % 对数正态参数估计
% f_logNormal = lognpdf(x_kde_wave, paramLogNormal(1), paramLogNormal(2));
% plot(x_kde_wave, f_logNormal, 'g--', 'LineWidth', 2, 'DisplayName', '对数正态分布');
% 
% legend('频率直方图', '核密度估计', '对数正态分布', 'Location', 'NorthEast');
% hold off;
%% 1. 风速频率分布直方图及拟合
figure;
hold on;
% 频率直方图
histogram(WindSpeed, 'Normalization', 'pdf', ...
    'FaceColor', [0.4, 0.7, 1], 'EdgeColor', 'black');
% 核密度估计
[f_kde, x_kde] = ksdensity(WindSpeed);
plot(x_kde, f_kde, 'r-', 'LineWidth', 1.8);
% 威布尔分布拟合
paramWeibull = wblfit(WindSpeed);
f_weibull = wblpdf(x_kde, paramWeibull(1), paramWeibull(2));
plot(x_kde, f_weibull, 'color', [1 0.5 0], 'LineStyle', '--', 'LineWidth', 1.8);
% 正态分布拟合
mu = mean(WindSpeed);
sigma = std(WindSpeed);
f_normal = normpdf(x_kde, mu, sigma);
plot(x_kde, f_normal, 'm-.', 'LineWidth', 1.8);
% 设置坐标轴标签
% xlabel('{\itv}/m{\cdot}s^{-1}', ...
%     'FontName', 'Times New Roman', ...  % 英文字母使用 Times New Roman
%     'FontSize', 8);
% ylabel('概率密度', ...
%     'FontName', '宋体', ...              % 中文使用宋体
%     'FontSize', 8);
% % 设置标题（中文）
% title('风速频率分布直方图及拟合', ...
%     'FontName', '宋体', ...
%     'FontSize', 8);
% % 设置图例（中文）
% legend({'频率直方图', '核密度估计', '威布尔分布', '正态分布'}, ...
%     'Location', 'northeast', ...
%     'FontName', '宋体', ...
%     'FontSize', 8);
% % 设置坐标轴刻度字体（新罗马）
% set(gca, ...
%     'FontName', 'Times New Roman', ...
%     'FontSize', 8);
% % 其他图形美化
% grid on;
% box on;
% hold off;
% 设置坐标轴标签（英文）
xlabel('{\itv}/m{\cdot}s^{-1}', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);
ylabel('Probability Density', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);

% 设置标题（英文）
% title('Wind Speed Histogram and Fitting Curves', ...
%     'FontName', 'Times New Roman', ...
%     'FontSize', 9);
% 设置图例（全部英文 + 新罗马字体）
legend({'Frequency Histogram', 'Kernel Density Estimation', ...
        'Weibull Distribution', 'Normal Distribution'}, ...
    'Location', 'northeast', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);
% 设置坐标轴刻度字体
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);
% 图形美化
grid on;
box on;
hold off;

%% 2. 浪高频率分布直方图及拟合
figure;
hold on;
% 直方图
histogram(WaveHeight, 'Normalization', 'pdf', ...
    'FaceColor', [0.4, 0.7, 1], 'EdgeColor', 'black');
% 核密度估计
[f_kde_wave, x_kde_wave] = ksdensity(WaveHeight);
plot(x_kde_wave, f_kde_wave, 'r-', 'LineWidth', 1.8);
% 对数正态拟合
paramLogNormal = lognfit(WaveHeight);
f_logNormal = lognpdf(x_kde_wave, paramLogNormal(1), paramLogNormal(2));
plot(x_kde_wave, f_logNormal, 'g--', 'LineWidth', 1.8);
% xlabel('{\ith}/m',...
%     'FontName', 'Times New Roman',...  % 变量斜体+单位用新罗马
%     'FontSize', 7.5);
% ylabel('概率密度',...
%     'FontName', 'SimSun',...           % 中文标签用宋体
%     'FontSize', 7.5);
% % 图标题
% title('浪高频率分布直方图及拟合',...
%     'FontName', 'SimSun-ExtB',...           % 中文标题用宋体
%     'FontSize', 7.5);
% % 图例（全部内容为中文）
% legend({'频率直方图', '核密度估计', '对数正态分布'},...
%     'Location', 'northeast',...
%     'FontName', 'SimSun-ExtB',...           % 图例中文用宋体
%     'FontSize', 7.5);
% % 统一设置坐标轴刻度字体（数字/英文用新罗马）
% set(gca,...
%     'FontName', 'Times New Roman',...
%     'FontSize', 7.5);
% grid on;
% box on;
% hold off;
% Set axis labels (italic variable + Roman unit)
xlabel('{\ith}/m', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);

ylabel('Probability Density', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);

% Set title (English)
% title('Wave Height Histogram and Fitting Curves', ...
%     'FontName', 'Times New Roman', ...
%     'FontSize', 9);

% Legend (all in English)
legend({'Frequency Histogram', 'Kernel Density Estimation', ...
        'Lognormal Distribution'}, ...
    'Location', 'northeast', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);

% Axis tick labels (numbers + text)
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 9);

% Plot beautification
grid on;
box on;
hold off;



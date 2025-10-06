oilp = readtable('MCOILBRENTEU.xlsx');
indprod = readtable('IPB50001N.xlsx');
gecon = readtable('GECON_indicator.xlsx');
funds = readtable('FEDFUNDS.xlsx');
cpi = readtable("CPIAUCSL.xlsx");


t_oil   = oilp{:,1};   y_oil   = oilp{:,2};
t_ip    = indprod{:,1}; y_ip   = indprod{:,2};
t_gecon = gecon{:,1};   y_gecon = gecon{:,2};
t_funds = funds{:,1};   y_funds = funds{:,2};
t_cpi   = cpi{:,1};     y_cpi   = cpi{:,2};

figure;
tiledlayout(3,2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t_oil, y_oil, 'LineWidth', 1.5);
title('Nominal Oil Price'); grid on;
nexttile;
plot(t_ip, y_ip, 'LineWidth', 1.5);
title('Industrial Production'); grid on;
nexttile;
plot(t_cpi, y_cpi, 'LineWidth', 1.5);
title('CPI (Price Index)'); grid on;
nexttile;
plot(t_funds, y_funds, 'LineWidth', 1.5);
title('Federal Funds Rate'); grid on;
nexttile;
plot(t_gecon, y_gecon, 'LineWidth', 1.5);
title('GECON Indicator'); grid on;
sgtitle('Raw Data Series');

real_oil = oilp{:,2} ./ cpi{:,2};
oilp.RealOil = real_oil;
%Log variable
oilp.Log = log(oilp{:,2}); 
oilp.Real = log(oilp{ :,3});  
indprod.Log = log(indprod{:,2});    
cpi.Log = log(cpi{:,2});        
%ADF test
real_oil_log = oilp.Real;
indprod_log  = indprod.Log;

real_oil_log = real_oil_log(~isnan(real_oil_log));
indprod_log  = indprod_log(~isnan(indprod_log));

[h_oil, p_oil, stat_oil, cValue_oil] = adftest(real_oil_log);
[h_ip, p_ip, stat_ip, cValue_ip] = adftest(indprod_log);

fprintf('ADF Test for Real Oil Price:\n');
if h_oil == 1
    fprintf('  Result: Stationary (rejects unit root)\n');
else
    fprintf('  Result: Nonstationary (fails to reject unit root)\n');
end
fprintf('  Test statistic: %.4f\n', stat_oil);
fprintf('  5%% critical value: %.4f\n', cValue_oil);
fprintf('  p-value: %.4f\n\n', p_oil);
fprintf('ADF Test for Industrial Production:\n');
if h_ip == 1
    fprintf('  Result: Stationary (rejects unit root)\n');
else
    fprintf('  Result: Nonstationary (fails to reject unit root)\n');
end
fprintf('  Test statistic: %.4f\n', stat_ip);
fprintf('  5%% critical value: %.4f\n', cValue_ip);
fprintf('  p-value: %.4f\n', p_ip);

%% HP filter sulla Industrial Production (log)
mask = ~isnan(indprod.Log);
t_ip_clean   = t_ip(mask);
ip_log_clean = indprod.Log(mask);

lambda_monthly = 129600;

[hp_trend, hp_cycle] = hpfilter_local(ip_log_clean, lambda_monthly);

indprod.HP_Trend = NaN(size(indprod.Log));
indprod.HP_Cycle = NaN(size(indprod.Log));
indprod.HP_Trend(mask) = hp_trend;
indprod.HP_Cycle(mask) = hp_cycle;

figure;
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(t_ip, indprod.Log, 'LineWidth', 1.2); hold on;
plot(t_ip_clean, hp_trend, 'LineWidth', 1.8);
title(sprintf('Industrial Production (log) e HP-Trend (\\lambda = %d)', lambda_monthly));
legend('Log IP','HP Trend','Location','best'); grid on; hold off;

nexttile;
plot(t_ip_clean, hp_cycle, 'LineWidth', 1.5); yline(0,'k-'); grid on;
title('Ciclo HP della Industrial Production'); xlabel('Data');

lambdas = [50000, 129600, 300000];
cycles = NaN(numel(ip_log_clean), numel(lambdas));
for j = 1:numel(lambdas)
    [tr_j, cy_j] = hpfilter_local(ip_log_clean, lambdas(j));
    cycles(:,j) = cy_j;
end

figure;
plot(t_ip_clean, cycles, 'LineWidth', 1.2); yline(0,'k-'); grid on;
legend(arrayfun(@(x) sprintf('\\lambda = %d',x), lambdas,'UniformOutput',false),'Location','best');
title('Ciclo HP per diverse scelte di \\lambda');

std_cycles = std(cycles,0,1);
corr_cycles = corr(cycles,'rows','pairwise');
fprintf('\n=== Sensibilità a lambda (HP ciclo su IP log) ===\n');
for j = 1:numel(lambdas)
    fprintf('lambda = %6d | std(cycle) = %.4f\n', lambdas(j), std_cycles(j));
end
fprintf('\nCorrelazioni tra cicli (ordine dei lambda):\n');
disp(corr_cycles);

function [trend, cycle] = hpfilter_local(y, lambda)
    y = y(:);
    T = numel(y);
    I = speye(T);
    e = ones(T,1);
    D = spdiags([e -2*e e], 0:2, T-2, T);  
    trend = (I + lambda*(D'*D)) \ y;
    cycle = y - trend;
end

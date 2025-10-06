oilp = readtable('MCOILBRENTEU.xlsx');
indprod = readtable('IPB50001N.xlsx');
gecon = readtable('GECON_indicator.xlsx');
funds = readtable('FEDFUNDS.xlsx');
cpi = readtable("CPIAUCSL.xlsx");
Plots
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
Real oil price
real_oil = oilp{:,2} ./ cpi{:,2};
oilp.RealOil = real_oil;
Log variable
oilp.Log = log(oilp{:,2}); 
oilp.Real = log(oilp{ :,3});  
indprod.Log = log(indprod{:,2});    
cpi.Log = log(cpi{:,2});        
ADF test
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

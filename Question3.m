% Question 3

TT_gecon = timetable(t_gecon, y_gecon, 'VariableNames', {'GECON'});
TT_ip    = timetable(t_ip   , log(y_ip) , 'VariableNames', {'IP_log'});         

try
    log_real_oil = oilp.Real; 
    TT_oil = timetable(t_oil, log_real_oil, 'VariableNames', {'Oil_logReal'});
catch
    real_oil = y_oil ./ y_cpi;
    TT_oil = timetable(t_oil, log(real_oil), 'VariableNames', {'Oil_logReal'});
end
TT_ff    = timetable(t_funds, y_funds, 'VariableNames', {'FFR'});                

if all(TT_gecon.GECON > 0, 'omitnan')
    TT_gecon.GECON = log(TT_gecon.GECON);
    gecon_name = 'GECON_log';
else
    gecon_name = 'GECON_level';
end
TT_gecon.Properties.VariableNames = {gecon_name};

TTall = synchronize(TT_gecon, TT_oil, TT_ip, TT_ff, 'intersection');

Y = TTall{:, [1 2 3 4]};     
dates = TTall.Properties.RowTimes;
varNames = {gecon_name, 'Oil_logReal', 'IP_log', 'FFR'};

mask = all(~isnan(Y),2);
Y = Y(mask,:);
dates = dates(mask);

maxLag = 12;
K = size(Y,2);
aicV = NaN(maxLag,1); bicV = NaN(maxLag,1); hqV = NaN(maxLag,1);

for p = 1:maxLag
    Mdl = varm(K, p);
    try
        [EstMdl, ~, logL, info] = estimate(Mdl, Y, 'Y0', Y(1:p,:)); 
        kparams = info.NumEstimatedParameters;
        T = size(Y,1);
        [aicV(p), bicV(p)] = aicbic(logL, kparams, T);
        hqV(p) = -2*logL + 2*kparams*log(log(T));
    catch
            continue
    end
end

[~, p_bic] = min(bicV);

fprintf('\n=== Lag-length selection (1..%d) ===\n', maxLag);
disp(table((1:maxLag)', aicV, bicV, hqV, 'VariableNames', {'p','AIC','BIC','HQ'}));
fprintf('-> Scelgo p = %d (BIC minimo)\n', p_bic);

p = p_bic;
Mdl = varm(K, p);
[EstMdl, ~, ~, info] = estimate(Mdl, Y, 'Y0', Y(1:p,:)); 

[~, ~, E] = infer(EstMdl, Y);   
SigmaU = cov(E,1);              

fprintf('\n=== Residual autocorrelation (Ljung–Box, up to 12 lags) ===\n');
for i = 1:K
    [h,pval] = lbqtest(E(:,i), 'Lags', 12);
    fprintf('Var %d (%s): h=%d, p=%.4f\n', i, varNames{i}, h, pval);
end

P = chol(SigmaU, 'lower');     

horizon = 36;  
IRF = irf(EstMdl, horizon);    

figure('Name','SVAR IRFs (Cholesky, ordering: GECON, Oil, IP, FFR)');
tiledlayout(K, K, 'TileSpacing','compact','Padding','compact');
for j = 1:K              
    for i = 1:K       
        nexttile;
        plot(0:horizon, squeeze(IRF(i,j,:)), 'LineWidth', 1.25); grid on;
        yline(0,'k-');
        title(sprintf('Resp %s a shock %s', varNames{i}, varNames{j}));
        xlabel('mesi'); 
    end
end

disp(EstMdl)
fprintf('\nOrdering (recursive):\n1) %s  2) %s  3) %s  4) %s\n', varNames{1},varNames{2},varNames{3},varNames{4});
fprintf('Lag scelto p = %d (BIC). Horizon IRF = %d mesi.\n', p, horizon);

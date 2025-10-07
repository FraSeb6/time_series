%% Q3 — SVAR (Monthly). Ordering: (1) Global activity, (2) Oil (real, log), (3) IP (log), (4) FFR (level)

clc; clear; close all;

%% 1) Read data
oilp   = readtable('MCOILBRENTEU.xlsx');
indprod= readtable('IPB50001N.xlsx');
gecon  = readtable('GECON_indicator.xlsx');
funds  = readtable('FEDFUNDS.xlsx');
cpi    = readtable('CPIAUCSL.xlsx');

get_ty = @(T) deal(local_to_datetime(T{:,1}), double(T{:,2}));
[t_oil,   y_oil]   = get_ty(oilp);
[t_ip,    y_ip]    = get_ty(indprod);
[t_gecon, y_gecon] = get_ty(gecon);
[t_funds, y_funds] = get_ty(funds);
[t_cpi,   y_cpi]   = get_ty(cpi);

%% 2) Transformations
real_oil     = y_oil ./ y_cpi;
log_real_oil = log(real_oil);
ip_log       = log(y_ip);
y_gecon = y_gecon(:);
if all(y_gecon(~isnan(y_gecon)) > 0)
    gecon_series = log(y_gecon); gecon_name = 'GECON_log';
else
    gecon_series = y_gecon;      gecon_name = 'GECON_level';
end
ffr_level = y_funds;

%% 3) Timetables & align (intersection). Order: GECON, OIL, IP, FFR
TT_gecon = timetable(t_gecon, gecon_series, 'VariableNames', {gecon_name});
TT_oil   = timetable(t_oil,   log_real_oil,  'VariableNames', {'Oil_logReal'});
TT_ip    = timetable(t_ip,    ip_log,        'VariableNames', {'IP_log'});
TT_ff    = timetable(t_funds, ffr_level,     'VariableNames', {'FFR'});
TTall    = synchronize(TT_gecon, TT_oil, TT_ip, TT_ff, 'intersection');
Y        = TTall{:,:};
varNames = TTall.Properties.VariableNames;

%% 4) Lag selection robusta (p=1..12) — AIC/BIC/HQ calcolati "a mano"
maxLag = 12; K = size(Y,2);
aicV = NaN(maxLag,1); bicV = NaN(maxLag,1); hqV = NaN(maxLag,1);
for p = 1:maxLag
    try
        Mdl_p = varm(K, p);
        Est_p = estimate(Mdl_p, Y, 'Y0', Y(1:p,:));
        E_p   = infer(Est_p, Y);
        Teff  = size(E_p,1);
        Sigma = (E_p' * E_p) / Teff;
        [U,pd] = chol(Sigma); if pd>0, continue; end
        logdetSigma = 2*sum(log(diag(U)));
        logL = -(Teff*K/2)*log(2*pi) - (Teff/2)*logdetSigma - (Teff*K/2);
        kparams = K*(K*p + 1) + K*(K+1)/2;
        aicV(p) = -2*logL + 2*kparams;
        bicV(p) = -2*logL + kparams*log(Teff);
        hqV(p)  = -2*logL + 2*kparams*log(log(Teff));
    catch
    end
end
valid = find(~isnan(bicV));
if isempty(valid), warning('Criteri non calcolabili: imposto p=1.'); p_bic = 1;
else, [~,ix] = min(bicV(valid)); p_bic = valid(ix); end
fprintf('\nLag selection (1..%d)\n', maxLag);
disp(table((1:maxLag)', aicV, bicV, hqV, 'VariableNames', {'p','AIC','BIC','HQ'}));
fprintf('Selected p = %d (BIC)\n', p_bic);

%% 5) Estimate VAR(p), stability, residual LB
p = p_bic; Mdl = varm(K, p);
[EstMdl,~,~,~] = estimate(Mdl, Y, 'Y0', Y(1:p,:));
[stableFlag, maxRoot] = local_isStable_var(EstMdl);
fprintf('\nStability (all eigenvalues < 1): %d | max |root|=%.6f\n', stableFlag, maxRoot);
E = infer(EstMdl, Y);
SigmaU = cov(E,1);
fprintf('\nResidual Ljung–Box (L=12)\n');
for i = 1:K
    [h,pval] = lbqtest(E(:,i), 'Lags', 12);
    fprintf('%s: h=%d, p=%.4f\n', varNames{i}, h, pval);
end

%% 6) IRF (Cholesky) — full grid
H = 36;
[IRF, tvec] = local_irf_chol(EstMdl, SigmaU, H);
figure('Name','SVAR IRFs (Cholesky ordering)'); tiledlayout(K,K,'TileSpacing','compact','Padding','compact');
for j = 1:K
    for i = 1:K
        nexttile; plot(tvec, squeeze(IRF(i,j,:)), 'LineWidth', 1.25);
        yline(0,'k-'); grid on;
        title(sprintf('Resp %s to shock %s', varNames{i}, varNames{j}), 'Interpreter','none'); xlabel('months');
    end
end

%% 7) IRF richieste: Global shock (1) e Oil shock (2)
shock = 1;
figure('Name','IRFs — Global activity shock'); tiledlayout(K,1,'TileSpacing','compact','Padding','compact');
for i = 1:K
    nexttile; plot(tvec, squeeze(IRF(i,shock,:)), 'LineWidth', 1.5);
    yline(0,'k-'); grid on; title(sprintf('%s shock → %s', varNames{shock}, varNames{i}), 'Interpreter','none');
end; xlabel('months');

shock = 2;
figure('Name','IRFs — Oil price shock'); tiledlayout(K,1,'TileSpacing','compact','Padding','compact');
for i = 1:K
    nexttile; plot(tvec, squeeze(IRF(i,shock,:)), 'LineWidth', 1.5);
    yline(0,'k-'); grid on; title(sprintf('%s shock → %s', varNames{shock}, varNames{i}), 'Interpreter','none');
end; xlabel('months');

%% 8) IRF — Monetary policy shock (FFR) → Industrial Production (IP_log)
shockIdx = find(strcmp(varNames,'FFR'), 1);     if isempty(shockIdx), shockIdx = 4; end
respIdx  = find(strcmp(varNames,'IP_log'), 1);  if isempty(respIdx),  respIdx  = 3; end
ip_irf = squeeze(IRF(respIdx, shockIdx, :));
figure('Name','IRF: Monetary policy shock on Industrial Production');
plot(tvec, ip_irf, 'LineWidth', 1.8); grid on; yline(0,'k-');
xlabel('months'); ylabel('log points'); title('FFR shock → IP');

%% 9) Summary
disp(EstMdl)
fprintf('\nOrdering:\n  1) %s  2) %s  3) %s  4) %s\n', varNames{1},varNames{2},varNames{3},varNames{4});
fprintf('Lag p = %d. IRF horizon = %d months.\n', p, H);

%% ===== Local functions =====
function dt = local_to_datetime(x)
    if isdatetime(x), dt = x; return; end
    if isstring(x) || iscellstr(x) || ischar(x)
        try, dt = datetime(x,'InputFormat','yyyy-MM-dd'); return; end
        dt = datetime(x); return
    end
    if isnumeric(x)
        try, dt = datetime(x,'ConvertFrom','excel'); return; end
        dt = datetime(1899,12,30) + days(x); return
    end
    error('Unrecognized date format.');
end

function [isStableFlag, maxRoot] = local_isStable_var(M)
    K = M.NumSeries; p = M.P;
    A = zeros(K*p);
    for i = 1:p
        Ai = M.AR{i}; if isempty(Ai), Ai = zeros(K); end
        A(1:K,(i-1)*K+1:i*K) = Ai;
    end
    if p > 1, A(K+1:end,1:K*(p-1)) = eye(K*(p-1)); end
    eigvals = eig(A); maxRoot = max(abs(eigvals)); isStableFlag = maxRoot < 1 - 1e-8;
end

function [IRF, t] = local_irf_chol(M, SigmaU, H)
    K = M.NumSeries; p = M.P;
    A = cell(p,1);
    for i = 1:p, Ai = M.AR{i}; if isempty(Ai), Ai = zeros(K); end, A{i} = Ai; end
    P = chol(SigmaU, 'lower');
    Theta = cell(H+1,1); Theta{1} = eye(K);
    for h = 1:H
        S = zeros(K);
        for i = 1:min(h,p), S = S + A{i} * Theta{h-i+1}; end
        Theta{h+1} = S;
    end
    IRF = zeros(K,K,H+1);
    for h = 0:H, IRF(:,:,h+1) = Theta{h+1} * P; end
    t = 0:H;
end

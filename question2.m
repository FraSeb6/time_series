%% Local: datetime parser
function dt = local_to_datetime(x)
    if isdatetime(x), dt = x; return; end
    if isstring(x) || iscellstr(x) || ischar(x)
        try, dt = datetime(x, 'InputFormat','yyyy-MM-dd'); return; end
        dt = datetime(x); return
    end
    if isnumeric(x)
        try, dt = datetime(x, 'ConvertFrom','excel'); return; end
        dt = datetime(1899,12,30) + days(x); return
    end
    error('Unrecognized date format.');
end

function [isStableFlag, maxRoot] = local_isStable_var(M)
    % Stability of VAR: all eigenvalues of the companion matrix inside unit circle
    K = M.NumSeries;
    p = M.P;
    % Build companion matrix
    A = zeros(K*p);
    for i = 1:p
        Ai = M.AR{i};
        if isempty(Ai), Ai = zeros(K); end
        A(1:K, (i-1)*K+1:i*K) = Ai;
    end
    if p > 1
        A(K+1:end, 1:K*(p-1)) = eye(K*(p-1));
    end
    eigvals = eig(A);
    maxRoot = max(abs(eigvals));
    isStableFlag = maxRoot < 1 - 1e-8;
end

function [IRF, t] = local_irf_chol(M, SigmaU, H)
    % IRF ortogonalizzate per VAR(p) stimato (M) con Cholesky su SigmaU
    % Output: IRF (K×K×(H+1)), t = 0:H
    K = M.NumSeries;
    p = M.P;

    % Matrici AR{1..p}
    A = cell(p,1);
    for i = 1:p
        Ai = M.AR{i};
        if isempty(Ai), Ai = zeros(K); end
        A{i} = Ai;
    end

    % Decomposizione Cholesky (impatto contemporaneo)
    P = chol(SigmaU, 'lower');

    % Ricorsione per le matrici di risposta ai movimenti ridotti-forma
    % Theta_0 = I_K ; Theta_h = sum_{i=1..min(h,p)} A_i * Theta_{h-i}
    Theta = cell(H+1,1);
    Theta{1} = eye(K);
    for h = 1:H
        S = zeros(K);
        for i = 1:min(h,p)
            S = S + A{i} * Theta{h - i + 1};
        end
        Theta{h+1} = S;
    end

    % IRF strutturali: Psi_h = Theta_h * P
    IRF = zeros(K, K, H+1);
    for h = 0:H
        IRF(:,:,h+1) = Theta{h+1} * P;
    end

    t = 0:H;
end



%% SVAR (Monthly) — Ordering: (1) Global activity, (2) Oil (real, log), (3) IP (log), (4) FFR (level)

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

y_gecon = y_gecon(:);  % assicurati sia vettore colonna
if all( y_gecon(~isnan(y_gecon)) > 0 )
    gecon_series = log(y_gecon); 
    gecon_name   = 'GECON_log';
else
    gecon_series = y_gecon;      
    gecon_name   = 'GECON_level';
end

ffr_level = y_funds;

%% 3) Timetables & align
TT_gecon = timetable(t_gecon, gecon_series, 'VariableNames', {gecon_name});
TT_oil   = timetable(t_oil,   log_real_oil,  'VariableNames', {'Oil_logReal'});
TT_ip    = timetable(t_ip,    ip_log,        'VariableNames', {'IP_log'});
TT_ff    = timetable(t_funds, ffr_level,     'VariableNames', {'FFR'});

TTall   = synchronize(TT_gecon, TT_oil, TT_ip, TT_ff, 'intersection');
Y       = TTall{:, :};
dates   = TTall.Properties.RowTimes; %#ok<NASGU>
varNames= TTall.Properties.VariableNames;

%% 4) Lag selection (p=1..12)
maxLag = 12; K = size(Y,2);
aicV = NaN(maxLag,1); bicV = NaN(maxLag,1); hqV = NaN(maxLag,1);

for p = 1:maxLag
    Mdl = varm(K, p);
    try
        [~,~,logL,info] = estimate(Mdl, Y, 'Y0', Y(1:p,:));
        kparams = info.NumEstimatedParameters; T = size(Y,1);
        [aicV(p), bicV(p)] = aicbic(logL, kparams, T);
        hqV(p) = -2*logL + 2*kparams*log(log(T));
    catch
    end
end
[~, p_bic] = min(bicV);
fprintf('\nLag selection (1..%d)\n', maxLag);
disp(table((1:maxLag)', aicV, bicV, hqV, 'VariableNames', {'p','AIC','BIC','HQ'}));
fprintf('Selected p = %d (BIC)\n', p_bic);

%% 5) Estimate VAR(p) & diagnostics
p = p_bic;
Mdl = varm(K, p);
[EstMdl,~,~,info] = estimate(Mdl, Y, 'Y0', Y(1:p,:)); %#ok<NASGU>
[isStableFlag, maxRoot] = local_isStable_var(EstMdl);
fprintf('\nStability (all eigenvalues < 1): %d  | max |root| = %.6f\n', isStableFlag, maxRoot);

E = infer(EstMdl, Y);    
SigmaU = cov(E,1);

fprintf('\nResidual Ljung–Box (L=12)\n');
for i = 1:K
    [h,pval] = lbqtest(E(:,i), 'Lags', 12);
    fprintf('%s: h=%d, p=%.4f\n', varNames{i}, h, pval);
end

%% 6) SVAR (Cholesky) & IRF — versione compatibile
horizon = 36;

% IRF ortogonalizzate: matrice K×K×(horizon+1), tempo t=0..H
[IRF, tvec] = local_irf_chol(EstMdl, SigmaU, horizon);

figure('Name','SVAR IRFs (Cholesky ordering)');
K = size(IRF,1);
tiledlayout(K, K, 'TileSpacing','compact','Padding','compact');
for j = 1:K
    for i = 1:K
        nexttile;
        plot(tvec, squeeze(IRF(i,j,:)), 'LineWidth', 1.25);
        yline(0,'k-'); grid on;
        title(sprintf('Resp %s to shock %s', varNames{i}, varNames{j}), 'Interpreter','none');
        xlabel('months');
    end
end


%% 7) Summary
disp(EstMdl)
fprintf('\nOrdering:\n  1) %s  2) %s  3) %s  4) %s\n', varNames{1},varNames{2},varNames{3},varNames{4});
fprintf('Lag p = %d. IRF horizon = %d months.\n', p, horizon);


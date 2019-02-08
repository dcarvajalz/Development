clc
clear all 
rng(1234)

%% SET UP 

% Set parameters
n = 1000;
T = 12*40; % Total periods: months
sigma2_e = 0.2;
sigma_e = sqrt(sigma2_e);
sigma2_u = 0.2; 
sigma_u = sqrt(sigma2_u);
beta = 0.99^(1/12); % Annual beta: 0.99
eta_1 = 1;
eta_2 = 2;
eta_4 = 4;

% Deterministic seasonal component
gm_m = [-0.147, -0.370, 0.141, 0.131, 0.090, 0.058, 0.036, 0.036, ...
    0.036, 0.002, -0.033, -0.082];
gm_h = [-0.293, -0.739, 0.282, 0.262, 0.180, 0.116, 0.072, 0.072, ...
    0.072, 0.004, -0.066, -0.164];
gm_l = [-0.073, -0.185, 0.071, 0.066, 0.045, 0.029, 0.018, 0.018, ... 
    0.018, 0.001, -0.017, -0.041];

% Stochastic seasonal component
sigma2_m_m = [0.085, 0.068, 0.290, 0.283, 0.273, 0.273, 0.239, ...
    0.205, 0.188, 0.188, 0.171, 0.137];
sigma2_m_h = [0.171, 0.137, 0.580, 0.567, 0.546, 0.546, 0.478, ...
    0.410, 0.376, 0.376, 0.341, 0.273];
sigma2_m_l = [0.043, 0.034, 0.145, 0.142, 0.137, 0.137, 0.119, ...
    0.102, 0.094, 0.094, 0.085, 0.068];

S_low = exp(kron(ones(n,40),gm_l));
S_mid = exp(kron(ones(n,40),gm_m));
S_high = exp(kron(ones(n,40),gm_h));

%% PART 1

% Create the consumption z
ln_u = mvnrnd(zeros(n,1),eye(n)*sigma2_u).'; % Nx1 random normal
z = exp(-sigma2_u/2)*exp(ln_u); 
Z = z * ones(1,T); % NxT permanent consumption

% Individuals errors
ln_e = zeros(n,T);
for i = 1:n
    for j = 0:39
        ln_e(i,(1+12*j):((j+1)*12)) = normrnd(0,sigma_e);
    end
end
E = exp(-sigma2_e/2)*exp(ln_e);

% Individual consumptions 
cl_s_r = Z.*S_low.*E; % Ind shock 
cm_s_r = Z.*S_mid.*E;
ch_s_r = Z.*S_high.*E;
c_r = Z.*E; 
cl_s = Z.*S_low; % No ind shock
cm_s = Z.*S_mid;
ch_s = Z.*S_high;
c = Z ; 

% Discounting matrix
beta_m = zeros(1,12);
beta_y = zeros(1,40);
for i = 1:12
    beta_m(1,i) = beta.^(i-1);
end
for i = 1:40 
    beta_y(1,i) = beta.^(12*i);
end
B = ones(n,1)*kron(beta_y,beta_m);


%% PART 1A

% Welfare gains without seasonal component (eta = 1)
wgl_1_s = zeros(n,1); % Here I store individuals "g"
wgm_1_s = zeros(n,1);
wgh_1_s = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum((B(i,:).*log(cl_s_r(i,:)*(1+x))).'-(B(i,:).*...
        log(c_r(i,:))).').');
    wgl_1_s(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum((B(i,:).*log(cm_s_r(i,:)*(1+x))).'-(B(i,:).*...
          log(c_r(i,:))) .').' );
    wgm_1_s(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum((B(i,:).*log(ch_s_r(i,:)*(1+x))).'-(B(i,:).*...
          log(c_r(i,:))).').');
    wgh_1_s(i,1) = fminbnd(fh,-2,5);
end

    % Results:
Results = [mean(wgl_1_s), mean(wgm_1_s), mean(wgh_1_s)];

disp('1A. Mean welfare gains of without seasonality (eta=1)')
disp(Results)


%% PART 1B

% Welfare gains without non-seasonal consumption risk (eta = 1)
wgl_1_r = zeros(n,1);
wgm_1_r = zeros(n,1);
wgh_1_r = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum((B(i,:).*log(cl_s_r(i,:)*(1+x))).'-(B(i,:).*...
        log(cl_s(i,:))).').');
    wgl_1_r(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum((B(i,:).*log(cm_s_r(i,:)*(1+x))).'-(B(i,:).*...
          log(cm_s(i,:))).').');
    wgm_1_r(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum((B(i,:).*log(ch_s_r(i,:)*(1+x))).'-(B(i,:).*...
          log(ch_s(i,:))).').');
    wgh_1_r(i,1) = fminbnd(fh,-2,5);
end

% Alternative: shut down seasonal risk, and measure the gains then (eta = 1)
wg_1_r_alt = zeros(n,1);
for i = 1:n
    f = @(x) abs(sum((B(i,:).*log(c_r(i,:)*(1+x))).'-(B(i,:)...
          .*log(c(i,:))).').');
wg_1_r_alt(i,1) = fminbnd(f,-2,5);
end

Results = [mean(wgl_1_r), mean(wgm_1_r), mean(wgh_1_r), mean(wg_1_r_alt); ...
std(wgl_1_r), std(wgm_1_r), std(wgh_1_r), std(wg_1_r_alt)];

% Results
disp('1B. Welfare gains of without the nonseasonal consumption risk (eta=1)')
disp(Results)

% Graph to compare distribution of welfare gains of without non-seasonal consumption risk
hold on
histogram(wg_1_r_alt);
xlabel('Individual welfare gain')
ylabel('Number of households')
legend('eta=1')
print('graph11B','-dpng')


%% PART 1D

% Welfare gains without seasonal component (eta = 2)
wgl_2_s = zeros(n,1);
wgm_2_s = zeros(n,1);
wgh_2_s = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum((B(i,:).*((cl_s_r(i,:)*(1+x))).^(1-eta_2)/...
        (1-eta_2)).'-(B(i,:).*((c_r(i,:)).^(1-eta_2)/(1-eta_2))).').');
    wgl_2_s(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum((B(i,:).*((cm_s_r(i,:)*(1+x))).^(1-eta_2)/...
        (1-eta_2)).'-(B(i,:).*((c_r(i,:)).^(1-eta_2)/(1-eta_2))).').');
    wgm_2_s(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum((B(i,:).*((ch_s_r(i,:)*(1+x))).^(1-eta_2)/...
        (1-eta_2)).'-(B(i,:).*((c_r(i,:)).^(1-eta_2)/(1-eta_2))).').');
    wgh_2_s(i,1) = fminbnd(fh,-2,5);
end

% Welfare gains without seasonal component (eta = 4)
wgl_4_s = zeros(n,1);
wgm_4_s = zeros(n,1);
wgh_4_s = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum((B(i,:).*((cl_s_r(i,:)*(1+x))).^(1-eta_4)/...
        (1-eta_4)).'-(B(i,:).*((c_r(i,:)).^(1-eta_4)/(1-eta_4))).').');
    wgl_4_s(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum((B(i,:).*((cm_s_r(i,:)*(1+x))).^(1-eta_4)/...
        (1-eta_4)).'-(B(i,:).*((c_r(i,:)).^(1-eta_4)/(1-eta_4))).').');
    wgm_4_s(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum((B(i,:).*((ch_s_r(i,:)*(1+x))).^(1-eta_4)/...
        (1-eta_4)).'-(B(i,:).*((c_r(i,:)).^(1-eta_4)/(1-eta_4))).').');
    wgh_4_s(i,1) = fminbnd(fh,-2,5);
end

% Results:
Results = [mean(wgl_2_s), mean(wgm_2_s), mean(wgh_2_s); ...
    mean(wgl_4_s), mean(wgm_4_s), mean(wgh_4_s)];
disp('1D. Mean welfare gains of without seasonality (eta=2 and eta=4)')
disp(Results)


% Welfare gains without non-seasonal consumption risk (eta = 2)
wgl_2_r = zeros(n,1);
wgm_2_r = zeros(n,1);
wgh_2_r = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum((B(i,:).*((cl_s_r(i,:)*(1+x))).^(1-eta_2)/...
        (1-eta_2)).'-(B(i,:).*((cl_s(i,:)).^(1-eta_2)/(1-eta_2))).').');
    wgl_2_r(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum((B(i,:).*((cm_s_r(i,:)*(1+x))).^(1-eta_2)/...
        (1-eta_2)).'-(B(i,:).*((cm_s(i,:)).^(1-eta_2)/(1-eta_2))).').');
    wgm_2_r(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum((B(i,:).*((ch_s_r(i,:)*(1+x))).^(1-eta_2)/...
        (1-eta_2)).'-(B(i,:).*((ch_s(i,:)).^(1-eta_2)/(1-eta_2))).').');
    wgh_2_r(i,1) = fminbnd(fh,-2,5);
end

% Alternative: shut down seasonal risk, and measure the gains then (eta = 2)
wg_2_r_alt = zeros(n,1);
for i = 1:n
    f = @(x) abs(sum((B(i,:).*((c_r(i,:)*(1+x))).^(1-eta_2)/...
        (1-eta_2)).'-(B(i,:).*((c(i,:)).^(1-eta_2)/(1-eta_2))).').');
    wg_2_r_alt(i,1) = fminbnd(f,-2,5);
end

% Welfare gains without non-seasonal consumption risk (eta = 4)
wgl_4_r = zeros(n,1);
wgm_4_r = zeros(n,1);
wgh_4_r = zeros(n,1);
for i = 1:n   
    fl = @(x) abs(sum((B(i,:).*((cl_s_r(i,:)*(1+x))).^(1-eta_4)/...
        (1-eta_4)).'-(B(i,:).*((cl_s(i,:)).^(1-eta_4)/(1-eta_4))).').');
    wgl_4_r(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum((B(i,:).*((cm_s_r(i,:)*(1+x))).^(1-eta_4)/...
        (1-eta_4)).'-(B(i,:).*((cm_s(i,:)).^(1-eta_4)/(1-eta_4))).').');
    wgm_4_r(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum((B(i,:).*((ch_s_r(i,:)*(1+x))).^(1-eta_4)/...
        (1-eta_4)).'-(B(i,:).*((ch_s(i,:)).^(1-eta_4)/(1-eta_4))).').');
    wgh_4_r(i,1) = fminbnd(fh,-2,5);
end

% Alternative: shut down seasonal risk, and measure the gains then (eta = 4)
wg_4_r_alt = zeros(n,1);
for i = 1:n
    f = @(x) abs(sum((B(i,:).*((c_r(i,:)*(1+x))).^(1-eta_4)/...
        (1-eta_4)).'-(B(i,:).*((c(i,:)).^(1-eta_4)/(1-eta_4))).').');
    wg_4_r_alt(i,1) = fminbnd(f,-2,5);
end

Results2 = [mean(wgl_2_r), mean(wgm_2_r), mean(wgh_2_r); ...
    std(wgl_2_r), std(wgm_2_r), std(wgh_2_r)];

Results4 = [mean(wgl_4_r), mean(wgm_4_r), mean(wgh_4_r); ...
    std(wgl_4_r), std(wgm_4_r), std(wgh_4_r)];

% Results (eta=2)
disp('1D Welfare gains of without non-seasonal consumption risk (eta=2)')
disp(Results2)


% Results (eta=4)
disp('1D Welfare gains of without non-seasonal consumption risk (eta=4)')
disp(Results4)


% Graph to compare distribution of welfare gains of without non-seasonal consumption risk
figure 
hold on
histogram(wg_1_r_alt,'FaceColor','b');
hold on
histogram(wg_2_r_alt,'FaceColor','y');
hold on
histogram(wg_4_r_alt,'FaceColor','r');
xlabel('Individual welfare gain')
ylabel('Number of households')
legend('eta=1','eta=2','eta=4')
print('graph11D','-dpng')


%% PART 2%%

% Seasonal components
sr_l = zeros(n,T);
for i = 1:n
    for j = 0:39
      sr_l(i,(1+j*12):(j+1)*12) = exp(-sigma2_m_l/2) .* ...
      exp(mvnrnd(zeros(12,1),ones(12,1)*sigma2_m_l.*eye(12)));
    end
end
sr_m = zeros(n,T);
for i = 1:n
    for j = 0:39
      sr_m(i,(1+j*12):(j+1)*12) = exp(-sigma2_m_m/2) .* ...
      exp( mvnrnd(zeros(12,1), ones(12,1) * sigma2_m_m .* eye(12) ) );
    end
end
sr_h = zeros(n,T);
for i = 1:n
    for j = 0:39
      sr_h(i,(1+j*12):(j+1)*12) = exp(-sigma2_m_h/2) .* ...
      exp( mvnrnd(zeros(12,1), ones(12,1) * sigma2_m_h .* eye(12) ) );
    end
end

c_sr_l_r = Z .* sr_l .* E; 
c_sr_m_r = Z .* sr_m .* E;
c_sr_h_r = Z .* sr_h .* E;
C_sm_r = Z .* S_mid .* E; 

c_sm_sr_l = Z .* S_mid .*  sr_l; 
c_sm_sr_m = Z .* S_mid .* sr_m;
c_sm_sr_h = Z .* S_mid .* sr_h;

c_sm_sr_l_r = Z .* S_mid .*  sr_l .* E; 
c_sm_sr_m_r = Z .* S_mid .* sr_m .* E;
c_sm_sr_h_r = Z .* S_mid .* sr_h .* E;


%% PART 2A

% Welfare gains without deterministic seasonal component (eta = 1)
wgl_1_s2 = zeros(n,1); 
wgm_1_s2 = zeros(n,1);
wgh_1_s2 = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_l_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_sr_l_r(i,:))) .').' );
    wgl_1_s2(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_m_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_sr_m_r(i,:))) .').' );
    wgm_1_s2(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_h_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_sr_h_r(i,:))) .').' );
    wgh_1_s2(i,1) = fminbnd(fh,-2,5);
end

% Welfare gains without stochastic seasonal component (eta = 1)
wgl_1_sr = zeros(n,1); 
wgm_1_sr = zeros(n,1);
wgh_1_sr = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_l_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(C_sm_r(i,:))) .').' );
    wgl_1_sr(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_m_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(C_sm_r(i,:))) .').' );
    wgm_1_sr(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_h_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(C_sm_r(i,:))) .').' );
    wgh_1_sr(i,1) = fminbnd(fh,-2,5);
end
  
% Welfare gains without both seasonal components (eta = 1)
wgl_1_s_sr = zeros(n,1); 
wgm_1_s_sr = zeros(n,1);
wgh_1_s_sr = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_l_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_r(i,:))) .').' );
    wgl_1_s_sr(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_m_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_r(i,:))) .').' );
    wgm_1_s_sr(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_h_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_r(i,:))) .').' );
    wgh_1_s_sr(i,1) = fminbnd(fh,-2,5);
end

% Results without deterministic component
Results = [mean(wgl_1_s2), mean(wgm_1_s2), mean(wgh_1_s2)];
disp('2A. Welfare gains of without deterministic seasonal component (eta=1)')
disp(Results)


% Results without stochastic component
Results = [mean(wgl_1_sr), mean(wgm_1_sr), mean(wgh_1_sr); ...
    std(wgl_1_sr), std(wgm_1_sr), std(wgh_1_sr)];
disp('2A. Welfare gains of without stochastic seasonal components (eta=1)')
disp(Results)

% Results without both stochastic and deterministic component
Results = [mean(wgl_1_s_sr), mean(wgm_1_s_sr), mean(wgh_1_s_sr); ...
    std(wgl_1_s_sr), std(wgm_1_s_sr), std(wgh_1_s_sr)];
disp('2A. Welfare gains of without both seasonality components (eta=1)')
disp(Results)


% Graph of the distribution of welfare gains of without the seasonal stochastic component
figure 
hold on
histogram(wgl_1_sr,'FaceColor','b');
hold on
histogram(wgm_1_sr,'FaceColor','y');
hold on
histogram(wgh_1_sr,'FaceColor','r');
xlabel('Individual welfare gain')
ylabel('Number of households')
legend('Low seasonality dispersion','Medium seasonality dispersion','High seasonality dispersion')
print('graph12A1','-dpng')

% Graph of the distribution of welfare gains of without both seasonal components
figure 
hold on
histogram(wgl_1_s_sr,'FaceColor','b');
hold on
histogram(wgm_1_s_sr,'FaceColor','y');
hold on
histogram(wgh_1_s_sr,'FaceColor','r');
xlabel('Individual welfare gain')
ylabel('Number of households')
legend('Low seasonality dispersion','Medium seasonality dispersion','High seasonality dispersion')
print('graph12A2','-dpng')

  
%% PART 2B

% Welfare gains without non-seasonal consumption risk (eta = 1)
wgl_1_r_q2 = zeros(n,1);
wgm_1_r_q2 = zeros(n,1);
wgh_1_r_q2 = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* log( c_sm_sr_l_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_sm_sr_l(i,:))) .').' );
    wgl_1_r_q2(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_m_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_sm_sr_m(i,:))) .').' );
    wgm_1_r_q2(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* log(c_sm_sr_h_r(i,:)*(1+x))) .' - ...
          (B(i,:) .* log(c_sm_sr_h(i,:))) .').' );
    wgh_1_r_q2(i,1) = fminbnd(fh,-2,5);
end

% Results
Results = [mean(wgl_1_r_q2), mean(wgm_1_r_q2), mean(wgh_1_r_q2); ...
    std(wgl_1_r_q2), std(wgm_1_r_q2), std(wgh_1_r_q2)];

disp('2B. Welfare gains of without non-seasonal consumption risk (eta=1)')
disp(Results)


% Graph to compare distribution of welfare gains of without non-seasonal consumption risk
figure
hold on
hist(wgl_1_r_q2);
xlabel('Individual welfare gain')
ylabel('Number of households')
legend('eta=1')
print('graph12B','-dpng')



%% PART 2D

% Welfare gains without deterministic seasonal component (eta = 2)
wgl_2_s2 = zeros(n,1); 
wgm_2_s2 = zeros(n,1);
wgh_2_s2 = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_sr_l_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgl_2_s2(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_sr_m_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgm_2_s2(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_sr_h_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgh_2_s2(i,1) = fminbnd(fh,-2,5);
end


% Welfare gains without deterministic seasonal component (eta = 4)
wgl_4_s2 = zeros(n,1); 
wgm_4_s2 = zeros(n,1);
wgh_4_s2 = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_sr_l_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgl_4_s2(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_sr_m_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgm_4_s2(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_sr_h_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgh_4_s2(i,1) = fminbnd(fh,-2,5);
end

% Results without deterministic component (eta = 2)
Results = [mean(wgl_2_s2), mean(wgm_2_s2), mean(wgh_2_s2)];
disp('2D. Welfare gains of without deterministic seasonal component (eta=2)')
disp(Results)


% Results without deterministic component (eta = 4)
Results = [mean(wgl_4_s2), mean(wgm_4_s2), mean(wgh_4_s2)];
disp('2D. Welfare gains of without deterministic seasonal component (eta=4)')
disp(Results)
disp('Rows: Means')
disp('Columns: low, mid and high seasonality')
disp(' ')
disp(' ')

% Welfare gains without stochastic seasonal component (eta = 2)
wgl_2_sr = zeros(n,1); 
wgm_2_sr = zeros(n,1);
wgh_2_sr = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((C_sm_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgl_2_sr(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((C_sm_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgm_2_sr(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((C_sm_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgh_2_sr(i,1) = fminbnd(fh,-2,5);
end


% Welfare gains without stochastic seasonal component (eta = 4)
wgl_4_sr = zeros(n,1); 
wgm_4_sr = zeros(n,1);
wgh_4_sr = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((C_sm_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgl_4_sr(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((C_sm_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgm_4_sr(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((C_sm_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgh_4_sr(i,1) = fminbnd(fh,-2,5);
end

% Results without stochastic component (eta=2)
Results = [mean(wgl_2_sr), mean(wgm_2_sr), mean(wgh_2_sr); ...
    std(wgl_2_sr), std(wgm_2_sr), std(wgh_2_sr)];
disp('2D. Welfare gains of without stochastic seasonal components (eta=2)')
disp(Results)


% Results without stochastic component (eta=4)
Results = [mean(wgl_4_sr), mean(wgm_4_sr), mean(wgh_4_sr); ...
    std(wgl_4_sr), std(wgm_4_sr), std(wgh_4_sr)];
disp('2D. Welfare gains of without stochastic seasonal components (eta=4)')
disp(Results)


% Welfare gains without both seasonal components (eta = 2)
wgl_2_s_sr = zeros(n,1); 
wgm_2_s_sr = zeros(n,1);
wgh_2_s_sr = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgl_2_s_sr(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgm_2_s_sr(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+x))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_r(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgh_2_s_sr(i,1) = fminbnd(fh,-2,5);
end

% Welfare gains without both seasonal components (eta = 4)
wgl_4_s_sr = zeros(n,1); 
wgm_4_s_sr = zeros(n,1);
wgh_4_s_sr = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgl_4_s_sr(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgm_4_s_sr(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_r(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgh_4_s_sr(i,1) = fminbnd(fh,-2,5);
end

% Results without both stochastic and deterministic component (eta=2):
Results = [mean(wgl_2_s_sr), mean(wgm_2_s_sr), mean(wgh_2_s_sr); ...
    std(wgl_2_s_sr), std(wgm_2_s_sr), std(wgh_2_s_sr)];
disp('2D. Welfare gains of without both seasonality components (eta=2)')
disp(Results)


% Results without both stochastic and deterministic component (eta=4):
Results = [mean(wgl_4_s_sr), mean(wgm_4_s_sr), mean(wgh_4_s_sr); ...
    std(wgl_4_s_sr), std(wgm_4_s_sr), std(wgh_4_s_sr)];
disp('2D. Welfare gains of without both seasonality components (eta=4)')
disp(Results)


% Welfare gains without non-seasonal consumption risk (eta = 2)
wgl_2_r_q2 = zeros(n,1);
wgm_2_r_q2 = zeros(n,1);
wgh_2_r_q2 = zeros(n,1);
for i = 1:n
    fl = @(gl) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+gl))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_sm_sr_l(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgl_2_r_q2(i,1) = fminbnd(fl,-2,5);
    fm = @(gm) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+gm))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_sm_sr_m(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgm_2_r_q2(i,1) = fminbnd(fm,-2,5);
    fh = @(gh) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+gh))).^(1-eta_2) / (1-eta_2)).' - ...
          (B(i,:) .* ((c_sm_sr_h(i,:)).^(1-eta_2) / (1-eta_2)) ) .').' );
    wgh_2_r_q2(i,1) = fminbnd(fh,-2,5);
end

% Welfare gains without non-seasonal consumption risk (eta = 4)
wgl_4_r_q2 = zeros(n,1);
wgm_4_r_q2 = zeros(n,1);
wgh_4_r_q2 = zeros(n,1);
for i = 1:n
    fl = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_l_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_sm_sr_l(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgl_4_r_q2(i,1) = fminbnd(fl,-2,5);
    fm = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_m_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_sm_sr_m(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgm_4_r_q2(i,1) = fminbnd(fm,-2,5);
    fh = @(x) abs(sum( (B(i,:) .* ((c_sm_sr_h_r(i,:)*(1+x))).^(1-eta_4) / (1-eta_4)).' - ...
          (B(i,:) .* ((c_sm_sr_h(i,:)).^(1-eta_4) / (1-eta_4)) ) .').' );
    wgh_4_r_q2(i,1) = fminbnd(fh,-2,5);
end

% Results
Results = [mean(wgl_2_r_q2), mean(wgm_2_r_q2), mean(wgh_2_r_q2); ...
    std(wgl_2_r_q2), std(wgm_2_r_q2), std(wgh_2_r_q2)];

disp('2D. Welfare gains of without non-seasonal consumption risk (eta=2)')
disp(Results)

% Results
Results = [mean(wgl_4_r_q2), mean(wgm_4_r_q2), mean(wgh_4_r_q2); ...
    std(wgl_4_r_q2), std(wgm_4_r_q2), std(wgh_4_r_q2)];

disp('2D. Welfare gains of without non-seasonal consumption risk (eta=4)')
disp(Results)

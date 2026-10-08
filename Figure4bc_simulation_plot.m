%% Recreate lineage diversity and survival panels
% Uses the full antigen-screening dynamics and individual-cell lineage
% bookkeeping described in the manuscript Methods.

clear; clc
output_dir = fullfile(fileparts(mfilename('fullpath')), 'Figure4_output');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

%% ---------------- User options ----------------
mutation_delay_days = 0;      % 0 for main Fig. 4; 5 for delayed-mutation analysis
use_antigen_decay   = false;  % main Fig. 4 uses constant antigen
save_results        = false;
save_figure         = false;
use_parallel        = false;  % set true if Parallel Computing Toolbox is available
num_workers         = [];     % [] = MATLAB default; or set, e.g., 5

s_values = 1:5;
Rep      = 15;
Tend_day = 20;

% Panel-b summaries are calculated over this interval. With a mutation
% delay, the interval automatically begins no earlier than the delay.
analysis_start_day = 7;

%% ---------------- Core model parameters ----------------
N_aff = 150;
N_T   = 3;
dt    = 3e-2;                 % hours

del_B  = 1/8;
del_T  = 1/18;
kappa  = 5;
lambda = 0.077;

B0          = 100;
T0_total    = 12;
iseed       = 15;
eps_aff     = 0.5;
ppos        = 0.05;
N_lineages  = 50;

alpha_A0 = 1000*del_B;
Ag_mid   = 25;
Ag_width = 5;

%% ---------------- Run the sensitivity sweep ----------------
nS = numel(s_values);
day_grid = 0:Tend_day;
nDay = numel(day_grid);

diversity       = nan(nS,Rep,nDay);
mean_mut_depth  = nan(nS,Rep,nDay);
mean_affinity_class = nan(nS,Rep,nDay);
affinity_shift  = nan(nS,Rep,nDay); % eps*(<n(t)>-<n(0)>)
survival        = nan(nS,Rep,nDay);
B_total         = nan(nS,Rep,nDay);
T_total         = nan(nS,Rep,nDay);
mutation_rate   = nan(nS,Rep);
diversity_panel = nan(nS,Rep);

nJobs = nS*Rep;
run_results = cell(nJobs,1);
fprintf('Running %d sensitivities x %d replicates (%d independent jobs)\n', ...
    nS,Rep,nJobs)

if use_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(num_workers)
            pool = parpool;
        else
            pool = parpool(num_workers);
        end
    end
    fprintf('Using %d parallel workers\n',pool.NumWorkers)

    % Flattening (s,replicate) into one parfor gives all workers useful
    % work even if some stochastic simulations finish before others.
    parfor job = 1:nJobs
        [is,rr] = ind2sub([nS,Rep],job);
        run_results{job} = run_one_job(s_values(is),day_grid, ...
            mutation_delay_days,use_antigen_decay,analysis_start_day, ...
            N_aff,N_T,dt,del_B,del_T,kappa,lambda,B0,T0_total,iseed, ...
            eps_aff,ppos,N_lineages,alpha_A0,Ag_mid,Ag_width);
    end
else
    for job = 1:nJobs
        [is,rr] = ind2sub([nS,Rep],job);
        run_results{job} = run_one_job(s_values(is),day_grid, ...
            mutation_delay_days,use_antigen_decay,analysis_start_day, ...
            N_aff,N_T,dt,del_B,del_T,kappa,lambda,B0,T0_total,iseed, ...
            eps_aff,ppos,N_lineages,alpha_A0,Ag_mid,Ag_width);
    end
end

% Assemble the compact worker results on the client. Keeping writes to the
% numeric output arrays outside parfor avoids sliced-variable restrictions.
for job = 1:nJobs
    [is,rr] = ind2sub([nS,Rep],job);
    run = run_results{job};
    diversity(is,rr,:)      = run.diversity;
    mean_mut_depth(is,rr,:) = run.mean_mut_depth;
    mean_affinity_class(is,rr,:) = run.mean_affinity_class;
    affinity_shift(is,rr,:) = run.affinity_shift;
    survival(is,rr,:)       = run.survival;
    B_total(is,rr,:)        = run.B_total;
    T_total(is,rr,:)        = run.T_total;
    mutation_rate(is,rr)    = run.mutation_rate;
    diversity_panel(is,rr)  = run.diversity_panel;
end
clear run_results

%% ---------------- Summary statistics ----------------
mut_mean = mean(mutation_rate,2,'omitnan');
mut_sem  = std(mutation_rate,0,2,'omitnan') ./ ...
           sqrt(sum(isfinite(mutation_rate),2));

div_mean = mean(diversity_panel,2,'omitnan');
div_sem  = std(diversity_panel,0,2,'omitnan') ./ ...
           sqrt(sum(isfinite(diversity_panel),2));

surv_mean = squeeze(mean(survival,2,'omitnan'));
surv_sem  = squeeze(std(survival,0,2,'omitnan')) ./ ...
            sqrt(squeeze(sum(isfinite(survival),2)));

aff_shift_mean = squeeze(mean(affinity_shift,2,'omitnan'));
aff_shift_sem  = squeeze(std(affinity_shift,0,2,'omitnan')) ./ ...
                 sqrt(squeeze(sum(isfinite(affinity_shift),2)));
mut_depth_mean = squeeze(mean(mean_mut_depth,2,'omitnan'));
mut_depth_sem  = squeeze(std(mean_mut_depth,0,2,'omitnan')) ./ ...
                 sqrt(squeeze(sum(isfinite(mean_mut_depth),2)));

%% ---------------- Plot panels b and c ----------------
green_map = [ ...
    0.10 0.29 0.17; ... % s = 1
    0.20 0.43 0.25; ...
    0.36 0.58 0.36; ...
    0.56 0.72 0.52; ...
    0.76 0.86 0.70];    % s = 5

figure('Color','w','Position',[100 100 820 340]);

% Panel b
subplot(1,2,1); hold on
for is = nS:-1:1
    errorbar(mut_mean(is),div_mean(is)/1e3,div_sem(is)/1e3,div_sem(is)/1e3, ...
        mut_sem(is),mut_sem(is),'o','Color',green_map(is,:), ...
        'MarkerFaceColor',green_map(is,:),'MarkerEdgeColor',[0.1 0.3 0.15], ...
        'LineWidth',0.8,'MarkerSize',6,'CapSize',6, ...
        'DisplayName',sprintf('s = %g',s_values(is)));
end
xlabel('Mutations/day')
ylabel('Diversity (10^3)')
box off
legend('Location','southwest','Box','off')
title('b','Units','normalized','Position',[-0.08 1.02], ...
    'FontWeight','bold','FontSize',13)

% Panel c: show replicate trajectories and mean +/- SEM for s=1 and s=5.
subplot(1,2,2); hold on
show_s = [nS 1];
for jj = 1:numel(show_s)
    is = show_s(jj);
    c = green_map(is,:);
    faint = 0.72 + 0.28*c;
    for rr = 1:Rep
        plot(day_grid,squeeze(survival(is,rr,:)),'-','Color',faint, ...
            'LineWidth',0.45,'HandleVisibility','off');
    end
    errorbar(day_grid,surv_mean(is,:),surv_sem(is,:),'o-', ...
        'Color',c,'MarkerFaceColor',c,'MarkerEdgeColor','k', ...
        'LineWidth',1.1,'MarkerSize',5,'CapSize',3, ...
        'DisplayName',sprintf('s = %g',s_values(is)));
end
xlabel('Days')
ylabel('Fraction surviving lineages')
xlim([0 Tend_day]); ylim([0 1.02])
box off
legend('Location','northeast','Box','off')
title('c','Units','normalized','Position',[-0.08 1.02], ...
    'FontWeight','bold','FontSize',13)

set(findall(gcf,'Type','axes'),'FontName','Arial','FontSize',10,'LineWidth',0.8)
fig_lineage = gcf;

% Select filenames automatically from the mutation-delay condition.
if mutation_delay_days > 0
    base_name = 'Figure4_bc_wDelay';
else
    base_name = 'Figure4_bc';
end

if save_figure
    savefig(fig_lineage,fullfile(output_dir,[base_name '.fig']));
    exportgraphics(fig_lineage,fullfile(output_dir,[base_name '.png']), ...
        'Resolution',300);
    exportgraphics(fig_lineage,fullfile(output_dir,[base_name '.pdf']), ...
        'ContentType','vector');
end

%% ---------------- Daily affinity and mutation figure ----------------
fig_timecourses = figure('Color','w','Position',[100 480 820 340]);

subplot(1,2,1); hold on
for is = nS:-1:1
    errorbar(day_grid,aff_shift_mean(is,:),aff_shift_sem(is,:),'o-', ...
        'Color',green_map(is,:),'MarkerFaceColor',green_map(is,:), ...
        'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',4,'CapSize',3, ...
        'DisplayName',sprintf('s = %g',s_values(is)));
end
xlabel('Days')
ylabel('\epsilon[\langle n(t)\rangle-\langle n(0)\rangle]')
xlim([day_grid(1) day_grid(end)])
box off
legend('Location','best','Box','off')
title('a','Units','normalized','Position',[-0.08 1.02], ...
    'FontWeight','bold','FontSize',13)

subplot(1,2,2); hold on
for is = nS:-1:1
    errorbar(day_grid,mut_depth_mean(is,:),mut_depth_sem(is,:),'o-', ...
        'Color',green_map(is,:),'MarkerFaceColor',green_map(is,:), ...
        'MarkerEdgeColor','k','LineWidth',1,'MarkerSize',4,'CapSize',3, ...
        'DisplayName',sprintf('s = %g',s_values(is)));
end
xlabel('Days')
ylabel('Average mutations per B cell')
xlim([day_grid(1) day_grid(end)])
box off
legend('Location','best','Box','off')
title('b','Units','normalized','Position',[-0.08 1.02], ...
    'FontWeight','bold','FontSize',13)
set(findall(fig_timecourses,'Type','axes'),'FontName','Arial', ...
    'FontSize',10,'LineWidth',0.8)

if save_figure
    savefig(fig_timecourses,fullfile(output_dir,[base_name '_timecourses.fig']));
    exportgraphics(fig_timecourses, ...
        fullfile(output_dir,[base_name '_timecourses.png']),'Resolution',300);
    exportgraphics(fig_timecourses, ...
        fullfile(output_dir,[base_name '_timecourses.pdf']), ...
        'ContentType','vector');
end

%% ---------------- Save outputs ----------------
settings = struct( ...
    'mutation_delay_days',mutation_delay_days, ...
    'use_antigen_decay',use_antigen_decay, ...
    'use_parallel',use_parallel,'num_workers',num_workers, ...
    's_values',s_values,'Rep',Rep,'Tend_day',Tend_day, ...
    'analysis_start_day',analysis_start_day,'dt_hours',dt);

if save_results
    save(fullfile(output_dir,[base_name '.mat']), ...
        'settings','day_grid', ...
        'diversity','mean_mut_depth','mean_affinity_class', ...
        'affinity_shift','survival','B_total','T_total', ...
        'mutation_rate','diversity_panel','mut_mean','mut_sem', ...
        'div_mean','div_sem','surv_mean','surv_sem', ...
        'aff_shift_mean','aff_shift_sem','mut_depth_mean','mut_depth_sem');
end

%% ---------------- Daily CSV exports ----------------
% Replicate-level tidy table: one row per sensitivity, replicate and day.
[is_grid,rep_grid,day_index] = ndgrid(1:nS,1:Rep,1:nDay);
daily_metrics = table( ...
    reshape(s_values(is_grid),[],1),rep_grid(:), ...
    reshape(day_grid(day_index),[],1),B_total(:),T_total(:), ...
    mean_affinity_class(:),affinity_shift(:),mean_mut_depth(:), ...
    diversity(:),survival(:), ...
    'VariableNames',{'sensitivity','replicate','day','B_population', ...
    'T_population','average_affinity_class','eps_delta_average_affinity', ...
    'average_mutations_per_B_cell','inverse_simpson_diversity', ...
    'fraction_surviving_lineages'});
writetable(daily_metrics,fullfile(output_dir,[base_name '_daily_metrics.csv']));

% Across-replicate daily mean and SEM for convenient plotting/statistics.
daily_summary = make_daily_summary(s_values,day_grid,B_total,T_total, ...
    mean_affinity_class,affinity_shift,mean_mut_depth,diversity,survival);
writetable(daily_summary,fullfile(output_dir,[base_name '_daily_summary.csv']));

%% ============================================================
% Local functions
% =============================================================
function run = run_one_job(s,day_grid,mutation_delay_days, ...
    use_antigen_decay,analysis_start_day,N_aff,N_T,dt,del_B,del_T, ...
    kappa,lambda,B0,T0_total,iseed,eps_aff,ppos,N_lineages, ...
    alpha_A0,Ag_mid,Ag_width)

run = simulate_one(s,day_grid,mutation_delay_days, ...
    use_antigen_decay,N_aff,N_T,dt,del_B,del_T,kappa,lambda, ...
    B0,T0_total,iseed,eps_aff,ppos,N_lineages, ...
    alpha_A0,Ag_mid,Ag_width);

run.mutation_rate = NaN;
run.diversity_panel = NaN;
t0 = max(analysis_start_day,mutation_delay_days);

keep = day_grid >= t0 & isfinite(run.mean_mut_depth);
if nnz(keep) >= 2
    p = polyfit(day_grid(keep),run.mean_mut_depth(keep),1);
    run.mutation_rate = p(1);
end

keepD = day_grid >= t0 & isfinite(run.diversity);
if any(keepD)
    run.diversity_panel = mean(run.diversity(keepD));
end
end

function run = simulate_one(s,day_grid,mutation_delay_days, ...
    use_antigen_decay,N_aff,N_T,dt,del_B,del_T,kappa,lambda, ...
    B0,T0_total,iseed,eps_aff,ppos,N_lineages, ...
    alpha_A0,Ag_mid,Ag_width)

nvec = (1:N_aff)';
M = round(day_grid(end)*24/dt);
save_step = round(24/dt);

% Two cells per founder for B0=100 and N_lineages=50.
Bcells = cell(B0,2);
for j = 1:B0
    Bcells(j,:) = {iseed,1 + mod(j-1,N_lineages)};
end
Btemp = zeros(N_aff,1);
Btemp(iseed) = B0;

% Use the same total initial TFH number as Fig3_new. Clones are identical
% within a sensitivity condition, so this is a scalar-s sweep.
Ttemp = zeros(N_T,1);
q = floor(T0_total/N_T);
Ttemp(:) = q;
Ttemp(1:mod(T0_total,N_T)) = Ttemp(1:mod(T0_total,N_T)) + 1;
s_vec = s*ones(N_T,1);

nmut_vec = zeros(500,1);
nmut_vec(1) = N_lineages;

nDay = numel(day_grid);
run.diversity     = nan(1,nDay);
run.mean_mut_depth = nan(1,nDay);
run.mean_affinity_class = nan(1,nDay);
run.affinity_shift = nan(1,nDay);
run.survival      = zeros(1,nDay);
run.B_total       = zeros(1,nDay);
run.T_total       = zeros(1,nDay);

run = record_stats(run,1,Bcells,Ttemp,N_lineages,iseed,eps_aff);

for it = 1:M
    if isempty(Bcells)
        % B-cell extinction: survival, diversity and population remain zero;
        % mutation depth is undefined after the population disappears.
        first_future = floor((it-1)/save_step) + 2;
        if first_future <= nDay
            run.diversity(first_future:end) = 0;
            run.survival(first_future:end) = 0;
            run.B_total(first_future:end) = 0;
            run.T_total(first_future:end) = 0;
        end
        break
    end

    day = (it-1)*dt/24;
    NB = sum(Btemp);
    nbar = sum(Btemp.*nvec)/NB;

    if use_antigen_decay
        Ag_frac = 1/(1 + exp((day-Ag_mid)/Ag_width));
    else
        Ag_frac = 1;
    end
    alpha_A = alpha_A0*Ag_frac;

    if sum(Ttemp) > 0
        s_av = sum(s_vec.*Ttemp)/sum(Ttemp);
    else
        % No B-cell births are possible without TFHs; this value is only
        % a harmless placeholder until the remaining B cells die.
        s_av = s;
    end
    f_vec = 1./(1 + exp(-eps_aff*(nvec-nbar)));
    psi_vec = 1./(1 + 1./(f_vec.*s_vec'));

    Bscr = Bscr_finder(Btemp,f_vec,alpha_A,del_B);
    Ctot = sum(Ttemp) + NB;

    dB_birth = dt*kappa*Bscr .* sum(psi_vec.*Ttemp'/Ctot,2) ./ Btemp;
    dB_birth(~isfinite(dB_birth)) = 0;
    dB_birth = min(max(dB_birth,0),1);
    dB_death = min(max(dt*del_B,0),1);

    dT_birth = dt*lambda*kappa*sum(psi_vec.*Bscr/Ctot,1)';
    dT_birth = min(max(dT_birth,0),1);
    dT_death = min(max(dt*del_T,0),1);

    mutation_on = day >= mutation_delay_days;
    [Bcells,nmut_vec] = bmut_delay(Bcells,dB_birth,eps_aff,s_av, ...
        nbar,nmut_vec,ppos,N_aff,mutation_on);
    Bcells = c_death(Bcells,dB_death);
    Btemp = Bcounter(Bcells,N_aff);

    for k = 1:N_T
        if Ttemp(k) > 0
            Ttemp(k) = Ttemp(k) + binornd(Ttemp(k),dT_birth(k)) ...
                - binornd(Ttemp(k),dT_death);
        end
    end
    Ttemp(Ttemp < 1) = 0;

    if mod(it,save_step) == 0
        iday = it/save_step + 1;
        run = record_stats(run,iday,Bcells,Ttemp,N_lineages,iseed,eps_aff);
    end
end
end

function run = record_stats(run,iday,Bcells,Ttemp,N_lineages,iseed,eps_aff)
N = size(Bcells,1);
run.B_total(iday) = N;
run.T_total(iday) = sum(Ttemp);
if N == 0
    run.diversity(iday) = 0;
    run.survival(iday) = 0;
    return
end

affinity_class = cell2mat(Bcells(:,1));
run.mean_affinity_class(iday) = mean(affinity_class);
run.affinity_shift(iday) = eps_aff*(mean(affinity_class)-iseed);

IDs = Bcells(:,2);
root_ID = cellfun(@(x) x(1),IDs);
founder_counts = histcounts(root_ID,0.5:(N_lineages+0.5));
run.survival(iday) = nnz(founder_counts)/N_lineages;

IDstr = cellfun(@(x) sprintf('%d_',x),IDs,'UniformOutput',false);
[~,~,ic] = unique(IDstr,'stable');
clone_sizes = accumarray(ic,1);
p = clone_sizes/sum(clone_sizes);
run.diversity(iday) = 1/sum(p.^2);

% Subtract one because the first ID element is the founder label, not a
% mutation event.
run.mean_mut_depth(iday) = mean(cellfun(@numel,IDs)-1);
end

function summary = make_daily_summary(s_values,day_grid,B_total,T_total, ...
    mean_affinity_class,affinity_shift,mean_mut_depth,diversity,survival)

nS = numel(s_values);
nDay = numel(day_grid);
[is_grid,day_index] = ndgrid(1:nS,1:nDay);

[B_mean,B_sem] = mean_sem(B_total);
[T_mean,T_sem] = mean_sem(T_total);
[aff_mean,aff_sem] = mean_sem(mean_affinity_class);
[shift_mean,shift_sem] = mean_sem(affinity_shift);
[mut_mean_daily,mut_sem_daily] = mean_sem(mean_mut_depth);
[div_mean_daily,div_sem_daily] = mean_sem(diversity);
[surv_mean_daily,surv_sem_daily] = mean_sem(survival);

summary = table(reshape(s_values(is_grid),[],1), ...
    reshape(day_grid(day_index),[],1), ...
    B_mean(:),B_sem(:),T_mean(:),T_sem(:), ...
    aff_mean(:),aff_sem(:),shift_mean(:),shift_sem(:), ...
    mut_mean_daily(:),mut_sem_daily(:), ...
    div_mean_daily(:),div_sem_daily(:), ...
    surv_mean_daily(:),surv_sem_daily(:), ...
    'VariableNames',{'sensitivity','day','B_population_mean', ...
    'B_population_sem','T_population_mean','T_population_sem', ...
    'average_affinity_class_mean','average_affinity_class_sem', ...
    'eps_delta_average_affinity_mean','eps_delta_average_affinity_sem', ...
    'average_mutations_per_B_cell_mean','average_mutations_per_B_cell_sem', ...
    'inverse_simpson_diversity_mean','inverse_simpson_diversity_sem', ...
    'fraction_surviving_lineages_mean','fraction_surviving_lineages_sem'});
end

function [m,se] = mean_sem(A)
m = squeeze(mean(A,2,'omitnan'));
se = squeeze(std(A,0,2,'omitnan')) ./ ...
     sqrt(squeeze(sum(isfinite(A),2)));
end

function Bscr = Bscr_finder(Btemp,f_vec,alpha_A,del_B)
Btemp = Btemp(:);
f_vec = f_vec(:);
A = alpha_A.*f_vec;
S = sum(f_vec.*Btemp);
Bscr0 = ones(size(Btemp));
tol = 1e-2;
maxit = 200;

for it = 1:maxit
    x = sum(f_vec.*Bscr0);
    Bscr = (Btemp.*A)./(A + del_B*(S-x));
    valid = Btemp > 0.5;
    if any(valid)
        eps_rel = mean(abs(Bscr(valid)-Bscr0(valid))./Btemp(valid));
    else
        eps_rel = 0;
    end
    if eps_rel <= tol && all(Bscr(Btemp>0) <= Btemp(Btemp>0)+1e-2)
        break
    end
    Bscr0 = Bscr;
end
end

function Bnew = c_death(Bcells,dB_death)
if isempty(Bcells)
    Bnew = Bcells;
else
    Bnew = Bcells(rand(size(Bcells,1),1) >= dB_death,:);
end
end

function [Bnew,nmut_vec] = bmut_delay(Bcells,dB_birth,eps_aff,s_av, ...
    nbar,nmut_vec,ppos,N_aff,mutation_on)
N = size(Bcells,1);
if N == 0
    Bnew = Bcells;
    return
end

n_i = cell2mat(Bcells(:,1));
birth = rand(N,1) < dB_birth(n_i);
ib = find(birth);
nb = numel(ib);
tail = N + (1:nb);

Bnew = cell(N+nb,2);
Bnew(1:N,:) = Bcells;
Bnew(tail,:) = Bcells(ib,:);

if ~mutation_on
    return
end

for k = 1:nb
    i = ib(k);
    n = Bcells{i,1};
    mut_ID = Bcells{i,2};
    psi = 1/(1 + (1 + exp(-eps_aff*(n-nbar)))/s_av);
    pmut = min(max(0.6-0.4*psi,0),1);
    pos = [i,tail(k)];

    for j = 1:2
        if rand < pmut
            [new_ID,nmut_vec] = mutgen(mut_ID,nmut_vec);
            new_n = n;
            if rand < 0.2
                if rand < ppos
                    new_n = n+1;
                else
                    new_n = n-1;
                end
            end
            new_n = min(max(new_n,1),N_aff);
            Bnew(pos(j),:) = {new_n,new_ID};
        end
    end
end
end

function [new_ID,nmut_vec] = mutgen(mut_ID,nmut_vec)
depth = numel(mut_ID);
if depth+1 > numel(nmut_vec)
    nmut_vec(end+500) = 0;
end
nmut_vec(depth+1) = nmut_vec(depth+1)+1;
new_ID = [mut_ID,nmut_vec(depth+1)];
end

function Bout = Bcounter(Bcells,N_aff)
if isempty(Bcells)
    Bout = zeros(N_aff,1);
    return
end
inds = min(max(cell2mat(Bcells(:,1)),1),N_aff);
Bout = accumarray(inds,1,[N_aff,1]);
end

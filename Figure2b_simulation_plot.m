% FIGURE2_PANEL_B  Reproduce Figure 2b from the stochastic GC model.
%
% Steady-state GC B-cell and TFH-cell populations are simulated for a
% 3-by-3 sweep of antigen availability and TFH sensitivity. Each plotted
% point is one stochastic replicate, summarized by its median population
% during the final steady-state window. The red dashed line is the model
% prediction B*/T* = delta_T/(lambda*delta_B).
%
% Requires Statistics and Machine Learning Toolbox (binornd).

clear; clc; close all

%% Parameters retained from simulation_code.m
N_aff = 150;
N_T = 3;
Tend = 14*24;                 % hours
dt = 3e-2;                    % hours
Rep = 3;

del_B0 = 1/8;
del_T0 = 1/18;
kappa = 5;
lambda = 0.077;

B0 = 100;
T0_total = 12;                % Figure 2 caption: 12 initial TFH cells
iseed = 15;
eps_aff = 0.5;
ppos = 0.05;
N_lineages = 50;

% Figure 2b sweep and visual encoding.
antigen_availability = [0.1 0.3 1];
tfh_sensitivity = [1 3 5];
marker_by_antigen = {'d','^','o'};      % low, medium, high antigen
marker_size = [5 7 9];
green = [0.08 0.28 0.12; 0.18 0.52 0.27; 0.60 0.82 0.61];

sample_every_days = 1;
steady_window_days = 3;
steps_per_sample = max(1,round(sample_every_days*24/dt));
M = round(Tend/dt) + 1;

% Fixed seed makes the panel exactly reproducible. Each condition/replicate
% nevertheless receives an independent random stream segment.
rng(2,'twister');

%% Run parameter sweep
nA = numel(antigen_availability);
nS = numel(tfh_sensitivity);
B_rep = nan(nA,nS,Rep);
T_rep = nan(nA,nS,Rep);

for ia = 1:nA
    alpha_A = 1000*del_B0*antigen_availability(ia);

    for is = 1:nS
        s_vec = tfh_sensitivity(is)*ones(N_T,1);

        for rr = 1:Rep
            fprintf('antigen %.1f, sensitivity %g, replicate %d/%d\n', ...
                antigen_availability(ia),tfh_sensitivity(is),rr,Rep);

            [B_daily,T_daily] = run_gc_replicate(N_aff,N_T,M,dt, ...
                del_B0,del_T0,kappa,lambda,B0,T0_total,iseed, ...
                eps_aff,ppos,N_lineages,alpha_A,s_vec,steps_per_sample);

            % Robust late-time estimate: median of daily counts over the
            % last three days (or all available samples after extinction).
            nlate = min(steady_window_days/sample_every_days+1, ...
                numel(B_daily));
            late = (numel(B_daily)-nlate+1):numel(B_daily);
            B_rep(ia,is,rr) = median(B_daily(late));
            T_rep(ia,is,rr) = median(T_daily(late));
        end
    end
end

B_ss = mean(B_rep,3,'omitnan');
T_ss = mean(T_rep,3,'omitnan');

%% Plot Figure 2b
fig = figure('Color','w','Units','inches','Position',[1 1 4.4 4.0]);
ax = axes(fig); hold(ax,'on'); box(ax,'on');

positive_T = T_rep(T_rep>0);
positive_B = B_rep(B_rep>0);
if isempty(positive_T) || isempty(positive_B)
    error('All simulations became extinct; no positive steady state to plot.')
end
xlo = min(1e1,10^floor(log10(min(positive_T))));
xhi = max(1e3,10^ceil(log10(max(positive_T))));
ylo = min(1e2,10^floor(log10(min(positive_B))));
yhi = max(1e4,10^ceil(log10(max(positive_B))));

rstar = del_T0/(lambda*del_B0);
xline = logspace(log10(xlo),log10(xhi),200);
plot(ax,xline,rstar*xline,'--','Color',[0.85 0.15 0.15], ...
    'LineWidth',1.25,'HandleVisibility','off');

for is = 1:nS
    for ia = 1:nA
        for rr = 1:Rep
            if T_rep(ia,is,rr)>0 && B_rep(ia,is,rr)>0
                plot(ax,T_rep(ia,is,rr),B_rep(ia,is,rr), ...
                    marker_by_antigen{ia},'LineStyle','none', ...
                    'MarkerSize',marker_size(ia), ...
                    'MarkerFaceColor',green(is,:), ...
                    'MarkerEdgeColor',[0.08 0.08 0.08], ...
                    'LineWidth',0.75,'HandleVisibility','off');
            end
        end
    end
end

set(ax,'XScale','log','YScale','log','FontName','Arial','FontSize',10, ...
    'TickDir','out','LineWidth',0.8);
xlabel(ax,'T_{FH} cells');
ylabel(ax,'GC B cells');
title(ax,'Model');
xlim(ax,[xlo xhi]);
ylim(ax,[ylo yhi]);
axis(ax,'square');

% Compact two-part legend reproducing the marker and color keys.
h = gobjects(nA+nS+1,1);
labels = cell(nA+nS+1,1);
for ia = 1:nA
    h(ia) = plot(ax,nan,nan,marker_by_antigen{ia},'LineStyle','none', ...
        'MarkerSize',marker_size(ia),'MarkerFaceColor',[0.72 0.72 0.72], ...
        'MarkerEdgeColor','k');
    labels{ia} = sprintf('Antigen %.1f',antigen_availability(ia));
end
for is = 1:nS
    h(nA+is) = plot(ax,nan,nan,'o','LineStyle','none','MarkerSize',6, ...
        'MarkerFaceColor',green(is,:),'MarkerEdgeColor','k');
    labels{nA+is} = sprintf('Sensitivity %g',tfh_sensitivity(is));
end
h(end) = plot(ax,nan,nan,'--','Color',[0.85 0.15 0.15],'LineWidth',1.25);
labels{end} = 'Analytic';
legend(ax,h,labels,'Location','northwest','Box','off','FontSize',8);

% Save both the plotted panel and all replicate-level values.
script_dir = fileparts(mfilename('fullpath'));
exportgraphics(fig,fullfile(script_dir,'Figure2b.png'),'Resolution',300);
save(fullfile(script_dir,'Figure2b_data.mat'),'B_rep','T_rep','B_ss','T_ss', ...
    'antigen_availability','tfh_sensitivity','rstar');

%% Local functions: stochastic dynamics from simulation_code.m
function [B_daily,T_daily] = run_gc_replicate(N_aff,N_T,M,dt, ...
    del_B0,del_T0,kappa,lambda,B0,T0_total,iseed,eps_aff,ppos, ...
    N_lineages,alpha_A,s_vec,steps_per_sample)

nvec = (1:N_aff)';
Bcells = cell(B0,2);
for j = 1:B0
    Bcells(j,:) = {iseed,1+mod(j,N_lineages)};
end
Btemp = Bcounter(Bcells,N_aff);

% Distribute the caption's 12 TFH cells evenly over the three clones.
Ttemp = floor(T0_total/N_T)*ones(N_T,1);
Ttemp(1:mod(T0_total,N_T)) = Ttemp(1:mod(T0_total,N_T))+1;

nmut_vec = zeros(500,1);
nmut_vec(1) = 1;
n_samples = floor((M-1)/steps_per_sample)+1;
B_daily = zeros(1,n_samples);
T_daily = zeros(1,n_samples);
B_daily(1) = sum(Btemp);
T_daily(1) = sum(Ttemp);
js = 1;

for it = 1:(M-1)
    if isempty(Bcells) || sum(Ttemp)==0
        % Extinction is absorbing; the remaining saved counts stay zero.
        break
    end

    nbar = sum(Btemp.*nvec)/sum(Btemp);
    s_av = sum(s_vec.*Ttemp)/sum(Ttemp);
    f_vec = 1./(1+exp(-eps_aff*(nvec-nbar)));
    psi_vec = 1./(1+1./(f_vec*s_vec'));

    Bscr = Bscr_finder(Btemp,f_vec,alpha_A,del_B0);
    Ctot = sum(Ttemp)+sum(Btemp);
    dB_birth = dt*kappa*Bscr.*sum(psi_vec.*Ttemp'/Ctot,2)./Btemp;
    dB_birth(isnan(dB_birth)) = 0;
    dB_death = dt*del_B0;
    dT_birth = dt*lambda*kappa*sum(psi_vec.*Bscr/Ctot,1)';
    dT_death = dt*del_T0;

    [Bcells,nmut_vec] = bmut(Bcells,dB_birth,eps_aff,s_av,nbar, ...
        nmut_vec,ppos,N_aff);
    Bcells = c_death(Bcells,dB_death);
    Btemp = Bcounter(Bcells,N_aff);

    Tnew = Ttemp;
    for k = 1:N_T
        Tnew(k) = Ttemp(k)+binornd(Ttemp(k),dT_birth(k)) ...
            -binornd(Ttemp(k),dT_death);
    end
    Tnew(Tnew<1) = 0;
    Ttemp = Tnew;

    if mod(it,steps_per_sample)==0
        js = js+1;
        B_daily(js) = sum(Btemp);
        T_daily(js) = sum(Ttemp);
    end
end
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
    denom_shift = del_B*(S-x);
    Bscr = (Btemp.*A)./(A+denom_shift);
    valid = Btemp>0.5;
    eps_rel = mean(abs(Bscr(valid)-Bscr0(valid))./Btemp(valid));
    if eps_rel<=tol && all(Bscr(Btemp>0)<=Btemp(Btemp>0)+1e-2)
        break
    end
    Bscr0 = Bscr;
end
end

function Bnew = c_death(Btemp_cell,dB_death)
N = size(Btemp_cell,1);
if N==0
    Bnew = Btemp_cell;
    return
end
Bnew = Btemp_cell(rand(N,1)>=dB_death,:);
end

function [Bnew,nmut_vec] = bmut(Btemp_cell,dB_birth,eps_aff, ...
    s_av,nbar,nmut_vec,ppos,N_aff)
N = size(Btemp_cell,1);
n_i = cell2mat(Btemp_cell(:,1));
birth = rand(N,1)<dB_birth(n_i);
Bnew = cell(N+sum(birth),2);
Bnew(1:N,:) = Btemp_cell;
ib = find(birth);
tail = N+(1:numel(ib));
Bnew(tail,1) = Btemp_cell(ib,1);
Bnew(tail,2) = Btemp_cell(ib,2);
for k = 1:numel(ib)
    i = ib(k);
    n = Btemp_cell{i,1};
    mut_ID = Btemp_cell{i,2};
    psi = 1/(1+(1+exp(-eps_aff*(n-nbar)))/s_av);
    pmut = 0.6-0.4*psi;
    pos = [i,tail(k)];
    for j = 1:2
        if rand<pmut
            [new_mut_ID,nmut_vec] = mutgen(mut_ID,nmut_vec);
            if rand>0.2
                new_n = n;
            elseif rand<ppos
                new_n = n+1;
            else
                new_n = n-1;
            end
            new_n = min(max(new_n,1),N_aff);
            Bnew(pos(j),:) = {new_n,new_mut_ID};
        end
    end
end
end

function [new_mut_ID,nmut_vec] = mutgen(mut_ID,nmut_vec)
n = length(mut_ID);
if n+1>numel(nmut_vec)
    nmut_vec(end+500) = 0;
end
nmut_vec(n+1) = nmut_vec(n+1)+1;
new_mut_ID = [mut_ID,nmut_vec(n+1)];
end

function Bout = Bcounter(Bcell,N_aff)
if isempty(Bcell)
    Bout = zeros(N_aff,1);
    return
end
inds = cell2mat(Bcell(:,1));
inds = min(max(inds,1),N_aff);
Bout = accumarray(inds,1,[N_aff,1]);
end

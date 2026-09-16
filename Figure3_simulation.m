% Figure 3 and Supplementary Figure 5: reciprocal TFH-B-cell simulation.
% Run this script in MATLAB R2025a. It contains all required local functions.
% Results are saved in a Figure3_output folder beside this script.
% The internal affinity index iseed=15 represents baseline affinity n=0.
% Set mutation_delay_day=5 only for the separate delayed-mutation analysis.

clear; clc; close all;

% ============================================================
%
% Conditions:
%
%   1) All TFH clones present from day 0
%        s = 1, 3, 5
%        initial TFH numbers = [4, 4, 4]
%
%   2) Only low-sensitivity TFH clone present
%        s = 1 only
%        initial TFH numbers = [12, 0, 0]
%
% Both conditions begin with 12 total TFH cells.
%
% Figure 3:
%   a) Antigen availability
%   b) TFH clone fractions: mean +/- SD shading
%   c) GC B-cell population: mean +/- SD shading
%
% Supplementary Figure 5:
%   a) Affinity shift epsilon*Delta n
%   b) B:TFH cell ratio
%
% Supplementary Figure 5:
%   - sampled every 2 days
%   - points show mean +/- SD
%   - no lines connecting points
%   - if B < 100, GC is treated as extinct
%
% ============================================================


%% ============================================================
% Parameters
% ============================================================

save_data = 0;

save_dir = fullfile(fileparts(mfilename('fullpath')), 'Figure3_output');
if ~isfolder(save_dir)
    mkdir(save_dir);
end


% ------------------------------------------------------------
% Simulation
% ------------------------------------------------------------

N_aff = 100;
N_T   = 3;

Tend = 40*24;
dt   = 3e-2;

M = round(Tend/dt) + 1;

Rep = 10;

t_days = (0:M-1)*dt/24;

nvec = (1:N_aff)';


% ------------------------------------------------------------
% Rates
% ------------------------------------------------------------

del_B0 = 1/8;
del_T0 = 1/18;

kappa  = 5;
lambda = 0.077;


% ------------------------------------------------------------
% Analytic B:TFH ratio
%
% r* = delta_T / (lambda * delta_B)
% ------------------------------------------------------------

ratio_analytic = ...
    del_T0 / (lambda * del_B0);


fprintf( ...
    'Analytic B:TFH ratio = %.3f\n', ...
    ratio_analytic);


% ------------------------------------------------------------
% GC extinction threshold
% ------------------------------------------------------------

B_extinct_threshold = 100;


% ------------------------------------------------------------
% Initial populations
% ------------------------------------------------------------

B0 = 100;

% All-TFH:
% [4,4,4]

T0_per_clone = 4;

% s=1-only:
% [12,0,0]

T0_total = 12;


% ------------------------------------------------------------
% TFH sensitivities
% ------------------------------------------------------------

s_vec = [1 3 5]';


% ------------------------------------------------------------
% Mutation delay
% ------------------------------------------------------------

mutation_delay_day = 0;  % 0 for the main figure; 5 for delayed-mutation analysis


% ------------------------------------------------------------
% Affinity / mutation
% ------------------------------------------------------------

iseed = 15;

eps  = 0.5;
ppos = 0.05;


% ------------------------------------------------------------
% Antigen
% ------------------------------------------------------------

alpha_A0 = ...
    1000 * del_B0;

Ag_mid   = 25;
Ag_width = 5;


% ------------------------------------------------------------
% Initial lineage labels
% ------------------------------------------------------------

N_lineages = 50;


%% ============================================================
% Storage
% ============================================================

% ------------------------------------------------------------
% All s = 1,3,5 condition
% ------------------------------------------------------------

Btot_all = ...
    nan(Rep,M);

Tnum_all = ...
    nan(Rep,N_T,M);

Ttot_all = ...
    nan(Rep,M);

Tfrac_all = ...
    nan(Rep,N_T,M);

aff_all = ...
    nan(Rep,M);


% ------------------------------------------------------------
% s = 1 only condition
% ------------------------------------------------------------

Btot_s1 = ...
    nan(Rep,M);

Ttot_s1 = ...
    nan(Rep,M);

aff_s1 = ...
    nan(Rep,M);


%% ============================================================
% Run both conditions
% ============================================================

total_runs = ...
    2*Rep;

run_count = ...
    0;

tic_total = tic;


for condition = 1:2

    if condition == 1

        condition_name = ...
            'All TFH clones: s = 1, 3, 5';

    else

        condition_name = ...
            'Only s = 1 TFH clone';

    end


    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('%s\n',condition_name);
    fprintf('====================================================\n');


    for rr = 1:Rep

        run_count = ...
            run_count + 1;


        fprintf('\n');

        fprintf( ...
            'Run %d/%d | replicate %d/%d\n', ...
            run_count,total_runs,rr,Rep);


        tic_run = tic;


        %% ----------------------------------------------------
        % Initialize B cells
        % -----------------------------------------------------

        Bcells = ...
            cell(B0,2);


        for j = 1:B0

            Bcells(j,:) = ...
                {iseed,1 + mod(j,N_lineages)};

        end


        Btemp = ...
            zeros(N_aff,1);

        Btemp(iseed) = ...
            B0;


        %% ----------------------------------------------------
        % Initialize TFH cells
        % -----------------------------------------------------

        if condition == 1

            % s = 1,3,5
            % four cells each

            Ttemp = ...
                T0_per_clone * ones(N_T,1);

        else

            % s = 1 only
            % all twelve TFH cells are s = 1

            Ttemp = ...
                zeros(N_T,1);

            Ttemp(1) = ...
                T0_total;

        end


        %% ----------------------------------------------------
        % Mutation bookkeeping
        % -----------------------------------------------------

        nmut_vec = ...
            zeros(500,1);

        nmut_vec(1) = ...
            1;


        %% ----------------------------------------------------
        % Initial state
        % -----------------------------------------------------

        if condition == 1

            Btot_all(rr,1) = ...
                sum(Btemp);


            Tnum_all(rr,:,1) = ...
                Ttemp;


            Ttot_all(rr,1) = ...
                sum(Ttemp);


            Tfrac_all(rr,:,1) = ...
                Ttemp / sum(Ttemp);


            aff_all(rr,1) = ...
                0;

        else

            Btot_s1(rr,1) = ...
                sum(Btemp);


            Ttot_s1(rr,1) = ...
                sum(Ttemp);


            aff_s1(rr,1) = ...
                0;

        end


        %% ====================================================
        % Main simulation
        % =====================================================

        progress_interval = ...
            round(5*24/dt);


        for it = 1:(M-1)

            day = ...
                t_days(it);


            % ------------------------------------------------
            % Extinction checks
            % ------------------------------------------------

            if isempty(Bcells) || sum(Ttemp) == 0

                fprintf( ...
                    '  extinct at day %.2f\n', ...
                    day);

                break

            end


            NB = ...
                sum(Btemp);


            if NB == 0

                fprintf( ...
                    '  B cells extinct at day %.2f\n', ...
                    day);

                break

            end


            % ------------------------------------------------
            % Current mean affinity
            % ------------------------------------------------

            nbar = ...
                sum(Btemp .* nvec) / NB;


            % ------------------------------------------------
            % Antigen availability
            % ------------------------------------------------

            Ag_frac = ...
                1 ./ ...
                (1 + exp((day-Ag_mid)/Ag_width));


            alpha_A = ...
                alpha_A0 * Ag_frac;


            % ------------------------------------------------
            % Mean TFH sensitivity
            % ------------------------------------------------

            s_av = ...
                sum(s_vec .* Ttemp) ...
                / max(1,sum(Ttemp));


            % ------------------------------------------------
            % Affinity-dependent interaction
            % ------------------------------------------------

            f_vec = ...
                1 ./ ...
                (1 + exp(-eps*(nvec-nbar)));


            psi_vec = ...
                1 ./ ...
                (1 + 1 ./ (f_vec .* s_vec'));


            % ------------------------------------------------
            % Antigen screening
            % ------------------------------------------------

            Bscr = ...
                Bscr_finder( ...
                Btemp, ...
                f_vec, ...
                alpha_A, ...
                del_B0);


            % ------------------------------------------------
            % Total interacting population
            % ------------------------------------------------

            Ctot = ...
                sum(Ttemp) + sum(Btemp);


            % ------------------------------------------------
            % B-cell birth probability
            % ------------------------------------------------

            dB_birth = ...
                dt * kappa .* ...
                Bscr .* ...
                sum(psi_vec .* Ttemp'/Ctot,2) ...
                ./ Btemp;


            dB_birth(~isfinite(dB_birth)) = ...
                0;


            dB_birth = ...
                min(max(dB_birth,0),1);


            % ------------------------------------------------
            % B-cell death probability
            % ------------------------------------------------

            dB_death = ...
                dt * del_B0;


            dB_death = ...
                min(max(dB_death,0),1);


            % ------------------------------------------------
            % TFH birth probability
            % ------------------------------------------------

            dT_birth = ...
                dt * lambda * kappa * ...
                sum( ...
                psi_vec .* Bscr / Ctot, ...
                1)';


            dT_birth = ...
                min(max(dT_birth,0),1);


            % ------------------------------------------------
            % TFH death probability
            % ------------------------------------------------

            dT_death = ...
                dt * del_T0;


            dT_death = ...
                min(max(dT_death,0),1);


            % ------------------------------------------------
            % Update B cells
            % ------------------------------------------------

            mutation_on = ...
                day >= mutation_delay_day;


            [Bcells,nmut_vec] = ...
                bmut_delay( ...
                Bcells, ...
                dB_birth, ...
                eps, ...
                s_av, ...
                nbar, ...
                nmut_vec, ...
                ppos, ...
                N_aff, ...
                mutation_on);


            Bcells = ...
                c_death( ...
                Bcells, ...
                dB_death);


            Btemp = ...
                Bcounter( ...
                Bcells, ...
                N_aff);


            % ------------------------------------------------
            % Mean affinity after B-cell update
            % ------------------------------------------------

            if sum(Btemp) > 0

                nbar_new = ...
                    sum(Btemp .* nvec) ...
                    / sum(Btemp);


                aff_shift = ...
                    eps * (nbar_new-iseed);

            else

                aff_shift = ...
                    NaN;

            end


            % ------------------------------------------------
            % Update TFH cells
            % ------------------------------------------------

            Tnew = ...
                Ttemp;


            for k = 1:N_T

                if Ttemp(k) > 0

                    births = ...
                        binornd( ...
                        Ttemp(k), ...
                        dT_birth(k));


                    deaths = ...
                        binornd( ...
                        Ttemp(k), ...
                        dT_death);


                    Tnew(k) = ...
                        Ttemp(k) ...
                        + births ...
                        - deaths;

                end

            end


            Tnew(Tnew < 1) = ...
                0;


            Ttemp = ...
                Tnew;


            % ------------------------------------------------
            % Save outputs
            % ------------------------------------------------

            if condition == 1

                Btot_all(rr,it+1) = ...
                    sum(Btemp);


                Tnum_all(rr,:,it+1) = ...
                    Ttemp;


                Ttot_all(rr,it+1) = ...
                    sum(Ttemp);


                if sum(Ttemp) > 0

                    Tfrac_all(rr,:,it+1) = ...
                        Ttemp / sum(Ttemp);

                end


                aff_all(rr,it+1) = ...
                    aff_shift;

            else

                Btot_s1(rr,it+1) = ...
                    sum(Btemp);


                Ttot_s1(rr,it+1) = ...
                    sum(Ttemp);


                aff_s1(rr,it+1) = ...
                    aff_shift;

            end


            % ------------------------------------------------
            % Progress
            % ------------------------------------------------

            if mod(it,progress_interval) == 0

                fprintf( ...
                    '  day %4.1f | B = %5d | T = %4d | mut=%d\n', ...
                    day, ...
                    sum(Btemp), ...
                    round(sum(Ttemp)), ...
                    mutation_on);

            end

        end


        %% ----------------------------------------------------
        % Timing
        % -----------------------------------------------------

        elapsed_run = ...
            toc(tic_run);


        elapsed_all = ...
            toc(tic_total);


        avg_run = ...
            elapsed_all/run_count;


        remaining = ...
            avg_run*(total_runs-run_count);


        fprintf( ...
            '  Run time: %.1f s\n', ...
            elapsed_run);


        fprintf( ...
            '  Overall progress: %.1f%%\n', ...
            100*run_count/total_runs);


        fprintf( ...
            '  Estimated remaining: %.1f min\n', ...
            remaining/60);

    end

end


%% ============================================================
% Supplementary Figure 5 quantities
%
% If B < 100, treat GC as extinct for Supplementary Figure 5.
% ============================================================

aff_fig2_all = ...
    aff_all;

aff_fig2_s1 = ...
    aff_s1;


ratio_all = ...
    Btot_all ./ Ttot_all;


ratio_s1 = ...
    Btot_s1 ./ Ttot_s1;


ratio_all(~isfinite(ratio_all)) = ...
    NaN;

ratio_s1(~isfinite(ratio_s1)) = ...
    NaN;


% ------------------------------------------------------------
% Extinction masks
% ------------------------------------------------------------

extinct_all = ...
    Btot_all < B_extinct_threshold;


extinct_s1 = ...
    Btot_s1 < B_extinct_threshold;


aff_fig2_all(extinct_all) = ...
    NaN;


aff_fig2_s1(extinct_s1) = ...
    NaN;


ratio_all(extinct_all) = ...
    NaN;


ratio_s1(extinct_s1) = ...
    NaN;


%% ============================================================
% Mean and SD trajectories
% ============================================================

% ------------------------------------------------------------
% GC B-cell population
% ------------------------------------------------------------

Bmean_all = ...
    mean(Btot_all,1,'omitnan');


Bmean_s1 = ...
    mean(Btot_s1,1,'omitnan');


Bsd_all = ...
    std(Btot_all,0,1,'omitnan');


Bsd_s1 = ...
    std(Btot_s1,0,1,'omitnan');


% ------------------------------------------------------------
% TFH fractions
% ------------------------------------------------------------

Tfrac_mean = ...
    squeeze( ...
    mean(Tfrac_all,1,'omitnan'));


Tfrac_sd = ...
    squeeze( ...
    std(Tfrac_all,0,1,'omitnan'));


% ------------------------------------------------------------
% Affinity
% ------------------------------------------------------------

aff_mean_all = ...
    mean(aff_fig2_all,1,'omitnan');


aff_mean_s1 = ...
    mean(aff_fig2_s1,1,'omitnan');


% ------------------------------------------------------------
% B:TFH ratio
% ------------------------------------------------------------

ratio_mean_all = ...
    mean(ratio_all,1,'omitnan');


ratio_mean_s1 = ...
    mean(ratio_s1,1,'omitnan');


%% ============================================================
% Antigen trajectory
% ============================================================

Ag_curve = ...
    1 ./ ...
    (1 + exp((t_days-Ag_mid)/Ag_width));


%% ============================================================
% Figure 3a sampling
%
% Every 6 hours
% ============================================================

plot_step = ...
    round(6/dt);


plot_idx = ...
    1:plot_step:M;


if plot_idx(end) ~= M

    plot_idx = ...
        [plot_idx M];

end


t_plot = ...
    t_days(plot_idx);


Ag_plot = ...
    Ag_curve(plot_idx);


%% ============================================================
% Supplementary Figure 5 sampling
%
% Every 2 days
% ============================================================

fig2_step = ...
    round(2*24/dt);


fig2_idx = ...
    1:fig2_step:M;


if fig2_idx(end) ~= M

    fig2_idx = ...
        [fig2_idx M];

end


t_fig2 = ...
    t_days(fig2_idx);


% ------------------------------------------------------------
% Affinity
% ------------------------------------------------------------

aff_fig2_all_sample = ...
    aff_fig2_all(:,fig2_idx);


aff_fig2_s1_sample = ...
    aff_fig2_s1(:,fig2_idx);


aff_fig2_mean_all = ...
    mean( ...
    aff_fig2_all_sample, ...
    1, ...
    'omitnan');


aff_fig2_mean_s1 = ...
    mean( ...
    aff_fig2_s1_sample, ...
    1, ...
    'omitnan');


aff_fig2_sd_all = ...
    std( ...
    aff_fig2_all_sample, ...
    0, ...
    1, ...
    'omitnan');


aff_fig2_sd_s1 = ...
    std( ...
    aff_fig2_s1_sample, ...
    0, ...
    1, ...
    'omitnan');


% ------------------------------------------------------------
% B:TFH ratio
% ------------------------------------------------------------

ratio_fig2_all_sample = ...
    ratio_all(:,fig2_idx);


ratio_fig2_s1_sample = ...
    ratio_s1(:,fig2_idx);


ratio_fig2_mean_all = ...
    mean( ...
    ratio_fig2_all_sample, ...
    1, ...
    'omitnan');


ratio_fig2_mean_s1 = ...
    mean( ...
    ratio_fig2_s1_sample, ...
    1, ...
    'omitnan');


ratio_fig2_sd_all = ...
    std( ...
    ratio_fig2_all_sample, ...
    0, ...
    1, ...
    'omitnan');


ratio_fig2_sd_s1 = ...
    std( ...
    ratio_fig2_s1_sample, ...
    0, ...
    1, ...
    'omitnan');


%% ============================================================
% Save sampling
%
% Every day
% ============================================================

save_step = ...
    round(24/dt);


save_idx = ...
    1:save_step:M;


if save_idx(end) ~= M

    save_idx = ...
        [save_idx M];

end


t_save = ...
    t_days(save_idx);


%% ============================================================
% Colors
% ============================================================

green_dark = ...
    [0.08 0.28 0.08];


green_mid = ...
    [0.25 0.58 0.25];


green_light = ...
    [0.58 0.78 0.50];


blue = ...
    [0.15 0.45 0.78];


black = ...
    [0 0 0];


red = ...
    [0.78 0.12 0.12];


%% ============================================================
% Plot formatting
% ============================================================

mean_line_width = ...
    2.5;


antigen_line_width = ...
    2.5;


errorbar_line_width = ...
    1.25;


marker_size = ...
    4;


errorbar_cap_size = ...
    5;


shade_alpha = ...
    0.20;


%% ============================================================
% FIGURE 3
%
% a) Antigen
% b) TFH fractions: mean +/- SD shading
% c) GC B-cell population: mean +/- SD shading
% ============================================================

figure( ...
    'Color','w', ...
    'Position',[200 80 620 900]);


tiledlayout( ...
    3,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ============================================================
% Figure 3a: Antigen
% ============================================================

nexttile;
hold on;


plot( ...
    t_plot, ...
    Ag_plot, ...
    '-', ...
    'Color',red, ...
    'LineWidth',antigen_line_width);


ylabel('Antigen (a.u.)');


xlim([0 Tend/24]);
ylim([0 1]);


box off;


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'XTickLabel',[]);


text( ...
    -0.10,1.05,'a', ...
    'Units','normalized', ...
    'FontSize',16, ...
    'FontWeight','bold');


%% ============================================================
% Figure 3b: TFH clone fractions
%
% Mean +/- SD shaded regions
% ============================================================

nexttile;
hold on;


% ------------------------------------------------------------
% s = 1 SD region
% ------------------------------------------------------------

upper_s1 = ...
    min( ...
    Tfrac_mean(1,:) + Tfrac_sd(1,:), ...
    1);


lower_s1 = ...
    max( ...
    Tfrac_mean(1,:) - Tfrac_sd(1,:), ...
    0);


fill( ...
    [t_days fliplr(t_days)], ...
    [upper_s1 fliplr(lower_s1)], ...
    green_dark, ...
    'FaceAlpha',shade_alpha, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');


% ------------------------------------------------------------
% s = 3 SD region
% ------------------------------------------------------------

upper_s3 = ...
    min( ...
    Tfrac_mean(2,:) + Tfrac_sd(2,:), ...
    1);


lower_s3 = ...
    max( ...
    Tfrac_mean(2,:) - Tfrac_sd(2,:), ...
    0);


fill( ...
    [t_days fliplr(t_days)], ...
    [upper_s3 fliplr(lower_s3)], ...
    green_mid, ...
    'FaceAlpha',shade_alpha, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');


% ------------------------------------------------------------
% s = 5 SD region
% ------------------------------------------------------------

upper_s5 = ...
    min( ...
    Tfrac_mean(3,:) + Tfrac_sd(3,:), ...
    1);


lower_s5 = ...
    max( ...
    Tfrac_mean(3,:) - Tfrac_sd(3,:), ...
    0);


fill( ...
    [t_days fliplr(t_days)], ...
    [upper_s5 fliplr(lower_s5)], ...
    green_light, ...
    'FaceAlpha',shade_alpha, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');


% ------------------------------------------------------------
% Mean trajectories
% ------------------------------------------------------------

h_s1 = ...
    plot( ...
    t_days, ...
    Tfrac_mean(1,:), ...
    '-', ...
    'Color',green_dark, ...
    'LineWidth',mean_line_width);


h_s3 = ...
    plot( ...
    t_days, ...
    Tfrac_mean(2,:), ...
    '-', ...
    'Color',green_mid, ...
    'LineWidth',mean_line_width);


h_s5 = ...
    plot( ...
    t_days, ...
    Tfrac_mean(3,:), ...
    '-', ...
    'Color',green_light, ...
    'LineWidth',mean_line_width);


ylabel('T_{FH} cell fraction');


xlim([0 Tend/24]);
ylim([0 1]);


legend( ...
    [h_s5 h_s3 h_s1], ...
    {'s = 5','s = 3','s = 1'}, ...
    'Location','east', ...
    'Box','off');


box off;


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'XTickLabel',[]);


text( ...
    -0.10,1.05,'b', ...
    'Units','normalized', ...
    'FontSize',16, ...
    'FontWeight','bold');


%% ============================================================
% Figure 3c: GC B-cell population
%
% Mean +/- SD shaded regions
% ============================================================

nexttile;
hold on;


% ------------------------------------------------------------
% All-TFH SD region
%
% Lower bound is clipped at 1 because this is plotted on
% a logarithmic axis.
% ------------------------------------------------------------

upper_B_all = ...
    Bmean_all + Bsd_all;


lower_B_all = ...
    max( ...
    Bmean_all - Bsd_all, ...
    1);


fill( ...
    [t_days fliplr(t_days)], ...
    [upper_B_all fliplr(lower_B_all)], ...
    blue, ...
    'FaceAlpha',shade_alpha, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');


% ------------------------------------------------------------
% s=1-only SD region
% ------------------------------------------------------------

upper_B_s1 = ...
    Bmean_s1 + Bsd_s1;


lower_B_s1 = ...
    max( ...
    Bmean_s1 - Bsd_s1, ...
    1);


fill( ...
    [t_days fliplr(t_days)], ...
    [upper_B_s1 fliplr(lower_B_s1)], ...
    black, ...
    'FaceAlpha',shade_alpha, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');


% ------------------------------------------------------------
% Mean trajectories
% ------------------------------------------------------------

h_all = ...
    plot( ...
    t_days, ...
    Bmean_all, ...
    '-', ...
    'Color',blue, ...
    'LineWidth',mean_line_width);


h_s1only = ...
    plot( ...
    t_days, ...
    Bmean_s1, ...
    '-', ...
    'Color',black, ...
    'LineWidth',mean_line_width);


set(gca, ...
    'YScale','log');


xlabel('Days');
ylabel('GC B cells');


xlim([0 Tend/24]);


legend( ...
    [h_all h_s1only], ...
    {'s = 1, 3, 5','s = 1 only'}, ...
    'Location','northeast', ...
    'Box','off');


box off;


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.2);


text( ...
    -0.10,1.05,'c', ...
    'Units','normalized', ...
    'FontSize',16, ...
    'FontWeight','bold');


if save_data == 1
    savefig(gcf, fullfile(save_dir, 'Figure3.fig'));
    exportgraphics(gcf, fullfile(save_dir, 'Figure3.pdf'), ...
        'ContentType', 'vector');
end

%% ============================================================
% SUPPLEMENTARY FIGURE 5
%
% Every 2 days
% Mean +/- SD
% Markers only; no connecting lines
% ============================================================

figure( ...
    'Color','w', ...
    'Position',[250 250 900 390]);


tiledlayout( ...
    1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% ============================================================
% Supplementary Figure 5a: Affinity shift
% ============================================================

nexttile;
hold on;


h_aff_all = ...
    errorbar( ...
    t_fig2, ...
    aff_fig2_mean_all, ...
    aff_fig2_sd_all, ...
    'o', ...
    'LineStyle','none', ...
    'Color',blue, ...
    'MarkerFaceColor',blue, ...
    'MarkerSize',marker_size, ...
    'LineWidth',errorbar_line_width, ...
    'CapSize',errorbar_cap_size);


h_aff_s1 = ...
    errorbar( ...
    t_fig2, ...
    aff_fig2_mean_s1, ...
    aff_fig2_sd_s1, ...
    'o', ...
    'LineStyle','none', ...
    'Color',black, ...
    'MarkerFaceColor',black, ...
    'MarkerSize',marker_size, ...
    'LineWidth',errorbar_line_width, ...
    'CapSize',errorbar_cap_size);


xlabel('Days');


ylabel( ...
    '$\varepsilon\,\Delta n$', ...
    'Interpreter','latex');


xlim([0 Tend/24]);


legend( ...
    [h_aff_all h_aff_s1], ...
    {'s = 1, 3, 5','s = 1 only'}, ...
    'Location','northwest', ...
    'Box','off');


box off;


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.2);


text( ...
    -0.15,1.05,'a', ...
    'Units','normalized', ...
    'FontSize',16, ...
    'FontWeight','bold');


%% ============================================================
% Supplementary Figure 5b: B:TFH ratio
% ============================================================

nexttile;
hold on;


% ------------------------------------------------------------
% Analytic prediction
% ------------------------------------------------------------

h_analytic = ...
    yline( ...
    ratio_analytic, ...
    ':', ...
    'Color',red, ...
    'LineWidth',2.5);


% ------------------------------------------------------------
% Simulation results
% ------------------------------------------------------------

h_ratio_all = ...
    errorbar( ...
    t_fig2, ...
    ratio_fig2_mean_all, ...
    ratio_fig2_sd_all, ...
    'o', ...
    'LineStyle','none', ...
    'Color',blue, ...
    'MarkerFaceColor',blue, ...
    'MarkerSize',marker_size, ...
    'LineWidth',errorbar_line_width, ...
    'CapSize',errorbar_cap_size);


h_ratio_s1 = ...
    errorbar( ...
    t_fig2, ...
    ratio_fig2_mean_s1, ...
    ratio_fig2_sd_s1, ...
    'o', ...
    'LineStyle','none', ...
    'Color',black, ...
    'MarkerFaceColor',black, ...
    'MarkerSize',marker_size, ...
    'LineWidth',errorbar_line_width, ...
    'CapSize',errorbar_cap_size);


xlabel('Days');
ylabel('B:T_{FH} ratio');


xlim([0 Tend/24]);


legend( ...
    [h_ratio_all h_ratio_s1 h_analytic], ...
    {'s = 1, 3, 5', ...
     's = 1 only', ...
     'Analytic prediction'}, ...
    'Location','northeast', ...
    'Box','off');


box off;


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.2);


text( ...
    -0.15,1.05,'b', ...
    'Units','normalized', ...
    'FontSize',16, ...
    'FontWeight','bold');


if save_data == 1
    savefig(gcf, fullfile(save_dir, 'Supplementary_Figure_5.fig'));
    exportgraphics(gcf, fullfile(save_dir, 'Supplementary_Figure_5.pdf'), ...
        'ContentType', 'vector');
end

%% ============================================================
% Save source data
% ============================================================

if save_data == 1

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('Saving source data\n');
    fprintf('====================================================\n');
    fprintf('Folder:\n%s\n\n',save_dir);


    %% --------------------------------------------------------
    % Daily sampled arrays
    % ---------------------------------------------------------

    Ag_save = ...
        Ag_curve(save_idx);


    Btot_all_save = ...
        Btot_all(:,save_idx);


    Btot_s1_save = ...
        Btot_s1(:,save_idx);


    Tfrac_all_save = ...
        Tfrac_all(:,:,save_idx);


    Ttot_all_save = ...
        Ttot_all(:,save_idx);


    Ttot_s1_save = ...
        Ttot_s1(:,save_idx);


    aff_all_save = ...
        aff_fig2_all(:,save_idx);


    aff_s1_save = ...
        aff_fig2_s1(:,save_idx);


    ratio_all_save = ...
        ratio_all(:,save_idx);


    ratio_s1_save = ...
        ratio_s1(:,save_idx);


    %% ========================================================
    % Common columns
    % =========================================================

    n_save = ...
        length(t_save);


    rep_col = ...
        repelem((1:Rep)',n_save);


    day_col = ...
        repmat(t_save(:),Rep,1);


    %% ========================================================
    % 1) GC B-cell population
    % =========================================================

    B_all_col = ...
        reshape( ...
        Btot_all_save.', ...
        [],1);


    B_s1_col = ...
        reshape( ...
        Btot_s1_save.', ...
        [],1);


    TB = table( ...
        rep_col, ...
        day_col, ...
        B_all_col, ...
        B_s1_col, ...
        'VariableNames', { ...
        'replicate', ...
        'day', ...
        'B_all_TFH', ...
        'B_s1_only'});


    writetable( ...
        TB, ...
        fullfile( ...
        save_dir, ...
        'Figure3_Bcell_population_raw.csv'));


    %% ========================================================
    % 2) TFH fractions
    % =========================================================

    TFH_s1_col = ...
        reshape( ...
        squeeze(Tfrac_all_save(:,1,:)).', ...
        [],1);


    TFH_s3_col = ...
        reshape( ...
        squeeze(Tfrac_all_save(:,2,:)).', ...
        [],1);


    TFH_s5_col = ...
        reshape( ...
        squeeze(Tfrac_all_save(:,3,:)).', ...
        [],1);


    TT = table( ...
        rep_col, ...
        day_col, ...
        TFH_s1_col, ...
        TFH_s3_col, ...
        TFH_s5_col, ...
        'VariableNames', { ...
        'replicate', ...
        'day', ...
        'TFH_fraction_s1', ...
        'TFH_fraction_s3', ...
        'TFH_fraction_s5'});


    writetable( ...
        TT, ...
        fullfile( ...
        save_dir, ...
        'Figure3_TFH_fractions_raw.csv'));


    %% ========================================================
    % 3) Total TFH populations
    % =========================================================

    T_all_col = ...
        reshape( ...
        Ttot_all_save.', ...
        [],1);


    T_s1_col = ...
        reshape( ...
        Ttot_s1_save.', ...
        [],1);


    TTtot = table( ...
        rep_col, ...
        day_col, ...
        T_all_col, ...
        T_s1_col, ...
        'VariableNames', { ...
        'replicate', ...
        'day', ...
        'TFH_total_all_TFH', ...
        'TFH_total_s1_only'});


    writetable( ...
        TTtot, ...
        fullfile( ...
        save_dir, ...
        'Figure3_TFH_total_raw.csv'));


    %% ========================================================
    % 4) Affinity
    % =========================================================

    aff_all_col = ...
        reshape( ...
        aff_all_save.', ...
        [],1);


    aff_s1_col = ...
        reshape( ...
        aff_s1_save.', ...
        [],1);


    TA = table( ...
        rep_col, ...
        day_col, ...
        aff_all_col, ...
        aff_s1_col, ...
        'VariableNames', { ...
        'replicate', ...
        'day', ...
        'epsilon_Delta_n_all_TFH', ...
        'epsilon_Delta_n_s1_only'});


    writetable( ...
        TA, ...
        fullfile( ...
        save_dir, ...
        'Figure3_affinity_raw.csv'));


    %% ========================================================
    % 5) B:TFH ratio
    % =========================================================

    ratio_all_col = ...
        reshape( ...
        ratio_all_save.', ...
        [],1);


    ratio_s1_col = ...
        reshape( ...
        ratio_s1_save.', ...
        [],1);


    TR = table( ...
        rep_col, ...
        day_col, ...
        ratio_all_col, ...
        ratio_s1_col, ...
        'VariableNames', { ...
        'replicate', ...
        'day', ...
        'B_to_TFH_all_TFH', ...
        'B_to_TFH_s1_only'});


    writetable( ...
        TR, ...
        fullfile( ...
        save_dir, ...
        'Figure3_B_to_TFH_ratio_raw.csv'));


    %% ========================================================
    % 6) Supplementary Figure 5 summary
    %
    % Sampled every two days.
    % Mean +/- SD.
    % =========================================================

    Tfig2 = table( ...
        t_fig2(:), ...
        aff_fig2_mean_all(:), ...
        aff_fig2_sd_all(:), ...
        aff_fig2_mean_s1(:), ...
        aff_fig2_sd_s1(:), ...
        ratio_fig2_mean_all(:), ...
        ratio_fig2_sd_all(:), ...
        ratio_fig2_mean_s1(:), ...
        ratio_fig2_sd_s1(:), ...
        repmat( ...
        ratio_analytic, ...
        length(t_fig2), ...
        1), ...
        'VariableNames', { ...
        'day', ...
        'affinity_mean_all_TFH', ...
        'affinity_SD_all_TFH', ...
        'affinity_mean_s1_only', ...
        'affinity_SD_s1_only', ...
        'ratio_mean_all_TFH', ...
        'ratio_SD_all_TFH', ...
        'ratio_mean_s1_only', ...
        'ratio_SD_s1_only', ...
        'ratio_analytic'});


    writetable( ...
        Tfig2, ...
        fullfile( ...
        save_dir, ...
        'Supplementary_Figure_5_2day_summary.csv'));


    %% ========================================================
    % 7) Antigen
    % =========================================================

    TAntigen = ...
        table( ...
        t_save(:), ...
        Ag_save(:), ...
        'VariableNames', { ...
        'day', ...
        'antigen'});


    writetable( ...
        TAntigen, ...
        fullfile( ...
        save_dir, ...
        'Figure3_antigen.csv'));


    %% ========================================================
    % Save workspace
    % =========================================================

    save( ...
        fullfile( ...
        save_dir, ...
        'Figure3_simulation_workspace.mat'));


    fprintf('\nSaved files:\n\n');


    dir( ...
        fullfile( ...
        save_dir, ...
        'Figure3_*'));

else

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('save_data = 0: no source data were saved.\n');
    fprintf('====================================================\n');

end


fprintf( ...
    '\nTotal runtime: %.1f min\n', ...
    toc(tic_total)/60);


%% ============================================================
% Functions
% ============================================================

function Bscr = Bscr_finder( ...
    Btemp, ...
    f_vec, ...
    alpha_A, ...
    del_B)

Btemp = ...
    Btemp(:);

f_vec = ...
    f_vec(:);


A = ...
    alpha_A .* f_vec;


S = ...
    sum(f_vec .* Btemp);


Bscr0 = ...
    ones(size(Btemp));


tol   = 1e-2;
maxit = 200;


for it = 1:maxit

    x = ...
        sum(f_vec .* Bscr0);


    denom_shift = ...
        del_B * (S-x);


    Bscr = ...
        (Btemp .* A) ./ ...
        (A + denom_shift);


    valid = ...
        Btemp > 0.5;


    if any(valid)

        eps_rel = ...
            mean( ...
            abs(Bscr(valid)-Bscr0(valid)) ...
            ./ Btemp(valid));


        if eps_rel <= tol

            if all( ...
                Bscr(Btemp>0) ...
                <= Btemp(Btemp>0)+1e-2)

                break

            end

        end

    end


    Bscr0 = ...
        Bscr;

end

end


function Bnew = ...
    c_death( ...
    Btemp_cell, ...
    dB_death)

N = ...
    size(Btemp_cell,1);


if N == 0

    Bnew = ...
        Btemp_cell;

    return

end


survive_mask = ...
    rand(N,1) >= dB_death;


Bnew = ...
    Btemp_cell(survive_mask,:);

end


function [Bnew,nmut_vec] = ...
    bmut_delay( ...
    Btemp_cell, ...
    dB_birth, ...
    eps, ...
    s_av, ...
    nbar, ...
    nmut_vec, ...
    ppos, ...
    N_aff, ...
    mutation_on)

N = ...
    size(Btemp_cell,1);


if N == 0

    Bnew = ...
        Btemp_cell;

    return

end


% ------------------------------------------------------------
% Determine which B cells divide
% ------------------------------------------------------------

n_i = ...
    cell2mat(Btemp_cell(:,1));


pb = ...
    dB_birth(n_i);


pb = ...
    min(max(pb,0),1);


birth = ...
    rand(N,1) < pb;


Nout = ...
    N + sum(birth);


% ------------------------------------------------------------
% Existing cells
% ------------------------------------------------------------

Bnew = ...
    cell(Nout,2);


Bnew(1:N,:) = ...
    Btemp_cell;


% ------------------------------------------------------------
% Daughter cells
% ------------------------------------------------------------

ib = ...
    find(birth);


nb = ...
    numel(ib);


tail = ...
    N + (1:nb);


Bnew(tail,1) = ...
    Btemp_cell(ib,1);


Bnew(tail,2) = ...
    Btemp_cell(ib,2);


% ------------------------------------------------------------
% Mutation off
% ------------------------------------------------------------

if ~mutation_on

    return

end


% ------------------------------------------------------------
% Mutation on
% ------------------------------------------------------------

for k = 1:nb

    i = ...
        ib(k);


    n = ...
        Btemp_cell{i,1};


    mut_ID = ...
        Btemp_cell{i,2};


    psi = ...
        1 / ...
        (1 + ...
        (1 + exp(-eps*(n-nbar))) / s_av);


    pmut = ...
        0.6 - 0.4*psi;


    pos = ...
        [i,tail(k)];


    for j = 1:2

        if rand < pmut

            [new_mut_ID,nmut_vec] = ...
                mutgen( ...
                mut_ID, ...
                nmut_vec);


            if rand > 0.2

                new_n = ...
                    n;

            else

                if rand < ppos

                    new_n = ...
                        n + 1;

                else

                    new_n = ...
                        n - 1;

                end

            end


            new_n = ...
                min(max(new_n,1),N_aff);


            Bnew(pos(j),:) = ...
                {new_n,new_mut_ID};

        end

    end

end

end


function [new_mut_ID,nmut_vec] = ...
    mutgen( ...
    mut_ID, ...
    nmut_vec)

n = ...
    length(mut_ID);


if n+1 > length(nmut_vec)

    nmut_vec(end+500) = ...
        0;

end


nmut_vec(n+1) = ...
    nmut_vec(n+1) + 1;


new_mut_ID = ...
    [mut_ID,nmut_vec(n+1)];

end


function Bout = ...
    Bcounter( ...
    Bcell, ...
    N_aff)

if isempty(Bcell)

    Bout = ...
        zeros(N_aff,1);

    return

end


inds = ...
    cell2mat(Bcell(:,1));


inds = ...
    min(max(inds,1),N_aff);


Bout = ...
    accumarray( ...
    inds, ...
    1, ...
    [N_aff,1]);

end
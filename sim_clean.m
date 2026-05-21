clear; clc

% =========================================================
% Parameters
% =========================================================

% Simulation
N_aff  = 150;          % number of affinity classes
N_T    = 3;            % number of TFH clones
Tend   = 35*24;        % duration (hours)
dt     = 3e-2;         % timestep (hours)
Rep    = 1;            % number of replicates

M      = round(Tend/dt) + 1;
t_days = (0:M-1) * dt / 24;
nvec   = (1:N_aff)';

% Rates
del_B0 = 1/8;          % B cell apoptosis rate
del_T0 = 1/18;         % TFH exit rate
kappa  = 5;            % max B cell proliferation rate
lambda = 0.077;        % TFH proliferation fraction

% Initial populations
B0 = 100;
T0 = 10;

% Affinity/mutation
iseed = 15;            % initial affinity class
eps   = 0.5;           % affinity spacing
ppos  = 0.05;          % probability mutation is beneficial

% Antigen and TFH sensitivity
alpha_A0 = 1000 * del_B0;
s_vec    = ones(N_T,1);

% Initial lineage labels
N_lineages = 50;

% Plotting
PrintFlag  = 1;
Uday       = 1;
ind_update = max(1, round(Uday*24/dt));
Mr         = round(M/ind_update) + 1;

% =========================================================
% Analytic steady-state estimates
% =========================================================

tau   = del_B0/(kappa*alpha_A0);
psi0  = 1/(1 + 2/max(s_vec));
rstar = del_T0/(lambda*del_B0);

Bstar = lambda*psi0/(lambda*del_B0*tau + del_T0*tau);
Tstar = lambda^2*psi0*del_B0*tau / ...
    ((del_T0*tau)*(lambda*del_B0*tau + del_T0*tau));

% =========================================================
% Run simulations
% =========================================================

out = cell(Rep,1);

for rr = 1:Rep
    fprintf('rep %d/%d...\n', rr, Rep)

    % -----------------------------
    % Initialize B cells
    % -----------------------------
    B = cell(M,1);
    B{1} = cell(B0,2);

    for j = 1:B0
        B{1}(j,:) = {iseed, 1 + mod(j,N_lineages)};
    end

    Bcc = zeros(N_aff, M);
    Bcc(iseed,1) = B0;

    % -----------------------------
    % Initialize TFH cells
    % -----------------------------
    T = zeros(N_T, M);
    T(1,1) = T0;

    % -----------------------------
    % Outputs
    % -----------------------------
    ccM          = 0;
    n_av         = zeros(M-1,1);
    nmut_vec     = zeros(500,1);
    nmut_vec(1)  = 1;

    shannon_node = zeros(Mr,1);
    mut_count    = zeros(Mr,1);
    lineage_count = zeros(Tend/24 + 1, N_lineages);

    % =====================================================
    % Main simulation loop
    % =====================================================
    for it = 1:(M-1)

        Bcells = B{it};
        Btemp  = Bcc(:,it);
        Ttemp  = T(:,it);

        if isempty(Bcells) || sum(Ttemp) == 0
            break
        end

        nbar     = sum(Btemp .* nvec) / sum(Btemp);
        n_av(it) = nbar;

        del_B = del_B0;
        del_T = del_T0;

        % -----------------------------
        % Lineage statistics
        % -----------------------------
        doStats = (~mod(it-1,ind_update) || it == M-1);

        if doStats
            ccM = ccM + 1;

            C = Bcells(:,2);
            root_ID = cellfun(@(x) x(1), C);

            lineage_count(ccM,:) = histcounts(root_ID, 0.5:(N_lineages+0.5));

            C_str = cellfun(@(x) sprintf('%g_', x), C, 'UniformOutput', false);
            [uStr, ~, ic] = unique(C_str, 'stable');
            clone_sizes = accumarray(ic, 1);

            shannon_node(ccM) = simpson_effnum(clone_sizes);

            nUnd = cellfun(@(s) sum(s == '_'), C_str);
            mut_count(ccM) = mean(nUnd);
        end

        % -----------------------------
        % Affinity-dependent rates
        % -----------------------------
        s_av = sum(s_vec .* Ttemp) / max(1,sum(Ttemp));

        f_vec   = 1 ./ (1 + exp(-eps*(nvec - nbar)));
        psi_vec = 1 ./ (1 + 1 ./ (f_vec .* s_vec'));

        Bscr = Bscr_finder(Btemp, f_vec, alpha_A0, del_B);
        Ctot = sum(Ttemp) + sum(Btemp);

        dB_birth = dt * kappa * Bscr .* sum(psi_vec .* Ttemp' / Ctot, 2) ./ Btemp;
        dB_birth(isnan(dB_birth)) = 0;

        dB_death = dt * del_B;

        dT_birth = dt * lambda * kappa * sum(psi_vec .* Bscr / Ctot, 1)';
        dT_death = dt * del_T;

        % -----------------------------
        % Update B cells
        % -----------------------------
        [B_birth_mut, nmut_vec] = bmut( ...
            Bcells, dB_birth, eps, s_av, nbar, nmut_vec, ppos, N_aff);

        B_next = c_death(B_birth_mut, dB_death);

        B{it+1}     = B_next;
        Bcc(:,it+1) = Bcounter(B_next, N_aff);

        % -----------------------------
        % Update TFH cells
        % -----------------------------
        Tnew = Ttemp;

        for k = 1:N_T
            Tnew(k) = Ttemp(k) ...
                + binornd(Ttemp(k), dT_birth(k)) ...
                - binornd(Ttemp(k), dT_death);
        end

        Tnew(Tnew < 1) = 0;
        T(:,it+1) = Tnew;

        % -----------------------------
        % Plot during simulation
        % -----------------------------
        if PrintFlag && doStats
            inds = 1:max(1,round(ind_update/3)):it;

            subplot(2,1,1)
            bar((1:N_aff)-iseed, Btemp, 'b')
            xlim([-10 20])
            ylim([0, Bstar/2])
            xlabel('Affinity class')
            ylabel('# B cells')

            subplot(2,1,2); cla
            plot(t_days(inds), sum(Bcc(:,inds),1), 'b.-'); hold on
            plot(t_days(inds), sum(T(:,inds),1), 'g.-')
            plot(xlim, [1 1]*Bstar, 'b-.')
            plot(xlim, [1 1]*Tstar, 'g-.')
            set(gca,'YScale','log')
            yticks([1e1 1e2 1e3 1e4])
            xlabel('Days')
            ylabel('# cells')
            xlim([0 t_days(end)])
            ylim([1e1 1e4])

            drawnow limitrate
        end
    end

    out{rr} = lineage_count;
end

% =========================================================
% Functions
% =========================================================

function Bscr = Bscr_finder(Btemp, f_vec, alpha_A, del_B)

Btemp = Btemp(:);
f_vec = f_vec(:);

A = alpha_A .* f_vec;
S = sum(f_vec .* Btemp);

Bscr0 = ones(size(Btemp));

tol   = 1e-2;
maxit = 200;

for it = 1:maxit
    x = sum(f_vec .* Bscr0);

    denom_shift = del_B * (S - x);
    Bscr = (Btemp .* A) ./ (A + denom_shift);

    valid = Btemp > 0.5;
    eps_rel = mean(abs(Bscr(valid) - Bscr0(valid)) ./ Btemp(valid));

    if eps_rel <= tol
        if all(Bscr(Btemp > 0) <= Btemp(Btemp > 0) + 1e-2)
            break
        end
    end

    Bscr0 = Bscr;
end

end

function Bnew = c_death(Btemp_cell, dB_death)

N = size(Btemp_cell,1);

if N == 0
    Bnew = Btemp_cell;
    return
end

survive_mask = rand(N,1) >= dB_death;
Bnew = Btemp_cell(survive_mask,:);

end

function [Bnew, nmut_vec] = bmut( ...
    Btemp_cell, dB_birth, eps, s_av, nbar, nmut_vec, ppos, N_aff)

N = size(Btemp_cell,1);
n_i = cell2mat(Btemp_cell(:,1));

pb = dB_birth(n_i);
birth = rand(N,1) < pb;

Nout = N + sum(birth);
Bnew = cell(Nout,2);

Bnew(1:N,:) = Btemp_cell;

ib = find(birth);
nb = numel(ib);

tail = N + (1:nb);
Bnew(tail,1) = Btemp_cell(ib,1);
Bnew(tail,2) = Btemp_cell(ib,2);

for k = 1:nb
    i = ib(k);

    n      = Btemp_cell{i,1};
    mut_ID = Btemp_cell{i,2};

    psi  = 1 / (1 + (1 + exp(-eps*(n-nbar))) / s_av);
    pmut = 0.6 - 0.4*psi;

    pos = [i, tail(k)];

    for j = 1:2
        if rand < pmut
            [new_mut_ID, nmut_vec] = mutgen(mut_ID, nmut_vec);

            if rand > 0.2
                new_n = n;
            else
                if rand < ppos
                    new_n = n + 1;
                else
                    new_n = n - 1;
                end
            end

            new_n = min(max(new_n,1),N_aff);
            Bnew(pos(j),:) = {new_n, new_mut_ID};
        end
    end

end

end

function [new_mut_ID, nmut_vec] = mutgen(mut_ID, nmut_vec)

n = length(mut_ID);
nmut_vec(n+1) = nmut_vec(n+1) + 1;
new_mut_ID = [mut_ID, nmut_vec(n+1)];

end

function Bout = Bcounter(Bcell, N_aff)

if isempty(Bcell)
    Bout = zeros(N_aff,1);
    return
end

inds = cell2mat(Bcell(:,1));
inds = min(max(inds,1),N_aff);

Bout = accumarray(inds, 1, [N_aff,1]);

end

function H = simpson_effnum(counts)

p = counts / sum(counts);
H = 1 / sum(p.^2);

end
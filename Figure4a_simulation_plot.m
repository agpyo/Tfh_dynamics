% SIMULATE_FIGURE4A_TREES
% Run the stochastic GC model and reproduce the two Figure 4a lineage-tree
% conditions: low TFH sensitivity (s = 1) and high sensitivity (s = 5).
%
% Figure conditions from the manuscript caption:
%   * trees collected on day 7
%   * normalized antigen availability alpha_A/(1000*delta_B) = 1
%   * 100 initial B cells distributed among 50 unique lineages
%   * all B cells begin in the same affinity class
%   * 12 initial TFH cells at the indicated sensitivity
%
% The stochastic birth, death, mutation, antigen-screening, and TFH update
% rules are retained from simulation_code.m. Mutation histories already
% encode the true parent-child paths, so the plotted genealogy is built
% directly from those paths. This avoids the very slow all-pairs distance
% matrix and neighbor-joining calculation.
%
% Requires Statistics and Machine Learning Toolbox for BINORND.

clear; clc; close all

%% Figure 4a parameters
N_aff = 150;
N_T = 3;
Tend = 7*24;                  % day-7 tree
dt = 3e-2;

del_B0 = 1/8;
del_T0 = 1/18;
kappa = 5;
lambda = 0.077;

B0 = 100;
T0_total = 12;
N_lineages = 50;
iseed = 15;                  % common starting affinity; plotted as n = 0
eps_aff = 0.5;
ppos = 0.05;

normalized_antigen = 1;
alpha_A0 = normalized_antigen*1000*del_B0;
sensitivity = [1 5];

% Fixed seed makes the two displayed examples reproducible. The two
% conditions use consecutive, independent segments of this random stream.
rng(4,'twister');

%% Simulate the two conditions
final_B = cell(1,numel(sensitivity));
final_T = cell(1,numel(sensitivity));

for ic = 1:numel(sensitivity)
    fprintf('Simulating day-7 tree for s = %g ...\n',sensitivity(ic));
    s_vec = sensitivity(ic)*ones(N_T,1);
    [final_B{ic},final_T{ic}] = simulate_condition( ...
        N_aff,N_T,Tend,dt,del_B0,del_T0,kappa,lambda,B0,T0_total, ...
        iseed,eps_aff,ppos,N_lineages,alpha_A0,s_vec);
    fprintf('  %d B cells and %d TFH cells survived.\n', ...
        size(final_B{ic},1),sum(final_T{ic}));
end

%% Plot the two exact lineage trees
fig = figure('Color','w','Units','inches','Position',[1 1 10 4.2]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
blue = [0.10 0.46 0.82];
tree_info = cell(size(sensitivity));
tree_axes = gobjects(size(sensitivity));

for ic = 1:numel(sensitivity)
    ax = nexttile(tl,ic);
    tree_axes(ic) = ax;
    tree_info{ic} = plot_encoded_genealogy(ax,final_B{ic}(:,2),blue);
    title(ax,sprintf('s = %g',sensitivity(ic)), ...
        'FontName','Arial','FontSize',11,'FontWeight','bold');
end

% Apply one coordinate system to both panels. One horizontal unit is one
% genotype position and one vertical unit is one accumulated mutation in
% either tree, so neither condition is stretched independently.
commonGenotypes = max(cellfun(@(x) x.genotypeCount,tree_info));
commonMutations = max(cellfun(@(x) x.maxMutations,tree_info));
xPad = max(1,0.02*commonGenotypes);
yPad = max(1,0.08*commonMutations);
commonXLim = [1-xPad,commonGenotypes+xPad];
commonYLim = [-commonMutations-yPad,yPad];

genotypeScale = choose_genotype_scale(commonGenotypes);
mutationScale = min(5,max(1,commonMutations));
for ic = 1:numel(tree_axes)
    xlim(tree_axes(ic),commonXLim);
    ylim(tree_axes(ic),commonYLim);
    add_scale_bars(tree_axes(ic),commonGenotypes,commonMutations, ...
        genotypeScale,mutationScale,xPad,yPad);
end

sgtitle(tl,'Day 7 B-cell phylogenetic trees', ...
    'FontName','Arial','FontSize',12,'FontWeight','bold');

script_dir = fileparts(mfilename('fullpath'));
exportgraphics(fig,fullfile(script_dir,'Figure4a_two_trees.png'), ...
    'Resolution',300);
save(fullfile(script_dir,'Figure4a_two_trees_data.mat'), ...
    'final_B','final_T','tree_info','sensitivity','normalized_antigen');

%% Stochastic GC simulation from simulation_code.m
function [Bcells,Ttemp] = simulate_condition(N_aff,N_T,Tend,dt, ...
    del_B0,del_T0,kappa,lambda,B0,T0_total,iseed,eps_aff,ppos, ...
    N_lineages,alpha_A0,s_vec)

M = round(Tend/dt)+1;
nvec = (1:N_aff)';

Bcells = cell(B0,2);
for j = 1:B0
    Bcells(j,:) = {iseed,1+mod(j,N_lineages)};
end
Btemp = Bcounter(Bcells,N_aff);

% Twelve total TFH cells, evenly distributed across the three clones.
Ttemp = floor(T0_total/N_T)*ones(N_T,1);
Ttemp(1:mod(T0_total,N_T)) = Ttemp(1:mod(T0_total,N_T))+1;

nmut_vec = zeros(500,1);
nmut_vec(1) = 1;

for it = 1:(M-1)
    if isempty(Bcells) || sum(Ttemp)==0
        warning('simulate_figure4a_trees:Extinction', ...
            'A population became extinct before day 7.');
        break
    end

    nbar = sum(Btemp.*nvec)/sum(Btemp);
    s_av = sum(s_vec.*Ttemp)/sum(Ttemp);

    f_vec = 1./(1+exp(-eps_aff*(nvec-nbar)));
    psi_vec = 1./(1+1./(f_vec*s_vec'));

    Bscr = Bscr_finder(Btemp,f_vec,alpha_A0,del_B0);
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
end
end

%% Fast exact tree construction and plotting
function info = plot_encoded_genealogy(ax,histories,color)
% One terminal leaf is drawn per unique mutation-defined genotype. The
% first history element is an initial-lineage label; subsequent elements
% are mutation IDs. A synthetic unmutated common ancestor joins the 50
% initial lineages at mutation depth zero.

if isempty(histories)
    axis(ax,'off');
    text(ax,0.5,0.5,'Population extinct','HorizontalAlignment','center');
    info = struct('genotypeCount',0,'nodeCount',0,'maxMutations',0);
    return
end

historyKeys = cellfun(@history_key,histories,'UniformOutput',false);
[~,uniqueIndex] = unique(historyKeys,'stable');
histories = histories(uniqueIndex);
nGenotypes = numel(histories);

% The maximum possible number of nodes is one plus the total number of
% entries across all histories. Preallocation avoids repeated growth.
maxNodes = 1+sum(cellfun(@numel,histories));
parent = zeros(maxNodes,1);
depth = zeros(maxNodes,1);
nodeKey = cell(maxNodes,1);
nodeKey{1} = 'UCA';
nNodes = 1;
terminalNode = zeros(nGenotypes,1);

keyToNode = containers.Map('KeyType','char','ValueType','double');
keyToNode('UCA') = 1;

for g = 1:nGenotypes
    h = histories{g};
    parentID = 1;
    for k = 1:numel(h)
        key = history_key(h(1:k));
        if isKey(keyToNode,key)
            nodeID = keyToNode(key);
        else
            nNodes = nNodes+1;
            nodeID = nNodes;
            keyToNode(key) = nodeID;
            nodeKey{nodeID} = key;
            parent(nodeID) = parentID;
            depth(nodeID) = k-1;  % initial lineage is depth zero
        end
        parentID = nodeID;
    end
    terminalNode(g) = parentID;
end

parent = parent(1:nNodes);
depth = depth(1:nNodes);
nodeKey = nodeKey(1:nNodes);

% Lexicographic terminal ordering keeps descendants of a common prefix
% adjacent. Accumulate descendant positions upward to place internal nodes.
terminalKeys = nodeKey(terminalNode);
[~,order] = sort(terminalKeys);
terminalY = zeros(nGenotypes,1);
terminalY(order) = 1:nGenotypes;

sumY = zeros(nNodes,1);
countY = zeros(nNodes,1);
for g = 1:nGenotypes
    nodeID = terminalNode(g);
    while nodeID~=0
        sumY(nodeID) = sumY(nodeID)+terminalY(g);
        countY(nodeID) = countY(nodeID)+1;
        nodeID = parent(nodeID);
    end
end
y = sumY./max(countY,1);

% Plot all rectangular branches in a single graphics object. Horizontal
% position represents genotype order; vertical depth is mutation count.
nEdges = nNodes-1;
X = nan(4*nEdges,1);
Y = nan(4*nEdges,1);
q = 0;
for child = 2:nNodes
    p = parent(child);
    idx = q+(1:4);
    X(idx) = [y(p);y(child);y(child);nan];
    Y(idx) = [-depth(p);-depth(p);-depth(child);nan];
    q = q+4;
end

plot(ax,X,Y,'Color',color,'LineWidth',0.45);
hold(ax,'on');

maxDepth = max(depth);
xPad = max(1,0.02*nGenotypes);
yPad = max(1,0.08*maxDepth);
xlim(ax,[1-xPad,nGenotypes+xPad]);
ylim(ax,[-maxDepth-yPad,yPad]);

% UCA marker and label.
plot(ax,y(1),0,'o','MarkerSize',4,'MarkerFaceColor','w', ...
    'MarkerEdgeColor',color,'LineWidth',0.8);
text(ax,y(1),yPad*0.35,'UCA','HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontName','Arial','FontSize',9);

axis(ax,'off');
box(ax,'off');

info = struct('genotypeCount',nGenotypes,'nodeCount',nNodes, ...
    'maxMutations',maxDepth);
end

function add_scale_bars(ax,nGenotypes,maxDepth,genotypeScale, ...
    mutationScale,xPad,yPad)
% Identical data-unit scale bars are added to both matched axes.
x0 = 1;
y0 = -maxDepth-yPad*0.25;
plot(ax,[x0 x0],[y0 y0+mutationScale],'k-','LineWidth',1.2);
text(ax,x0+xPad*0.25,y0+mutationScale/2, ...
    sprintf('%d mutations',mutationScale),'FontName','Arial', ...
    'FontSize',8,'VerticalAlignment','middle');

xBar = max(1,nGenotypes-genotypeScale);
yBar = -maxDepth-yPad*0.62;
plot(ax,[xBar xBar+genotypeScale],[yBar yBar],'k-','LineWidth',1.2);
text(ax,xBar+genotypeScale/2,yBar-yPad*0.12, ...
    sprintf('%d genotypes',genotypeScale),'FontName','Arial', ...
    'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','top');
end

function scale = choose_genotype_scale(n)
if n>=1000
    scale = 1000;
elseif n>=500
    scale = 500;
elseif n>=200
    scale = 200;
elseif n>=100
    scale = 100;
elseif n>=50
    scale = 50;
else
    scale = max(1,10*floor(n/20));
end
end

function key = history_key(history)
key = sprintf('%.17g,',history);
end

%% Original stochastic helper functions
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
    if eps_rel<=tol
        if all(Bscr(Btemp>0)<=Btemp(Btemp>0)+1e-2)
            break
        end
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

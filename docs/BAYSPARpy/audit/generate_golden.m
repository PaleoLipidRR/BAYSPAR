function generate_golden()
% GENERATE_GOLDEN  Produce the MATLAB reference values the Python port is tested against.
%
% Run from the BAYSPAR repository root. The '>>' below is the MATLAB prompt --
% these are typed inside MATLAB, not in a shell.
%
%     >> cd /path/to/BAYSPAR
%     >> addpath('docs/BAYSPARpy/audit')
%     >> generate_golden
%
% Or, entirely from a shell, in the repository root:
%
%     matlab -batch "addpath('docs/BAYSPARpy/audit'); generate_golden"
%
% GNU Octave works too, and needs no MATLAB licence -- this script was developed
% against Octave 8.4 with the statistics package (apt install octave
% octave-statistics), which is where prctile comes from:
%
%     octave --no-gui --quiet --eval "pkg load statistics; addpath('docs/BAYSPARpy/audit'); generate_golden"
%
% MATLAB needs the Statistics Toolbox for the same function.
%
% Writes docs/BAYSPARpy/audit/golden_matlab.mat and prints the same values in
% the layout of reference_trace.txt, so the two can be diffed by eye.
%
% Everything computed here is DETERMINISTIC -- it all happens upstream of the
% first randn -- so MATLAB and Python must agree to the digit. The stochastic
% comparison is a separate job (see HANDOFF.md).

if ~exist('ModelOutput', 'dir')
    error('run this from the BAYSPAR repository root (no ModelOutput/ directory here)');
end

G = struct();
G.generated = datestr(now, 'yyyy-mm-dd HH:MM:SS');
G.matlab_version = version();

%% A. Thinning index -------------------------------------------------------
fprintf('\n== A. Thinning index ==\n');
nsamps_list = [10 999 1000 1001 20000];
G.thinning = struct();
for k = 1:numel(nsamps_list)
    n = nsamps_list(k);
    ind = round(linspace(1, 20000, n));
    G.thinning.(sprintf('n%d', n)) = ind;
    fprintf('  Nsamps=%5d: first 5 [%s]  last [%d]  unique=%d\n', ...
        n, num2str(ind(1:min(5,end)), '%d '), ind(end), numel(unique(ind)));
end

%% B. Chordal distances ----------------------------------------------------
fprintf('\n== B. Chordal distance (km) ==\n');
pairs = [   0.0000    0.0000    0.0000    1.0000
            0.0000    0.0000    1.0000    0.0000
          -17.6635    9.1660  -17.5000    9.5000
          100.0000  -45.0000 -100.0000   45.0000
          180.0000    0.0000 -180.0000    0.0000];
G.distances = zeros(size(pairs, 1), 1);
for k = 1:size(pairs, 1)
    G.distances(k) = EarthChordDistances_2(pairs(k, 1:2), pairs(k, 3:4));
    fprintf('  (%9.4f, %8.4f) -> (%9.4f, %8.4f)  =  %.9f\n', pairs(k, :), G.distances(k));
end
G.distance_pairs = pairs;

%% C. Prior mean and grid cell, per demo series ---------------------------
fprintf('\n== C. Prior mean and grid cell ==\n');
load('ModelOutput/tex_testdata.mat', 'castaneda2010', 'lopes_santos2010', 'shevenell2011');
series = {castaneda2010, lopes_santos2010, shevenell2011};
names  = {'castaneda2010', 'lopes_santos2010', 'shevenell2011'};
runs   = {'SST', 'subT'};
grid_half_space = 10; min_num = 1; max_dist = 500;
G.prior = struct();
for r = 1:2
    runname = runs{r};
    if strcmp(runname, 'SST')
        load('ModelOutput/obsSST', 'locs_st_obs', 'st_obs_ave_vec');
    else
        load('ModelOutput/obssubT', 'locs_st_obs', 'st_obs_ave_vec');
    end
    load(['ModelOutput/Output_SpatAg_', runname, '/params_standard'], 'Locs_Comp');
    fprintf('  --- runname = ''%s'' ---\n', runname);
    for s = 1:3
        st = series{s};
        lon = st.lon; lat = st.lat;
        dists = EarthChordDistances_2([lon, lat], locs_st_obs);
        [vals, inds] = sort(dists);
        num_below = find(vals < max_dist, 1, 'last');
        if num_below > min_num
            pm = mean(st_obs_ave_vec(inds(1:num_below)));
        else
            pm = mean(st_obs_ave_vec(inds(1:min_num)));
            num_below = 0;
        end
        cell_ix = find(abs(Locs_Comp(:,1)-lon) <= grid_half_space & ...
                       abs(Locs_Comp(:,2)-lat) <= grid_half_space);
        key = sprintf('%s_%s', runname, names{s});
        G.prior.(key) = struct('lon', lon, 'lat', lat, 'n', numel(st.tex86), ...
                               'prior_mean', pm, 'n_below', num_below, ...
                               'grid_loc', Locs_Comp(cell_ix, :));
        fprintf('  %-18s N=%4d  site=(%9.4f, %8.4f)  PriorMean=%.10f  from %d obs  GridLoc=[%s]\n', ...
            names{s}, numel(st.tex86), lon, lat, pm, num_below, ...
            num2str(Locs_Comp(cell_ix, :), '%d '));
    end
end

%% D. Analogue selection, Wilson Lake -------------------------------------
fprintf('\n== D. Analogue selection (wilsonlake, SST) ==\n');
load('ModelOutput/wilsonlake', 'wilsonlake');
load('ModelOutput/Data_Input_SpatAg_SST', 'Data_Input');
dats = wilsonlake.tex86(:);
search_tol = std(dats) * 2;
N_bg = length(Data_Input.Locs(:,1));
spatialMean = NaN(N_bg, 1);
for i = 1:N_bg
    spatialMean(i) = mean(Data_Input.Obs_Stack(Data_Input.Inds_Stack == i));
end
inder_g = spatialMean >= (mean(dats) - search_tol) & spatialMean <= (mean(dats) + search_tol);
G.analog = struct('n', numel(dats), 'mean_tex', mean(dats), 'std_tex', std(dats), ...
                  'search_tol', search_tol, 'cell_means', spatialMean, ...
                  'selected', find(inder_g), 'locs', Data_Input.Locs(inder_g, :));
fprintf('  N=%d  mean=%.10f  std=%.10f  search_tol=%.10f\n', ...
    numel(dats), mean(dats), std(dats), search_tol);
fprintf('  cell means span %.6f to %.6f; %d analogues selected\n', ...
    min(spatialMean), max(spatialMean), sum(inder_g));

%% E. Per-draw posterior, before noise ------------------------------------
fprintf('\n== E. post_mean / post_sig, lopes_santos2010, subT, prior_std=6 ==\n');
load('ModelOutput/Output_SpatAg_subT/params_standard', ...
     'alpha_samples_comp', 'beta_samples_comp', 'tau2_samples', 'Locs_Comp');
dats = lopes_santos2010.tex86(:);
lon = lopes_santos2010.lon; lat = lopes_santos2010.lat;
prior_mean = G.prior.subT_lopes_santos2010.prior_mean;
prior_std = 6;
ind_s = round(linspace(1, length(tau2_samples), 1000));
cell_ix = find(abs(Locs_Comp(:,1)-lon) <= grid_half_space & ...
               abs(Locs_Comp(:,2)-lat) <= grid_half_space);
a = alpha_samples_comp(cell_ix, ind_s);
b = beta_samples_comp(cell_ix, ind_s);
t2 = tau2_samples(ind_s);
pinv_cov = prior_std^-2;
den = pinv_cov + b.^2 ./ t2;
num = pinv_cov * prior_mean + (b ./ t2) .* (dats(1) - a);
G.posterior = struct('dats1', dats(1), 'prior_mean', prior_mean, 'prior_std', prior_std, ...
                     'alpha', a(1:5), 'beta', b(1:5), 'tau2', t2(1:5), ...
                     'post_mean', num(1:5) ./ den(1:5), 'post_sig', sqrt(1 ./ den(1:5)));
fprintf('  dats(1)=%.10f  PriorMean=%.10f\n', dats(1), prior_mean);
for j = 1:5
    fprintf('  draw %5d  alpha %.10f  beta %.10f  tau2 %.10f  post_mean %.10f  post_sig %.10f\n', ...
        ind_s(j), a(j), b(j), t2(j), num(j)/den(j), sqrt(1/den(j)));
end

%% F. Percentile convention ------------------------------------------------
fprintf('\n== F. prctile convention ==\n');
G.prctile_1to10 = prctile(1:10, [5 50 95]);
fprintf('  prctile(1:10, [5 50 95]) = [%s]   %% Python: method="hazen"\n', ...
    num2str(G.prctile_1to10, '%.6f '));

%% Save --------------------------------------------------------------------
out = fullfile('docs', 'BAYSPARpy', 'audit', 'golden_matlab.mat');
% -v7 matters: Octave's default save format is its own ASCII text, which
% scipy.io.loadmat cannot read (it reports "Unknown mat file type"). Octave wants
% the option before the filename, MATLAB after it.
if exist('OCTAVE_VERSION', 'builtin')
    save('-v7', out, '-struct', 'G');
else
    save(out, '-struct', 'G', '-v7');
end
fprintf('\nwrote %s\n', out);
fprintf('Now run, from the repository root:\n');
fprintf('    cd BAYSPARpy && python -m pytest tests/test_golden.py -v\n\n');
end

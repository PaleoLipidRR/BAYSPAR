function verify_original()
% VERIFY_ORIGINAL  Run the ORIGINAL MATLAB functions and save their output.
%
% generate_golden.m recomputes the intermediates inline, which checks the
% arithmetic but not the reading of the functions themselves. This script calls
% bayspar_tex.m, bayspar_tex_analog.m and TEX_forward.m unmodified, so the Python
% is compared against the reference actually executing.
%
% From the repository root:
%
%     matlab -batch "addpath('docs/BAYSPARpy/audit'); verify_original"
%
% or, under Octave 8+ with the statistics package:
%
%     octave --no-gui --quiet --eval "pkg load statistics; addpath('docs/BAYSPARpy/audit'); verify_original"
%
% OCTAVE LIMITATION, in the reference and not in the port: bayspar_tex.m line 66
% is `if runname=="SST"`. MATLAB reads "SST" as a string and compares as a whole;
% Octave reads it as a char array and compares elementwise, so a 4-character
% runname ('subT') raises "nonconformant arguments". This script therefore runs
% the standard-mode case as 'SST' only under Octave. Under MATLAB it does both --
% pass 'both' as the argument to force it.
%
% Writes docs/BAYSPARpy/audit/original_matlab.mat.

if ~exist('ModelOutput', 'dir')
    error('run this from the BAYSPAR repository root');
end
is_octave = exist('OCTAVE_VERSION', 'builtin') > 0;

V = struct();
V.is_octave = is_octave;
if is_octave
    V.engine = ['GNU Octave ', OCTAVE_VERSION];
else
    V.engine = ['MATLAB ', version()];
end
V.generated = datestr(now, 'yyyy-mm-dd HH:MM:SS');
V.n_draws = 1000;

%% Standard mode, through bayspar_tex.m itself
load('ModelOutput/tex_testdata.mat', 'lopes_santos2010');
s = lopes_santos2010;
V.standard_input = struct('tex86', s.tex86(:), 'lon', s.lon, 'lat', s.lat, ...
                          'prior_std', 6, 'runname', 'SST');
O = bayspar_tex(s.tex86, s.lon, s.lat, 6, 'SST', V.n_draws, 1);
V.standard = struct('PriorMean', O.PriorMean, 'GridLoc', O.GridLoc, ...
                    'SiteLoc', O.SiteLoc, 'Preds', O.Preds, 'PredsEns', O.PredsEns);
fprintf('bayspar_tex   : PriorMean %.10f  GridLoc [%s]  Preds %dx%d\n', ...
    O.PriorMean, num2str(O.GridLoc, '%d '), size(O.Preds));

%% Analogue mode, through bayspar_tex_analog.m itself
load('ModelOutput/wilsonlake', 'wilsonlake');
dats = wilsonlake.tex86(:);
tol = std(dats) * 2;
V.analog_input = struct('tex86', dats, 'prior_mean', 30, 'prior_std', 20, ...
                        'search_tol', tol, 'runname', 'SST');
A = bayspar_tex_analog(dats, 30, 20, tol, 'SST', V.n_draws, 1);
V.analog = struct('AnLocs', A.AnLocs, 'Preds', A.Preds, ...
                  'PredsEns_size', size(A.PredsEns), ...
                  'PriorMean', A.PriorMean, 'PriorStd', A.PriorStd);
fprintf('bayspar_tex_analog: %d analogues  Preds %dx%d  PredsEns %dx%d\n', ...
    size(A.AnLocs, 1), size(A.Preds), size(A.PredsEns));

%% Forward model, through TEX_forward.m itself
t_in = [22; 25; 28];
V.forward_input = struct('lat', s.lat, 'lon', s.lon, 't', t_in, 'runname', 'SST');
F = TEX_forward(repmat(s.lat, 3, 1), repmat(s.lon, 3, 1), t_in, 'SST');
V.forward = F;
fprintf('TEX_forward   : %dx%d, mean per row [%s]\n', size(F), num2str(mean(F, 2)', '%.4f '));

out = fullfile('docs', 'BAYSPARpy', 'audit', 'original_matlab.mat');
if is_octave
    save('-v7', out, '-struct', 'V');
else
    save(out, '-struct', 'V', '-v7');
end
fprintf('\nwrote %s\n', out);
end

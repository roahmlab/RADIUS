% clear
addpath(genpath('/data')); 
addpath(genpath('/simulator'));
%% User Inputs
input_dir = '/data';
output_dir = '/data';

is_left_turn = true; % set it to false to preprocess FRS of highway scenarios

if ~is_left_turn
    input_fname = 'CUDA_Highway_frs.mat';
    cpp_processed_output_fname = 'cpp_processed_CUDA_Highway_frs.frs.bin';
    input_fpath = string(fullfile(input_dir, input_fname));
    cpp_processed_fpath = string(fullfile(output_dir, cpp_processed_output_fname));
    process_file(input_fpath, cpp_processed_fpath, is_left_turn);
else
    is_left_turn = true;
    input_fname = 'CUDA_LeftTurn_frs.mat';
    cpp_processed_output_fname = "cpp_processed_CUDA_LeftTurn_frs.frs.bin";
    input_fpath = string(fullfile(input_dir, input_fname))
    cpp_processed_fpath = string(fullfile(output_dir, cpp_processed_output_fname))
    process_file(input_fpath, cpp_processed_fpath, is_left_turn)
end

function process_file(input_fpath, cpp_processed_fpath, is_left_turn)
%% Load Full FRS If it is not loaded
% clear frs_full loaded_fpath
if (~exist('loaded_fpath', 'var')) || (~exist('frs_full', 'var')) || any(loaded_fpath ~= input_fpath)
    if ~is_left_turn
        frs_full = load(input_fpath);
        loaded_fpath = input_fpath;
    else
        frs_in = load(input_fpath);
        frs_in.LeftTurnFRS.manu_type = 3;
        frs_in.LeftTurnFRS.cuda_FRS.manu_type = 3;
        loaded_fpath = input_fpath;
        lan_cells = cell(1, 1);
        lan_cells{1,1} = frs_in.LeftTurnFRS;
        empty_map =  containers.Map(["lan", "Au", "dir", "lantb", "autb", "dirtb"], {{}, {}, {}, {}, {}, {}});
        mega_1_map =  containers.Map(["lan", "Au", "dir", "lantb", "autb", "dirtb"], {{}, {}, {}, {}, {}, {}});
        mega_1_map('lan') = lan_cells;
        mega_cells = {empty_map, mega_1_map};
        frs_full = struct('M_mega', []);
        frs_full.M_mega = mega_cells;
        
    end
end

frs_empty_removed = frs_full;
while true
    key_all_empty = @(mega, key) (isKey(mega, key) && all(cellfun(@isempty, mega(key)), 'all'));
    frses_empty = @(mega) key_all_empty(mega, 'Au') && key_all_empty(mega, 'dir') && key_all_empty(mega, 'lan');
    mega0 = frs_empty_removed.M_mega{1};
    if frses_empty(mega0)
        frs_empty_removed.M_mega(1) = [];
    else
        break;
    end
end

vehrs_keys = {'Au', 'dir', 'lan'};
to_zono_mat = @(zono) zono.Z;

for i = 1:length(frs_empty_removed.M_mega)
    frs_empty_removed.M_mega{i} = map_to_struct(frs_empty_removed.M_mega{i});

    for vehrs_key = vehrs_keys
        vehrs_key = vehrs_key{1};
        if ~isfield(frs_empty_removed.M_mega{i}, vehrs_key)
            continue;
        end
        vehrs_cells = frs_empty_removed.M_mega{i}.(vehrs_key);
        if isempty(vehrs_cells)
            continue;
        end
        i
        vehrs_cells = cellfun(@(x) check_min_time(x, true), vehrs_cells, 'UniformOutput', false)
        cellfun(@(x) check_min_time(x, false), vehrs_cells, 'UniformOutput', false)
        frs_empty_removed.M_mega{i}.(vehrs_key) = cellfun(@(x) fix_manu(x), vehrs_cells, 'UniformOutput', false);
    end
end

tic
MATLAB_PROCESS_FRS(frs_empty_removed, cpp_processed_fpath)
toc
end
%% Helper Functions

function manu_frses = fix_manu(manu_frses)
if nargin < 1 || isempty(manu_frses)
    manu_frses = [];
    return;
end
    fix_fcn = @(zono) zono.Z;
    manu_frses.vehRS_save = cellfun(fix_fcn, manu_frses.vehRS_save, 'UniformOutput', false);
end

function ret = map_to_struct(map)
    ret = struct();
    for key = map.keys()
        key = key{1};
        ret.(key) = map(key);
    end
end

function single_frs = check_min_time(single_frs, fix_errors)
    if isempty(single_frs)
        return;
    end
    vehrs = single_frs.vehRS_save;
    for i = 1:length(vehrs)
        if fix_errors
            [vals, idxs] = find(vehrs{i}.Z(end,2:end) ~= 0);
            % TODO disable this
            for j = 2:length(idxs)
                vehrs_i = vehrs{i};
                Zmat = vehrs_i.Z;
                Zmat(end,1+idxs(j)) = 0;
                vehrs_i = zonotope(Zmat);
                vehrs{i} = vehrs_i;
            end
            [vals, idxs] = find(vehrs{i}.Z(end,2:end) ~= 0);
            assert(length(idxs) == 1)
            t_gen = (double(int64((vehrs{i}.Z(end, idxs(1)+1) * 2) * 100)) / 100) / 2;

            if t_gen > 0.08
                t_gen = 0.08;
            end

            vehrs_i_new = vehrs{i}.Z;
            vehrs_i_new(end, idxs(1)+1) = t_gen;
            vehrs{i} = zonotope(vehrs_i_new);
            single_frs.vehRS_save{i} = vehrs{i};
        end

        % TODO SHOULD JUST KEEP THIS
        zono_mat = vehrs{i}.Z;
        num_time_gens = sum(zono_mat(end,2:end) > 0, 'all');
        assert(num_time_gens == 1, "num time gens != 1");
        time_sum_sec = 2*sum(abs(zono_mat(end, 2:end)));
        assert(time_sum_sec - 0.01 > -1.0e-6);
        assert(time_sum_sec - 0.16 <  1.0e-6);
        % TODO should check center value for divisibility by dt 
        % between mu sigmas, but currently it should be rounded
        assert(mod(int64(time_sum_sec * 1000), 10) == 0);
    end
end

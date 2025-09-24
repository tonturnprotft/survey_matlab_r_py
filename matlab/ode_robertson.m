% ode_robertson.m — Solve Robertson ODE with ode15s and write one CSV row.
% Usage (from MATLAB -batch):
%   addpath('/Users/tree/Simu/matlab'); ode_robertson(1e5, 1e-6, 1e-8, 'ode15s', 'results_ode_robertson_matlab.csv', 'ref_robertson_t1e5.csv')
function ode_robertson(t_end, rtol, atol, solver, outcsv, refcsv)
    if nargin < 1, t_end = 1e5; end
    if nargin < 2, rtol = 1e-6; end
    if nargin < 3, atol = 1e-8; end
    if nargin < 4 || isempty(solver), solver = 'ode15s'; end
    if nargin < 5 || isempty(outcsv), outcsv = 'results_ode_robertson_matlab.csv'; end
    if nargin < 6, refcsv = ''; end

    y0 = [1;0;0];

    % sample times: 0 plus logspace(1e-6 to t_end)
    n = 500;
    if t_end <= 1e-6
        t_eval = linspace(0, t_end, n);
    else
        t_eval = [0, logspace(-6, log10(t_end), n-1)];
    end

    opts = odeset('RelTol', rtol, 'AbsTol', atol);
    f = @(t,y)[-0.04*y(1) + 1.0e4*y(2)*y(3); 0.04*y(1) - 1.0e4*y(2)*y(3) - 3.0e7*(y(2)^2); 3.0e7*(y(2)^2)];
    cpu0 = cputime;
    tic;
    switch lower(solver)
        case 'ode15s'
            sol = ode15s(f, [0 t_end], y0, opts);
        case 'ode23s'
            sol = ode23s(f, [0 t_end], y0, opts);
        case 'ode23t'
            sol = ode23t(f, [0 t_end], y0, opts);
        otherwise
            sol = ode15s(f, [0 t_end], y0, opts);
    end
    
    wall = toc;
    cpu = cputime - cpu0;

    % Interpolate onto t_eval using solver's continuous output
    y_interp = deval(sol, t_eval).';   % N x 3

    rmse = ""; linf = "";
    if ~isempty(refcsv) && isfile(refcsv)
        ref = readmatrix(refcsv);
        % Expect columns: t,y1,y2,y3 with header row; handle header
        if any(isnan(ref(1,:)))
            ref = ref(2:end,:);
        end
        tref = ref(:,1);
        y1r = interp1(tref, ref(:,2), t_eval, 'linear', 'extrap')';
        y2r = interp1(tref, ref(:,3), t_eval, 'linear', 'extrap')';
        y3r = interp1(tref, ref(:,4), t_eval, 'linear', 'extrap')';
        diff = [y_interp(:,1)-y1r, y_interp(:,2)-y2r, y_interp(:,3)-y3r];
        rmse = sqrt(mean(sum(diff.^2,2)));
        linf = max(max(abs(diff)));
    end

    row = struct();
    row.family = "ODE";
    row.benchmark_id = "ode_robertson";
    row.platform = "matlab";
    row.tool = "ode15s";
    row.version = version;
    row.solver = string(solver);
    row.step_mode = "adaptive";
    row.abs_tol = atol;
    row.rel_tol = rtol;
    row.dt_init = "";
    row.dt_max = "";
    row.seed = "";
    row.params_json = jsonencode(struct("t_end",t_end,"rtol",rtol,"atol",atol,"solver",solver));
    row.scale_param = "t_end";
    row.scale_value = t_end;
    row.repeat_idx = 0;
    row.run_id = "matlab_" + string(floor(posixtime(datetime('now'))));
    row.start_ts = "";
    row.end_ts = "";
    row.wall_time_s = wall;
    row.cpu_time_s = cpu;
    row.peak_mem_mb = "";
    row.n_steps = length(t_eval)-1;
    row.n_events = "";
    row.rmse = rmse;
    row.linf = linf;
    row.event_time_err_mean = "";
    row.event_time_err_max = "";
    row.throughput = "";
    row.L_timeavg = "";
    row.status = "OK";
    row.error_message = "";
    row.host = string(system_dependent('getos'));

    field_order = ["family","benchmark_id","platform","tool","version","solver","step_mode", ...
        "abs_tol","rel_tol","dt_init","dt_max","seed","params_json","scale_param","scale_value", ...
        "repeat_idx","run_id","start_ts","end_ts","wall_time_s","cpu_time_s","peak_mem_mb", ...
        "n_steps","n_events","rmse","linf","event_time_err_mean","event_time_err_max", ...
        "throughput","L_timeavg","status","error_message","host"];

    T = struct2table(orderfields(row, field_order));
    if isfile(outcsv)
        writetable(T, outcsv, 'WriteMode','append');
    else
        writetable(T, outcsv);
    end
    fprintf("status=OK n_steps=%d wall=%.3fs rmse=%s linf=%s\n", length(t_eval)-1, wall, string(rmse), string(linf));
end
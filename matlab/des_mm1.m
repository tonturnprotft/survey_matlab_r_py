% des_mm1.m — Plain MATLAB M/M/1 DES producing one CSV row compatible with results_schema.csv.
% Usage: In MATLAB command window:
%   des_mm1(0.9, 1.0, 100000, 42, 'results_des_mm1_matlab.csv')
function des_mm1(lambda, mu, horizon, seed, outcsv)
    if nargin < 1, lambda = 0.9; end
    if nargin < 2, mu = 1.0; end
    if nargin < 3, horizon = 100000; end
    if nargin < 4, seed = 42; end
    if nargin < 5, outcsv = 'results_des_mm1_matlab.csv'; end

    rng(seed,'twister');

    t0cpu = cputime;
    t0 = tic;

    % Event times
    t = 0.0;
    t_next_arr = -log(1 - rand) / lambda;
    t_next_dep = inf;

    n_system = 0;
    arrivals = 0;
    departures = 0;
    events_processed = 0;

    % time-average N(t)
    last_t = 0.0;
    area = 0.0;

    while true
        t_next = min(t_next_arr, t_next_dep);
        if t_next > horizon
            % close integral to horizon
            area = area + n_system * (horizon - last_t);
            t = horizon;
            break;
        end
        t = t_next;
        area = area + n_system * (t - last_t);
        last_t = t;

        if t_next_arr <= t_next_dep
            % arrival
            arrivals = arrivals + 1;
            events_processed = events_processed + 1;
            if n_system == 0
                t_next_dep = t + (-log(1 - rand) / mu);
            end
            n_system = n_system + 1;
            t_next_arr = t + (-log(1 - rand) / lambda);
        else
            % departure
            n_system = n_system - 1;
            departures = departures + 1;
            events_processed = events_processed + 1;
            if n_system > 0
                t_next_dep = t + (-log(1 - rand) / mu);
            else
                t_next_dep = inf;
            end
        end
    end

    wall = toc(t0);
    cpu = cputime - t0cpu;

    throughput = departures / horizon;
    L_timeavg = area / horizon;
    rho = lambda / mu;

    % Build row struct following results_schema.csv
    row = struct();
    row.family = "DES";
    row.benchmark_id = "des_mm1";
    row.platform = "matlab";
    row.tool = "plain_des";
    row.version = version;
    row.solver = "event_loop";
    row.step_mode = "event";
    row.abs_tol = "";
    row.rel_tol = "";
    row.dt_init = "";
    row.dt_max = "";
    row.seed = seed;
    row.params_json = jsonencode(struct("lambda",lambda,"mu",mu,"horizon",horizon));
    row.scale_param = "rho";
    row.scale_value = rho;
    row.repeat_idx = 0;
    row.run_id = "matlab_" + string(floor(posixtime(datetime('now'))));
    row.start_ts = "";
    row.end_ts = "";
    row.wall_time_s = wall;
    row.cpu_time_s = cpu;
    row.peak_mem_mb = "";
    row.n_steps = events_processed;
    row.n_events = departures;
    row.rmse = "";
    row.linf = "";
    row.event_time_err_mean = "";
    row.event_time_err_max = "";
    row.throughput = throughput;
    row.L_timeavg = L_timeavg;
    row.status = "OK";
    row.error_message = "";
    row.host = string(getenv('COMPUTERNAME'));
    if row.host == "", row.host = system_dependent('getos'); end

    % Convert to table and append to CSV
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

    fprintf("rho=%.3f departures=%d L_timeavg≈%.3f thr≈%.4f wall=%.3fs\n", rho, departures, L_timeavg, throughput, wall);
end
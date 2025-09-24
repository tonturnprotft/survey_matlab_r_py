function hybrid_bounce(t_end, rtol, atol, e, h0, v0, outpath, refcsv)
% Usage: hybrid_bounce(10,1e-6,1e-8,0.9,10,0,'results_hybrid_bounce_matlab.csv','ref_bounce.csv')
if nargin<8, refcsv=''; end
opts = odeset('RelTol',rtol,'AbsTol',atol,'Events',@eventsfun);
G=9.81; Y0=[h0; v0]; T=0; Y=Y0.';
t0 = tic;
while T(end) < t_end
    [t,y,te,ye] = ode45(@(t,y)[y(2); -G],[T(end) t_end],Y0,opts);
    T=[T; t(2:end)]; Y=[Y; y(2:end,:)]; %#ok<AGROW>
    if isempty(te), break; end
    Y0 = [0; -e*ye(end,2)];
end
grid = linspace(0,t_end,1001);
h = interp1(T, Y(:,1), grid, 'linear','extrap');
v = interp1(T, Y(:,2), grid, 'linear','extrap');
rmse=NaN; linf=NaN;
if ~isempty(refcsv) && isfile(refcsv)
    ref = readtable(refcsv);
    [rmse, linf] = piecewise_aligned_rmse(T, Y, ref.t, ref.h, ref.v, t_end);
end
wall = toc(t0);
row.platform="matlab"; row.family="HYBRID"; row.benchmark_id="hybrid_bounce";
row.tool="ode45"; row.solver="events"; row.status="OK";
row.wall_time_s=wall; row.rmse=rmse; row.linf=linf; row.n_steps=size(Y,1);
row.params_json=jsonencode(struct('t_end',t_end,'rtol',rtol,'atol',atol,'e',e,'h0',h0,'v0',v0));
Trow=struct2table(row);
if isfile(outpath), writetable(Trow,outpath,'WriteMode','Append','QuoteStrings',true); else, writetable(Trow,outpath,'QuoteStrings',true); end
fprintf('status=OK steps=%d wall=%.3fs rmse=%g linf=%g\n', size(Y,1), wall, rmse, linf);
    function [value,isterminal,direction]=eventsfun(t,y), value=y(1); isterminal=1; direction=-1; end
end

function ev = detect_events(t, h)
% Detect impact times from dense reference by height zero-crossings
    ev = [];
    n = numel(t);
    if n < 2, return; end
    for i = 1:(n-1)
        h0 = h(i); h1 = h(i+1);
        if ( (h0 > 0 && h1 <= 0) || (h0 >= 0 && h1 < 0) )
            alpha = h0 / (h0 - h1);
            tz = t(i) + alpha * (t(i+1) - t(i));
            ev(end+1) = tz; %#ok<AGROW>
        end
    end
    if ~isempty(ev)
        ev = unique(round(ev, 12));
    end
end

function [rmse, linf] = piecewise_aligned_rmse(T, Y, t_ref, h_ref, v_ref, t_end)
% Compare candidate (T,Y) vs reference between impacts using NRMSE.
% Normalization: RMSE(h)/range(h_ref) and RMSE(v)/max(|v_ref|), averaged.
    ev = detect_events(t_ref, h_ref);
    ev = ev(ev > 0 & ev < t_end);
    bounds = [0.0, ev(:).', t_end];
    sse_h = 0.0; sse_v = 0.0; count = 0; linf = 0.0;
    if numel(bounds) < 2, rmse = NaN; return; end

    % reference scales
    h_scale = max(1e-12, max(h_ref) - min(h_ref));
    v_scale = max(1e-12, max(abs(v_ref)));

    for k = 1:numel(bounds)-1
        a = bounds(k); b = bounds(k+1);
        if b <= a, continue; end
        eps = min(1e-6, 1e-6*t_end);
        aa = a + eps; bb = b - eps;
        if bb <= aa, aa = a; bb = b; end
        grid = linspace(aa, bb, 201);
        hc = interp1(T, Y(:,1), grid, 'linear', 'extrap');
        vc = interp1(T, Y(:,2), grid, 'linear', 'extrap');
        hrefg = interp1(t_ref, h_ref, grid, 'linear', 'extrap');
        vrefg = interp1(t_ref, v_ref, grid, 'linear', 'extrap');
        dh = hc - hrefg; dv = vc - vrefg;
        sse_h = sse_h + sum(dh.^2);
        sse_v = sse_v + sum(dv.^2);
        count = count + numel(dh);
        linf = max(linf, max([max(abs(dh))/h_scale, max(abs(dv))/v_scale]));
    end
    if count == 0
        rmse = NaN;
    else
        rmse_h = sqrt(sse_h / count);
        rmse_v = sqrt(sse_v / count);
        rmse = 0.5 * (rmse_h / h_scale + rmse_v / v_scale); % NRMSE
    end
end
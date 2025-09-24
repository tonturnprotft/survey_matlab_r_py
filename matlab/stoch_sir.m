function stoch_sir(beta, gamma, N, I0, t_end, n_rep, seed, outpath, refcsv)
% Usage: stoch_sir(0.3/1000,0.1,1000,10,160,10,42,'results_ssa_sir_matlab.csv','ref_sir_mean.csv')
    if nargin<9, refcsv=''; end
    S0 = N - I0; rng(seed); 
    grid = 0:1:t_end; 
    acc = zeros(numel(grid),3); 
    steps_total=0;
    
    %timing
    t0 = tic;
    for r=1:n_rep
        [T,X] = sir_once(beta,gamma,S0,I0,0,t_end);
        steps_total = steps_total + (numel(T)-1);
        acc = acc + series_to_grid(T,X,grid);
    end
    wall = toc(t0);
    mean_series = acc/n_rep; 
    rmse=NaN; linf=NaN;
    if ~isempty(refcsv) && isfile(refcsv)
        ref = readtable(refcsv); 
        refm = [ref.S_mean, ref.I_mean, ref.R_mean];
        D = mean_series - refm; 
        rmse = sqrt(mean(D(:).^2)); 
        linf = max(abs(D(:)));
    end
    row.platform="matlab"; row.family="STOCH"; row.benchmark_id="ssa_sir";
    row.tool="event_loop"; row.solver="gillespie"; row.status="OK";
    row.wall_time_s= wall; row.rmse=rmse; row.linf=linf; row.n_steps=steps_total;
    row.params_json=jsonencode(struct('beta',beta,'gamma',gamma,'N',N,'I0',I0,'t_end',t_end,'n_rep',n_rep,'seed',seed));
    T = struct2table(row);
    if isfile(outpath), writetable(T,outpath,'WriteMode','Append','QuoteStrings',true); else, writetable(T,outpath,'QuoteStrings',true); end
    fprintf('status=OK n_steps=%d rmse=%g linf=%g\n', steps_total, rmse, linf);

function [times,states]=sir_once(beta,gamma,S,I,R,t_end)
    t=0; times=0; states=[S,I,R];
    while t<t_end && I>0
        a1=beta*S*I; a2=gamma*I; a0=a1+a2; if a0<=0, break; end
        tau=-log(rand)/a0; t=t+tau;
        if rand*a0 < a1 && S>0, S=S-1; I=I+1; else, I=I-1; R=R+1; end
        times(end+1)=min(t,t_end); states(end+1,:)=[S,I,R]; %#ok<AGROW>
    end
    if times(end)<t_end, times(end+1)=t_end; states(end+1,:)=[S,I,R]; end
    end

function Y=series_to_grid(T,X,grid)
    S=interp1(T,X(:,1),grid,'previous','extrap'); I=interp1(T,X(:,2),grid,'previous','extrap'); R=interp1(T,X(:,3),grid,'previous','extrap'); Y=[S(:),I(:),R(:)];
    end
end
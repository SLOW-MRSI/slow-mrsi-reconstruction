function runParallel = setupTsEPSIParallel(useParallel)
%SETUPTSEPSIPARALLEL Always use a MATLAB thread pool when possible.

runParallel = false;

if ~useParallel || ~license('test','Distrib_Computing_Toolbox')
    return
end

try
    pool = gcp('nocreate');

    % If an old process pool is open, close it first.
    if ~isempty(pool) && ~isa(pool, 'parallel.ThreadPool')
        delete(pool);
        pool = [];
    end

    % Start thread pool if no pool is active.
    if isempty(pool)
        pool = parpool('threads');
    end

    runParallel = ~isempty(pool);

    if runParallel
        fprintf('Parallel processing enabled with THREAD pool: %d workers.\n', pool.NumWorkers);
    end

catch ME
    warning('MATLAB:RunTsEPSIParallelSetup', ...
        'Thread-pool setup failed (%s). Running serial loops.', ME.message);
end
end
function iters = iters_from_maxfes(MaxFEs, N)
% init evaluates N solutions once, so subtract that.
if MaxFEs <= N
    iters = 0;
else
    iters = floor((MaxFEs - N) / N);
end
end

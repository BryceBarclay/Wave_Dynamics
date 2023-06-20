function samples = unifrndintervals(intervals,N)
% sample N values uniformly from a collection of disjoint intervals

widths = intervals(:,2) - intervals(:,1);
cumwidth = cumsum(widths);
samples = zeros(N,1);
for i = 1:N
    random = unifrnd(0,cumwidth(end));
    intchoice = sum(cumwidth < random) + 1;
    samples(i) = random - cumwidth(intchoice) + intervals(intchoice,2);
end

end
function [q_pdf, q_range] = pdf_var(q, res, scaleflag, varargin)
% Returns the pdf and range vectors of quantity q.
% Inputs:
% q: vector containing observations of a variable
% res: number of bins to use in the pdf
% scaleflag: use 'true' for log-scaled bins, 'false' for linear bins
% qlims: (optional) specify [qmin qmax]

if isempty(q)
    q_pdf = nan;
    q_range = nan;
    return
end

% define upper and lower bin bounds
if ~isempty(varargin)
    qlims = varargin{1};
else
    qlims = [min(q), max(q)];
end

% define bin edges
if scaleflag
    q_range = real(logspace(log10(qlims(1)),log10(qlims(2)),res+1));
else
    q_range = linspace(qlims(1),qlims(2),res+1);
end

% histogram counts
q_pdf = histcounts(q(:),q_range);

% bin centers
q_range = q_range(2:end) - diff(q_range)/2;

% normalize so integral = 1
q_pdf = q_pdf/(sum(q_pdf)*diff(q_range(1:2)));

end
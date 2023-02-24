function [outputsearchvals]=ReduceSearchRangeAroundBest(inputsearchvals,bestval,minval,maxval)
% helper file for fitting functions: reduces search range to a smaller
% subset of values with a finer spacing around the best value of the
% previous range. Spacing will be half previous spacing, number of values
% will be the same (so half the total range)

% inputsearch values is a vector of values
% bestval is the value to search more finely around (single value)
% minval and maxval are extremes of the new range allowed. Range will be adjusted to stay within bounds
% outputsearchvalues is the new set of values to test (with finer spacing & smaller range)

if nargin<3
    minval=-Inf;maxval=Inf;
end
npts=length(inputsearchvals);
newdiff=mean(diff(inputsearchvals))/2;
if mod(npts,2) %odd
    outputsearchvals=(-(npts-1)/2*newdiff+bestval):newdiff:((npts+1)/2*newdiff+bestval);
else % even
    outputsearchvals=(-npts/2*newdiff+bestval):newdiff:(npts/2*newdiff+bestval);
end

if outputsearchvals(1)<minval
    outputsearchvals=outputsearchvals-outputsearchvals(1)+minval;
end
if outputsearchvals(end)>maxval
    outputsearchvals=outputsearchvals-outputsearchvals(end)+maxval;
end
end

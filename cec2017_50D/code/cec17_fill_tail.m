function curve = cec17_fill_tail(curve)
% Fill unused curve tail with the latest best-so-far objective value.
last = find(~isnan(curve), 1, 'last');

if isempty(last)
    curve(:) = inf;
else
    curve(last + 1:end) = curve(last);
end
end

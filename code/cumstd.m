function cstd = cumstd(y)

cstd = zeros(size(y));
for j = 1:length(y)
    cstd(j) = std(y(1:j));
end
cstd = cstd./std(y);
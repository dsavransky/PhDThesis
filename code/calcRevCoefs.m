function xkmat = calcRevCoefs(n)
% calcRevCoefs - Find series reversion coefficient sets
%
% INPUTS
%  calcRevCoefs(n) returns all sets of natural numbers such that
%  sum(i*x_i) = n. n must be greater than 1. 
% 
% OUTPUT
%  xkmat    N x n matrix with each row representing a unique set,
%           such that xkmat(j,:)*(1:n).' = n

%check input
if ~isa(n,'numeric') || floor(n) ~= n || n < 2
    error('calcRevCoefs:inputError',...
        'n must be a natural number greater than 1');
end

%generate original pairs (excluding (0,n))
matr = genBitPairs(n);

%expand
for k = 3:n
    %there will be this many new subsets:
    newCombs = floor(matr(:,end)/2); 
    oDim = size(matr); %orig dimensions
    newDim = oDim + [sum(newCombs),1]; %new dimensions
    matr(newDim(1),newDim(2)) = 0; %resize matrix
    
    %only expand rows with new combinations
    for j = find(newCombs).' 
        if matr(j,oDim(2)) > 1
            matr(oDim(1)+1+sum(newCombs(1:j-1)):...
                oDim(1)+sum(newCombs(1:j)),:) = ...
                [repmat(matr(j,1:oDim(2)-1),newCombs(j),1),...
                genBitPairs(matr(j,oDim(2)))];
        end
    end
    %keep only unique combinations
    [~,i] = unique(sort(matr,2),'rows','first');
    matr = matr(i,:);
end
y = sort(matr,2);

%xkmat will be the same size as y
xkmat = zeros(size(y));
%for all non-zero entries in y, add 1 to xkmat 
%in the same row as y, and column given by the y entry 
%(i.e., if y has an i anywhere in row j, add 1
%to xkmat at index (i-1)*size(y,1) + j
inds = (y-1)*size(y,1);
good = inds >= 0;
inds = inds + repmat((1:size(y,1)).',1,n);
inds = inds(good);
uinds = unique(inds);
%returns number of all unique indices:
xkmat(uinds) = histc(inds,uinds); 

%pad output with initial row and package
xkmat = [zeros(1,n);xkmat]; 
xkmat(1,n) = 1;
end

function bitpairs = genBitPairs(n)
%generate 2 x floor(n/2) matrix containing all
%integer pairs summing to n excluding (0,n)

if n <= 1
    bitpairs = [];
    return
end

bitpairs = zeros(floor(n/2),2); 
bitpairs(:,1) = 1:floor(n/2); 
bitpairs(:,2) = n - bitpairs(:,1);

end
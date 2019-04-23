function sdists = genCostFunTerm1(sinds, u_ts)
%  genCostFunTerm1(sinds, u_ts) generates the first cost 
%  function term for available target indices (s_inds) and
%  target unit look vectors (u_ts)

%combinations of visits
combs = nchoosek(sinds,2);

%angle between each pair of target stars
angdists = acos(dot(u_ts(combs(:,1),:),...
    u_ts(combs(:,2),:),2)); 

%properly ordered cost function term
sdists = zeros(length(sinds));
sdists((tril(ones(length(sinds))) - ...
    eye(length(sinds)) == 1)) = angdists;
sdists = sdists+sdists.';
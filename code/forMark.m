dr = '~/Documents/NRO/';

[allres,res,data] = analyzeDRM_v3([dr,'nro_small_kepler_mod_case2b.mat']);

inds = [];
for j = 1:length(res.universes)
        inds = [inds;res.universes{j}.planInds(logical(res.universes{j}.observed))];
end

xs = 1:length(data.stars.NAME);
ns = hist(inds,xs);

[nss,is] = sort(ns,2,'descend');
xss = xs(is);

names = data.stars.NAME(xss(1:50));

fid = fopen([dr,'top_targs.txt'],'w');
for j =1:50, fprintf(fid,'%s\n',names{j}); end
fclose(fid)
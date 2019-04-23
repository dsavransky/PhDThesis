flds = {'HIP','RADEG','DECDEG','DIST','PARX','PMRA','PMDEC','V','V_I','J__MAG','K__MAG','H__MAG','L'};%,'SPEC'};
typs = {'%10u','%12.10f','%12.10f','%12.10f','%12.10f','%12.10f','%12.10f','%12.10f','%12.10f','%12.10f','%12.10f','%12.10f','%s'};

fid = fopen('TPF_TLDB.txt','w');
for k=1:length(dat.S.HIP)
    for j = 1:length(flds)
        fprintf([typs{j},'\t'],dat.S.(flds{j})(k));
    end
    fprintf(fid,'\n');
end
fclose(fid)



out = zeros(length(dat.S.HIP),length(flds));
for j = 1:length(flds)
    out(:,j) = dat.S.(flds{j});
end
save 'TPF_TLDB.txt' out -ASCII

specs = dat.S.SPEC;
fid = fopen('TPF_TLDB_specs.txt','w');
for j = 1:length(specs)
    fprintf(fid,'%s\n',specs{j});
end
fclose(fid)
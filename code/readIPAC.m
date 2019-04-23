function [headers,units,data] = readIPAC(fname)

if ~exist(fname,'file')
    error('readIPAC:fileNotFound',[fname,' does not exist.']);
end

fid = fopen(fname);
raw = textscan(fid, '%s', 'delimiter', '\n');
commEnd = max(strmatch('\',raw{1})); %last comment line

headers = raw{1}{commEnd+1};
headers = textscan(headers, '%s', 'delimiter', ',');
headers = headers{1};

fstring = repmat('%s ',1,length(headers));

units = raw{1}{commEnd+2};
units = textscan(units, '%s', 'delimiter', ',');
units = units{1};

frewind(fid);
data = textscan(fid,fstring,'delimiter', ',','headerlines',commEnd+2);

data = reshape(cat(1,data{:}),length(data{1}),length(data));
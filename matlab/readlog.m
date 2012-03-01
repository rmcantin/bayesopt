%
%
%

fid = fopen('micha_bopt.log')
data = textscan(fid, '%f|%f|%f|[%d](%f,%f,%f,%f,%f)');
fclose(fid);
function val = read_directional_spectra_bulkstat(csvfile,directional_methods)


opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = directional_methods;
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
valTab = readtable(csvfile, opts);

for i = 1:length(opts.VariableNames)
    eval(['val(i) = valTab.',char(opts.VariableNames(i)),';'])
end

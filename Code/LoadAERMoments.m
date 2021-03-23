%% LOAD PUBLIC DATA FOR REPLICATION - AER moments

%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\asus\Desktop\Giulia\UBC\Year2\567 - Empirical IO\AFT2017-Replication\Code\Original files\AER1_moments_2007.out
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2021/03/22 15:57:05

%% Initialize variables.
filename = 'C:\Users\asus\Desktop\Giulia\UBC\Year2\567 - Empirical IO\AFT2017-Replication\Code\Original files\AER1_moments_2007.out';
delimiter = ',';

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,4,5,6,7,8,9,10,11,12]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,3,4,5,6,7,8,9,10,11,12]);
rawStringColumns = string(raw(:, 2));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
AER1moments2007 = table;
AER1moments2007.code_TD = cell2mat(rawNumericColumns(:, 1));
AER1moments2007.iso3 = rawStringColumns(:, 1);
AER1moments2007.fe = cell2mat(rawNumericColumns(:, 2));
AER1moments2007.distw = cell2mat(rawNumericColumns(:, 3));
AER1moments2007.contig = cell2mat(rawNumericColumns(:, 4));
AER1moments2007.comlang_off = cell2mat(rawNumericColumns(:, 5));
AER1moments2007.control_of_corruption = cell2mat(rawNumericColumns(:, 6));
AER1moments2007.gdp = cell2mat(rawNumericColumns(:, 7));
AER1moments2007.rule_of_law = cell2mat(rawNumericColumns(:, 8));
AER1moments2007.exp_cost = cell2mat(rawNumericColumns(:, 9));
AER1moments2007.inputs_rounded = cell2mat(rawNumericColumns(:, 10));
AER1moments2007.firm_count_rounded = cell2mat(rawNumericColumns(:, 11));

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R;

%% Delete first row (NaN)
AER1moments2007(1,:) = [];

%% Save to excel format
filename = 'AER1moments2007.csv';
writetable(AER1moments2007,filename,'Delimiter',',')

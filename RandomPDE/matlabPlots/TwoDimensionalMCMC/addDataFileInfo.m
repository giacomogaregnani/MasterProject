function addDataFileInfo(filename, text)

fileID = fopen(filename, 'a');
 
for i = 1 : length(text)
    fprintf(fileID, ['\n', text{i}, '\n']);
end
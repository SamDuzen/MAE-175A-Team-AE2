function [dummy,column]=readecp(filename,var1,var2,var3,var4,var5,var6,var7,var8,var9)
%READECP - READ raw data format file generated by ECP
% [data,datainfo]=readecp(filename,var1,var2,...)
%will read a ASCII file generated by ECP using the Export Raw Data option.
%Output is data matrix and datainfo that contains information on the names 
%of the columns of the data. Additional column names var1, var2 etc. can be
%specified to only select certain data columns from the ASCII file generated 
%by ECP using the Export Raw Data option. 

% Written by R.A. de Callafon, Dept. of MAE, UCSD (2006-2025)
% Report errors in this software to <callafon@ucsd.edu>


if nargin<1,
    error('filename required')
end

% reset empty column names
for k=nargin:9,
    eval(['var ' num2str(k) '=[];']);
end

if exist(filename)~=2,
    error(['cannot find ''' filename '''']);
end        

fid=fopen(filename);
if fid<0,
    error(['cannot open ''' filename '''']);
end        
    
% read it binary
%filedata=fread(fid);
% unfortunately we have to do it in chuncks of 'Nmax' characters and save
% information linewise to avoid a maxsize error warning in case 
% of a student version of Matlab. A bit more work, but possible:

% to be on the save side:
Nmax=15000;

% first we have to determine the longest line in the file so
% we know how large the matrix will be to store the data in file
temp='get started!';
maxlinelength=0;
lines=0;
while ~isempty(temp),
    temp=fread(fid,Nmax);
    full_length_temp=length(temp);
    % to handle the char(10) that appears after the char(13) in a DOS ASCII file
    if temp(1)==10,
        temp=temp(2:end);
    end
    % find the char(13)'s to determine the largest linelength
    temp_chr13index = find(temp==char(13));
    if isempty(temp_chr13index),
        % no char(13) found, line must be crazy long
        maxlinelength=Nmax;
        % stop iteration to read file
        full_length_temp=0;
    else
        n=length(temp_chr13index);
        maxlinelength=max([maxlinelength;temp_chr13index-[0;temp_chr13index(1:n-1)]-1]);
        lines=lines+n;
    end
    if full_length_temp<Nmax,
        % there is nothing more to read
        temp=[];
    else
        % reset file pointer on last char(13)+1 for next read operation on file
        fseek(fid,temp_chr13index(n)-length(temp),0);
    end    
end    

if maxlinelength>=Nmax,
    % something is wrong, give error message
    error(['lines too long: file ''' filename ''' is missing EOL characters or not a valid ECP raw data format file'])
end    
if lines>Nmax,
    % something is wrong, give error message
    error(['too many lines: file ''' filename ''' is too large or not a valid ECP raw data format file'])
end    
fclose(fid);

% size of variable geing created:
%maxlinelength, lines


% reopen and read the file, but now we will actually save 
% the information linewise. We don't do any error check anymore
% to speed up the process

fid=fopen(filename);

% create initial data matrix to be filled
filedata=zeros(lines,maxlinelength);
temp='get started!';
linecounter=0;
while ~isempty(temp),
    temp=fread(fid,Nmax);
    full_length_temp=length(temp);
    % find all char(13) to store all lines
    temp_chr13index = find(temp==char(13));
    n=length(temp_chr13index);
    % add a 0 index the following for loop to work
    temp_chr13index=[0;temp_chr13index];
    for k=1:n,
        linedata=temp(temp_chr13index(k)+1:temp_chr13index(k+1)-1);
        if ~isempty(linedata),
            % to handle the char(10) that appears after the char(13) in a DOS ASCII file        
            if linedata(1)==10,
                linedata=linedata(2:end);
            end
        end
        linecounter=linecounter+1;
        % assign the line (without chr(13) or chr(10)) to the filedata variable
        filedata(linecounter,:)=[linedata' zeros(1,maxlinelength-length(linedata))];
    end
    % remove 0 index again
    temp_chr13index=temp_chr13index(2:end);

    if full_length_temp<Nmax,
        % assign the last line if there was more behind last char(13)
        if temp_chr13index(n)<length(temp),
            linedata=temp(temp_chr13index(n)+1:full_length_temp);
            % to handle the char(10) that appears after the char(13) in a DOS ASCII file
            if linedata(1)==10,
                linedata=linedata(2:end);
            end
            linecounter=linecounter+1;
            % assign the line (without chr(13) or chr(10)) to the filedata variable
            filedata(linecounter,:)=[linedata' zeros(1,maxlinelength-length(linedata))];

            % there is nothing more to read
            temp=[];
        end            
    else
        % reset file pointer on last char(13)+1 for next read operation of file
        fseek(fid,temp_chr13index(n)-length(temp),0);
    end
end    
fclose(fid);

% find begin and end of data by looking for `[' and a `]'
[data_begin1,data_begin2]=find(filedata==91);
[data_end1,data_end2]=find(filedata==93);
if (isempty(data_begin1))|(isempty(data_end1)),
    % something is wrong, give error message
    error(['missing ''['' or '']'': file ' filename ' is not a valid ECP raw data format file'])
end
if ((size(data_begin1)>1)|(size(data_end1)>1)),
    % something is wrong, give error message
    error(['too many ''['' or '']'': file ' filename ' is not a valid ECP raw data format file'])
end

% find the line that identifies what was measured by looking for the word 'Sample'
[sample_index1,sample_index2]=find(filedata==83);
line_begin=[];
for k=1:length(sample_index1),
    % for each letter 'S' we found, see if the word `Sample' can be constructed
    if sample_index2(k)<(maxlinelength-5),
        if strcmp(char(filedata(sample_index1(k),sample_index2(k):sample_index2(k)+5)),'Sample'),
            line_begin=[line_begin sample_index2(k)];
        end
    end
end            
if isempty(line_begin)
    % something is wrong, give error message
    error(['missing ''Sample'' codeword: file ' filename ' is not a valid ECP raw data format file'])
end
if length(line_begin)>1
    % something is wrong, give error message
    error(['too many ''Sample'' codewords: file ' filename ' is not a valid ECP raw data format file'])
end

% find the first char(0) in the row  of filedata 
% where we found the code word 'Sample'
line_end=find(filedata(sample_index1(1),:)==0);
if isempty(line_end),
    line_end=maxlinelength+1;
end

% get data info
data_info=char(filedata(sample_index1(1),line_begin:line_end-1));

% scan data info for information
% delete 'Sample'
data_info=data_info(7:end);
% eliminate trailing spaces
index=find(data_info==' ');
while index(1)==1,
    data_info=data_info(2:end);
    index=find(data_info==' ');
end
% next item should be 'Time'
if strcmp(data_info(1:4),'Time')~=1
    % something is wrong, give error message
    error(['missing ''Time'' codeword in data info string: file ' filename ' is not a valid ECP raw data format file'])
end
% delete 'Time'
data_info=data_info(5:end);
% eliminate trailing spaces
index=find(data_info==' ');
while index(1)==1,
    data_info=data_info(2:end);
    index=find(data_info==' ');
end

% Now find information on data columns by looking for two spaces
column_counter=1;
column{column_counter} = 'Sample';
column_counter=2;
column{column_counter} = 'Time';
addional_column=1;
while addional_column==1,
    addional_column=0;    
    index=find(data_info==' ');
    for k=1:length(index),
        if (index(k)<length(data_info))&(addional_column==0),
            if data_info(index(k)+1)==' ', 
                % there is an additional data column
                addional_column=1;
                column_end=index(k);
            end
        end
    end
    column_counter=column_counter+1;    
    if addional_column==0,
        % apparently no additional data columns, so select remaining data info
        column{column_counter} = data_info(1:end);
    else
        % apparently an additional data columns, so select current info
        column{column_counter} = data_info(1:column_end-1);
        % and delete this from data_info
        data_info=data_info(column_end:end);
        % eliminate trailing spaces
        index=find(data_info==' ');
        while index(1)==1,
            data_info=data_info(2:end);
            index=find(data_info==' ');
        end
    end
end    

disp(['- Found and read ' num2str(max(size(column))) ' data columns from ''' filename ''':']);
disp(column)

% now data to the dummy variable by evaluation
% assign initial size to dummy
dummy=zeros(data_end1(1)-data_begin1(1)+1,size(column,2));
% find chr(0) we added to make filedata a square matrix, as we need to
% eliminate there before evaluation
chr0_index=find(filedata(data_begin1(1),1:maxlinelength)==0);
% first row
eval(['dummy(1,:) = ' char(filedata(data_begin1(1),data_begin2(1):chr0_index-1)) '];' ]);
for k=2:data_end1(1)-data_begin1(1),
    % subsequent rows
    chr0_index=find(filedata(data_begin1(1)+k-1,1:maxlinelength)==0);    
    eval(['dummy(' num2str(k) ',:) = [' char(filedata(data_begin1(1)+k-1,1:chr0_index-1)) '];' ]);
end
% last row
k=data_end1(1)-data_begin1(1)+1;
chr0_index=find(filedata(data_begin1(1)+k-1,1:maxlinelength)==0);    
t=['dummy(' num2str(k) ',:) = [' char(filedata(data_begin1(1)+k-1,1:chr0_index-1)) ';' ];
eval(['dummy(' num2str(k) ',:) = [' char(filedata(data_begin1(1)+k-1,1:chr0_index-1)) ';' ]);

% now we can clear filedata
clear filedata
    
if nargin>1,
    columnindex=[];
    % select specified columns
    columnselection=['- Selecting data columns: '];    
    for m=1:nargin-2,
        columnselection=[columnselection '''' eval(['var' num2str(m)]) ''' and '];
    end
    m=nargin-1;
    columnselection=[columnselection '''' eval(['var' num2str(m)]) ''''];
    disp(columnselection)
    for m=1:nargin-1,
        eval(['columnname=var' num2str(m) ';']);
        columnnamematch=0;    
        for k=1:max(size(column)),
            if strcmp(column{k},columnname),
                columnindex=[columnindex k];
                columnnamematch=1;
            end
        end
        if columnnamematch==0,
            error(['column data ''' columnname ''' was not found in ''' filename '''!'])
        end
    end
    dummy=dummy(:,columnindex);
    column=column(columnindex);
end
    

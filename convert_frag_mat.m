close all
clearvars


%%% Converting a fragment file (ProbHAP) to sparse matrix



% N number of reads
% l haplotype length




File='/SinaMc/University/WUR/WURcode/wur_code1/convertingfragtomatrix/24_1/fragFiltered';
a1=readtable(File,'Delimiter','\t');%,'HeaderLines',0
a=table2cell(a1);
fragment_cell=a;
N=size(a1,1);
l=175;
R=sparse(N,l);
    

for i=1:N
    row=zeros(1,l);
    %string = a{i}; %1 chr22_SPH2_1940 3 100 C==
    
    
    
    
    %%%%%%% for haplogenerator 
    string = a{N-i+1};
    
    
    
    
    
    
    
    ii=1;
    String_BlockNumber=[]; %for two digits BlockNumber
    while isspace(string(ii))~=1
        String_BlockNumber=[String_BlockNumber,string(ii)];
        ii=ii+1;
    end
    BlockNumber=str2double(String_BlockNumber);
    StartingPoint=zeros(1,BlockNumber);
    ii=ii+1;
    while isspace(string(ii))~=1 % for going affter 'Chr2_...'
        ii=ii+1;
    end
    for j=1:BlockNumber
        ii=ii+1;
        String_StartingPoint=[]; %for two digits starting point
        while isspace(string(ii))~=1
            String_StartingPoint=[String_StartingPoint,string(ii)];
            ii=ii+1;
        end
        StartingPoint(j)=str2double(String_StartingPoint);
        site=StartingPoint(j);
        ii=ii+1;
        while isspace(string(ii))~=1  %
            row(site)=2*str2double(string(ii))-1;   % put the reads  in 'row'  for each block
            site=site+1;
            ii=ii+1;
        end 
    end
    R(i,:)=row;
end


clearvars -except R  fragment_cell

save('24_1/R24_1.mat','-v7.3')

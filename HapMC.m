clearvars
load('/SinaMc/University/WUR/WURcode/Data24_1/R24_1.mat')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Haplotyping using HapSVT


% Input: read matrix in .mat format
% output: a text file containg haplotypes and corresponding variant index


% the longest block is considered 2000


%Sina Majidian Dec 2018
%Iran University of Science and Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% removing those caloumn which no  read covered it.
% maybe it covered at first to be used in freebayes, maybe some of reads covering less than two snp
%are removed
nonzero_col_idx=find(sum(abs(R))~=0);%index of those columns with at least one elemnt
R_in=R(:,nonzero_col_idx);
R_in=R_in(1:50,1:50);

fileID = fopen('out/hap_opt24_1.txt','w'); % The output file

block_idx=0;columnNumber_blocks_accm=0;rowNumber_blocks_acc=0; % just for starting
idxi_vector=[];idxj_vector=[];block_num_col=[]; start=1;block_idx=0;
rowNumber_blocks=[];columnNumber_blocks=[];

while columnNumber_blocks_accm(end)<size(R_in,2)-1  % run for each block until last column
    block_idx=block_idx+1
    
    if size(R_in,2)<columnNumber_blocks_accm(end)+2000
        % the start column is the last value of accumaltion all haploblock length +1
        % the end column is .+2000
        R_sliced2000=R_in(rowNumber_blocks_acc(end)+1:end,columnNumber_blocks_accm(end)+1:end);
    else
        R_sliced2000=R_in(rowNumber_blocks_acc(end)+1:end,columnNumber_blocks_accm(end)+1:columnNumber_blocks_accm(end)+2000);
    end
    
    if start==1  % for the first block
        start=0;rowNumber_blocks_acc=[];columnNumber_blocks_accm=[];
    end
    
    %Exctracting each block of read matrix which all reads have overlap
    [rowNumber_block,columnNumber_block,R_block]=first_block_extractor(R_sliced2000);%block_length is the haploBlock length
    
    rowNumber_blocks=[rowNumber_blocks, rowNumber_block];columnNumber_blocks=[columnNumber_blocks, columnNumber_block];
    rowNumber_blocks_acc=cumsum(rowNumber_blocks);columnNumber_blocks_accm=cumsum(columnNumber_blocks);
    if  size(R_block,2)>1
        
        if block_idx==1 %for reporting the index of estimated allele of variant to not lossing
            %the index in overall, we remove some columns
            indces_block=nonzero_col_idx(1:columnNumber_blocks_accm(1));
        else
            indces_block=nonzero_col_idx(columnNumber_blocks_accm(block_idx-1)+1:columnNumber_blocks_accm(block_idx));
        end
        %processing of read matrix block
        nonzeor_idx_row=find(sum(abs(R_block'))>1); % those rows with at least two nonzero
        R_used=R_block(nonzeor_idx_row,:);
        %%%%% haplotpying  using HapOPT
        [X,S_opt,Y,dist] = OptSpace(R_used,3,500,.00001); %(R,[],500,.001)  matrix, rank,number iter, toleranc
        X_opt=X*S_opt*Y';
        %%% SVT
        
        % % % %%% Hapsvt
        % % Omega=find(R_used);
        % % Rt=(R_used+1)/2;
        % % mu_final=.01;%100*max(N,l);
        % % [U,S,V,numiter]= FPC(size(R_used),Omega,Rt(Omega),mu_final);  %FPC(n,Omega,b,mu_final,maxiter,tol)
        % % X1=U*S*V';[u_sv,S_sv,v_sv]=svds(X1,3); X_svt=u_sv*S_sv*v_sv';
        % % A=X_svt';
        A=X_opt';
        [~,colind] = rref(A);
        Xsub = A(:, colind);
        H=Xsub'>0;
        H_with_ind=[indces_block',H'];
        
        fprintf(fileID,'Block len=%d\n',columnNumber_block);
        fprintf(fileID,'%d\t%d\t%d\t%d\n',H_with_ind');
    end
end
fclose(fileID);


% % %norm(R_used(Omega)- X_svt(Omega))
% % R_logc=R>0;
% % a1=sum(R_logc);
% % a2=sum(abs(R));
% % a2./a1

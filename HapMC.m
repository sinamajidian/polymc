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


% R1=full(R(1:200,1:400));
% % % 
% % % nonzero_col_idx=find(sum(abs(R))~=0);%index of those columns with at least one elemnt 
% % % %zero_idx_col=find(sum(abs(R))==0);
% % % R_in=R(:,nonzero_col_idx);


R_in=R(1:200,1:400);
fileID = fopen('out/hap_opt24_1.txt','w'); % The output file

block_idx=0;columnNumber_block=0; % just for starting
idxi_vector=[];idxj_vector=[];block_num_col=[]; start=1;


while columnNumber_blocks_accm(end)<size(R_in,2)-1  % run for each block until last column
    if start==1  % for the first block
        start=0;block_idx=0;rowNumber_blocks_acc=[];columnNumber_blocks_accm=[];
    end
    block_idx=block_idx+1
    
    
    if size(R_in,2)<haploblock_length_accm(end)+2000
        % the start column is the last value of accumaltion all haploblock length +1
        % the end column is .+2000
        R_sliced2000=R_in(rowNumber_blocks_acc(end)+1:end,columnNumber_blocks_accm(end)+1:end);
    else
        R_sliced2000=R_in(rowNumber_blocks_acc(end)+1:end,columnNumber_blocks_accm(end)+1:columnNumber_blocks_accm(end)+2000);
    end
    %Exctracting each block of read matrix which all reads have overlap    
    [rowNumber_block,columnNumber_block,R_block]=first_block_extractor(R_sliced2000);%block_length is the haploBlock length

    rowNumber_blocks=[rowNumber_blocks, rowNumber_block];columnNumbers_block=[columnNumber_blocks, columnNumber_block];
    rowNumber_blocks_acc=cumsum(rowNumber_blocks);columnNumber_blocks_accm=cumsum(columnNumber_blocks);
    if  size(R_block,2)>1
        if block_idx==1
            indces_block=nonzero_col_idx(1:columnNumber_blocks_accm(1));
        else
            indces_block=nonzero_col_idx(columnNumber_blocks_accm(block_idx-1)+1:columnNumber_blocks_accm(block_idx));
        end
        
        
        
        
        R2=R_block;
        %%%%%% haplotpying
        nonzeor_idx_row=find(sum(abs(R2'))>1); % those rows with at least two nonzero
        R_used=R2(nonzeor_idx_row,:);
        % %%%%% Hap opt
        [X,S_opt,Y,dist] = OptSpace(R_used,1,500,.00001); %(R,[],500,.001)  matrix, rank,number iter, toleranc
        X_opt=X*S_opt*Y';
        A=X_opt';
        [~,colind] = rref(A);
        Xsub = A(:, colind);
        h=Xsub'>0;
        
        
        H=[inces_block',h'];
        %h_all=zeros(1,size(R1,2))+NaN;
        %h_all(nonzeor_idx_col)=h
        
        % dlmwrite('file_OPT_H_21_8_a',H,'\t')
        fprintf(fileID,'Block len=%d\n',block_num_col);
        fprintf(fileID,'%d\t%d\n',H');
        
    end
end



fclose(fileID);









%
% R=full(R);
% R1=R(1:20,1:20);
% %R2=R1;
% % find(R2~=0)
% % Rtst=R2;
% % Rtst(find(R2==0))=NaN;
% % Rtst=.5*(Rtst+1);
% nonzeor_idx_col=find(sum(abs(R1))~=0); % column with all elemnt zero
% R2=R1(:,nonzeor_idx_col);
% nonzeor_idx_row=find(sum(abs(R2'))>1); % those rows with at least one nonzero
% R_used=R2(nonzeor_idx_row,:);
%
%
%
% %%%%% Hap opt
% [X,S_opt,Y,dist] = OptSpace(R_used,3,500,.00001); %(R,[],500,.001)  matrix, rank,number iter, toleranc
% X_opt=X*S_opt*Y';
% % X_tr=X_opt;
% % [Q_qr, R_qr, E_qr] = qr(X_tr,0);
% % Xsub=X_tr(sort(E_qr(1:3)),:);
% %
% A=X_opt';
% [~,colind] = rref(A);
% Xsub = A(:, colind);
% H=Xsub>0
% dlmwrite('H_OPT24_1_1',H,'\t')
%
%
%
%
%
% %%%%
% %
% % nonzeor_idx_col=find(sum(abs(R1))~=0); % column with all elemnt zero
% % R2=R1(:,nonzeor_idx_col);
% % nonzeor_idx_row=find(sum(abs(R2'))>1); % those rows with at least one nonzero
% % R_used=R2(nonzeor_idx_row,:);
%
% % % %%% Hapsvt
% % Omega=find(R_used);
% % Rt=(R_used+1)/2;
% % mu_final=.01;%100*max(N,l);
% % [U,S,V,numiter]= FPC(size(R_used),Omega,Rt(Omega),mu_final);  %FPC(n,Omega,b,mu_final,maxiter,tol)
% % X1=U*S*V';[u_sv,S_sv,v_sv]=svds(X1,3); X_svt=u_sv*S_sv*v_sv';
% % A=X_svt';
% % [~,colind] = rref(A);
% % Xsub = A(:, colind);
% % H=Xsub>0;
% % dlmwrite('H_SVT21_9_a.txt',H,'\t')
% %
% % %norm(R_used(Omega)- X_svt(Omega))
% % R_logc=R>0;
% % a1=sum(R_logc);
% % a2=sum(abs(R));
% % a2./a1
%

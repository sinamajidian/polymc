clearvars
%warning off
load('/SinaMc/University/WUR/WURcode/wur_code1/convertingfragtomatrix/Data24_1/R24_1.mat')

R=full(R);
R1=R(1:20,1:20);
%R2=R1;
% find(R2~=0)
% Rtst=R2;
% Rtst(find(R2==0))=NaN;
% Rtst=.5*(Rtst+1);
nonzeor_idx_col=find(sum(abs(R1))~=0); % column with all elemnt zero
R2=R1(:,nonzeor_idx_col);
nonzeor_idx_row=find(sum(abs(R2'))>1); % those rows with at least one nonzero
R_used=R2(nonzeor_idx_row,:);



%%%%% Hap opt
[X,S_opt,Y,dist] = OptSpace(R_used,3,500,.00001); %(R,[],500,.001)  matrix, rank,number iter, toleranc
X_opt=X*S_opt*Y';
% X_tr=X_opt;
% [Q_qr, R_qr, E_qr] = qr(X_tr,0);
% Xsub=X_tr(sort(E_qr(1:3)),:);
% 
A=X_opt';
[~,colind] = rref(A);
Xsub = A(:, colind);
H=Xsub>0
dlmwrite('H_OPT24_1_1',H,'\t')





%%%%

nonzeor_idx_col=find(sum(abs(R1))~=0); % column with all elemnt zero
R2=R1(:,nonzeor_idx_col);
nonzeor_idx_row=find(sum(abs(R2'))>1); % those rows with at least one nonzero
R_used=R2(nonzeor_idx_row,:);

% % %%% Hapsvt
% Omega=find(R_used);
% Rt=(R_used+1)/2;
% mu_final=.01;%100*max(N,l);
% [U,S,V,numiter]= FPC(size(R_used),Omega,Rt(Omega),mu_final);  %FPC(n,Omega,b,mu_final,maxiter,tol)
% X1=U*S*V';[u_sv,S_sv,v_sv]=svds(X1,3); X_svt=u_sv*S_sv*v_sv';
% A=X_svt';
% [~,colind] = rref(A);
% Xsub = A(:, colind);
% H=Xsub>0;
% dlmwrite('H_SVT21_9_a.txt',H,'\t')
% 
% %norm(R_used(Omega)- X_svt(Omega))
% R_logc=R>0;
% a1=sum(R_logc);
% a2=sum(abs(R));
% a2./a1


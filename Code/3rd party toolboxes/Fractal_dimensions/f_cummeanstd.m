%F_CUMMEANSTD "Cumulative" mean and standard deviation
%
%[m_t,s_t]=f_cummeanstd(y)
%
% y   is an MxN matrix.
% m_t is the MxN matrix of the means. The element m_t(j,k) is the 
%    mean of the first k elements of the jth row of y, that is: 
%    m_t(j,k)=mean(y(j,1:k));
% s_t is the MxN matrix of the standard deviations. The element 
%    s_t(j,k) is the standard deviation of the first k elements 
%    of the jth row of y, that is: 
%    s_t(j,k)=std(y(j,1:k));
%
% Note: s_t(j,k) and std(y(j,1:k)) could be not exactly equal 
%      beacuse of precision problems.
%
%See Also: MEAN, STD
%
%Gabriele Bulian, ruga@libero.it (30/05/2002)
function [m_t,s_t]=f_cummeanstd(y)
[R,Npti] = size(y);
z=y;
n=[1:1:Npti]; %Note that n(k)=k;
N=repmat(n,R,1);
m_t=cumsum(z,2)./N; %MEAN...
%m_t(j,k)=1/N(j,k)*sum{i=1,k}(z(j,i))=1/k*sum{i=1,k}(z(j,i))
%so m_t(j,k) is the mean of the firt k elements of the jth row of z(=y)
N(:,1)=NaN; %avoid division by zero for v_t(:,1)
v_t=(cumsum(z.^2,2)-N.*m_t.^2)./(N-1) ; %VARIANCE...
%v(j,k)=Var(z(j,1:k))=1/(N(j,k)-1)*sum{i=1,k}((z(j,i)-m_t(j,k))^2)=...
%1/(N(j,k)-1)*(sum{i=1,k}(z(j,i)^2)-N(,jk)*m_t(j,k)^2) 
%is the variance of the first k elements of the jth row.
s_t=sqrt(v_t);

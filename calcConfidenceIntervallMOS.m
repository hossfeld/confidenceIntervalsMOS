% CALCCONFIDENCEINTERVALLMOS  Computes confidence intervals for QoE data on a discrete 5-point rating scale
% The confidence intervals for the Mean Opinion Scores are computed.
% 
% Different approaches are used to derive the confidence intervals for MOS
% values. The first three approaches 1-3 estimate the confidence interval for
% mean values.
%
% The approaches 4-6 consider the fact that user ratings follow a discrete,
% truncated random variable for a certain TC. Then, the binomial distribution
% may be used as an upper bound (for the emerging variance due to the binomial
% sum variance inequality). Accordingly, CI estimators for binomial proportions
% are then utilized.
%
% 1. CI for Mean: Bootstrap algorithm
% 2. CI for Mean: Normal approximation for the MOS
% 3. CI for Mean: Student-T approximation for the MOS
% 4. CI for Bino: Wald interval employing normal approximation
% 5. CI for Bino: Wilson score interval with continuity correction
% 6. CI for Bino: Clopper-Pearson (using beta distribution)
% 7. Simultaneous CI: S. Jin, Computing Exact Confidence Coefficients of Simultaneous  Confidence, Intervals for Multinomial Proportions and their Functions. Department of Statistics, Uppsala University, 2013.
% 8. Jeffrey's interval: H. Jeffreys, The theory of probability. OUP Oxford, 1998.
%
% Input data: 
%    alpha  - significance level, default value is 0.05
%    y  - user ratings y of size k x n with k test conditions and n subjects
%
% Output data: The return value is a struct containing the following data.
%      str - name of the confidence interval estimators
%      shortstr - abbreviation (for plots)
%      mos - mean opinion score averaged over the n subjects for all k test conditions
%      ciUpper - the upper confidence intervall bound
%      ciLower - the lower confidence intervall bound
%      ciWidth - the size of the confidence interval; note that the CI is not symmetric around the MOS
%
% See also BOOTCI
function calc=calcConfidenceIntervallMOS(alpha,y)
if nargin<1,alpha=0.05;end
% If no data is given, then random user ratings are generated
if nargin<2    
    k=101; % number of test conditions
    n=20;  % number of subjects    
    
    y=zeros(k,n);
    p=zeros(k,1);
    % generate random user ratings
    for i=1:k
        p(i)=(i-1)/(k-1);
        y(i,:)=binornd(4,p(i),n,1)+1;
    end
    mos=mean(y,2);    
end
ci.alpha=alpha;
%%
str={'bootstrap','normal','student-t','Wald','Wilson','Clopper-Pearson','sim. CI','Jeffreys'};
algs=8;
res = zeros(size(y,1),algs); % confidence interval width for the different approaches
up = zeros(size(y,1),algs); % upper confidence interval for the different approaches
low = zeros(size(y,1),algs); % lower confidence interval for the different approaches
z=norminv(1-ci.alpha/2,0,1);
z2=z.^2;
for i=1:size(y,1)        
    yi=y(i,:);
    n=length(yi);
    
    % CI for Mean: Bootstrap algorithm
    bootstrap = bootci(2000,{@mean,yi})';
    res(i,1)=diff(bootstrap);
    low(i,1)=bootstrap(1);    up(i,1)=bootstrap(2);
    
    % CI for Mean: Normal approximation for the MOS
    res(i,2) = z.*std(yi)/sqrt(n) *2;
    low(i,2)=mean(yi)-res(i,2)/2;    up(i,2)=mean(yi)+res(i,2)/2;
    
    % CI for Mean: Student-T approximation for the MOS
    res(i,3)= tinv(1-ci.alpha/2  ,n-1).*std(yi)/sqrt(n) *2;
    low(i,3)=mean(yi)-res(i,3)/2;    up(i,3)=mean(yi)+res(i,3)/2;
    
    % CI for Bino: Wald interval employing normal approximation
    mu=mean(yi)-1; m=4;
    p=mu/m;
    s = sqrt(p*(1-p)/n);
    p1 = p-z*s;p2=p+z*s;
    low(i,4)=p1*m+1; up(i,4)=p2*m+1;
    res(i,4)=(up(i,4)-low(i,4));    
    
    % CI for Bino: Wilson score interval with continuity correction
    low(i,5)=max(1, (2*p*n*m + z2 - (1+z*sqrt(z2-1./(n*m) + 4*n*m*p*(1-p)+(4*p-2)) )  )/(2*(n*m+z2)) *m+1 );
    up(i,5)=min(m+1, (2*p*n*m + z2 + (1+z*sqrt(z2-1./(n*m) + 4*n*m*p*(1-p)+(4*p-2)) )  )/(2*(n*m+z2)) *m+1 );
    res(i,5)=(up(i,5)-low(i,5));
    
    % CI for Bino: Clopper-Pearson (using beta distribution)
    numSuc=sum(yi-1);trials=n*m;
    low(i,6) = max(1, betainv(ci.alpha/2,numSuc,trials-numSuc+1) *m + 1);
    up(i,6) = min(m+1, betainv(1-ci.alpha/2,numSuc+1,trials-numSuc) *m + 1);
    res(i,6)=(up(i,6)-low(i,6));
    
    % simultaneous confidence intervalls
    di=1:5;
    xi = hist(yi,di);
    a = chi2inv(1-ci.alpha,length(di)-1);  % a = chi^2(1-alpha,d-1)
    a = chi2inv(1-ci.alpha/length(di),1); % bonferroni
    scim = sum((1:5).*xi/n);
    low(i,7)=scim-sqrt(a/n.*(sum((1:5).^2.*xi/n)- (sum((1:5).*xi/n)).^2 ));
    up(i,7)=scim+sqrt(a/n.*(sum((1:5).^2.*xi/n)- (sum((1:5).*xi/n)).^2 ));
    res(i,7)=(up(i,7)-low(i,7));
    
    % Jeffreys interval
    if numSuc==0
        low(i,8)=1;
    else
        low(i,8)=betainv(ci.alpha/2,numSuc+1/2,trials-numSuc+1/2)*m+1;
    end
    if numSuc==trials
        up(i,8)=5;
    else
        up(i,8)=betainv(1-ci.alpha/2,numSuc+1/2,trials-numSuc+1/2)*m+1;
    end
    res(i,8)=(up(i,8)-low(i,8));
    
end

%% Return results
% str={'bootstrap','normal','student-t','Wald','Wilson','Clopper-Pearson','sim. CI'};
calc.str=str;
calc.shortstr={'boot.','norm.','stud.','Wald','Wils.','C-P','sim.CI','Jeff.'};
calc.mos=mos;
calc.ciUpper=up;
calc.ciLower=low;
calc.ciWidth=res;
calc.alpha=ci.alpha;

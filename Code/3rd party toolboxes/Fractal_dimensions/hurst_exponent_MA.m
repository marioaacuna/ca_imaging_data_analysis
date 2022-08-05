function h = hurst_exponent_MA(data0)

d = 8;

if max(size(d)) == 1
    % For scalar d set dmin=d and find the 'optimal' vector d
    dmin = d;
    % Find such a natural number OptN that possesses the largest number of 
    % divisors among all natural numbers in the interval [0.99*N,N] 
    N = length(data0); 
    N0 = floor(0.99*N);
    dv = zeros(N-N0+1,1);
    for i = N0:N
        dv(i-N0+1) = length(divisors(i,dmin));
    end
    OptN = N0 + max(find(max(dv) == dv)) - 1;
    % Use the first OptN values of x for further analysis
    data0 = data0(1:OptN);
    % Find the divisors of x
    d = divisors(OptN,dmin);
else
    OptN = length(x);
end

N = length(d);
yvals = zeros(N,1);

% Calculate empirical R/S
for i=1:N
   yvals(i) = RScalc(data0,d(i));
end
xvals = d;

% convert the numbers to log
logx = log10(xvals);
logy = log10(yvals);

P = polyfit(logx, logy,1);
h = P (1); % Hurst Exponent
return
% plotting
figure
plot(logx, logy, 'ro-')
hold on
f1 = polyval(P, logx);
plot(logx, f1)
xlabel('log10(scale)')
ylabel('log10(R/S)')
text (1,1,['H = ', num2str(h)])




function d = divisors(n,n0)
% Find all divisors of the natural number N greater or equal to N0
i = n0:floor(n/2);
d = find((n./i) == floor(n./i))' + n0 - 1;

function rs = RScalc(Z,n)
% Calculate (R/S)_n for given n
m = length(Z)/n;
Y = reshape(Z,n,m);
E = mean(Y);
S = std(Y);
for i=1:m
    Y(:,i) = Y(:,i) - E(i);
end
Y = cumsum(Y, 1);
% Find the ranges of cummulative series
MM = max(Y) - min(Y);
% Rescale the ranges by the standard deviations
CS = MM./S;
rs = mean(CS);
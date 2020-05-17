%% Part I

% part d - estimating distribution of X
%where X = 0 wp 1/3
%      X = 1 wp 1/2
%      X = 2 wp 1/6
Nvec=[100,1000,10e5,10e7];
 figure;
 X = [0:2];
 Uni = [1/3,1/2,1/6];
for w = (1:4)
    numberOfPoints=Nvec(w);
        Fx = rand(1,numberOfPoints); %generating iid sequence of uniform RV's
        Map = zeros(1,length(numberOfPoints));  
    for p=1:length(Fx) % performing transforation to map distrbution to same dist as given RV
        if(Fx(p)>.3333) && (Fx(p)<.8333)
            Map(p)=1;
        elseif(Fx(p)>.8333)
            Map(p) =2;
        end
    end

        bins = [0:.1:2];
        [yvalues,xvalues]=hist(Map,bins);
        yvalues = yvalues/(numberOfPoints); % we normalize to produce the probability densities
        subplot(2,2,w);
        hold on;
        a = stem(X,Uni,'r'); %using stem function to approximate deltas spikes for our PDF
        bar(xvalues,yvalues,'c');
        uistack(a,'top')
        xlabel('x');
        str = sprintf('Number of Points = %d',Nvec(w) );
        title(str);
end

% Part e expontential
figure;
lambda = [1,3,10];  
N = 1000000;
bins = [0:0.01:5];
X = rand(N,1);
for w=1:length(lambda)
    E = -log(X)/lambda(w); %taking the transformation Y = -log(X)/lambda takes a unifrom dist and
    [yvalues,xvalues]=hist(E, bins);% transforms it to fit our exponential dist (proof in project)
    subplot(3,1,w)
    yvalues = yvalues/(N*0.01); %normalizing
    bar(xvalues,yvalues);
    z = [0:0.01:5];
    realPdf = inline('lambda*exp(-lambda.*z)','z','lambda'); %compare to theoretical values
    hold on;
    plot(z,realPdf(z,lambda(w)),'red');
    xlabel('x');
    str = sprintf('Lambda= %d', lambda(w));
    title(str);
end

%Part f z1
N = 1000000;
U = rand(N,1);
Z1= zeros(1,length(U));
bins = [-5:0.01:5];
for w=1:length(U)
    if U(w)> 1/2
        Z1(w)=-log(2*(1-U(w))); %applying the given transformation for each range
    else
        Z1(w)=log(2*U(w));
    end
end
    [yvalues,xvalues]=hist(Z1, bins);
    yvalues = yvalues/(N*0.01);
    figure;
    bar(xvalues,yvalues);
    z1 = [-5:0.01:0];
    z2 = [0:0.01:5];
    realPdf1 = inline('(exp(z1))/2','z1'); %comparing to real pdf, see written section for derivaiton.
    realPdf2 = inline('(exp(-z2))/2','z2');
    hold on;
    plot(z1,realPdf1(z1),'red');
    plot(z2,realPdf2(z2),'red');
    xlabel('Value');
    str = sprintf('Z1 Pdf');
    title(str);

figure;
N = 1000000;
U = rand(N,1);
Z3=sqrt(-2*log(1-U)); %applying transformation
bins = [0:0.01:5];
[yvalues,xvalues]=hist(Z3, bins); %monte carlo
yvalues = yvalues/(N*0.01);
bar(xvalues,yvalues);
z = [0:0.01:5];
realPdf = inline('z.*exp((-z.^2)/2)','z');
hold on;
plot(z,realPdf(z),'red'); %plotting the theoretically obtained pdf to compare, see written report for derivation
xlabel('Value');
str = sprintf('Z2 Pdf');
title(str);

%% Part II

%part a


figure;
index=[1 2 3 6 7 200];
bins = [-2:0.01:2];

for i =1:length(index)
    Ck=zeros(N,1);
    for j=1:index(i)
        X = rand(N,1)-0.5; % generate two uniform distributions
        Y= rand(N,1)-0.5;
        Z=X+Y; % use of addition convolution property seen in class to obtain triangle pdf
        Ck=Z+Ck; %calculate partial sum for a given index
       
    end
    Ck=(1/sqrt(index(i)))*Ck; %normalize
    subplot(3,2,i);
    [yvalues,xvalues]=hist(Ck, bins);
    yvalues = yvalues/(N*0.01);
    bar(xvalues,yvalues); %graphs become more gaussian for higher index
    xlabel('x');
    str = sprintf('Partial sum for C%d',index(i));
    title(str);
end

%part b
% code is similar to part a, without the use of the convolution property.
figure;
index=[1 2 3 6 7 200];
bins = [-2:0.01:2];
for i =1:length(index)
    Ck=zeros(N,1);
    for j=1:index(i)
        X = rand(N,1)-0.5;
        Ck=X+Ck;
       
    end
    Ck=(1/sqrt(index(i)))*Ck;
    subplot(3,2,i);
    [yvalues,xvalues]=hist(Ck, bins);
    yvalues = yvalues/(N*0.01);
    bar(xvalues,yvalues);
    xlabel('x');
    str = sprintf('Partial sum for C%d',index(i));
    title(str);
end 
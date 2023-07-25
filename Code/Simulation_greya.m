% greya_e.m
% last modified 9/9/22
% host-greya moth model
% latest model version with multiple preference rounds
% output is a figure and a csv file

clear;

%er = .2 %rand(1) %.14
%s = .1 %rand(1)%.3 %.1
%r = .1 %.5*rand(1)%.5 %0.25
%b = .2 %rand(1) %.3
%T=501;
% define x, insect genotype freq, & q, plant genotype freq,
% and D, linkage disequilibrium, & ICs
% random ICs
%x0 = [rand(1),rand(1),rand(1),rand(1)]; y0 = rand(1);
%x0 = [.01,.05,.05,1]; y0 = .7;

% regimes
% insects do not switch - no switch1 regime
%er = .2
%s = .1
%r = .1
%b = .2
%T=501; x0 = [1,.05,.05,.01];
%y0 = .7; %(plants nor insects switch)
%y0 = .3; %(plants switch)


%no switch2 regime
%er = .2
%s = .1
%r = .1
%b = .2
%T=501;
%x0 = [.01,.05,.05,1]; y0 = .7;

%switch1 regime
%er = .2; s = .1; r = .1; b = .2; T=201; x0 = [1,.05,.05,.01]; y0 = .2;

%switch2 regime
%er = .2; s = .1; r = .1; b = .2; T=501; x0 = [.01,.05,.05,1]; y0 = .8;

%switch3 regime
er = .2; s = .1; r = .1; b = .2; T=101; x0 = [1,.05,.05,.01]; y0 = .2;

%special cases
%er = 0; s = .1; r = .1; b = .2; T=501; x0 = [rand(1),rand(1),rand(1),rand(1)]; y0 = rand(1);
%er = 1; s = .1; r = .1; b = .2; T=501; x0 = [.5,.5,.5,.3]; y0 = .7;



x(1,:) = x0/sum(x0);
y(1,:) = [y0, 1-y0];
D(1,1) = x0(1)*x0(4)-x0(2)*x0(3);

% define matrices
% plant offspring prob given parent is m and n = P(m,n->l) = S(l,m,n)
S = zeros(2,2,2);
S(1,:,:) = [1 .5; .5 0];
S(2,:,:) = [0 .5; .5 1];
% insect recombination, P(i,j->k) = R(k,i,j)
R = zeros(4,4,4);
R(1,:,:) = [1 .5 .5 (1-r)/2; .5 0 r/2 0; .5 r/2 0 0; (1-r)/2 0 0 0];
R(2,:,:) = [0 .5 0 r/2; .5 1 (1-r)/2 .5; 0 (1-r)/2 0 0; r/2 .5 0 0];
R(3,:,:) = [0 0 .5 r/2; 0 0 (1-r)/2 0; .5 (1-r)/2 1 .5; r/2 0 .5 0];
R(4,:,:) = [0 0 0 (1-r)/2; 0 0 r/2 .5; 0 r/2 0 .5; (1-r)/2 .5 .5 1];
% local adaptation, S(n on k) = sel(n,k)
W = [1 1-s; 1 1-s; 1-s 1; 1-s 1];
% define the preference matrix P(i->k) = pref(i,k) each gen
pref = [(1+er)/2 (1-er)/2; (1-er)/2 (1+er)/2; (1+er)/2 (1-er)/2; (1-er)/2 (1+er)/2];

for a = 1:(T-1)
    %_____________________________________________________________________
    %normalize each pref matrix so that sumP(k) = sumoveri(pref(i,k)*x(a,i))

    %Conditional probability that given plant m, insect i will visit it is
    %P(i,m)
    for m = 1:2
        if (x(a,:)*pref(:,m))==0
            P(:,m) = zeros(4,1);
        else
            P(:,m) = pref(:,m)/(x(a,:)*pref(:,m));
        end
        xplant(:,m) = P(:,m).*transpose(x(a,:)); %freq of insect i on plant m
    end

    %____________________________________________________________________
    % Let cm be the overall proportion of adult insect found on plant type
    % m.  This will help us later on to calculate larval frequencies and
    % then to normalize after larval selection

    c = zeros(1,2);
    for m=1:2
        c(m) = y(a,m)*(x(a,:)*pref(:,m))/(x(a,:)*pref*(y(a,:))');
    end

    %____________________________________________________________________
    % define the mating probabilities/encounters for insects

    % insects, we assume random mating on host
    % M(m,i,j) = mating freq of insects i and j on plant m
    M = zeros(2,4,4);
    for m = 1:2
        M(m,:,:) = xplant(:,m)*transpose(xplant(:,m));
    end

    %___________________________________________________________________
    % calculation of larval proportions on each plant so that we can
    % calculate population after selection on each plant

    % E(k,i) is the proportion of eggs with genotype k, carried by mother
    % of genotype i, so E is 4x4
    E = zeros(4,4);
    for m=1:2
        for k=1:4
            for i = 1:4
                E(k,i) = E(k,i)+ c(m)*sum(M(m,i,:).*R(k,i,:));
            end
        end
    end

    % Now we calculate the proportions of eggs laid on each plant, e(k,n)
    % then the amount after selection es(k,n) & normalized
    e = zeros(4,2);
    es = zeros(4,2);
    for n = 1:2
        e(:,n) = (E*pref(:,n))/sum(E*pref(:,n));
        % insects after selection on plant n surviving to adulthood
        es(:,n) = (e(:,n).*W(:,n))/(transpose(e(:,n))*W(:,n));
    end

    % final adult insect frequecies in population
    x(a+1,:) = c*transpose(es);

    %_______________________________________________________________
    % Calculation of plant offspring - there is no selection, so
    % normalization is not necessary

    %Conditional probability that given plant m, insect i will visit it is
    %P(i,m)
    Q = zeros(4,2);
    for i = 1:4
        if pref(i,:)*transpose(y(a,:))== 0
            Q(i,:) = zeros(1,2);
        else
            Q(i,:) = pref(i,:)/(pref(i,:)*transpose(y(a,:))); %=pi/sumoverm of pi*x
        end
    end

    % The new offspring sprouting in the next generation
    yn = zeros(1,2);
    for i = 1:4
        for m = 1:2
            for n = 1:2
                yn = yn + x(a,i)*Q(i,m)*y(a,m)*(S(:,m,n))'*Q(i,n)*y(a,n);
            end
        end
    end


    % now incorporate the past generation to find the new generation plant
    % composition
    y(a+1,:) = (1-b)*y(a,:) + b*yn;


    %____________________________________________________________________
    % Just to check behavior in the middle of code
    %if a == 3
    %    P
    %    Q
    %    c
    %    x(1:a,:)
    %    y(1:a,:)
    %    M
    %end

    %D(a+1,1) = x(a+1,1)*x(a+1,4)-x(a+1,2)*x(a+1,3);
end


t = (0:1:(T-1));
p1 = x(:,1) + x(:,2);
p2 = x(:,1) + x(:,3);
finalplantfreq = y(T,:);

%plotting
plot(t, p1, 'r', 'LineWidth', 1, t, p2, 'b','LineWidth', 1, t, y(:,1),'k','LineWidth', 1)
%plot(t, p1, 'r', t, p2, 'b', t, y(:,1),'k', t, D, 'm','LineWidth',2)
legend('A at Local Adaptation Locus','B at Preference Locus','C Plant Type')
legend('boxoff')
%title('Figure 3a: Change in allele frequencies over time')
xlabel('Generation')
ylabel('Allele Frequency in Total Population')
axis([0 T -.05 1.05])
set(gca,'YTick',0:0.2:1)

%export csv
greya = [t.',p1,p2,y(:,1)];
csvwrite ("greya_switch.csv",greya)

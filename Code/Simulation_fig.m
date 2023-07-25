%fig2pop_multfigures.m
% last modified 1/7/09
% this is to run simulations of the plant-insect model
% corrected version which incorporate's SG feedback
% with matrix approach from scratch
% uses preferences as (1+e)/2 , (1-e)/2
% eliminates a plant visit, so
% that it mates and lays its eggs on the same plant.

clear;

% define x, insect genotype freq, & q, plant genotype freq,
% and D, linkage disequilibrium, & ICs
s = .6;
r = .2;
b = .02;
T=2001; %51;
nrun = 0;
plotcsv = []; %to store the vectors for export
%er=.5;
erarray = [0,.25,.5,.75,1];
% random ICs/param
for ics = 1:3
    if ics == 1
        xplant0 = [0.1,.1;0.5,0.1;0.1,0.5;1,1]; y0 = .1; z0 = .7;
    elseif ics == 2
        %host switch
        xplant0 = [0.853,0.6221;0.351,0.5132;0.4018,0.0760;0.2399,0.1233]; y0 = .1839; z0 = .24;
    elseif ics == 3
        %xplant0 = rand(4,2)
        %y0 = rand(1)
        %z0 = rand(1)
        xplant0 = [0.7060,0.0971;0.0318,0.8235;0.2769,0.6948;0.0462,0.3171];y0=0.9502;z0=0.0344;

    end
    for pstep = 0:4 % vary this
        er = erarray(pstep+1);
        nrun = nrun+1;
        %s = sarray(pstep+1);
        %r = rarray(pstep+1);
        %b = barray(pstep+1);
%random ICs/param
% x is the proportion of females flying around in the air
% xplant is the proportions of insect gentotypes within each plant, so
% that the proportions add to 1 on each plant
%er =rand(1); s = rand(1); r = .5*rand(1); b = rand(1); T=5001;
%xplant0 = [rand(1),rand(1);rand(1),rand(1);rand(1),rand(1);rand(1),rand(1)]; y0 = rand(1); z0 = rand(1);

xplant=zeros(4,2);
for m = 1:2
    xplant(:,m) = xplant0(:,m)/sum(xplant0(:,m));
end

y(1,:) = [y0, 1-y0];
z(1,:) = [z0, 1-z0];
x(1,:) = xplant*transpose(y(1,:));
D(1,1) = x(1,1)*x(1,4)-x(1,2)*x(1,3);

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
% local adaptation, W(n on k) = W(n,k)
W = [1 1-s; 1 1-s; 1-s 1; 1-s 1];
% define the preference matrix P(i->k) = pref(i,k) each gen
pref = [(1+er)/2 (1-er)/2; (1-er)/2 (1+er)/2; (1+er)/2 (1-er)/2; (1-er)/2 (1+er)/2];

for a = 1:(T-1)
    %____________________________________________________________________
    % define which subpopulation of plants is in male phase and which are
    % in female phase
    if mod(a,2) == 1
        malefig = y(a,:); %frequency of the interfloral and male fig
        femalefig = z(a,:); %frequency of the female fig
    else
        malefig = z(a,:);
        femalefig = y(a,:);
    end

    %INSECT LIFE CYCLE

    %_____________________________________________________________________
    % define the mating probabilities/encounters for insects
    % we assume random mating on host
    % M(m,i,j) = mating freq of insects i and j on plant m
    M = zeros(2,4,4);
    for m = 1:2
        M(m,:,:) = xplant(:,m)*transpose(xplant(:,m));
    end

    %_____________________________________________________________________
    % frequency of females of genotype i have eggs k in the entire
    % population emerging - that is because each fig type contributes
    % females in propotion to its frequency
    f = zeros(4,4); %f(i,k) = frequency of female i in air carrying egg k
    F = zeros (1,4); % F(i) = frequency of female i in air
    for i = 1:4
        for k = 1:4
            for m = 1:2
                for j = 1:4
                    f(i,k) = f(i,k) + malefig(m)*M(m,i,j)*R(k,i,j);
                end
            end
            F(i) = F(i)+f(i,k);
        end
    end

    % get info for what the frequency is of females i from plant m
    Fp = zeros (4,2); %Fp(i,m) = freq of female i emerging from m
    for i = 1:4
        for m = 1:2
            for k = 1:4
                for j = 1:4
                    Fp(i,m) = Fp(i,m) + malefig(m)*M(m,i,j)*R(i,j,k);
                end
            end
        end
    end

    %_____________________________________________________________________
    %Conditional probability that given plant n, insect i will visit it is
    %P(i,n)
    P = ones(4,2);
    for n = 1:2
        % if (x(a,:)*pref(:,n))==0
        %    P(:,n) = zeros(4,1);
        %else
             P(:,n) = pref(:,n)/(F*pref(:,n));
        %end
    end

    %___________________________________________________________________
    % calculation of larval proportions on each plant so that we can
    % calculate population after selection on each plant

    % e(i,n) is the proportion of eggs with genotype i on plant n, so e is
    % 4x2

    e = zeros(4,2); %freq of eggs
    for n=1:2
        for k = 1:4
            e(k,n) = e(k,n)+ transpose(f(:,k))*P(:,n);
        end
    end

    %________________________________________________________________
    % SELECTION
    % Now we calculate the proportions of eggs laid on each plant, e(k,n)
    % then the amount after selection xplant(k,n) & normalized

    %es = zeros(4,2);
    for n = 1:2
        % insects after selection on plant n surviving to adulthood
        xplant(:,n) = (e(:,n).*W(:,n))/(transpose(e(:,n))*W(:,n));
    end

    %____________________________________________________________
    % xplant give the number on each plant, but all together, the total
    % insect frequencies are:
    x(a+1,:) = xplant*transpose(femalefig);

  %PLANT LIFE CYCLE

    %_______________________________________________________________
    % Calculation of plant offspring - there is no selection, so
    % normalization is not necessary

    % Q = Conditional probability that given insect i, it will visit plant n
    Q = zeros(4,2);
    for i = 1:4
        if pref(i,:)*transpose(femalefig)== 0
            Q(i,:) = zeros(1,2);
        else
            Q(i,:) = pref(i,:)/(pref(i,:)*transpose(femalefig)); %=pi/sumoverm of pi*x
        end
    end

    % The new offspring sprouting in the next generation
    fign = zeros(1,2);
    for i = 1:4
        for m = 1:2
            for n = 1:2
                fign = fign + Fp(i,m)*(S(:,m,n))'*Q(i,n)*femalefig(n);
            end
        end
    end
    fign = fign/sum(fign);

    % now incorporate the past generation to find the new generation plant
    % composition

    if mod(a,2) == 1
        z(a+1,:) = femalefig; %frequency of the interfloral and male fig
        y(a+1,:) = (1-b)*malefig + b*fign; %frequency of the female fig
    else
        y(a+1,:) = femalefig;
        z(a+1,:) = (1-b)*malefig + b*fign;
    end

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

end

t = (0:1:(T-1));
p1 = x(:,1) + x(:,2);
p2 = x(:,1) + x(:,3);
plotcsv=vertcat(plotcsv, [t',p1,p2,y(:,1),z(:,1),repelem(ics,T,1),repelem(er,T,1)]);

finalplantfreq = [y(T,:);z(T,:)];
subplot(6,3,3*pstep+ics)
plot(t, p1, 'r', t, p2, 'b', t, y(:,1),'k', t,z(:,1),'g','LineWidth',2)
axis([0 T -.05 1.05])
if pstep < 4
    set(gca,'xtick',[])
end
if ics > 1
    set(gca,'yTick',[])
else
    set(gca,'YTick',0:0.5:1)
end
    end
end
xlabel('Pollinator Generations')
ylabel('Allele Frequency in Total Population')
legend('A at Local Adaptation Locus','B at Preference Locus','C Plant Type Odd Generations', 'C Plant Type Even Generations')

csvwrite('multifig_er.csv',plotcsv)

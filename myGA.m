function [fxbest,xbest] = myGA(funchandle, xrange)

xrange = abs(xrange);
lb = (-1.*xrange)./2;
ub = (xrange./2);
n = 2000;

samples = (ub-lb).*rand(n,max(size(xrange))) + lb;

fit_samples = funchandle(samples);
fit_samples = max(fit_samples) - fit_samples;

fnorm = fit_samples / sum(fit_samples);
cum_fnorm = cumsum(fnorm);

n_parents = 1000;
repro = linspace(0,1-(1/n_parents),n_parents)';
repro = repro + rand(1);
repro = [repro(repro<1);repro(repro>=1)-1];
repro = sort(repro,'ascend');
repro = repmat(repro,1,n);
repro = repro';
cum_fnorm = repmat(cum_fnorm,1,n_parents);
h = (cum_fnorm-repro);
h(h<0) = nan;
[~,idx] = min(h);
init = samples(idx,:);

iter = 1;

hh=1;
count = 50;
while  hh > 0.0001 || (iter < 500 && count > 1)
    count = count - 1;
    fit_init = funchandle(init);
    fit_init = max(fit_init) - fit_init;
    fit_init_norm = fit_init / sum(fit_init);
    cum_fit_init_norm = cumsum(fit_init_norm);
    cum_fit_init_norm = repmat(cum_fit_init_norm,1,n_parents);
    
    repro = linspace(0,1-(1/n_parents),n_parents)';
    repro = repro + rand(1);
    repro = [repro(repro<1);repro(repro>=1)-1];
    repro = repmat(repro,1,n_parents);
    repro = repro';

    h = (cum_fit_init_norm-repro);
    h(h<0) = nan;
    [~,idx] = min(h);
    
    init = init(idx,:);
    
    if iter>1
        tremp = init;
        idxx = randi(n_parents,n_parents - floor(elit_percent*n_parents),1);
        tremp2 = tremp(idxx,:);
        init = [tremp2;elit_parents];
    end

    initbin = numtobin(init);

    initbin2 = join(initbin,"");

    pair_idx = randperm(n_parents,n_parents)';

    paired_parents = initbin2(pair_idx,:);

    crossed_parents = paired_parents;
    len = size(char(crossed_parents(1,1)),2);
    for i = 1:n_parents/2
        probnum = rand(1);
        if probnum >= 0.25
            mutpoint1 = randi([1,len-1],1);
            mutpoint2 = randi([mutpoint1+1,len],1);
            string1 = char(paired_parents(2*i-1,1));
            temp1 = string1(mutpoint1:mutpoint2);
            string2 = char(paired_parents(2*i,1));
            temp2 = string2(mutpoint1:mutpoint2);
            string1(mutpoint1:mutpoint2) = temp2;
            string2(mutpoint1:mutpoint2) = temp1;            
            crossed_parents(2*i-1,1) = cellstr(string1);
            crossed_parents(2*i,1) = cellstr(string2);
        else
            crossed_parents(2*i-1,1) = paired_parents(2*i-1,1);
            crossed_parents(2*i,1) = paired_parents(2*i,1);
        end
    end

    mutated_parents=crossed_parents;

    tot_size = n_parents*len;
    prob_mut = 0.02;
    n_mutations = prob_mut * tot_size;
    nos = randi(tot_size-1,n_mutations,1);
    rno_cno = [floor(nos/len)+1 rem(nos,len)+1];
    mutated_parents = char(mutated_parents);
    for i = 1:n_mutations
        mutated_parents(rno_cno(i,1),rno_cno(i,2))=num2str(abs(str2num(mutated_parents(rno_cno(i,1),rno_cno(i,2)))-1));
    end

    ff = mutated_parents;
    ind_len = len/size(xrange,2);
    for i = 1:size(xrange,2)      
        rr(:,i) = cellstr(ff(:,ind_len*(i-1)+1:ind_len*i));
    end
    mutated_parents = rr;

    numeric_parents = bintonum(mutated_parents);   

    fit_mut_parents = funchandle(numeric_parents);
    fit_mut_parents = max(fit_mut_parents) - fit_mut_parents;
    fit_norm_mut_parents = fit_mut_parents/sum(fit_mut_parents);
    [~,idx] = max(fit_mut_parents);
    fit_best(iter,1) = funchandle(numeric_parents(idx,:));
    if iter>5
        h1 = abs(fit_best(iter,1) - fit_best(iter-5,1));
        h2 = abs(fit_best(iter,1) - fit_best(iter-4,1));
        h3 = abs(fit_best(iter,1) - fit_best(iter-3,1));
        h4 = abs(fit_best(iter,1) - fit_best(iter-2,1));
        h5 = abs(fit_best(iter,1) - fit_best(iter-1,1));
        hh = max([h1,h2,h3,h4,h5]);
    end
        
    iter = iter+1;

    elit_percent = 0.3;
    [sorted_parents,sort_idx] = sort(fit_norm_mut_parents,'descend');
    elit_parents = numeric_parents(sort_idx(1:floor(elit_percent*n_parents)),:);
    
end

[~,idx] = max(fit_mut_parents);
fxbest = funchandle(numeric_parents(idx,:));
format long
xbest = numeric_parents(idx,:);


function[gray] = numtobin(num)

% Dynamic Range
maxval = max(max(abs(num)));
bits = nextpow2(floor(maxval));
if bits<2
    bits = 2;
end
% Vectorizing the Function
[m,n] = size(num);

sign = num<0;
sign = reshape(sign,m*n,1);
ipp = abs(num); 
int = floor(ipp);
float = 100*(ipp - int);

intbin = dec2bin(int,bits);
floatbin = dec2bin(float,7);

gray = cellstr([num2str(sign) bintogray(intbin) bintogray(floatbin)]);
gray = reshape(gray,m,n);

function g = bintogray(b)
g(:,1) = b(:,1);
for i = 2 : size(b,2)
    x = xor(str2num(b(:,i-1)),str2num(b(:,i)));
    g(:,i) = num2str(x);
end

function[num] = bintonum(gray)

[m,n] = size(gray);

gray = char(reshape(gray,m*n,1));
[~,bits] = size(gray);
bits = bits-8;

intbin = graytobin(gray(:,2:bits+1));
floatbin = graytobin(gray(:,bits+2:end));
sign = str2num(gray(:,1));
sign(sign == 1) = -1;
sign(sign == 0) = 1;

intnum = bin2dec(intbin);
floatnum = bin2dec(floatbin);
num = sign.*(intnum + floatnum/100);
num = reshape(num,m,n);

function b = graytobin(g)
b(:,1) = g(:,1);
for i = 2 : size(g,2)
    x = xor(str2num(b(:,i-1)), str2num(g(:,i)));
    b(:,i) = num2str(x);
end

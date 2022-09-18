clc
clear all
close all

% Opponent strategy
opstrat = strcat(dec2bin(randi(65535),16),num2str(randi([0 1])),num2str(dec2bin(randi(15),4)));

% Player strategy
n_parents = 100;
a = dec2bin(randi(65535,n_parents,1),16);
b = num2str(randi([0 1],n_parents,1));
c = num2str(dec2bin(randi(15,n_parents,1),4));
pstrat(:,1) = cellstr([a,b,c]);
             
iter = 1;
while iter < 100
    n_crimes = 10;

    sent = fitness(n_parents,pstrat,opstrat,n_crimes);

    fit = n_crimes*5 - sent;
    fit_norm = fit/sum(fit);
    fit_cum = cumsum(fit_norm);

    fit_cum = repmat(fit_cum,1,n_parents);

    % Selection
    repro = linspace(0,1-(1/n_parents),n_parents)';
    repro = repro + rand(1);
    repro = [repro(repro<1);repro(repro>=1)-1];
    repro = repmat(repro,1,n_parents);
    repro = repro';

    h = (fit_cum-repro);
    h(h<0) = nan;
    [~,idx] = min(h);

    init = pstrat(idx,:);
    
    if iter>1
        tremp = init;
        idxx = randi(n_parents,n_parents - floor(elit_percent*n_parents),1);
        tremp2 = tremp(idxx,:);
        elit_parents = cellstr(elit_parents);
        init = [tremp2;elit_parents];
    end
    
    pstrat = init;
    
    initbin2 = init;
    
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
    
    numeric_parents = mutated_parents;
    if iter > 1
        numeric_parents(randi(size(numeric_parents,1)),:) = best_parent;
    end

    fit_mut_parents = fitness(n_parents,numeric_parents,opstrat,n_crimes);
    fit_mut_parents = n_crimes*5 - fit_mut_parents;
    fit_norm_mut_parents = fit_mut_parents/sum(fit_mut_parents);
    [~,idx] = max(fit_mut_parents);
    best_parent = numeric_parents(idx,:);
    fit_best(iter,1) = n_crimes*5 - fit_mut_parents(idx,1);
    
        
    iter = iter+1;
    
    elit_percent = 0.3;
    [sorted_parents,sort_idx] = sort(fit_norm_mut_parents,'descend');
    elit_parents = numeric_parents(sort_idx(1:floor(elit_percent*n_parents)),:);
    
end

format short

disp("Player strategy should be")
disp(best_parent)

disp(" ")
disp("Robinrounds for 10 crimes")
disp("Always testifying opponent")
[psent,opsent] = robinround(best_parent,"000000000000000000000",10);
disp("Player sentence is:")
disp(psent)
disp("Opponent sentence is:")
disp(opsent)

disp(" ")
disp("Always silent opponent")
[psent,opsent] = robinround(best_parent,"111111111111111111111",10);
disp("Player sentence is:")
disp(psent)
disp("Opponent sentence is:")
disp(opsent)

disp(" ")
disp("Tit for tat opponent")
[psent,opsent] = robinround(best_parent,"100110000111100001111",10);
disp("Player sentence is:")
disp(psent)
disp("Opponent sentence is:")
disp(opsent)

disp(" ")
disp("100 Random opponents")
for i = 1:100
    opstrat = strcat(dec2bin(randi(65535),16),num2str(randi([0 1])),num2str(dec2bin(randi(15),4)));
    [psent(i,1),opsent(i,1)] = robinround(best_parent,char(opstrat),10);
end
disp("Player sentence is (Average):")
disp(mean(psent))
disp("Opponent sentence is (Average):")
disp(mean(opsent))


function [sent] = fitness(n_parents,pstrat,opstrat,n_crimes)

decision   = [0 0 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
              0 0 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
              1 0 0 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
              1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];

for i = 1:n_parents
    temp_pstrat = char(pstrat(i,:));
    temp_ostrat = char(opstrat);
    count = 0;
    for j = 1:n_crimes
        if j == 1
            ps(j) = temp_pstrat(j);
            ops(j) = temp_ostrat(j);
                if str2num(ps(j)) == 1 && str2num(ops(j)) == 1
                    count = count + 1;
                elseif str2num(ps(j)) == 1 && str2num(ops(j)) == 0
                    count = count + 5;

                elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 1
                    count = count;

                elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 0
                    count = count + 3;
                end
        elseif j == 2
            h = [str2num(ps(j-1));str2num(ps(j-1));str2num(ops(j-1));
                str2num(ops(j-1))];
            [~, idx]=ismember((decision(:,1:4))',h','rows');
            [~,idxx] = max(idx); 
            ps(j) = (temp_pstrat(idxx+1));
            ops(j) = (temp_ostrat(idxx+1));
                if str2num(ps(j)) == 1 && str2num(ops(j)) == 1
                    count = count + 1;
                elseif str2num(ps(j)) == 1 && str2num(ops(j)) == 0
                    count = count + 5;

                elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 1
                    count = count;

                elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 0
                    count = count + 3;
                end
        else
            h = [str2num(ps(j-2));str2num(ps(j-1));str2num(ops(j-2));
                str2num(ops(j-1))];
            [~, idx]=ismember((decision(:,5:20))',h','rows');
            [~,idxx] = max(idx);
            ps(j) = (temp_pstrat(idxx+5));
            ops(j) = (temp_ostrat(idxx+5));
                if str2num(ps(j)) == 1 && str2num(ops(j)) == 1
                    count = count + 1;
                elseif str2num(ps(j)) == 1 && str2num(ops(j)) == 0
                    count = count + 5;

                elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 1
                    count = count;

                elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 0
                    count = count + 3;
                end
        end
    end
    sent(i,1) = count;
end
end

function [psent,opsent] = robinround(pstrat,opstrat,n_crimes)

temp_pstrat = char(pstrat);
temp_ostrat = char(opstrat);

decision   = [0 0 1 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1;
              0 0 1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1;
              1 0 0 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1;
              1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];

p_count = 0; 
op_count = 0;

for j = 1:n_crimes
    if j == 1
        ps(j) = temp_pstrat(j);
        ops(j) = temp_ostrat(j);
            if str2num(ps(j)) == 1 && str2num(ops(j)) == 1
                p_count = p_count + 1;
                op_count = op_count + 1;
                
            elseif str2num(ps(j)) == 1 && str2num(ops(j)) == 0
                p_count = p_count + 5;
                op_count = op_count;

            elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 1
                p_count = p_count;
                op_count = op_count + 5;

            elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 0
                p_count = p_count + 3;
                op_count = op_count + 3;
            end
    elseif j == 2
        h = [str2num(ps(j-1));str2num(ps(j-1));str2num(ops(j-1));
            str2num(ops(j-1))];
        [~, idx]=ismember((decision(:,1:4))',h','rows');
        [~,idxx]= max(idx);
        ps(j) = (temp_pstrat(idxx+1));
        ops(j) = (temp_ostrat(idxx+1));
            if str2num(ps(j)) == 1 && str2num(ops(j)) == 1
                p_count = p_count + 1;
                op_count = op_count + 1;
                
            elseif str2num(ps(j)) == 1 && str2num(ops(j)) == 0
                p_count = p_count + 5;
                op_count = op_count;

            elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 1
                p_count = p_count;
                op_count = op_count + 5;

            elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 0
                p_count = p_count + 3;
                op_count = op_count + 3;
            end

    else
        h = [str2num(ps(j-2));str2num(ps(j-1));str2num(ops(j-2));
            str2num(ops(j-1))];
        [~, idx]=ismember((decision(:,5:20))',h','rows');
        [~,idxx] = max(idx);
        ps(j) = (temp_pstrat(idxx+5));
        ops(j) = (temp_ostrat(idxx+5));
            if str2num(ps(j)) == 1 && str2num(ops(j)) == 1
                p_count = p_count + 1;
                op_count = op_count + 1;
                
            elseif str2num(ps(j)) == 1 && str2num(ops(j)) == 0
                p_count = p_count + 5;
                op_count = op_count;

            elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 1
                p_count = p_count;
                op_count = op_count + 5;

            elseif str2num(ps(j)) == 0 && str2num(ops(j)) == 0
                p_count = p_count + 3;
                op_count = op_count + 3;
            end

    end
end

psent = p_count;
opsent = op_count;

end
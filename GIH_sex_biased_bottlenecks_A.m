% MATLAB code simulates the effects of sex-biased bottlenecks
% plausible estimates of theta_w and pi statistics generated
% code written by Joe Lachance (jlachance6@biology.gatech.edu)
% last updated 9/26/17

% number of times the simulation is run from start
runs=20;

% number of times data is sampled at the end of each rung
number_of_sample_runs=50;

% defines maximum sizes population size and sets bottleneck sizes
max_males=10000;
max_females=10000;
GIH_male_bottleneck=300;
GIH_female_bottleneck=200;

% initializes data matrices
% each row is a different run
% each column is a different sample run
theta_A_matrix=zeros(runs,number_of_sample_runs);
pi_A_matrix=zeros(runs,number_of_sample_runs);

% number of generations to run simulation
max_generation=2400;
generation=1;

% defines growth rate per generation
population_growth_rate=zeros(1,max_generation);
for generation=1:2320
    population_growth_rate(generation)=0.005;
end
for generation=2320:1:2400
    population_growth_rate(generation)=0.1;
end

% defines population size vectors
number_of_males=zeros(1,max_generation);
number_of_females=zeros(1,max_generation);
number_of_males(1)=max_males;
number_of_females(1)=max_females;
number_of_males(2)=GIH_male_bottleneck;
number_of_females(2)=GIH_female_bottleneck;
for generation=3:1:max_generation
    number_of_males(generation)=round(max_males/(1+((max_males-((GIH_male_bottleneck+GIH_female_bottleneck)/2))/((GIH_male_bottleneck+GIH_female_bottleneck)/2))*exp(-population_growth_rate(generation)*(generation-2))));
    number_of_females(generation)=round(max_females/(1+((max_females-((GIH_male_bottleneck+GIH_female_bottleneck)/2))/((GIH_male_bottleneck+GIH_female_bottleneck)/2))*exp(-population_growth_rate(generation)*(generation-2))));
end

for i=1:1:runs
    % sets numbers of sites, mutation rate, and initial theta
    total_number_of_sites=10^6;
    mutation_rate=1.2*10^-8;
    theta_autosome_initial=0.001;
    autosomal_population_state=zeros(1,(2*max_males+2*max_females));
    autosomal_population_state_nextgen=zeros(1,(2*max_males+2*max_females));
     
    % generates initial autosome population state vector
    % equilibrium SFS assumed with probability proportional to 1/(allele freq)
    initial_probs=zeros(1,(2*max_males+2*max_females));
    initial_polymorphic_sites=round(total_number_of_sites*theta_autosome_initial*harmonic(2*max_males+2*max_females-1));
    for i_SFS=2:1:(2*max_males+2*max_females)
        initial_probs(i_SFS)=1/(i_SFS-1);
    end
    polymorphic_probability_sum=sum(initial_probs(2:2*max_males+2*max_females));
    initial_probs=initial_probs/polymorphic_probability_sum;
    for i_SFS=2:1:(2*max_males+2*max_females)
        autosomal_population_state(i_SFS)=bnldev(initial_polymorphic_sites,initial_probs(i_SFS));
    end
    % leftover copies go to monomorphic class
    autosomal_population_state(1)=total_number_of_sites-sum(autosomal_population_state(2:2*max_males+2*max_females));
        
    for generation=1:1:max_generation-1
        % mutation step (monomorphic to singleton)
        % note: mutation rates are assumed are the same for X and A
        % this means that empirical data needs to use corrected theta and pi
        new_mutations=round(autosomal_population_state(1)*(2*number_of_males(generation)+2*number_of_females(generation))*mutation_rate);
        autosomal_population_state(1)=autosomal_population_state(1)-new_mutations;
        autosomal_population_state(2)=autosomal_population_state(2)+new_mutations;
            
        % autosomal genetic drift step
        % speeds things up for monomorphic sites
        autosomal_population_state_nextgen(1)=autosomal_population_state(1);
        % goes down the population state vector
        for j=2:1:(2*max_males+2*max_females)
            % binomial sampling k times
            % where k is the number of sites with j-1 alleles
            for k=1:1:autosomal_population_state(j)
                copies_plus_one=bnldev((2*number_of_males(generation+1)+2*number_of_females(generation+1)),(j-1)/((2*number_of_males(generation)+2*number_of_females(generation))))+1;
                % fixation events become monomorphic sites
                if copies_plus_one==(2*number_of_males(generation+1)+2*number_of_females(generation+1))+1
                    autosomal_population_state_nextgen(1)=autosomal_population_state_nextgen(1)+1;
                else
                    % increments next generation's population state vector
                    autosomal_population_state_nextgen(copies_plus_one)=autosomal_population_state_nextgen(copies_plus_one)+1;
                end
            end
        end
       
        % shows which run and generation it is on
        [i, generation]
        
        % updates population state vectors
        autosomal_population_state=autosomal_population_state_nextgen;
        % resets population states for the next generation
        autosomal_population_state_nextgen=zeros(1,(2*max_males+2*max_females));
    end
    
    for sample_runs=1:1:number_of_sample_runs
        % generates sample of autosome data
        autosomal_sample_size=8;
        autosomal_sample_state=zeros(1,autosomal_sample_size);
        % speeds things up for monomorphic sites
        autosomal_sample_state(1)=autosomal_population_state(1);
        % goes down the population state vector
        for j=2:1:(2*max_males+2*max_females)
            for k=1:1:autosomal_population_state(j)
                % binomial sampling k times, where k is the number of sites
                % with j-1 allele
                copies_plus_one=bnldev(autosomal_sample_size,(j-1)/((2*number_of_males(generation+1)+2*number_of_females(generation+1))))+1;
                % deals with monomorphic sites
                if copies_plus_one==autosomal_sample_size+1
                    autosomal_sample_state(1)= autosomal_sample_state(1)+1;
                else
                    autosomal_sample_state(copies_plus_one)= autosomal_sample_state(copies_plus_one)+1;
                end
            end
        end
        
        theta_A_matrix(i,sample_runs)=(total_number_of_sites-autosomal_sample_state(1))/(harmonic(autosomal_sample_size-1)*total_number_of_sites);
        
        for sample_alleles=1:1:autosomal_sample_size
            p=(sample_alleles-1)/(autosomal_sample_size);
            q=(autosomal_sample_size+1-sample_alleles)/(autosomal_sample_size-1);
            pi_A_matrix(i,sample_runs)=pi_A_matrix(i,sample_runs)+autosomal_sample_state(sample_alleles)*2*p*q/total_number_of_sites;
        end
        
    end
 
end


% changes formatting of significant digits
format short G

% generates population genetic statistics for each run
theta_A_mean=mean(theta_A_matrix(:));
theta_A_95CI=prctile(theta_A_matrix(:), [2.5 50 97.5]);
pi_A_mean=mean(pi_A_matrix(:));
pi_A_95CI=prctile(pi_A_matrix(:), [2.5 50 97.5]);


% output file summarizes statistics
% row names: theta_w, pi
% column names: mean, 2.5%ile, median, 97.5%ile
clean_output_matrix=[theta_A_mean theta_A_95CI; 
    pi_A_mean pi_A_95CI];
output_matrix=[theta_A_mean theta_A_95CI; 
    pi_A_mean pi_A_95CI];
output_matrix=dataset({output_matrix 'mean','percentile_2pt5','median','percentile_97pt5'}, ...
              'obsnames', {'theta_w','pi'});
output_filename=['GIH_A_output_' int2str(2*GIH_male_bottleneck+2*GIH_female_bottleneck) 'autosomes_' int2str(runs) 'runs.txt']
save (output_filename,'clean_output_matrix', '-ascii','-tabs');
output_filename
output_matrix




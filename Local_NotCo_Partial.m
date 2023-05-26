%Ultimatum game simulations with local, independent mutations, and partial observability
%Akdeniz & van Veelen (2023)

%%Population parameters
s = 100; %number of agents
w = 0.1; %intensity of selection
u = 0.01; %mutation probability
p = 0.1; %probability of observing other's mao
scale = 5; %division scale of mutation
x = 1; %size of the pie to be split
t_max = (2^11)*(10^6); %maximum number of time periods
t_timesteps = 10^6; %determines how often a snapshot is taken
t_snapshots = t_max/t_timesteps; %number of times a snapshot will be taken
summary = zeros(18,1); %final average offers and acceptance thresholds

average_proposal_over_time_not_observed = zeros(t_max,1);
average_proposal_over_time_observed = zeros(t_max,1);
average_threshold_over_time = zeros(t_max,1);

%stopping rule
epsilon = 0.005; %when to stop
m_final_prop_not_observed = zeros(s,1);
m_final_prop_observed = zeros(s,1);
m_final_mao = zeros(s,1);

rng(40); %setting the seed

    counter1 = 0; %number of times mutation occurred at reproduction
    counter2 = 0; %number of times mutation occurred due to fixation
    m_snapshot1 = zeros(s,t_snapshots); %matrix to record snapshots of proposals in time
    m_snapshot2 = zeros(s,t_snapshots); %matrix to record snapshots of acceptance thresholds in time
    stop = 0;
    
    %%Random initialization of the population
    m = rand(s,3); %random initial str, first 2 proposer str. and the last mao
    m_new = zeros(s,3);
    m_total = m; %initialization of the total matrix
    
    average_proposal_over_time_not_observed(1,1) = sum(m(:,1))/s;
    average_proposal_over_time_observed(1,1) = sum(m(:,2))/s;
    average_threshold_over_time(1,1) = sum(m(:,3))/s;

    t = 2;
    %%Game play
    while t <=t_max %as long as maximum time step is not reached
        t;
        payoff = w*zeros(s,1); %payoff vector in actuall payoffs
        probs = rand(s,s); %matrix of observation probabilities for each interaction
        for i=1:s %game play for each player i
            for j=1:s %for each opponent j
                if j~=i
                    if (probs(i,j) > p) && (m(i,1) >= m(j,3)) %i is a proposer, doesn't observe and her offer is accepted
                        payoff(i,1) = payoff(i,1) + w*x - w*m(i,1); %payoff increase for i, 1ccounted for intensity
                        payoff(j,1) = payoff(j,1) + w*m(i,1); %payoff increase for j, 1ccounted for intensity
                    elseif (probs(i,j) <= p) && (m(i,2) >= m(j,3)) %i is a proposer, observes and she matches
                        payoff(i,1) = payoff(i,1) + w*x - w*m(j,3); %payoff increase for i, 1ccounted for intensity
                        payoff(j,1) = payoff(j,1) + w*m(j,3); %payoff increase for j, 1ccounted for intensity
                    end
                    if (probs(j,i) > p) && (m(j,1) >= m(i,3)) %i is a responder, isn't observed and she accepts
                        payoff(i,1) = payoff(i,1) + w*m(j,1); %payoff increase for i, 1ccounted for intensity
                        payoff(j,1) = payoff(j,1) + w*x - w*m(j,1); %payoff increase for j, 1ccounted for intensity
                    elseif (probs(j,i) <= p) && (m(j,2) >= m(i,3)) %i is a responder, is observed and j matches
                        payoff(i,1) = payoff(i,1) + w*m(i,3); %payoff increase for i, 1ccounted for intensity
                        payoff(j,1) = payoff(j,1) + w*x - w*m(i,3); %payoff increase for j, 1ccounted for intensity
                    end
                end
            end
        end

        payoff = payoff./(2*(s-1)); %averaging the payoffs over all interactions
        payoff_exp = exp(payoff); %exponential transformation of the payoff vector
        payoff_proportional = payoff_exp/sum(payoff_exp); %payoff vector with proportional payoffs

          %%%%Wright-Fisher version
        r = rand(s,1); %random vector to determine who are reproducing
        mut = rand(s,3); %random vector to determine if a mutation occurs at spot k
        for k=1:s %for each k, to determine who is reproducing to the location m(k,:)
        	index = 0; %to determine the location of the reproducing agent
            while r(k,1) > 0     %While loop to pick an agent with the given probabilities
                index  = index + 1;
                r(k,1) = r(k,1) - payoff_proportional(index);    
            end
            
            if mut(k,1) < u %mutation occurs at proposals w/o observing
                z1 = (rand - 0.5)/(scale); %random number to determine the new mutation
                counter1 = counter1 + 1; %number of times mutation occurred
                m_new(k,1) = m(index,1) + z1; 
                if mut(k,2) < u %mutation occurs at proposals w/ observing
                    z2 = (rand - 0.5)/(scale); %random number to determine the new mutation
                    counter1 = counter1 + 1; %number of times mutation occurred
                    m_new(k,2) = m(index,2) + z2;
                    if mut(k,3) < u %mutation occurs at mao's
                        z3 = (rand - 0.5)/(scale); %random number to determine the new mutation
                        counter1 = counter1 + 1; %number of times mutation occurred
                        m_new(k,3) = m(index,3) + z3;
                    else %no mutation occurs at mao's
                        m_new(k,3) = m(index,3);
                    end
                else %no mutation at proposals w/ observing
                    m_new(k,2) = m(index,2); 
                    if mut(k,3) < u %mutation occurs at mao's
                        z3 = (rand - 0.5)/(scale); %random number to determine the new mutation
                        counter1 = counter1 + 1; %number of times mutation occurred
                        m_new(k,3) = m(index,3) + z3;
                    else %no mutation occurs at mao's
                        m_new(k,3) = m(index,3);
                    end
                end
             else %no mutation occurs at proposals w/o observing
                m_new(k,1) = m(index,1); 
                if mut(k,2) < u %mutation occurs at proposals w/ observing
                    z2 = (rand - 0.5)/(scale); %random number to determine the new mutation
                    counter1 = counter1 + 1; %number of times mutation occurred
                    m_new(k,2) = m(index,2) + z2;
                    if mut(k,3) < u %mutation occurs at mao's
                        z3 = (rand - 0.5)/(scale); %random number to determine the new mutation
                        counter1 = counter1 + 1; %number of times mutation occurred
                        m_new(k,3) = m(index,3) + z3;
                    else %no mutation occurs at mao's
                        m_new(k,3) = m(index,3);
                    end
                else %no mutation at proposals w/ observing
                    m_new(k,2) = m(index,2); 
                    if mut(k,3) < u %mutation occurs at mao's
                        z3 = (rand - 0.5)/(scale); %random number to determine the new mutation
                        counter1 = counter1 + 1; %number of times mutation occurred
                        m_new(k,3) = m(index,3) + z3;
                    else %no mutation occurs at mao's
                        m_new(k,3) = m(index,3);
                    end
                end
            end
            
            if m_new(k,1) > 1
                m_new(k,1) = 1;
            elseif m_new(k,1) < 0
                m_new(k,1) = 0;
            end 
            if m_new(k,2) > 1
                m_new(k,2) = 1;
            elseif m_new(k,2) < 0
                m_new(k,2) = 0;
            end  
            if m_new(k,3) > 1
                m_new(k,3) = 1;
            elseif m_new(k,3) < 0
                m_new(k,3) = 0;
            end
            
        end
        m(:,:) = m_new(:,:); %updating to the old population state to the new population

        if rem(t,t_timesteps) == 0
            m_snapshot1(:,t/t_timesteps) = m(:,1); %recording proposals present at point t
            m_snapshot2(:,t/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
        end
        
        average_proposal_over_time_not_observed(t,1) = sum(m(:,1))/s;
        average_proposal_over_time_observed(t,1) = sum(m(:,2))/s;
        average_threshold_over_time(t,1) = sum(m(:,3))/s;

        m_total = m_total + m; %adding the current state of the population
        t_reached = t;
        
          for j=4:11
              if t == t_timesteps*(2^j)
                  if abs(sum(average_proposal_over_time_not_observed(1:(t/2), 1)) - sum(average_proposal_over_time_not_observed((t/2+1):t, 1)))/(t/2) <= epsilon
                      if abs(sum(average_threshold_over_time(1:(t/2), 1)) - sum(average_threshold_over_time((t/2+1):t, 1)))/(t/2) <= epsilon
                          if abs(sum(average_proposal_over_time_observed(1:(t/2), 1)) - sum(average_proposal_over_time_observed((t/2+1):t, 1)))/(t/2) <= epsilon
                            stop = j;
                          end
                      end
                  end
              end
          end

        if max(m(:,1)) == min(m(:,1)) && max(m(:,2)) == min(m(:,2)) && max(m(:,3)) == min(m(:,3))
            number_of_reproductions = geornd(u); %drawing a random number of reproductions until next mutation
            number_of_time_periods = fix(number_of_reproductions/(3*s)); %number of generations passed until next mutation
            mut_location = rem(number_of_reproductions,3*s); %location of the mutation
            if t + number_of_time_periods <= t_max
                for b=1:number_of_time_periods
                    average_proposal_over_time_not_observed(t_reached+b,1) = sum(m(:,1))/s;
                    average_proposal_over_time_observed(t_reached+b,1) = sum(m(:,2))/s;
                    average_threshold_over_time(t_reached+b,1) = sum(m(:,3))/s;
                    if rem(t_reached+b,t_timesteps) == 0
                        m_snapshot1(:,(t_reached+b)/t_timesteps) = m(:,1); %recording proposals present at point t
                        m_snapshot2(:,(t_reached+b)/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
                    end
                      for j=4:11
                          if t_reached+b == t_timesteps*(2^j)
                              if abs(sum(average_proposal_over_time_not_observed(1:((t_reached+b)/2), 1)) - sum(average_proposal_over_time_not_observed(((t_reached+b)/2+1):(t_reached+b), 1)))/((t_reached+b)/2) <= epsilon
                                  if abs(sum(average_threshold_over_time(1:((t_reached+b)/2), 1)) - sum(average_threshold_over_time(((t_reached+b)/2+1):(t_reached+b), 1)))/((t_reached+b)/2) <= epsilon
                                      if abs(sum(average_proposal_over_time_observed(1:((t_reached+b)/2), 1)) - sum(average_proposal_over_time_observed(((t_reached+b)/2+1):(t_reached+b), 1)))/((t_reached+b)/2) <= epsilon
                                        stop = j;
                                      end
                                  end
                              end
                          end
                      end
                end
                counter2 = counter2 + 1; %number of times mutation occurred
                new_z = (rand - 0.5)/(scale); %random number to determine the new mutation
                if mut_location == 0
                    m_total = m_total + m*(number_of_time_periods-1); %accounting for periods in between where nothing happens
                        m(s,3) = m(s,3) + new_z; %the mutant replaces the mao of the last individual
                        if m(s,3) > 1
                            m(s,3) = 1;
                        elseif m(s,3) < 0
                            m(s,3) = 0;
                        end
                    t = t + number_of_time_periods; %fastforwarding to the new time period in which mutation occurred
                    t_reached = t;
                else %mut_location > 0
                    m_total = m_total + m*(number_of_time_periods); %accounting for periods in between where nothing happens
                    if rem(mut_location,3) == 1 %mutation in proposals w/o observing
                        mut_individual = fix(mut_location/3)+1;
                        m(mut_individual,1) = m(mut_individual,1) + new_z; %the mutant replaces the strategy at location (k,1)
                        if m(mut_individual,1) > 1
                           m(mut_individual,1) = 1;
                        elseif m(mut_individual,1) < 0
                           m(mut_individual,1) = 0;
                        end
                    elseif rem(mut_location,3) == 2 %mutation in proposals w/ observing
                        mut_individual = fix(mut_location/3)+1;
                        m(mut_individual,2) = m(mut_individual,2) + new_z; %the mutant replaces the strategy at location (k,2)
                        if m(mut_individual,2) > 1
                           m(mut_individual,2) = 1;
                        elseif m(mut_individual,2) < 0
                           m(mut_individual,2) = 0;
                        end
                    else %mutation in mao's
                        mut_individual = fix(mut_location/3);
                        m(mut_individual,3) = m(mut_individual,3) + new_z; %the mutant replaces the strategy at location (k,3)
                        if m(mut_individual,3) > 1
                            m(mut_individual,3) = 1;
                        elseif m(mut_individual,3) < 0
                            m(mut_individual,3) = 0;
                        end
                    end
                    mut_index = 0; %check for additional mutations within the same generation
                    while mut_index == 0
                        number_of_reproductions2 = geornd(u); %drawing a random number of reproductions until next mutation
                        if mut_location + number_of_reproductions2 <= 3*s
                            counter2 = counter2 + 1; %number of times mutation occurred
                            mut_location = mut_location + number_of_reproductions2;
                            new_z = (rand - 0.5)/(scale); %random number to determine the new mutation
                            if rem(mut_location,3) == 1 %mutation in proposals w/o observing
                                mut_individual = fix(mut_location/3)+1;
                                m(mut_individual,1) = m(mut_individual,1) + new_z; %the mutant replaces the strategy at location (k,1)
                                if m(mut_individual,1) > 1
                                   m(mut_individual,1) = 1;
                                elseif m(mut_individual,1) < 0
                                   m(mut_individual,1) = 0;
                                end
                            elseif rem(mut_location,3) == 2 %mutation in proposals w/ observing
                                mut_individual = fix(mut_location/3)+1;
                                m(mut_individual,2) = m(mut_individual,2) + new_z; %the mutant replaces the strategy at location (k,2)
                                if m(mut_individual,2) > 1
                                   m(mut_individual,2) = 1;
                                elseif m(mut_individual,2) < 0
                                   m(mut_individual,2) = 0;
                                end
                            else %mutation in mao's
                                mut_individual = fix(mut_location/3);
                                m(mut_individual,3) = m(mut_individual,3) + new_z; %the mutant replaces the strategy at location (k,3)
                                if m(mut_individual,3) > 1
                                    m(mut_individual,3) = 1;
                                elseif m(mut_individual,3) < 0
                                    m(mut_individual,3) = 0;
                                end
                            end    
                        else
                            mut_index = 1; %next mutation goes out of the current generation
                        end
                    end
                    t = t + number_of_time_periods + 1; %fastforwarding to the new time period in which mutation occurred
                    t_reached = t;
                end

                m_total = m_total + m; %adding the current state of the population

                average_proposal_over_time_not_observed(t,1) = sum(m(:,1))/s;
                average_proposal_over_time_observed(t,1) = sum(m(:,2))/s;
                average_threshold_over_time(t,1) = sum(m(:,3))/s;
                if rem(t,t_timesteps) == 0
                    m_snapshot1(:,t/t_timesteps) = m(:,1); %recording proposals present at point t
                    m_snapshot2(:,t/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
                end
                
            else %mutation happens after the max time period is reached
                for b=1:(t_max-t_reached)
                    average_proposal_over_time_not_observed(t_reached+b,1) = sum(m(:,1))/s;
                    average_proposal_over_time_observed(t_reached+b,1) = sum(m(:,2))/s;
                    average_threshold_over_time(t_reached+b,1) = sum(m(:,3))/s;
                    if rem(t_reached+b,t_timesteps) == 0
                        m_snapshot1(:,(t_reached+b)/t_timesteps) = m(:,1); %recording proposals present at point t
                        m_snapshot2(:,(t_reached+b)/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
                    end
                end
                m_total = m_total + m*(t_max - t);
                t_reached = t_max;
                m_final_prop_not_observed(:,1) = m(:,1); %proposals at the last point
                m_final_prop_observed(:,1) = m(:,2); %proposals at the last point
                m_final_mao(:,1) = m(:,3); %mao's at the last point
                t = t_max + 1;
            end
        end

        t = t + 1; %moving to the next period
         if stop ~= 0
             t = t_max + 1;
         end
    end
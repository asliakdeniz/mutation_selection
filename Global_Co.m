%Ultimatum game simulations with global mutations
%Akdeniz & van Veelen (2023)

%%Population parameters
s = 10000; %number of agents
w = 1; %intensity of selection
u = 0.125; %mutation probability
x = 1; %size of the pie to be split
t_max = 2*(10^4); %2*(10^7); %maximum number of time periods
t_timesteps = 2*(10^3); %2*(10^6); %determines how often a snapshot is taken
t_snapshots = t_max/t_timesteps; %number of times a snapshot will be taken
summary = zeros(15,1); %final average offers and acceptance thresholds
dist = zeros(100,2); %matrix to keep track of how common a certain strategy is
dist_threshold = 10^2; %when to start keeping track of the strategy distributions
bins = 100;
bin_size = 1/bins;
check = 0;

average_proposal_over_time = zeros(t_max,1);
average_threshold_over_time = zeros(t_max,1);

%stopping rule
epsilon = 0.005; %when to stop
m_final_prop = zeros(s,1);
m_final_mao = zeros(s,1);

rng(40); %setting the seed

    counter1 = 0; %number of times mutation occurred at reproduction
    counter2 = 0; %number of times mutation occurred due to fixation
    snapshot1 = zeros(s,t_snapshots); %matrix to record snapshots of proposals in time
    snapshot2 = zeros(s,t_snapshots); %matrix to record snapshots of acceptance thresholds in time
    stop = 0;
    
    %%Random initialization of the population
    m = rand(s,2); %random draw of strategies
    m_new = zeros(s,2);
    m_total = m; %initialization of the total matrix
    
    average_proposal_over_time(1,1) = sum(m(:,1))/s;
    average_threshold_over_time(1,1) = sum(m(:,2))/s;

    t = 2;
    %%Game play
    while t <=t_max %as long as maximum time step is not reached
        t;
        payoff = w*zeros(s,1); %payoff vector in actuall payoffs
        payoff_responder = zeros(s,1); 
        
        for i=1:s %game play for each player i
            for j=1:s %for each opponent j
                if j~=i
                    if m(i,1) >= m(j,2) %i is a proposer, and her offer is accepted
                        payoff(i,1) = payoff(i,1) + w*x - w*m(i,1); %payoff increase for i, accounted for intensity
                        payoff(j,1) = payoff(j,1) + w*m(i,1); %payoff increase for j, accounted for intensity
                    end
                    if m(j,1) >= m(i,2) %i is a responder, and she accepts the offer
                        payoff(i,1) = payoff(i,1) + w*m(j,1); %payoff increase for i, accounted for intensity
                        payoff(j,1) = payoff(j,1) + w*x - w*m(j,1); %payoff increase for j, accounted for intensity
                    end
                end
                if m(j,1) >= m(i,2) %i is a responder, and she accepts the offer
                    payoff_responder(i,1) = payoff_responder(i,1) + w*m(j,1);
               end
            end
        end

        payoff = payoff./(2*(s-1)); %averaging the payoffs over all interactions
        payoff_responder = payoff_responder./(s-1);
        payoff_exp = exp(payoff); %exponential transformation of the payoff vector
        payoff_proportional = payoff_exp/sum(payoff_exp); %payoff vector with proportional payoffs
        
        for f=1:s %to check if payoff is always higher for lower than higher mao
            for g=1:s
                if m(f,2) < m(g,2) && payoff_responder(f,1) < payoff_responder(g,1)
                    check = check + 1;
                end
            end
        end

          %%%%Wright-Fisher version
        r = rand(s,1); %random vector to determine who are reproducing
        mut = rand(s,1); %random vector to determine if a mutation occurs at spot k
        for k=1:s %for each k, to determine who is reproducing to the location m(k,:)
        	index = 0; %to determine the location of the reproducing agent
            while r(k,1) > 0     %While loop to pick an agent with the given probabilities
                index  = index + 1;
                r(k,1) = r(k,1) - payoff_proportional(index);    
            end

            if mut(k,1) < u %mutation occurs
                z = rand(1,2); %random number to determine the new
                m_new(k,:) = z; %the mutant replaces the strategy at location (k,1)
                counter1 = counter1 + 1; %number of times mutation occurred
             else %mut(k,1) >= u, mutation does not occur
                 m_new(k,:) = m(index,:); %the reproducing agent replaces the strategy at location (k,1)
             end
        end
        m(:,:) = m_new(:,:); %updating to the old population state to the new population

         if rem(t,t_timesteps) == 0
             snapshot1(:,t/t_timesteps) = m(:,1); %recording proposals present at point t
             snapshot2(:,t/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
         end
        
        average_proposal_over_time(t,1) = sum(m(:,1))/s;
        average_threshold_over_time(t,1) = sum(m(:,2))/s;

        m_total = m_total + m; %adding the current state of the population
        t_reached = t;
        
          for j=4:9
              if t == t_timesteps*(2^j)
                  if abs(sum(average_proposal_over_time(1:(t/2),1)) - sum(average_proposal_over_time((t/2+1):t,1)))/(t/2) <= epsilon
                      if abs(sum(average_threshold_over_time(1:(t/2),1)) - sum(average_threshold_over_time((t/2+1):t,1)))/(t/2) <= epsilon
                          stop = j;
                      end
                  end
              end
          end
        if t > dist_threshold
            for e=1:s
                interval1 = fix(m(e,1)/(bin_size)) + 1;
                if interval1 > bins
                    interval1 = bins;
                end
                interval2 = fix(m(e,2)/(bin_size)) + 1;
                if interval2 > bins
                    interval2 = bins;
                end
                dist(interval1,1) = dist(interval1,1) + 1;
                dist(interval2,2) = dist(interval2,2) + 1;
            end
        end
        
        if max(m(:,1)) == min(m(:,1)) && max(m(:,2)) == min(m(:,2))
            number_of_reproductions = geornd(u); %drawing a random number of reproductions until next mutation
            number_of_time_periods = fix(number_of_reproductions/s); %number of generations passed until next mutation
            mut_location = rem(number_of_reproductions,s); %location of the mutation
            if t + number_of_time_periods <= t_max
                for b=1:number_of_time_periods
                    average_proposal_over_time(t_reached+b,1) = sum(m(:,1))/s;
                    average_threshold_over_time(t_reached+b,1) = sum(m(:,2))/s;
                     if rem(t_reached+b,t_timesteps) == 0
                         snapshot1(:,(t_reached+b)/t_timesteps) = m(:,1); %recording proposals present at point t
                         snapshot2(:,(t_reached+b)/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
                     end
                      for j=4:9
                          if t_reached+b == t_timesteps*(2^j)
                              if abs(sum(average_proposal_over_time(1:((t_reached+b)/2),1)) - sum(average_proposal_over_time(((t_reached+b)/2+1):(t_reached+b),1)))/((t_reached+b)/2) <= epsilon
                                  if abs(sum(average_threshold_over_time(1:((t_reached+b)/2),1)) - sum(average_threshold_over_time(((t_reached+b)/2+1):(t_reached+b),1)))/((t_reached+b)/2) <= epsilon
                                      stop = j;
                                  end
                              end
                          end
                      end
                end
                
                counter2 = counter2 + 1; %number of times mutation occurred
                new_z = rand(1,2); %random number to determine the new mutation
                if mut_location == 0
                    if t + number_of_time_periods > dist_threshold
                        for e=1:s
                            interval1 = fix(m(e,1)/(bin_size)) + 1;
                            if interval1 > bins
                                interval1 = bins;
                            end
                            interval2 = fix(m(e,2)/(bin_size)) + 1;
                            if interval2 > bins
                                interval2 = bins;
                            end
                            if t <= dist_threshold
                                dist(interval1,1) = dist(interval1,1) + t + number_of_time_periods - dist_threshold - 1;
                                dist(interval2,2) = dist(interval2,2) + t + number_of_time_periods - dist_threshold - 1;
                            else
                                dist(interval1,1) = dist(interval1,1) + number_of_time_periods - 1;
                                dist(interval2,2) = dist(interval2,2) + number_of_time_periods - 1;
                            end
                        end
                    end
                    m_total = m_total + m*(number_of_time_periods-1); %accounting for periods in between where nothing happens
                    m(s,:) = new_z; %the mutant replaces the dying agent's strategy
                    t = t + number_of_time_periods; %fastforwarding to the new time period in which mutation occurred
                    t_reached = t;
                else %mut_location > 0
                    if t + number_of_time_periods > dist_threshold
                        for e=1:s
                            interval1 = fix(m(e,1)/(bin_size)) + 1;
                            if interval1 > bins
                                interval1 = bins;
                            end
                            interval2 = fix(m(e,2)/(bin_size)) + 1;
                            if interval2 > bins
                                interval2 = bins;
                            end
                            if t <= dist_threshold
                                dist(interval1,1) = dist(interval1,1) + t + number_of_time_periods - dist_threshold;
                                dist(interval2,2) = dist(interval2,2) + t + number_of_time_periods - dist_threshold;
                            else
                                dist(interval1,1) = dist(interval1,1) + number_of_time_periods;
                                dist(interval2,2) = dist(interval2,2) + number_of_time_periods;
                            end
                        end
                    end
                    m_total = m_total + m*(number_of_time_periods); %accounting for periods in between where nothing happens
                    m(mut_location,:) = new_z; %the mutant replaces the dying agent's strategy
                    mut_index = 0; %check for additional mutations within the same generation
                    while mut_index == 0
                        number_of_reproductions2 = geornd(u); %drawing a random number of reproductions until next mutation
                        if mut_location + number_of_reproductions2 <= s
                            counter2 = counter2 + 1; %number of times mutation occurred
                            mut_location = mut_location + number_of_reproductions2;
                            new_z = rand(1,2); %random number to determine the new mutation
                            m(mut_location,:) = new_z; %the mutant replaces the dying agent's strategy
                        else
                            mut_index = 1; %next mutation goes out of the current generation
                        end
                    end
                    t = t + number_of_time_periods + 1; %fastforwarding to the new time period in which mutation occurred
                    t_reached = t;
                end

                m_total = m_total + m; %adding the current state of the population

                average_proposal_over_time(t,1) = sum(m(:,1))/s;
                average_threshold_over_time(t,1) = sum(m(:,2))/s;
                 if rem(t,t_timesteps) == 0
                     snapshot1(:,t/t_timesteps) = m(:,1); %recording proposals present at point t
                     snapshot2(:,t/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
                 end
                if t > dist_threshold
                    for e=1:s
                        interval1 = fix(m(e,1)/(bin_size)) + 1;
                        if interval1 > bins
                            interval1 = bins;
                        end
                        interval2 = fix(m(e,2)/(bin_size)) + 1;
                        if interval2 > bins
                            interval2 = bins;
                        end
                        dist(interval1,1) = dist(interval1,1) + 1;
                        dist(interval2,2) = dist(interval2,2) + 1;
                    end
                end
                
            else %mutation happens after the max time period is reached
                for b=1:(t_max-t_reached)
                    average_proposal_over_time(t_reached+b,1) = sum(m(:,1))/s;
                    average_threshold_over_time(t_reached+b,1) = sum(m(:,2))/s;
                     if rem(t_reached+b,t_timesteps) == 0
                         snapshot1(:,(t_reached+b)/t_timesteps) = m(:,1); %recording proposals present at point t
                         snapshot2(:,(t_reached+b)/t_timesteps) = m(:,2); %recording acceptance thresholds present at point t
                     end
                end
                if t + number_of_time_periods > dist_threshold
                    for e=1:s
                        interval1 = fix(m(e,1)/(bin_size)) + 1;
                        if interval1 > bins
                            interval1 = bins;
                        end
                        interval2 = fix(m(e,2)/(bin_size)) + 1;
                        if interval2 > bins
                            interval2 = bins;
                        end
                        dist(interval1,1) = dist(interval1,1) + t_max - t;
                        dist(interval2,2) = dist(interval2,2) + t_max - t;
                    end
                end
                m_total = m_total + m*(t_max - t);
                t_reached = t_max;
                m_final_prop(:,1) = m(:,1); %proposals at the last point
                m_final_mao(:,1) = m(:,2); %mao's at the last point
                t = t_max + 1;
            end
        end

        t = t + 1; %moving to the next period
         if stop ~= 0
             t = t_max + 1;
         end
    end
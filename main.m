%% The script corresponds to the Algorithm 5 with r = ( 1 1 1 1)
clc;
clear ;
close all ;
%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the proposed ACOSQ is an iterative algorithm, the distortion
% value for every iteration is provided in this file for given a epsilon and delta.
FileID = fopen ('Results.txt' , 'a') ;
alpha = 200000 ;



%% Number of Quantization level
numLevel = 16 ;

%% Channel's cross-over probability epsilon
epsilon = unique ([10^-5 : 2 * 10^-5 : 10^-4 , 10 ^ -4 : 10^-4 : 10^-3 10^-3 0.005 0.01 0.05 0.1]);

% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use the increase-decrease method.
SIZE = length(epsilon) ;
noise = [1 : SIZE , SIZE : -1 : 1 , 1 : SIZE] ;

%% Distortion parameter
D_4_prime = zeros(length(noise) , 1) ;
D_4 = zeros(length(noise) , 1) ;

%% SDR parameters
SDR_4 = zeros(length(noise) , 1) ;
SDR_3 = zeros(length(noise) , 1) ;
SDR_2 = zeros(length(noise) , 1) ;
SDR_1 = zeros(length(noise) , 1) ;
Final_SDR_4 = zeros(4 , length(noise)) ;
Final_SDR_3 = zeros(36 , length(noise)) ;
Final_SDR_2 = zeros(576 , length(noise)) ;

%% Noise correlation
% The variable delta determines the amount of noise correlation.
%% Noise correlation
for delta = [ 0 5 10]
    % Set the channel's cross-over probability
    for k = 1 : length(noise)
        counter_1 = 0 ;
        counter_2 = 0 ;
        counter_3 = 0 ;
        
        i = noise(k);
        
        Pr_1 = [1 - epsilon(i) , epsilon(i) ;
            epsilon(i) , 1 - epsilon(i)] ;
        
        Pr_z = [(1 - epsilon(i) + delta) / (1 + delta)  , epsilon(i) / (1 + delta) ;
            (1 - epsilon(i)) / (1 + delta)  , (epsilon(i) + delta) / (1 + delta)] ;
        
        % Find the channels transition distribution for a given number of
        % quantization levels i.e. numLevel
        Pr = Channel_with_Memory(numLevel , epsilon(i) , delta) ;
        
        % Set up the parameters for computing the integrals. Note that in
        % this script all integrals are computed numerically using the
        % Riemann summation.
        delta_u = 8 / 2^ 11 ;
        T_1(: , 1) = -4 : delta_u : 4 ;
        u = T_1(: , 1) ;
        % Compute the source pdf. We herein consider a zero-mean
        % unit-variance Gaussian source distribution.
        
        f =  1 ./ (sqrt (2 .* pi)) .* exp (-u .^ 2 ./ 2) ;
        f = f ./ (sum(f) .* delta_u) ;
        Bit_index_2 = [2 2 ; 2 3 ; 2 4 ; 3 2 ; 3 3 ; 3 4 ; 4 2 ; 4 3 ; 4 4]' ;
        Bit_index_3 = [3 3 3 3 ; 3 3 3 4 ; 3 3 4 3 ; 3 3 4 4 ; 3 4 3 3 ; 3 4 3 4 ; 3 4 4 3 ; 3 4 4 4 ; 4 3 3 3 ; 4 3 3 4 ; 4 3 4 3 ; 4 3 4 4 ; 4 4 3 3 ; 4 4 3 4 ; 4 4 4 3 ; 4 4 4 4]' ;
        
        % As noted in the Thesis, the ultimate codebook otabined in the
        % last step of the ACOSQ decribed in Section 4.1 is used as the
        % initial state of the proposed ACOSQ with identical noise
        % correlation and the smallest cross-over probability.
        if (k == 1)
            LOAD  = ['ACOSQ_1_1_1_delta_' num2str(delta)] ;
            load (LOAD) ;
            codebook_1 = hold_codebook_4 ;
        else
            % We slightly increase the channel's cross-over probability,
            % setting the codebook from the system with small epsilon as
            % the initial state of the system with new epsilon.
            load codebook_1
        end
        % The first step of the proposed ACOSQ. Design a 4 bit COSQ
        [SDR_1(k) , ~ , T_1 , codebook_1] = ACOSQ_step_1(f , Pr , numLevel , T_1 , codebook_1 , delta_u , 1 : 16) ;
        
        % save the codebook to initialize the system with the next value of
        % epsilon.
        save('codebook_1'  , 'codebook_1') ;
        
        % save the partition set and codebook for computing the experimental results.
        Data = ['T\T_1_k_' num2str(k) '_delta_' num2str(delta)] ;
        save(Data , 'T_1' , 'codebook_1') ;
        
        % Exhaustively search for the single bit to transmit over the channel.
        for bit_index_1 = 1 : 4
            counter_1 = counter_1 + 1 ;
            
            % Scramble the partition indexes based on the bit chosen for
            % transmision such that the bit chosen for transmision is
            % always located first in the channel input sequence.
            [codebook_2 , ini_T_1] = scramble_labels_step_1(f , Pr , T_1 , numLevel , bit_index_1) ;
            
            % The second step of the proposed ACOSQ design. In this step,
            % the remaining three bits i.e., bits 2 3 4 are generated adaptive to
            % the received bit y_1 corresponding to the previous
            % transmision in step one.
            [SDR_2(k) , T_2 , codebook_2] ...
                = ACOSQ_step_2(f , Pr , Pr_z , codebook_2 , ini_T_1 , numLevel , delta_u) ;
            Final_SDR_2(counter_1 , k) = SDR_2(k) ;
            
            % save the partition set and codebook for the computing the experimental results.
            Data = ['T\T_2_k_' num2str(k) '_counter_1_' num2str(counter_1) '_delta_' num2str(delta)] ;
            save(Data , 'T_2' , 'codebook_2') ;
            
            % As mentioned earlier, this script implements the proposed
            % ACOSQ with r = (1 1 1 1). In the second step, a 4 bit COSQ
            % designed such that the remaining 3 bits are generated
            % adaptive to the y_1 received over the feedback link where y_1
            % is the channel output corresponding to the single bit
            % transmitted in the first step. We herein for every value of y_1 exhaustively search
            % for the best 1 bit out of the generated 3-tuple for transmision.
            
            for index_2 = 1 : 9
                counter_2 = counter_2 + 1 ;
                bit_index_2 = Bit_index_2(: , index_2) ;
                
                % Scramble the quantization cell indexes such that the
                % single bit chosen for transmision is always located in
                % the second place in the channel input
                % sequence.
                [codebook_3 , init_T_2] = scramble_labels_step_2(f , Pr , T_2 , numLevel , bit_index_2) ;
                
                % In the third step of the proposed ACOSQ, a 4 bit COSQ is
                % designed such that the remaining two bits i.e. the third
                % and the fourth bit are generated adaptive to the received
                % sequence y_1y_2, where y_2 is the channel output
                % corresponding to the transmitted sequence in the second
                % step. 
                [SDR_3(k) , T_3 , codebook_3] = ACOSQ_step_3(f , Pr , Pr_z , codebook_3 , init_T_2 , numLevel , delta_u) ;
                Final_SDR_3(counter_2 , k) = SDR_3(k) ;
                
                Data = ['T\T_3_k_' num2str(k) '_counter_2_' num2str(counter_2) '_delta_' num2str(delta)] ;
                save(Data , 'T_3' , 'codebook_3') ;
                
                % In the third step we design a 4 bit COSQ such that the
                % last two bits are generated adaptive to the received
                % sequence y_1y_2. For every value of y_1y_2 we
                % exhaustively search for the best bit out of the 2-tuple
                % generated in the third step for transmision. 
                for index_3 = 1 : 16
                    bit_index_3 = Bit_index_3(: , index_3) ;
                    counter_3 = counter_3 + 1 ;
                    
                    % Scramble the quantization cell indexes such that the
                    % single bit chosen for transmision is always located in
                    % the third place in the channel input
                    % sequence.
                    [codebook_4 , T_3] = scramble_labels_step_3(f , Pr , T_3 , numLevel , bit_index_3) ;
                    
                    % In the fourth step of the proposed ACOSQ, a 4 bit COSQ is
                    % designed such that the remaining one bits i.e. the
                    % the fourth bit is generated adaptive to the received
                    % sequence y_1y_2y_3, where y_3 is the channel output
                    % corresponding to the transmitted sequence in the
                    % third step. 
                    [SDR_4(k) , T_4 , codebook_4 , D_4(k)] = ACOSQ_step_4(f , Pr , Pr_z , codebook_4 , T_3 , numLevel , delta_u) ;
                    
                    Data = ['T\T_4_k_' num2str(k) '_counter_3_' num2str(counter_3) '_delta_' num2str(delta)] ;
                    save(Data , 'T_4' , 'codebook_4') ;
                    Final_SDR_4(counter_3 , k) = SDR_4(k) ;
                end
            end
        end
    end
    myFinal_SDR_2 = max(Final_SDR_2 , [] , 1) ;
    
    myFinal_SDR_3 = max(Final_SDR_3 , [] , 1) ;
    
    myFinal_SDR_4 = max(Final_SDR_4 , [] , 1) ;
    
    % Pick the best SDR value after the end of the so-called
    % increase-decrease method. 
    % The variable final_SDR_2 provides the SDR value corresponding to the Proposed ACOSQ with r = (1 3). 
    final_SDR_2 = zeros(SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var  = myFinal_SDR_2(index) ;
        hold_var = hold_var (:) ;
        [final_SDR_2(i)  , index] = max(hold_var) ;
        fprintf (FileID , 'i = %d\n' , i) ;
        fprintf (FileID , '\nfinal_SDR_2 = %f' , final_SDR_2(i)) ;
        fprintf (FileID , '\nindex = %d\n' , index) ;
    end
    % The variable final_SDR_3 provides the SDR value corresponding to the Proposed ACOSQ with r = (1 1 2). 
    final_SDR_3 = zeros(SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var  = myFinal_SDR_3(index) ;
        hold_var = hold_var (:) ;
        [final_SDR_3(i)  , index] = max(hold_var) ;
        fprintf (FileID , 'i = %d\n' , i) ;
        fprintf (FileID , '\nfinal_SDR_3 = %f' , final_SDR_3(i)) ;
        fprintf (FileID , '\nindex = %d\n' , index) ;
    end
    % The variable final_SDR_4 provides the SDR value corresponding to the Proposed ACOSQ with r = (1 1 1 1). 
    final_SDR_4 = zeros(SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var  = myFinal_SDR_4(index) ;
        hold_var = hold_var (:) ;
        [final_SDR_4(i)  , index] = max(hold_var) ;
        fprintf (FileID , 'i = %d\n' , i) ;
        fprintf (FileID , '\nfinal_SDR_4 = %f' , final_SDR_4(i)) ;
        fprintf (FileID , '\nindex = %d\n' , index) ;
    end
    % Save the best SDR values for every channel parameters.
    Data = ['Proposed_ACOSQ_delta_' num2str(delta)] ;
    save(Data , 'final_SDR_4' , 'myFinal_SDR_4' , 'Final_SDR_4' , 'final_SDR_3' , 'myFinal_SDR_3' , 'Final_SDR_3' , 'final_SDR_2' , 'myFinal_SDR_2' , 'Final_SDR_2')
end

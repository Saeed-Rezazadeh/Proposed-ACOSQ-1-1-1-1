function [SDR_2 , T_2 , codebook] = ACOSQ_step_2(f , Pr , Pr_z , codebook , T , numLevel , delta)
FileID = fopen ('Results.txt' , 'a') ;
D_2 = [2 1] ;
Threshold = 0.001 ;
while ((D_2(1) - D_2(2)) / D_2(2)) > Threshold /4
    D_2(1) = D_2(2) ;
    
    %% Optimal Partitions
    T_u = zeros(length(T) , 2) ;
    for u_index = 1 : length(T)
        summation = 0 ;
        u = T(u_index , 1) ;
        hold_x = T(u_index , 2) ;
        
        binary_x = de2bi(hold_x - 1 , log2(numLevel) , 'left-msb') ;
        
        d = zeros(2 , 8) ;
        
        x_1 = binary_x(1) + 1;
        
        for y_1 = 1 : 2
            for x_2 = 1 : 2
                for x_3 = 1 : 2
                    for x_4 = 1 : 2
                        x_prime = (x_2 - 1) * 4 + (x_3 - 1) * 2 + x_4 ;
                        for y_2 = 1 : 2
                            for y_3 = 1 : 2
                                for y_4 = 1 : 2
                                    y = (y_1 - 1) * 8 + (y_2 - 1) * 4 + (y_3 - 1) * 2 + y_4 ;
                                    
                                    summation = summation + Pr_z(xor(x_1 - 1 , y_1 - 1) + 1 , xor(x_2 - 1 , y_2 - 1) + 1) * ...
                                        Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1) * ...
                                        Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1) * ...
                                        (u - codebook(y)) ^ 2 ;
                                end
                            end
                        end
                        d(y_1 , x_prime) = summation ;
                        summation = 0 ;
                    end
                end
            end
        end
        [~ , partition_index] = min(d , [] , 2) ;
        T_u(u_index , : ) = (x_1 - 1) * 8 + partition_index' ;
    end
    T_2 = cat(2 , T(: , 1) , T_u) ;
    %% Optimal Centroids
    for y_1 = 1 : 2
        for y_prime = 1 : 8
            y = (y_1 - 1) * 8 + y_prime ;
            
            numerator = 0 ;
            denominator = 0 ;
            for x = 1 : numLevel
                
                u_index = find (T_2 (: , 1 + y_1) == x) ;
                u = T_2(u_index , 1) ;
                
                numerator = numerator + Pr(x , y) * sum(u .* f(u_index)) ;
                denominator = denominator + Pr(x , y) * sum(f(u_index)) ;
            end
            codebook(y) = numerator / denominator ;
        end
    end
    %% Distortion
    [D_2(2)] = distortion_2(f , Pr, T_2 , codebook , numLevel , delta) ;
    
    fprintf (FileID , 'Overall D_2 = %f\n' ,D_2(2)) ;
end
fprintf (FileID , 'Overall D_2 = %f\n' ,D_2(2)) ;
SDR_2 = 10 * log10(1 / D_2(2)) ;
fprintf (FileID , 'Overall SDR_2 = %f\n' , SDR_2 ) ;

fprintf (FileID , '=================\n') ;
fclose (FileID) ;
end
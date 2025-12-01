function [Mat] = stitch(X, Nelx,Nely, N)

    C = cell(Nelx, Nely);   % grid of matrices
    
    e = 1;
    for j = 1:Nely
        for i = 1:Nelx
    
            %% To account for repeated rows
            if(j==1)
                y_start = 1;
            else
                y_start = 2;
            end
                
            if(i==1)
                x_start = 1;
                x_end   = N+1;
            elseif(i==N+1)
                x_start = 2;
                x_end   = N;
            else
                x_start = 2;
                x_end   = N+1;
            end

            x = reshape(X(:,e,:), N+1, N+1)';
            x = x(y_start:end,x_start:x_end);

            C{j,i} = x;
            e = e + 1;
        end
    end

    Mat = cell2mat(C);

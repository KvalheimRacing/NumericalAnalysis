% Function that calculates the approximate value of lambda 1 for a companion matrix, given a polynomial P as a row vector
function root = sdroot(p)


    format shortG                       % Pretty format
    maxIterations = 100;                % Numer of iterations
    tol = 1.e-6;                        % Errortolerance
    u = 0;                              % Initializing value for estimation of the dominant lambda

    % Check for leading coefficient == 1
    if ~(p(:,1) == 1)
        fprintf("Leading coefficient does not equal 1, try again")
        exit
    end

    % Makeing Companion matrix C out of p
    n = length(p)-1;                    % Matrix dimension for nxn C
    C = zeros(n-1,1);                   % Zeros above a0
    C(n,1:n) = -fliplr(p(1,2:n+1));     % Adding input vector
    C(1:n-1,2:n) = eye(n-1)             % Adding identety part

    % Making initial iteration vector
    x = zeros(n,1);
    x(n) = 1;


    % Approximating lambda 1
    for i=1:maxIterations


        tmpX = u;                       % Storing previous lambda
        x = C*x;                        % Iterate closer to correct eigenvalues

        [u,nr]=max(abs(x));             % get the dominant eigenvalue
        x = (1/u)*x;                    % Calculate new eigenvector

        % Break out of the loop if the error between the previous and the current lambda is less then the tolerance
        if max(abs(tmpX - u)) < tol
            fprintf('The dominant root of C is %d, achieved after %d iterations', u, i)
            break
        end

    end


    if max(abs(tmpX - u)) > tol
        fprintf(['The desired errortolerance was not achieved in the sequence of 100 iterations, ',...
                 'and C does not have a strictly dominant eigenvalue\n'])
    end

    % Returns approximate value of the dominant lambda for the companion matrix C
    root = max(abs(C*x));


end

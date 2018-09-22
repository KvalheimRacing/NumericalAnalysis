% Function for factorization of a Vandermonde Matrix and approximating of a dataset with a m-1 degree polynomial
function [] = NumericalAnalysis(m)

    % As we try to fit a higher degree of polynomial to our datasets, the factorization becomes more and more ill conditioned. I set the "advisable limit" at m=7
    if m > 7
        warning('With m = %d, you are clearly overfitting dramaticly. This leads to very ill conditioned matrices, and is not advisable!\n', m)
    end

    close all;                                                      % Close all previous figures
    n = 30;                                                         % n observations of data
    A = ones(n,m);                                                  % Initializing A
    x = linspace(-2,2,n);                                           % Where -2 is the start value and 2 is the stop value
    epsilon = 1;                                                    % Noise
    rng(1);                                                         % Initializing the random number generator using a seed of 1 for generating repeatable numbers
    r = rand(1,n) * epsilon;                                        % Generate the repetable random numbers with some noise
    y_1 = (x.*(cos(r+0.5*x.^3)+sin(0.5*x.^3))).';                   % Dataset 1 transposed
    y_2 = (4*x.^5 - 5*x.^4 - 20*x.^3 + 10*x.^2 + 40*x + 10 + r).';  % Dataset 2 transposed

    % Making A a Vandermonde Matrix
    for j=2:m
        for i=1:n

            A(i,j) = x(1,i)^(j-1);

        end
    end



    % ------------------------------------------------------------------------------------------------------------------------%
    %                                                   QR Factorization                                                      %
    % ------------------------------------------------------------------------------------------------------------------------%


    tic % starts a stopwatch timer to measure performance

    % QR decomposition of the Vandermonde matrix
    [Q, R1] = qr(A,0); % When using the qr function with (A, 0), it returns the R1 matrix (originally R from QR is R = [R1 ; 0])

    elapsedTime = toc; % the value of the internal timer at the execution of the toc command
    fprintf('\nElapsed time for QR factiorization is--------: %f\n', elapsedTime)
    tic % restarts timer

    c_1 = backSubstitution2(R1,(Q')*y_1);

    dataset1window = figure('Name','Approximation of dataset 1 using QR factorization and least squares','NumberTitle','off');
    movegui(dataset1window,'north')
    plot(x,A*c_1)   % Where A*c_1 is the approximation of the first dataset with a polynomial of degree m-1
    hold on         % And also
    plot(x,y_1,'o') % The first dataset
    title(['QR - Dataset 1 - Polynomial of degree (m-1) = ' num2str(m-1)])


    c_2 = backSubstitution2(R1,(Q')*y_2);

    dataset2window = figure('Name','Approximation of dataset 2 using QR factorization and least squares','NumberTitle','off');
    movegui(dataset2window,'northeast')
    plot(x,A*c_2)   % Where A*c_2 is the approximation of the second dataset with a polynomial of degree m-1
    hold on         % And also
    plot(x,y_2,'o') % The second dataset
    title(['QR - Dataset 2 - Polynomial of degree (m-1) = ' num2str(m-1)])

    elapsedTime = toc; % the value of the internal timer at the execution of the toc command
    fprintf('Elapsed time using QR factiorization is------: %f\n', elapsedTime)
    tic % restarts timer



    % ------------------------------------------------------------------------------------------------------------------------%
    %                                                Cholesky Factorization                                                   %
    % ------------------------------------------------------------------------------------------------------------------------%


    % Convert the Vandermonde matrix into a square matrix
    B = (A.')*A;

    % Cholesky factorizatio, returns R in the factorization (A.')A = R(R.')
    R = Cholesky(B);

    elapsedTime = toc; % the value of the internal timer at the execution of the toc command
    fprintf('Elapsed time for Cholesky factiorization is--: %f\n', elapsedTime)
    tic % restarts timer

    S = chol(B);       % Not used, this is just to check the time used on Matlabs internal Cholesky factorization

    elapsedTime = toc; % the value of the internal timer at the execution of the toc command
    fprintf('Elapsed time for Matlabs internal chol(B) is-: %f\n', elapsedTime)
    tic % restarts timer

    c1 = backSubstitution2(R.', forwardSubstitution(R,(A.')*y_1));

    dataset3window = figure('Name','Approximation of dataset 1 using Cholesky factorization and least squares','NumberTitle','off');
    movegui(dataset3window,'south')
    plot(x,A*c1)    % Where A*c1 is the approximation of the first dataset with a polynomial of degree m-1
    hold on         % And also
    plot(x,y_1,'o') % The first dataset
    title(['Cholesky - Dataset 1 - Polynomial of degree (m-1) = ' num2str(m-1)])

    c2 = backSubstitution2(R.', forwardSubstitution(R,(A.')*y_2));

    dataset4window = figure('Name','Approximation of dataset 2 using Cholesky factorization and least squares','NumberTitle','off');
    movegui(dataset4window,'southeast')
    plot(x,A*c2)    % Where A*c2 is the approximation of the second dataset with a polynomial of degree m-1
    hold on         % And also
    plot(x,y_2,'o') % The second dataset
    title(['Cholesky - Dataset 2 - Polynomial of degree (m-1) = ' num2str(m-1)])

    elapsedTime = toc; % the value of the internal timer at the execution of the toc command
    fprintf('Elapsed time using Cholesky factiorization is: %f\n\n', elapsedTime)

    % ------------------------------------------------------------------------------------------------------------------------%
    %                              Discussion: Conditioning of QR vs Conditioning of Cholesky                                 %
    % ------------------------------------------------------------------------------------------------------------------------%


    fprintf('Condition number of A is---------------------: %.3g\n', cond(A))   % Using pre defined matlab function for calculating the condition number.
    fprintf('Condition number of B in Cholesky is---------: %.6g\n', norm(B)*norm(inv(B)))   % Using the definition of condition number.
    fprintf('Condition number of R in Cholesky is---------: %.3g\n', norm(R)*norm(inv(R)))   % Using the definition of condition number.
    fprintf('Condition number of R1 in QR is -------------: %.3g\n', norm(R1)*norm(inv(R1))) % This will be the same as for the whole R matrix
    fprintf('Condition number of Q  in QR is -------------: %.3g\n\n', cond(Q)) % This is obviously 1 since Q*Q^T = I

    % We see that when we compare Matlabs own QR and Cholesky algorihms, Cholesky is always twice (or more) as fast. This corresponds to the litterature
    % ref; http://www.math.ualberta.ca/ijnam/Volume-13-2016/No-1-16/2016-01-07.pdf
    % Further, we see that both methods has the same condition number for the "R" matrix, meaning in QR, the R matrix has the same condition number as the
    % R matrix in Cholesky, where R = LD^(1/2). However since the Q matrix in QR consists of an orthonormal set of vectors, Q^TQ = I, therfore the condition
    % number of Q is 1. Since Cholesky factorizes the matrix B into R*R^T we have two potentially ill conditioned matrices, where in QR, we only have one.
    % This difference will likely not have much to say for low condition numbers, but for high condition numbers and big matrices, that would lead us to prefere
    % the QR factorization, which the litterature also agrees with.
    % Another point, which is very obvoius is that in order to use Cholesky, we need to make A into a symmetric matrix B. It's pretty straight forward to see
    % that the Condition number og this (B) matrix is very bad. Therfore, for non symmetric matrices, I believe the best choice is QR factorization. Wee also see
    % that the implementation of the factorization in order to solve the least squares problem is quicker with QR.
    % But if A initially was a symmetric matrix, and we were trying to optimize for speed, I beleve Cholesky is the best choice (for low condition numbers).


end




% Finds x, where Ux = b, for An Upper Triangular Matrix U and right hand side b
function x = backSubstitution2(u,b)

    [n,m] = size(u);
    x = zeros(n,1);

    % Check for singular values
    if (abs(u(n,m)) > 1e-12)
        x(n) = b(n)/u(n,m);
    else
        disp('U is singular')
        return
    end

    for k=(n-1):-1:1

        x(k) = (b(k)-u(k,k+1:n)*x(k+1:n))/u(k,k);

    end

end




% Cholesky decomposition of a square, positive matrix: A = LD(L^T). Returns R = L*D^(1/2)
function [LD_squareRoot] = Cholesky(A)


    [n,m] = size(A);
    L = zeros(n,n);
    D  = L;
    B = A;

    if n == m
        for k=1:n

            L(:,k) = B(:,k)/B(k,k);
            D(k,k) = B(k,k);
            B = B - D(k,k)*L(:,k)*(L(:,k))';

        end
    else
        fprintf('Input is not a squate matrix')
    end

    fprintf('Condition number of L in Cholesky is---------: %.3g\n', cond(L))
    fprintf('Condition number of D^(1/2) in Cholesky is---: %.3g\n', cond((D^(1/2))))

    % This additional multiplication is most likely one of the biggest contributors to the difference
    % in time used on this cholesky factorization, contra matlabs internal chol(B)
    LD_squareRoot = L*D^(1/2);


end




% Returns the solution x, where A is a square lower triangular matrix in the equation Ax = b
function x = forwardSubstitution(A,b)

    [n,m] = size(A);
    if n ~= m
        error('Input is not a square matrix');
    end

    if size(b,1) ~= n
        error('Input dimensions do not match')
    end

    x = zeros(n,1);

    for k=1:n
        if abs(A(k,k)) > 1e-12

            tmp = 0;
            if k > 1
                for j=1:k-1
                    tmp = tmp + A(k,j)*x(j);
                end
            end

            x(k) = (b(k)-tmp)/A(k,k);

        else
            error('Input Singular')
            return;
        end
    end
end

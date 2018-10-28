fprintf('\n--------------------------------------TASK 4---------------------------------------\n')

p1 = [1, 3, -1, -3, -1, 1];
p2 = [1, -1, 4, -4];

companionMatrices(p1)
companionMatrices(p2)

roots_of_p1 = roots(p1)
roots_of_p2 = roots(p2)

fprintf(['When comparing the result from the power method with matlabs root function,\nwe see ',...
'that the method finds the dominant root of p1, but not of p2.\nThis is expected as ',...
'p2 has two complex conjugate poles that have the biggest absolute value of all roots.\n'])


fprintf('\n--------------------------------------TASK 5---------------------------------------\n')

Pascal = [ 1 1 1  1  1   1
           1 2 3  4  5   6
           1 3 6  10 15  21
           1 4 10 20 35  56
           1 5 15 35 70  126
           1 6 21 56 126 252 ];

p3 = poly(Pascal);

companionMatrices(p3)

l = sprintf('\n%d', eig(Pascal));
fprintf('Eigenvalues of Pascal matrix using eig function =\n%s\n', l);

fprintf(['\nWhen comparing the result from the power method with matlabs eig function,\nwe see ',...
'that the method finds the dominant root of the Pascal matrix, thus the answer is reasonable.\n'])


fprintf('\n--------------------------------------TASK 6---------------------------------------\n')

Pascal_Eigen_Values_On_Diagonal_After_Using_Iterations_With_QR = Pascal;

for i=1:9

    [Q,R] = qr(Pascal_Eigen_Values_On_Diagonal_After_Using_Iterations_With_QR);
    Pascal_Eigen_Values_On_Diagonal_After_Using_Iterations_With_QR = R*Q;
    if i==9 Pascal_Eigen_Values_On_Diagonal_After_Using_Iterations_With_QR, end

end

fprintf(['When comparing the result from iterations of QR factorization with matlabs eig function,\nwe see ',...
'that the iterative approach finds the eigenvalues of the Pascal matrix, located on the diagonal.\n'])

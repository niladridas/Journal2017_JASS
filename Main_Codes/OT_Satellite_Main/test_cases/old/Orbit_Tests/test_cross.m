function C = test_cross(A,B)
    C = [A(2,1)*B(3,1)-A(3,1)*B(2,1);
         A(3,1)*B(1,1)-A(1,1)*B(3,1);
         A(1,1)*B(2,1)-A(2,1)*B(1,1)];
end
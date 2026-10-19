function Mattrix = Constell_Mat(K,N)
    if K>1
        tempmat1 = kron(ones(N,1),Constell_Mat(K-1,N));
        tempmat2 = repmat([1:N],N^(K-1),1);
        Mattrix = [tempmat2(:),tempmat1];
    else
        Mattrix = [1:N]';
end
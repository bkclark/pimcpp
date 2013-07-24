N = 41;

A = rand(N,N);
Acopy = A;
Ainv = A^(-1);

mydet = 1.0;

ipiv = [1:N];
for k = [1:N]
   [abspivot,kb] = max([zeros(k-1,1);abs(A(k:N,k))]);
   pivot = A(kb,k);
   [A(k,:),A(kb,:)] = swap(A(k,:),A(kb,:));
   [ipiv(k),ipiv(kb)] = swap(ipiv(k),ipiv(kb));
   A(:,k) = -[A(1:k-1,k);0;A(k+1:N,k)]/pivot;
   A = A + A(:,k)*A(k,:);
   A(k,:) = [A(k,1:k-1),1,A(k,k+1:N)]/pivot;
   mydet = mydet * pivot;
   if (kb != k)
     mydet = -mydet;
   end
end

A(:,ipiv) = A;
(A*Acopy);
mydet
det(Acopy)


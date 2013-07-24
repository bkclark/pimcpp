N = 128;
BS = 32;

A = rand(N,N);
Acopy = A;
NB = N/BS;

mydet = eye(BS);

z = zeros(BS,BS);
one = eye(BS);

dets = zeros(1,NB);
mask = ones(1,NB);

ipiv = [1:NB];

dopivot = 0

for k = 1:NB
   x = (k-1)*BS + 1;
   y = k*BS;

   for i = 1:NB
     dets(i) = mask(i)*det(A((i-1)*BS+1:i*BS,x:y));
   end
   [abspivot, kb] = max(abs(dets));
   
   xb = (kb-1)*BS + 1;
   yb = (kb*BS);

   if (dopivot == 1)
     [A(x:y,:),A(xb:yb,:)] = swap(A(x:y,:),A(xb:yb,:));
     [ipiv(k),ipiv(kb)] = swap(ipiv(k), ipiv(kb));
   end

   pivot = A(x:y,x:y);
   mydet = pivot * mydet;
   pivotInv = pivot^(-1);
   A(:,x:y) = -[A(1:x-1,x:y);z;A(y+1:N,x:y)] *pivotInv;
   A      = A + A(:,x:y)*A(x:y,:);
   A(x:y,:) = pivotInv*[A(x:y,1:x-1),one,A(x:y,y+1:N)];
   mask(k) = 0.0;
end

% Reverse block permutations
rows = reshape([1:N],N/NB,NB);
rows = rows(:,ipiv);
rows = reshape(rows,N,1);
A(:,rows) = A;

nopivot = sqrt(norm(A*Acopy-eye(N),'fro')/(N*N))


A = Acopy;

dopivot = 1
ipiv = [1:NB];
mask = ones(1,NB);
for k = 1:NB
   x = (k-1)*BS + 1;
   y = k*BS;

   for i = 1:NB
     dets(i) = mask(i)*det(A((i-1)*BS+1:i*BS,x:y));
   end
   [abspivot, kb] = max(abs(dets));
   
   xb = (kb-1)*BS + 1;
   yb = (kb*BS);

   if (dopivot == 1)
     [A(x:y,:),A(xb:yb,:)] = swap(A(x:y,:),A(xb:yb,:));
     [ipiv(k),ipiv(kb)] = swap(ipiv(k), ipiv(kb));
   end

   pivot = A(x:y,x:y);
   mydet = pivot * mydet;
   pivotInv = pivot^(-1);
   A(:,x:y) = -[A(1:x-1,x:y);z;A(y+1:N,x:y)] *pivotInv;
   A      = A + A(:,x:y)*A(x:y,:);
   A(x:y,:) = pivotInv*[A(x:y,1:x-1),one,A(x:y,y+1:N)];
   mask(k) = 0.0;
end

% Reverse block permutations
rows = reshape([1:N],N/NB,NB);
rows = rows(:,ipiv);
rows = reshape(rows,N,1);
A(:,rows) = A;

withpivot = sqrt(norm(A*Acopy-eye(N),'fro')/(N*N))

norm(A*Acopy-eye(N))/norm(Acopy^(-1)*Acopy - eye(N))

nopivot/pivot


(det(mydet)-det(Acopy))/det(Acopy)
det(Acopy)

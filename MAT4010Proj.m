function [] = MAT4010Proj()
  %error tolerance
  epsilon = 10^-2;
  %n values for each iteration, to be called by index
  n = [5,10,15,20,25,30];
  %w for SOR method
  w = 3/2;
  %computation time t
  t = 0;
  %vector of average computation times for each method and each n
  avgt = zeros(6,3);
  %vector of residual vectors, only for when n = 30. Number of columns are 1 
  %for now, as number of columns is equal to iteration number k
  r(:,1) = zeros(30,1);
  %vector of infinity norms of the residuals for each method. Number of columns
  %is equal to iteration number 1 until k is established
  rnorm(:,1) = zeros(1,3);
  %vector of error vectors. Number of columns is equal to 1 until k is
  %established
  err(:,1) = zeros(30,1);
  %vector of infinity norms of error
  errnorm(:,1) = zeros(1,3);
  %vector to hold values in x vector when k = 10,
  %only for when x is a 30x1 vector
  kten = zeros(30,1);
  %vector to hold values in x vector when k is at its last iteration,
  %only for when x is a 30x1 vector
  kfinal = zeros(30,3);
  %vector to hold number of iterations required before system is solved
  kiterates = zeros(6,3);
  
  for nindex = 1:6
    %initialize nx1 x, y, and z vectors to hold initial guesses
    %between -20 and 20. Vectors will be nx1, and n will be
    %chosen by index
    x = 40.*rand(n(nindex),1) - 20;
    y = x;
    z = x;
    
    %initialize nx1 vector b and nxn matrix A
    b = ones(n(nindex),1);
    A = full(gallery('tridiag',n(nindex),1,-2,1));
    %assign true solutions of SNLE to xTrue for error analysis
    xTrue = A \ b;
    
    %jacobi method
    for count = 1:5
      %printf('jacobi')
      k = 1;
      while ( norm((A*x(:,k) - b), inf ) > epsilon )
        tic;
        for i = 1:n(nindex)
          sum = 0;
          for j = 1:n(nindex)
            if (j ~= i)
            sum = sum + A(i,j).*x(j,k);
            end
          end
          x(i,k+1) = (1/A(i,i))*(b(i)-sum);
        end
          
        if (n(nindex) == 30)
          %error vector
          err(:,1) = (x(:,k)-xTrue);
          %error norm
          errnorm(k,1) = norm(err(:,1), inf);
          %residual vector
          r(:,1) = A*err;
          %residual norm
          rnorm(k,1) = norm(r, inf);
          %hold vector x10 in a vector
          if (k == 10)
            kten(:,1) = x(:,k);
          end
        end
        %kiterates(nindex, 1) = k;
        k = k + 1;
      end
      kiterates(nindex,1) = k;
      %hold vector xk in a vector
      if (size(x(:,nindex)) == 30)
        kfinal(:,1) = x(:,k-1)
      end
      t = toc;
      avgt(nindex,1) = (avgt(nindex,1) + t) / 5;
    end
    
    %gauss-seidel method
    for count = 1:5
      %printf('gauss-seidel');
      k = 1;
      while (norm((A*y(:,k) - b), inf) > epsilon)
        tic;
        k = k + 1;
        for i = 1:n(nindex)
          sigma = 0;
          
          for j = i+1:n(nindex)
            sigma = sigma + A(i,j).*y(j,k-1);
          end
          
          for j = 1:i-1
            sigma = sigma + A(i,j).*y(j,k);
          end
          
          y(i,k) = (1/A(i,i))*(b(i)-sigma);
        end
        
        if (n(nindex) == 30)
          %error vector
          err(:,2) = (y(:,k-1)-xTrue);
          %error norm
          errnorm(k-1,2) = norm(err(:,2), inf);
          %residual vector
          r(:,2) = A*err(:,2);
          %residual norm
          rnorm(k-1,2) = norm(r(:,2), inf);
          %hold vector y10 in a vector
          if (k == 10)
            kten(:,2) = y(:,k);
          end
        end
        %kiterates(nindex, 2) = k;
        %k = k + 1;
      end
      kiterates(nindex, 2) = k;
      %hold vector yk in a vector
      if (size(y(:,nindex)) == 30)
        kfinal(:,2) = y(:,k);
      end
      t = toc;
      avgt(nindex,2) = (avgt(nindex,2) + t) / 5;
    end
    
    %SOR method
    for count = 1:5
      %printf('sor');
      k = 1;
      while (norm((A*z(:,k) - b), inf) > epsilon)
        k = k + 1;
        for i = 1:n(nindex)
          sigma = 0;
          
        for j = i+1:n(nindex)
          sigma = sigma + A(i,j).*z(j,k-1);
        end
        
        for j = 1:i-1
          sigma = sigma + A(i,j).*z(j,k);
        end
        
        z(i,k) = (1-w).*z(i,k-1) + (w/A(i,i))*(b(i)-sigma);
        end
        
        if (n(nindex) == 30)
          %error vector
          err(:,3) = (z(:,k-1)-xTrue);
          %error norm
          errnorm(k-1,3) = norm(err(:,3));
          %residual vector
          r(:,3) = A*err(:,3);
          %residual norm
          rnorm(k-1,3) = norm(r(:,3));
          %hold vector z10 in a vector
          if (k == 10)
            kten(:,3) = z(:,k);
          end
        end
        %kiterates(nindex,3) = k;
        %k = k + 1;
      end
      kiterate(nindex,3) = k;
      %hold vector zk in a vector
      if (size(z(:,nindex)) == 30)
      kfinal(:,3) = z(:,k);
      end
      t = toc;
      avgt(nindex,3) = (avgt(nindex,3) + t) / 5;
    end
    printf('\n');
  end
  
  %Plot of infinity norm of the error versus k
  figure(1)
  clf;
  plot(1:size(errnorm(:,1),1) , errnorm(1:size(errnorm(:,1),1)) , 'k');
  hold on;
  plot(1:size(errnorm(:,2),2) , errnorm(1:size(errnorm(:,2),2)) , 'b');
  plot(1:size(errnorm(:,3),3) , errnorm(1:size(errnorm(:,3),3)) , 'r');
  hold off;
  xlabel('k');
  ylabel('error norm');
  title('Infinity norm of the error versus k');
  
  %Plot of infinity norm of the residual versus k
  figure(2)
  clf;
  plot(1:size(rnorm(:,1),1) , rnorm(1:size(rnorm(:,1),1)) , 'k');
  hold on;
  plot(1:size(rnorm(:,2),2) , rnorm(1:size(rnorm(:,2),2)) , 'b');
  plot(1:size(rnorm(:,3),3) , rnorm(1:size(rnorm(:,3),3)) , 'r');
  hold off;
  xlabel('k');
  ylabel('residual norm');
  title('Infinity norm of the residual versus k');
  
  %Plot of the 10th x vector versus index i
  figure(3)
  clf;
  plot(1:30 , kten(1:30,1) , 'k');
  hold on;
  plot(1:30 , kten(1:30,2) , 'b');
  plot(1:30 , kten(1:30,3) , 'r');
  hold off;
  xlabel('i');
  ylabel('x');
  title('10th vector of x values versus i');
  
  %Plot of the final x vector versus index i
  figure(4)
  clf;
  plot(1:30 , kfinal(1:30,1) , 'k');
  hold on;
  plot(1:30 , kfinal(1:30,2) , 'b');
  plot(1:30 , kfinal(1:30,3) , 'r');
  hold off;
  xlabel('i');
  ylabel('x');
  title('final x vector values versus index i');
  
  %Plot of the average time taken for each method versus n
  figure(5)
  clf;
  plot(1:6 , avgt(1:6,1) , 'k');
  hold on;
  plot(1:6 , avgt(1:6,2) , 'b');
  plot(1:6 , avgt(1:6,3) , 'r');
  hold off;
  xlabel('n');
  ylabel('average computation time');
  title('average time taken to compute results versus n');
  
  %Plot of the number of iterations k required for each method versus n
  figure(6)
  clf;
  plot(1:6 , kiterates(1:6,1) , 'k');
  hold on;
  plot(1:6 , kiterates(1:6,2) , 'b');
  plot(1:6 , kiterates(1:6,3) , 'r');
  hold off;
  xlabel('n');
  ylabel('iterations k');
  title('number of iterations k required for each method versus n');
  
end
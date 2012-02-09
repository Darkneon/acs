function n_o = LEN(mat)
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The wisdom of stack-overflow %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % http://stackoverflow.com/questions/7209521/vector-norm-of-an-array-of-vectors-in-matlab
    % twoNorm = sqrt(sum(abs(M).^2,1)); %# The two-norm of each column
    
    n_o = sqrt(sum(abs(mat).^2,1));

    %%%%%%%%%%%%%%%%
    % inefficient  %
    %%%%%%%%%%%%%%%%
    
    %disp(mat);
    %[r,c]=size(mat)
    %for k=1:c
    %    n_o(:,k)=norm(mat(:,k),2);
    %end
    %disp(n_o);
    
end
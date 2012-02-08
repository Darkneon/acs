function n_o = LEN(mat)

[r,c]=size(mat)
for k=1:c
    n_o(k)=norm(mat(:,k));
end
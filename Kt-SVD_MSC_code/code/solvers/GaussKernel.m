function  KX =GaussKernel(Xv)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
     dis=0;
     [m,n]=size(Xv);
     KX=zeros(n,n);
     for i=1:n
         KX(i,i)=1;
         for j=i+1:n
             disij=norm(Xv(:,i)-Xv(:,j),2);
             KX(i,j)=-disij^2;
             KX(j,i)=-disij^2;
             dis=dis+disij;
         end
     end
     %dis
     
     sigma=dis/n*(n-1);
     %sigma
     KX(:,:)=exp(KX/(2*sigma^2));
     %KX(:,:)=exp(KX/sigma);
end


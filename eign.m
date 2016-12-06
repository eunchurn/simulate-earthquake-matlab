function [Lambda,Phi,Psi]=eign(A,B)
%     eign: Solve the generalized eigenvalue problem.
%           Order & normalize the eigenvectors.
%
%     [Lambda,Phi,Psi]=eign(A,B)
%
%              A x = [Lambda] B x
%     with the special (`structural') normalizations:     
%              Phi(i)'*B*Phi(i) = 1
%              Psi(i)'*B*psi(j) = kronecker delta(i,j)
%     where    
%              Phi(i) -- i-th right eigenvector
%              Psi(i) -- i-th left  eigenvector
%                 '   -- transpose
%
%     Caution: real mode -- imag(Lambda) < tolerance
%              tolerance = 1E-13
%
%     reference: Junkins & Kim, Dyn. & Ctrl of Structures, Ch. 2,4
%
%
%     programmed by  Youdan Kim
%                    Dept. of Aerospace Engineering
%                    Texas A&M University
%
%     revised date : May  24, 1989
%                    APR. 25, 1991
%                    Dec. 16, 1991
n=max(size(A)); Lambda=zeros(n,1); Phi=zeros(n); Psi=zeros(n);
tolerance=1.e-13;
%
%    Solve Left and Right Eigenvalue Problem
%
[VR,DR]=eig(A,B); [VL,DL]=eig(A',B');
kr=zeros(n,1); kc=kr; er=kr; ec=kr;
%
%    Sort Right Eigenvectors and Eigenvalues
%
indr=0; indc=0;
for i=1:n;
     if abs(imag(DR(i,i))) <= tolerance;
          indr=indr+1; kr(indr)=i; er(indr)=DR(i,i);
     elseif imag(DR(i,i)) > tolerance;
          indc=indc+1; kc(indc)=i; ec(indc)=DR(i,i);
     end
end
er=real(er(1:indr)); ec=ec(1:indc);
ind=1; [lr,krn]=sort(er); [lc,kcn]=sort(imag(ec));
for i=1:indr+indc;
     if i <= indr;
          Phi(:,i)=real(VR(:,kr(krn(i))));
          Lambda(i)=real(DR(kr(krn(i)),kr(krn(i))));
          ind=ind+1;
     else
          ii=i-indr;
          Phi(:,ind)=VR(:,kc(kcn(ii)));
          Phi(:,ind+1)=conj(Phi(:,ind));
          Lambda(ind)=DR(kc(kcn(ii)),kc(kcn(ii)));
          Lambda(ind+1)=conj(Lambda(ind));
          ind=ind+2;
     end
end
%
%     Sort Left Eigenvectors
%
indr=0; indc=0;
for i=1:n;
     if abs(imag(DL(i,i))) <= tolerance;
          indr=indr+1; kr(indr)=i; er(indr)=DL(i,i);
     elseif imag(DL(i,i)) > tolerance;
          indc=indc+1; kc(indc)=i; ec(indc)=DL(i,i);
     end
end
er=real(er(1:indr)); ec=ec(1:indc);
ind=1; [lr,krn]=sort(er); [lc,kcn]=sort(imag(ec));
for i=1:indr+indc;
     if i <= indr;
          Psi(:,i)=real(VL(:,kr(krn(i))));
          ind=ind+1;
     else
          ii=i-indr;
          Psi(:,ind)=VL(:,kc(kcn(ii)));
          Psi(:,ind+1)=conj(Psi(:,ind));
          ind=ind+2;
     end
end
%
%     Normalize Right and Left Eigenvectors
%
for i=1:n;
     xi=Phi(:,i);
     yi=Psi(:,i);
     sc1=conj(xi')*B*xi;       Phi(:,i)=Phi(:,i)/sqrt(sc1);
     sc2=conj(yi')*B*Phi(:,i); Psi(:,i)=Psi(:,i)/sc2;
end

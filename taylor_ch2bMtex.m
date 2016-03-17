function [  OUT MT ROT]=taylor_ch2bMtex(Eij,Eul)
%% T.H. SIMM 2016- To compute Taylor model using Bunge, H. J. (1970).
%% Some applications of the Taylor theory of polycrystal plasticity. 
%% Kristall Und Technik, 5(1), 145–175. http://doi.org/10.1002/crat.19700050112
%%
%%for use with cubic and stress tensor input from Mtex that have been
%%rotated to the crystal symmetry e.g. as below where M is the strain tensor e.g for tensile shown....
% % % Eul=niceEuler(inv(ori));
% % % M = zeros(3);M(1,1) = 1;M(2,2)=-.5;M(3,3)=-.5;sigma001 = tensor(M,'name','stress');
% % % sigmaCS = rotate(sigma001,inv(ori));
% % % for ii=1:length(sigmaCS)
% % %     eijn=double(sigmaCS(ii));
% % %     [bn(:,ii) ang_(ii)]=taylor_ch2bMtex(eijn,Eul(ii,:));
% % % end

%% speeds up computation by reducing available options to account for
%% a matrix having no inverse i.e. the code is run for all combinations and
%% if the determinant is non-zero then this is assigned as a valid
%% combination and saved
load('fcc12sy','sy_b')
%% reduce size of deformation tensor so the strain is small here= 0.1% strain
step=1e-3;%
Eij=Eij*step;
%% %define slip systems plane normals, pl, Burgers vectors, d, here for FCC
pl(:,:,1)=[1 1 1;1 1 1;1 1 1];
pl(:,:,2)=[-1 1 1;-1 1 1;-1 1 1];
pl(:,:,3)=[1 -1 1;1 -1 1;1 -1 1];
pl(:,:,4)=[1 1 -1;1 1 -1;1 1 -1];
d(:,:,1)=[1 -1 0; 1 0 -1; 0 1 -1];
d(:,:,2)=[1 1 0; 1 0 1; 0 1 -1];
d(:,:,3)=[1 1 0; 1 0 -1; 0 1 1];
d(:,:,4)=[1 -1 0; 1 0 1; 0 1 1];

%% normalise vectors
d(:,:,1:4)=d(:,:,1:4)*sqrt(0.5);
pl=pl*1/sqrt(3);
%% define CRSS, if more than one we may define a different value for
%% different systems- or if twinning we may want a directional CRSS
CRSS=ones(1,12);

%% create deformation tensor in unit cell symmetry --> removed as found by
%% mtex

% % % %%          create matrices- applied deformation tensor
% % % Eij_m=zeros(3,3);
% % % for i=1:3
% % %     for j=1:3
% % %         for n=1:3
% % %             for m=1:3
% % %                 Eij_m(i,j)=Eij_m(i,j)+g(i,m)*g(j,n)*t(m,n);
% % %             end
% % %         end
% % %     end
% % % end
% % % Eij=[Eij_m(1,:),Eij_m(2,:),Eij_m(3,:)];

%%          create matrices (deformation tensors) of slip systems
h=1;
eijn_mat=zeros(3,3,size(d,3)*3);
for k=1:size(d,3)
    for m=1:3
        eijn_mat(:,:,h)=pl(m,:,k)'*d(m,:,k);
        h=h+1;
    end
end
%% Create the symmetric and anti-symmetric components of the deformation
%% tensor - eijnc is pure strain and used to match the imposed deformation
%% and calculate slip systems
%% and eijnc_II is pure rotation used to calculate rotation changes
eijnc=zeros(3,3,size(eijn_mat,3));
eijnc_II=zeros(3,3,size(eijn_mat,3));
for n=1:size(eijn_mat,3)
    for i=1:3
        for j=1:3
            eijnc(i,j,n)=0.5*(eijn_mat(i,j,n)+eijn_mat(j,i,n));
            eijnc_II(i,j,n)=0.5*(eijn_mat(i,j,n)-eijn_mat(j,i,n));
        end
    end
end
%% turn matrices into vectors 1 x 9
eijn=zeros(size(eijn_mat,3),9);
eijnII=zeros(size(eijn_mat,3),9);
for n=1:size(eijn_mat,3)
    eijn(n,:)=[eijnc(:,1,n);eijnc(:,2,n);eijnc(:,3,n)];
    eijnII(n,:)=[eijnc_II(:,1,n);eijnc_II(:,2,n);eijnc_II(:,3,n)];
end

%%          only need 5 of the def tensor components due to symmetry and
%%          since exx+eyy+ezz=1, so reduce to 5 components

n1=1;n2=2;n3=3;n4=5;n5=6;
eijn_ck=[eijn(:,n1),eijn(:,n2),eijn(:,n3),eijn(:,n4),eijn(:,n5)];%%%%the independent ij values
Eij_ck=[Eij(n1),Eij(n2),Eij(n3),Eij(n4),Eij(n5)];

%%  Processing part
% create variable bn which represents the activity of the 5 slip systems
% for all possible slip system combinations where 5 of the 12 slip systems are active and the rest have no activity 
bn2(:,:)=zeros(5,size(sy_b,2));

h=1;
% loop through each slip system combination
      for n=1:size(sy_b,2)
                sy=sy_b(:,n);
                fsy=find(sy>0);%find the active slip systems for the combination
                % create a matrix eijn_ck2- each row represents a different
                % slip system and each column the independent components of
                % their strain tensor
                eijn_ck2=[eijn_ck(fsy(1),:);eijn_ck(fsy(2),:);eijn_ck(fsy(3),:);eijn_ck(fsy(4),:);eijn_ck(fsy(5),:)];
                if det(eijn_ck2)>eps(1)%%%a matrix A has no inverse if it is singular or detA=0 so we want to ignore these eps is a small number
                    %% the main solver part here we are effectively finding the slip activity 
                    %% of different systems (bn) that when multiplied by their strain tensors (eijn_ck2)
                    %% give the applied strain tensor Eij_ck
                    bn2(:,h)=eijn_ck2'\Eij_ck';
                    % and the active slip systems of the h-th element of bn2 
                    sy_b(:,h)=sy_b(:,n);
                     if sum(abs(bn2(:,h)))~=0 % another check for no solutions else next variable prints over this one
                         h=h+1;
                     end
                end
      end
      
 sy_b=sy_b(:,1:(h-1));
%% Find bn combination with minimum work i.e. sum(bn) is minimised
  
   siz_bn=zeros(1,h-1);
   for m=1:(h-1)%%%add up all bn for all non inf combos
       ssss=find(sy_b(:,m)~=0);%fthe active slip systems
       M_add=0;
       % find Taylor factor for each combination -> siz_bn
       for ff=1:5           
            M_add=M_add+abs(bn2(ff,m))*CRSS(ssss(ff));
       end
       siz_bn(m)=M_add;
   end
   %% 
    bnmin_1st=find(siz_bn==min((siz_bn)),1) ; %% finds position of slip sys with minimum taylor factor
    siz_bn=abs(siz_bn(bnmin_1st)-siz_bn);
    %%allow for loss of accuracy from g i.e. +-90 Euler 3 should give same
    %%results (but different order of bn)
    bnmin=find(siz_bn<10e-014); %% find position of combinations with minimum TaylorFactor
    lenbnmin=length(bnmin)     ;         
    clear siz_bn
    %define variables
    bn=zeros(5,lenbnmin);slipcomb=zeros(size(eijn_mat,3),lenbnmin);slipcomb_no=zeros(5,lenbnmin);
    bn_all=zeros(size(eijn_mat,3),lenbnmin);
    
    for x=1:lenbnmin
        bn(:,x)=bn2(:,bnmin(x)) ;   %%%the bn of the leat work (min Taylor factor sums)
        %identify which slip systems are active in the low M options
        slipcomb(:,x)=(sy_b(:,bnmin(x))) ;
        slipcomb_no(:,x)=find(slipcomb(:,x)>0);%%%and their slip systems- the 5 of the 12
        % find slip activity of active slip systems
        for n=1:5
            bn_all(slipcomb_no(n,x),x)=bn(n,x);
        end
    end
    clear bn2 bn sy_b
    %% if there are multiple solutions we take the average of these - other
    %% methods to determine the solution exist
    if lenbnmin>1;
        bn_mean=mean(bn_all');
    else
        bn_mean=bn_all';
    end
%% Rotation The crystal lattice rotation is then the sum of the rotations
%% for each slip system weighted by their activity bn
    
        rk = -bn_mean*eijnII;
        rkM = [rk(1) , rk(2), rk(3);rk(4) , rk(5), rk(6);rk(7) , rk(8), rk(9)];
        Eul=Eul*pi/180;phi1=Eul(1);PHI=Eul(2);phi2=Eul(3);
        d_Phi1=sin(phi2)/sin(PHI);
% the rotation vector        
rk2=[rk(1:3);rk(4:6);rk(7:9)];
r3=rk2(1,2)-rk2(2,1); 
r2=rk2(1,3)-rk2(3,1);
r1=rk2(2,3)-rk2(3,2);
r1=-r1;
r2;
r3=-r3;
r=[r1 r2 r3];
d_phi1=(sin(phi2)/sin(PHI))*r1 +(cos(phi2)/sin(PHI))*r2; 
d_PHI=cos(phi2)*r1-sin(phi2)*r2;
d_phi2=-(sin(phi2)/tan(PHI))*r1 -(cos(phi2)/tan(PHI))*r2 + r3  ;

d_phi1=d_phi1*180/pi;if isnan(d_phi1)==1;d_phi1=0;end%isnan checks for NaN values from sin(PHI)/tan(PHI)=0
d_PHI=d_PHI*180/pi;if isnan(d_PHI)==1;d_PHI=0;end
d_phi2=d_phi2*180/pi;if isnan(d_phi2)==1;d_phi2=0;end

newEul=Eul*180/pi +[d_phi1 , d_PHI, d_phi2];
newEul=newEul';

%%          OUT        
%The Taylor factor
MT=sum(abs(bn_mean))/step;

ROT=newEul;
%the absolute value of the activity of the different slip systems- adjusted
%so that sum is equal to Taylor factor
OUT=(abs([bn_mean]'/step));


end
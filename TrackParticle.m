function TrackParticle(varargin)
% rajouter les temsp de pause dans les stats


    if nargin==0
        Prop.TUnit='s';
        Prop.XUnit='µm';
        Prop.VMin=.05;%Vmin pour etre considéré en mouvement
        Prop.TMean=3;% ODD :Durée de moyenne glissante pour la vitesse
        Prop.TMin=3;% Durée min d'une période de mouvement
        Prop.ux=0.13;% Unite en [micron]/pxl
        Prop.ut=.3;% Unite en [min]/image

        Prop.CurDir=strcat(cd,'\condition1_Control Impwt');
        Prop.LDirection=cat(2,[1,-1,-1],[1,-1,-1,1,1,1],[-1,-1],1,[-1,1,1,1],[1,1],-1);
                
        % Prop.CurDir=strcat(cd,'\condition2_Imp DQ');
        % Prop.LDirection=cat(2,[-1,1,-1,1,1],[-1,1,1],[-1,-1,1],[1,-1,1,-1,-1,-1],[-1,1,1,1,1,1],[1,1,1,1]);
        % 
        % Prop.CurDir=strcat(cd,'\Trial2_Control_ImpBwt');
        % Prop.LDirection=cat(2,[-1,1],[-1,1],[-1,1],-1,[-1,1,1,1,1],[1,1,1,-1],[1,1,-1],[-1,1],-1);
        % 
        % Prop.CurDir=strcat(cd,'\trial2_ImpBDQ');
        % Prop.LDirection=cat(2,1,[-1,-1,1,-1],[1,-1],[-1,-1],[-1,-1],[-1,1,-1,1,-1,1,1]);      

        Prop.LFile=GetXLSList(Prop.CurDir);
        
        Prop.Ratio=0.0;% Ratio de déplacement pour fixer le mode
        Prop.ListMode={'Stable','Forward','Backward','Oscillator'};
        Prop.iDendrite=1;
        Prop.MotileMap=[.5,.5,.5;0,0,1;1,0,0;1,0,1];

%For testing
       Prop.LFile=Prop.LFile(2+(1:2));
       Prop.LDirection=Prop.LDirection(1:2);
        
        for i1=1:length(Prop.LFile)
            Dendrite(i1)=CreateDendrite(Prop,fullfile(Prop.CurDir,Prop.LFile{i1}),Prop.LDirection(i1));%#ok
        end
        IO.Prop=Prop;
        IO.Graph=Graphs('Initialize',Prop);
        IO.Dendrite=Dendrite;
        
        guidata(gcf,IO);
        
        Graphs('Dendrite',Prop,IO.Graph,Dendrite(1));
        Graphs('MotilePhase',Prop,IO.Graph,Dendrite(1));
        [~,OutXLS]=Graphs('FullData',Prop,IO.Graph,Dendrite);

% print it properly in A4
        [~,FigName]=fileparts(Prop.CurDir);
        set(gcf,'Name',FigName);
        set(gcf,'Units','centimeters');
        % set(gcf,'Position',[1.5,1.5,29.7,21]);
        set(gcf,'Position',[1.5,1.5,24,18]);
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        IO.Graph.TabGpe.SelectedTab=IO.Graph.Tab{1}.hTab;
        print(gcf,FigName,'-dpdf','-r0')
        
%Export to XLS
        XLSName=[FigName,'.xls'];
        if exist(XLSName,'file')
            delete(XLSName);
        end
        warning off
        for i1=1:length(OutXLS)
            if ~isempty(OutXLS(i1).Data)
                class(OutXLS(i1).Data)
                disp(OutXLS(i1))
                if isa(OutXLS(i1).Data,'cell')
                    writecell(OutXLS(i1).Data,XLSName,'Sheet',OutXLS(i1).Name);
                elseif isa(OutXLS(i1).Data,'double')
                    Graph.HFig=figure(9);
                    writematrix(OutXLS(i1).Data,XLSName,'Sheet',OutXLS(i1).Name);
                end
            end
        end
        warning on
    else
        feval(str2func(varargin{1}),varargin{2:nargin});
        return
    end
end

function TmpDendrite=CreateDendrite(varargin)
%    IN : PRop,FileName,Direction
    Prop=varargin{1};
    if nargin==1
        TmpDendrite=ImportXLS(Prop);
    else
        FullName=varargin{2};
        TmpDendrite=ImportXLS(Prop,FullName);
        TmpDendrite.Direction=varargin{3};
        k=length(FullName);
        ShowName=FullName(max(1,k-40):k);
    end

    disp([ShowName,'>GetProjection...']);
    [TmpDendrite.Projection,TmpDendrite.Box]=GetProjection(TmpDendrite.Track,Prop);
    disp([ShowName,'>GetSkel...']);
    [TmpDendrite.CurvDist,TmpDendrite.LineSkel]=GetSkel(TmpDendrite.Projection,TmpDendrite.Direction);
    disp([ShowName,'>GetCurvDist...']);
    TmpDendrite.Track=GetCurvDist(TmpDendrite.Track,TmpDendrite.CurvDist,Prop,TmpDendrite.Box);
    disp([ShowName,'>Analyze...']);
    TmpDendrite=Analyze(TmpDendrite,Prop);
    disp([ShowName,'>StatPerMobileSection...']);
    TmpDendrite.Motile=StatPerMobileSection(TmpDendrite.Track,Prop);
end

%==========================================================================
function Dendrite=ImportXLS(varargin)
%verifie l'ordre et la bonne existence des differents parametres d'entree:
% fichier de donnée, tableau de valeur, numero du groupe
    Prop=varargin{1};
    if nargin==1
        MonPath=Prop.CurDir;
        if ~exist(MonPath,'dir')
            MonPath='';
        end
        [FileName,PathName]=uigetfile('*.xls;*.xlsx','Import File from Excel',MonPath);
        [~,~,Extension]=fileparts(fullfile(PathName,FileName));
        if FileName==0
            return
        end
    else
        [PathName,FileName,Extension]=fileparts(fullfile(varargin{2}));
    end
    if strcmp(Extension,'.xls') || strcmp(Extension,'.xlsx')
        Dendrite.FileName=FileName;
        Dendrite.PathName=PathName;
        [~,~,k]=xlsread(fullfile(PathName,FileName),'Sheet1');
        f=find(strcmp(k(1,:),'Slice n°'));
        if ~isempty(f)
            DataTemp.T=cell2mat(k(2:size(k,1),f));
        end
        f=find(strcmp(k(1,:),'X'));
        if ~isempty(f)
            DataTemp.X=cell2mat(k(2:size(k,1),f));
        end
        f=find(strcmp(k(1,:),'Y'));
        if ~isempty(f)
            DataTemp.Y=cell2mat(k(2:size(k,1),f));
        end
        f=find(strcmp(k(1,:),'Z'));
        if ~isempty(f)
            DataTemp.Z=cell2mat(k(2:size(k,1),f));
        end
    end
    LIsNan=[0;find(isnan(DataTemp.T));length(DataTemp.T)+1];
    Debut=LIsNan(1:length(LIsNan)-1)+1;
    Fin=LIsNan(2:length(LIsNan))-1;
    k=Fin-Debut>1;
    Debut=Debut(k);
    Fin=Fin(k);
    for i1=1:length(Debut)
        Dendrite.Track(i1).TN=DataTemp.T(Debut(i1):Fin(i1));
        Dendrite.Track(i1).X=DataTemp.X(Debut(i1):Fin(i1));
        Dendrite.Track(i1).Y=DataTemp.Y(Debut(i1):Fin(i1));
        if isfield(DataTemp,'Z')
            Dendrite.Track(i1).Z=DataTemp.Z(Debut(i1):Fin(i1));
        end
    end
end

%==========================================================================
function [Projection,Box]=GetProjection(Track,Prop)
    % On construit l'image de projection des tracks de la dendrite
    n=length(Track);
    Box=round([min(cat(1,Track(:).X)/Prop.ux),...
              max(cat(1,Track(:).X)/Prop.ux),...
              min(cat(1,Track(:).Y)/Prop.ux),...
              max(cat(1,Track(:).Y)/Prop.ux)]);
          
    Projection=zeros(Box(4)-Box(3)+3,Box(2)-Box(1)+3);
    Sze=size(Projection);
    % On construit l'image de projection de la dendrite
    for i1=1:n
        x=Track(i1).X/Prop.ux-Box(1)+2;
        y=Track(i1).Y/Prop.ux-Box(3)+2;
        s=[0;cumsum(abs(diff(x+1i*y)))];
        Nz=[true;abs(diff(x+1i*y))>0];   
    % on inscrit les segments
        if max(s)>0
            k=unique(round(interp1(s(Nz),x(Nz)+1i*y(Nz),0:max(s))));   
        else
            k=x(1)+1i*y(1);
        end
        L=sub2ind(Sze,imag(k),real(k));
        Projection(L)=Projection(L)+1;
    % on inscrit les points redondants pour renforcer leur importance par la
    % suite
        k2=round(x(~Nz)+1i*y(~Nz));
        L2=sub2ind(Sze,imag(k2),real(k2));
        Projection(L2)=Projection(L2)+1;
    end
    
% On relis les blocs si besoin
    test=true;
    while test
    % On recherche les 2 parties les plus proches
        BWL=bwlabel(Projection>0);
        MBWL=max(BWL(:));
        Score=[0,0,inf,0,0];%[Label1, Label2, Distance, Index1,Index2]
        for L1=1:MBWL
            MD=bwdist(BWL==L1);
            f=find((BWL>0) & (BWL~=L1));
            [d,k]=min(MD(f));
            ind2=f(k);
            if Score(3)>d
                MD2=bwdist(BWL==BWL(ind2));
                f=find(BWL==L1);
                [~,k]=min(MD2(f));
                ind1=f(k);
                Score=[L1,BWL(ind2),d,ind1,ind2];
            end        
        end

    % on inscrit le segment intermédiaire
        if ~isinf(Score(3))
            [x,y]=ind2sub(Sze,Score(4:5));
            s=[0;abs(diff(x+1i*y))];
            IntermP=unique(round(interp1(s,x+1i*y,0:max(s)))).';
            IntermP=setdiff(IntermP,x+1i*y);
            L=sub2ind(Sze,real(IntermP),imag(IntermP));
            Projection(L)=1;
        else
            test=false;
        end
    end
end

%==========================================================================
function [CurvDist,LineSkel]=GetSkel(Projection,Direction)
% On vire les points qui ne sont pas indispensable en commençant par les
% moins intenses qui ne sont pas des extremites pour garder un unique objet
% Le skel est CurvDist>=0
% TODO : améliorer le skel on fittant un peu la courbe
    Sze=size(Projection);
    Skel=bwmorph(Projection>0,'skel');
    InvDist=.5./(max(.5,bwdist(Skel)));
    Tips=(conv2(double(Skel),ones(3,3),'same').*double(Skel))==2;
    if sum(Tips(:))<2
        Tips=Tips | (conv2(double(Skel),ones(3,3),'same').*double(Skel))==3;
        disp('No Extremal Tips')
    end
    [X,Y]=find(Tips);
    Z=X+1i*Y;
    n=length(Z);
    M=abs(repmat(Z,1,n)-repmat(Z,1,n).');
    [~,k]=max(sum(M));
    IndExtr(1)=sub2ind(Sze,X(k),Y(k));
    [~,k]=max(M(k,:));
    IndExtr(2)=sub2ind(Sze,X(k),Y(k));
    if Direction==-1
        IndExtr=fliplr(IndExtr);
    end
    
    LPoint=find(Projection>0);
    LPoint=setdiff(LPoint,IndExtr);
    LPoint(:,2)=Projection(LPoint)+InvDist(LPoint);%Le Invdist est pour secreger les points de même fréquence
    LPoint=sortrows(LPoint,-2);
    Skeleton=Projection;
    for i1=length(LPoint):-1:1
        M2=Skeleton>0;
        M2(LPoint(i1,1))=false;
        k=bwlabel(M2);
        if k(IndExtr(1))==k(IndExtr(2))
            Skeleton(LPoint(i1,1))=0;
        end
    end
    Skeleton=Skeleton>0;   

% On calcule la distance curviligne
    CurvDist=zeros(Sze);
    CurvDist(IndExtr(1))=1;
    test=true;
    kern=[1.4,1,1.4;1,0,1;1.4,1,1.4];
    d=0;
    while test
        M=(conv2(double(CurvDist>0),kern,'same').*Skeleton>0) & (CurvDist==0);
        d=d+max(M(:));
        CurvDist=CurvDist+(M>0)*d;
        test=any(Skeleton(:) & ~(CurvDist(:)>0));
    end
    CurvDist=CurvDist-1;
    
    k=find(CurvDist>=0);
    k(:,2)=CurvDist(k);
    k=sortrows(k,2);
    [x,y]=ind2sub(Sze,k(:,1));
    LineSkel=[x,y,k(:,2)];
end

%==========================================================================
function Track=GetCurvDist(Track,CurvDist,Prop,Box)
% On écrit la distance curviligne pour chaque point [µm]
    n=length(Track);
    Sze=size(CurvDist);
    LSkel=find(CurvDist>=0);
    for i1=1:n
        Track(i1).T=Track(i1).TN*Prop.ut;
        Track(i1).XN=Track(i1).X/Prop.ux;
        Track(i1).YN=Track(i1).Y/Prop.ux;
        
        x=round(Track(i1).XN-Box(1)+2);
        y=round(Track(i1).YN-Box(3)+2);
        m=length(x);
        s=zeros(m,1);
        for i2=1:m
            M=false(Sze);
            M(y(i2),x(i2))=1;
            D=bwdist(M);
            [~,ind]=min(D(LSkel));
            s(i2)=CurvDist(LSkel(ind));
        end
%         s=s(1)+round(cumsum(.6*sin((1:length(s))/10)).');
        Track(i1).SN=s;
        Track(i1).S=s*Prop.ux;
    end
end

%==========================================================================
function Dendrite=Analyze(Dendrite,Prop)
% Analyze current dendrite
    n=length(Dendrite.Track);
    v=cell(n,1);
    v2=cell(n,1);
    for i1=1:n
        v{i1}=diff(Dendrite.Track(i1).S./Prop.ut);        
        kn=length(v{i1,1});
        v2{i1,1}=interp1((-.5:kn+.5),[v{i1}(1);v{i1};v{i1}(kn)],0:kn,'linear','extrap').';
    end
    Dendrite.Track=SetMode(Dendrite.Track,Prop);
%    [hy,LVX]=hist(cell2mat(v),-1:0.1:1);
    edge=-1:0.1:1;
    hy=histcounts(cell2mat(v),edge);
    LVX=0.5*(edge(1:end-1)+edge(2:end));
    LVY=(hy)/sum(hy)*100;
    Classe=2-uint8(LVX<=-Prop.VMin)+uint8(LVX>=Prop.VMin);
    Dendrite.HistoV.X=LVX;
    Dendrite.HistoV.Y=LVY;
    Dendrite.HistoV.Classe=Classe;
    Dendrite.HistoV.Color=Prop.MotileMap([2,1,3],:);
    
    Distrib=zeros(1,length(Prop.ListMode));
    ListPerMode=cell(size(Distrib));
    for i1=1:length(Prop.ListMode)
        k=~cellfun(@isempty,strfind({Dendrite.Track(:).Mode},Prop.ListMode{i1}));
        Distrib(i1)=sum(k);
        ListPerMode{i1}=find(k);
    end
    Dendrite.Distrib=Distrib;
    Dendrite.ListPerMode=ListPerMode;
    
%=======================
    function Track=SetMode(Track,Prop)
        for i0=1:length(Track)
            [Track(i0).Filter,Track(i0).V]=FilterShortMov(Track(i0).S,Prop,false);
             
            if sum(double(Track(i0).Filter==3))/length(Track(i0).Filter)>Prop.Ratio
                if sum(double(Track(i0).Filter==1))/length(Track(i0).Filter)>Prop.Ratio
                    Track(i0).Mode=Prop.ListMode{4};
                else
                    Track(i0).Mode=Prop.ListMode{2};
                end
            else
                if sum(double(Track(i0).Filter==1))/length(Track(i0).Filter)>Prop.Ratio
                    Track(i0).Mode=Prop.ListMode{3};
                else
                    Track(i0).Mode=Prop.ListMode{1};
                end
            end
            
% get pauses
            Lf=length(Track(i0).Filter);
            BWL=bwlabel(Track(i0).Filter==2);
            if BWL(1)>0
                BWL(BWL==BWL(1))=0;
            end
            if BWL(Lf)>0
                BWL(BWL==BWL(Lf))=0;
            end
            k=regionprops(BWL>0,'Area');
            Track(i0).ListPause=cat(1,k(:).Area);
            Track(i0).ListPauseN=Track(i0).ListPause/Prop.ut;       
        end
    end
%=======================
    function [filter,V]=FilterShortMov(S,Prop,Show)
% on classe en 3 classes [Back,Stable,Forward]
        [FFor,V]=GetSteps(S,Prop,1);
        [FBack,~]=GetSteps(S,Prop,-1);
        filter=2-FBack+FFor;
        if Show
            figure(8)
            clf
            subplot 221
                plot(S,'-ok')
                title('Track')
            subplot 222
                plot(diff(S),'-dk');
                hold on
                plot(SlidingMean(diff(S),Prop.TMean),'-og');
                plot([1,length(S)],Prop.VMin*[1,1],':m');
                plot([1,length(S)],-Prop.VMin*[1,1],':m');
                plot([1,length(S)],[0,0],':m');
                hold off
                legend({'DS=Diff(S)','SlidingMean(DS)','Speed Thresh','Speed Thresh','no speed'})
            subplot 223
                plot(FFor+6,'-go')
                hold on
%                plot(SFor+4,':gs')
                plot(FBack-6,'-ro')
%                plot(SBack-4,':rs')
                plot(filter-2,'-md')
                hold off
            subplot 224
                plot(FFor+2,'-go')
                hold on
                plot(FBack-3,'-ro')
                plot(filter==2,'-md')
                hold off
            figure(9)
        end
    end
%=======================
   function [Select,V]=GetSteps(S,Prop,Orient)
        V=MeanSpeed(S,Prop.TMean)./Prop.ut;
        n=length(S);

        if Orient==1
            Select=V>=Prop.VMin;
        elseif Orient==-1
            Select=V<=-Prop.VMin;
        end        
        s0=Select;
        Select=[0;Select;0];
        kplus=find(diff(double(Select))==1);
        kmoins=find(diff(double(Select))==-1);
        for j1=1:length(kplus)
            kstop=find(kmoins>kplus(j1),1);
            if ~isempty(kstop)
                if (kmoins(kstop)-kplus(j1))<Prop.TMin
                    Select(kplus(j1)+1:kmoins(kstop))=0;
                end
            end
        end
        Select=Select(2:n+1);
        % figure(10)
        % subplot 411
        % plot(S,'bo-');
        % xlim([1 length(S)])
        % subplot 412
        % plot(V,'bo-');
        % xlim([1 length(S)])
        % subplot 413
        % plot(s0,'bo-');
        % xlim([1 length(S)])
        % subplot 414
        % plot(Select,'bo-');
        % xlim([1 length(S)])
   end
%=======================
    function m=SlidingMean(x,L)
        x2=[repmat(x(1),floor((L-1)/2),1);x;repmat(x(1),L-1-floor((L-1)/2),1)];
        m=conv2(x2,ones(L,1)/L,'valid');
    end
%=======================
    function v=MeanSpeed(y,L)
        n=length(y);
        delta=floor((L-1)/2);
        y2=[repmat(y(1),delta,1);y;repmat(y(length(y)),L-1-delta,1)];       
        if L>1
            kernel=zeros(1,L);
            kernel(1)=-1/(2*delta);
            kernel(L)=1/(2*delta);
            v=conv(y2,-kernel,'valid');
        else
            v=diff([y2(1);y2]);
        end
    end
end

%==========================================================================
function Motile=StatPerMobileSection(Track,Prop)
    n=length(Track);
    Motile.Forw.V=[];
    Motile.Forw.T=[];
    Motile.Forw.L=[];
    Motile.Back=Motile.Forw;
    for i1=1:n
        iTrack=Track(i1);
        k=[2;iTrack.Filter;2];
        
        F1=find(diff(k==3)==1);
        F2=find(diff(k==3)==-1);
        for i2=1:length(F1)
            dt=(F2(i2)-F1(i2)+1)*Prop.ut;
            Motile.Forw.V=cat(1,Motile.Forw.V,iTrack.V(F1(i2):F2(i2)-1));
            Motile.Forw.L=cat(1,Motile.Forw.L,cumsum(abs(diff(iTrack.S(F1(i2):F2(i2)-1)))));
            Motile.Forw.T=cat(1,Motile.Forw.T,dt);
        end
        % figure(10)
        % subplot 411
        % plot(iTrack.Filter,'b-o')
        % xlim([1,length(iTrack.Filter)])
        % subplot 412
        % plot(diff(iTrack.S),'b-x')
        % hold on
        % plot([1 length(iTrack.S)-1;1 length(iTrack.S)-1].',[.05,.05;-.05,-.05].','k:')
        % hold off
        % xlim([1,length(iTrack.Filter)])
        % subplot 413
        % plot(diff(iTrack.SN),'b-o')
        % xlim([1,length(iTrack.Filter)])
        % subplot 414
        % plot(iTrack.Filter)
        % xlim([1,length(iTrack.Filter)])
        % pause(.1)
        
        F1=find(diff(k==1)==1);
        F2=find(diff(k==1)==-1);
        for i2=1:length(F1)
            dt=(F2(i2)-F1(i2)+1)*Prop.ut;
            Motile.Back.V=cat(1,Motile.Back.V,iTrack.V(F1(i2):F2(i2)-1));
            Motile.Back.L=cat(1,Motile.Back.L,cumsum(abs(diff(iTrack.S(F1(i2):F2(i2)-1)))));
            Motile.Back.T=cat(1,Motile.Back.T,dt);
        end
    end
    
end

%==========================================================================
function UpdateExample(iCurve)%#ok
    IO=guidata(gcf);
    Graphs('DendriteExample',IO.Prop,IO.Graph,IO.Dendrite(IO.Prop.iDendrite),iCurve);    
    IO.Graph.TabGpe.SelectedTab=IO.Graph.Tab{2}.hTab;
    set(gcf,'Name',IO.Dendrite(IO.Prop.iDendrite).FileName);
end


function ListXLS=GetXLSList(mydir)
    L=dir(mydir);
    ListXLS={};
    for i1=3:length(L)
        [~,~,ext]=fileparts(fullfile(mydir,L(i1).name));
        if strcmp(ext,'.xls') || strcmp(ext,'.xlsx')
            ListXLS=cat(1,ListXLS,L(i1).name);
        end
    end
end
 


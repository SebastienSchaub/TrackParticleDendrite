function  [Graph,OutXLS]=Graphs(varargin)
    if nargin==0
       TrackParticle;
       return 
    end
    Mode=varargin{1};
    Prop=varargin{2};
    switch Mode
        case 'Initialize'
            Graph.HFig=figure(9);
            clf;drawnow;
            Graph.TabGpe= uitabgroup(Graph.HFig,'Position',[0,0,1,1]);
            Graph.Tab{1}.hTab= uitab(Graph.TabGpe,'Title','Full Data','BackgroundColor',[.7,.8,.8]);
            Graph.Tab{1}.hAxe(1)=axes('Parent',Graph.Tab{1}.hTab,'OuterPosition',[0,.5,.5,.5]);
            Graph.Tab{1}.hAxe(2)=axes('Parent',Graph.Tab{1}.hTab,'OuterPosition',[0,0,.5,.5]);
            Graph.Tab{1}.hAxe(3)=axes('Parent',Graph.Tab{1}.hTab,'OuterPosition',[.5,.75,.5,.25]);
            Graph.Tab{1}.hAxe(4)=axes('Parent',Graph.Tab{1}.hTab,'OuterPosition',[.5,.5,.5,.25]);
            Graph.Tab{1}.hAxe(5)=axes('Parent',Graph.Tab{1}.hTab,'OuterPosition',[.5,.25,.5,.25]);
            Graph.Tab{1}.hAxe(6)=axes('Parent',Graph.Tab{1}.hTab,'OuterPosition',[.5,.0,.5,.25]);
            Graph.Tab{2}.hTab= uitab(Graph.TabGpe,'Title','Full Dendrite','BackgroundColor',[.8,.7,.8]);
            Graph.Tab{2}.hAxe(1)=axes('Parent',Graph.Tab{2}.hTab,'OuterPosition',[0,.5,.5,.5]);
            Graph.Tab{2}.hAxe(2)=axes('Parent',Graph.Tab{2}.hTab,'OuterPosition',[0,0,.5,.5]);
            Graph.Tab{2}.hAxe(3)=axes('Parent',Graph.Tab{2}.hTab,'OuterPosition',[.5,.75,.5,.25]);
            Graph.Tab{2}.hAxe(4)=axes('Parent',Graph.Tab{2}.hTab,'OuterPosition',[.5,.5,.5,.25]);
            Graph.Tab{2}.hAxe(5)=axes('Parent',Graph.Tab{2}.hTab,'OuterPosition',[.5,.25,.5,.25]);
            Graph.Tab{2}.hAxe(6)=axes('Parent',Graph.Tab{2}.hTab,'OuterPosition',[.5,.0,.5,.25]);

            Graph.Tab{3}.hTab= uitab(Graph.TabGpe,'Title',Prop.ListMode{1},'BackgroundColor',[.8,.8,.8]);
            Graph.Tab{3}.hAxe(1)=axes('Parent',Graph.Tab{3}.hTab,'OuterPosition',[0,0,1,1]);

            Graph.Tab{4}.hTab= uitab(Graph.TabGpe,'Title',Prop.ListMode{2},'BackgroundColor',[.8,.7,.7]);
            Graph.Tab{4}.hAxe(1)=axes('Parent',Graph.Tab{4}.hTab,'OuterPosition',[0,0,1,1]);

            Graph.Tab{5}.hTab= uitab(Graph.TabGpe,'Title',Prop.ListMode{3},'BackgroundColor',[.7,.8,.7]);
            Graph.Tab{5}.hAxe(1)=axes('Parent',Graph.Tab{5}.hTab,'OuterPosition',[0,0,1,1]);

            Graph.Tab{6}.hTab= uitab(Graph.TabGpe,'Title',Prop.ListMode{4},'BackgroundColor',[.8,.8,.7]);
            Graph.Tab{6}.hAxe(1)=axes('Parent',Graph.Tab{6}.hTab,'OuterPosition',[0,0,1,1]);

            
%==========================================================================
        case 'Dendrite'
            if nargin==2
                IO=guidata(gcf);
                Graph=IO.Graph;
                Prop=IO.Prop;
                IO.Prop.iDendrite=varargin{2};
                Dendrite=IO.Dendrite(IO.Prop.iDendrite);
                guidata(gcf,IO);
            else
                Graph=varargin{3};
                Dendrite=varargin{4};
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image de la dendrite 
            set(gcf,'CurrentAxe',Graph.Tab{2}.hAxe(1))
                s=Prop.ux;
                x0=(Dendrite.Box(1)-2)*s;
                y0=(Dendrite.Box(3)-2)*s;
                y=(Dendrite.Box(3):Dendrite.Box(4))*s;
                x=(Dendrite.Box(1):Dendrite.Box(2))*s;
                imagesc(x,y,Dendrite.Projection);
                hold on
                plot(Dendrite.LineSkel(:,2)*s+x0,Dendrite.LineSkel(:,1)*s+y0,'r.-')
% plot de la flêche
                n=size(Dendrite.LineSkel,1);
                v=.5*diff(Dendrite.LineSkel(n+[-2,0],2))+1i*diff(Dendrite.LineSkel(n+[-2,0],1));
                v=v./abs(v)*2;
                p=Dendrite.LineSkel(n,2)+1i*Dendrite.LineSkel(n,1);
                v2=p+cumsum(v.*[0,.5*exp(1i*pi/2),exp(-1i*pi/6),exp(1i*-5*pi/6),.5*exp(1i*pi/2)]);
                plot(v2*s+x0+1i*y0,'r-','Linewidth',2)

                Graphs('DendriteExample',Prop,Graph,Dendrite,1);
                hold off
                axis image;
                axis ij
                set(gca,'Color',[0,0,0])
                xlabel(Prop.XUnit);
                ylabel(Prop.XUnit);
                
                colormap(gca,gray(255))
                title('Projection of dendrite')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogramme des différents modes PAR DENDRITE
            set(gcf,'CurrentAxe',Graph.Tab{2}.hAxe(2))
            M=Dendrite.Distrib./sum(Dendrite.Distrib)*100;
            for i1=1:4
                bar(i1,M(i1),'FaceColor',Prop.MotileMap(i1,:))
                hold on
            end
            hold off
            set(gca,'XTick',1:4,'XTickLabel',Prop.ListMode)
            ylabel('% vesicle')
            title(['Distribution of Modes % ',mat2str(Dendrite.Distrib)])
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GRAPH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogramme de "Curvilign Speed" i.e. des biais de vitesse
            set(gcf,'CurrentAxe',Graph.Tab{2}.hAxe(3))
            for i1=1:3
                k=Dendrite.HistoV.Classe==i1;
                bar(Dendrite.HistoV.X,Dendrite.HistoV.Y.*k,...
                    'FaceColor',Dendrite.HistoV.Color(i1,:))
                xlim([min(Dendrite.HistoV.X) max(Dendrite.HistoV.X)])
                hold on
            end
            hold off
            title('Curvilign Speed')
            xlabel(['Speed ',Prop.XUnit,'/',Prop.TUnit])
            ylabel('% vesicle')

%Tracing individual vesicle path
            for i1=1:length(Prop.ListMode)
                ListPerMode=Dendrite.ListPerMode{i1};
                col=hsv(length(ListPerMode));
                set(gcf,'CurrentAxe',Graph.Tab{2+i1}.hAxe(1));
                delete(allchild(Graph.Tab{2+i1}.hAxe(1)));
                mx=max(cat(1,Dendrite.Track(ListPerMode(:)).X));
                for i2=1:length(ListPerMode)
                    s=Dendrite.Track(ListPerMode(i2)).S;
                    t=(1:length(s))*Prop.ut;
                    f1=Dendrite.Track(ListPerMode(i2)).Filter==1;
                    f3=Dendrite.Track(ListPerMode(i2)).Filter==3;
                    txt=['T',num2str(ListPerMode(i2))];
                    plot(t,s,'-','Color',col(i2,:),'Tag',txt,'ButtonDownFcn',['TrackParticle(''UpdateExample'',',num2str(ListPerMode(i2)),')']);            
                    hold on
                    text(t(length(t))+.03*mx,s(length(s)),txt,'Color',col(i2,:))
                    plot(t(f1),s(f1),'s','Color',col(i2,:),'ButtonDownFcn',['TrackParticle(''UpdateExample'',',num2str(ListPerMode(i2)),')']);
                    plot(t(f3),s(f3),'d','Color',col(i2,:),'ButtonDownFcn',['TrackParticle(''UpdateExample'',',num2str(ListPerMode(i2)),')']);           
                end
                hold off
                title(Prop.ListMode{i1})
                xlabel(['Time [',Prop.TUnit,']'])
                ylabel (['Curviligne distance [',Prop.XUnit,']'])                
                set(gca,'Color',.1*[1,1,1],'YLim',[0,Dendrite.LineSkel(end,3)*Prop.ux])                
                Graph.TabGpe.SelectedTab=Graph.Tab{2}.hTab;
            end
%==========================================================================
        case 'DendriteExample'
            Graph=varargin{3};
            iDendrite=varargin{4};
            iCurve=varargin{5};
            set(gcf,'CurrentAxe',Graph.Tab{2}.hAxe(1))
            delete(findobj(Graph.Tab{2}.hAxe(1),'Tag','Example'));
            hold on
%             plot(iDendrite.Track(iCurve).X/Prop.ux-iDendrite.Box(1)+2,iDendrite.Track(iCurve).Y/Prop.ux-iDendrite.Box(3)+2,'g.-','Tag','Example');
            plot(iDendrite.Track(iCurve).X,iDendrite.Track(iCurve).Y,'g.-','Tag','Example');
            hold off
        
%==========================================================================
        case 'MotilePhase'
            nargin
            if nargin==5
                OutXLS=varargin{5};
            end
            Graph=varargin{3};
            Dendrite=varargin{4};
            
            
            nd=length(Dendrite);
            LVB=[];%Vitesse Backward
            LTB=[];%Temps Backward
            LLB=[];%Longueur Backward
            LVF=[];
            LTF=[];
            LLF=[];
            for i1=1:nd
                LVB=cat(1,LVB,-Dendrite(i1).Motile.Back.V);
                LTB=cat(1,LTB,Dendrite(i1).Motile.Back.T);
                LLB=cat(1,LLB,Dendrite(i1).Motile.Back.L);
                LVF=cat(1,LVF,Dendrite(i1).Motile.Forw.V);
                LTF=cat(1,LTF,Dendrite(i1).Motile.Forw.T);
                LLF=cat(1,LLF,Dendrite(i1).Motile.Forw.L);
            end
            iTab=1+uint8(nd==1);
            
% Speed histogram
            [~,edge]=histcounts(cat(1,LVF,LVB),10);
            xV=0.5*(edge(1:end-1)+edge(2:end));
            hV=[histcounts(LVF,edge);histcounts(LVB,edge)];
            mV=[mean(LVF),mean(LVB)];

% Length histogram
            [~,edge]=histcounts(cat(1,LLF,LLB),10);
            xL=0.5*(edge(1:end-1)+edge(2:end));
            hL=[histcounts(LLF,edge);histcounts(LLB,edge)];
            mL=[mean(LLF),mean(LLB)];

% Duration histogram
            [~,edge]=histcounts(cat(1,LTF,LTB),10);
            xT=0.5*(edge(1:end-1)+edge(2:end));
            hT=[histcounts(LTF,edge);histcounts(LTB,edge)];
            mT=[mean(LTF),mean(LTB)];

            set(gcf,'CurrentAxe',Graph.Tab{iTab}.hAxe(4))
                h=bar(xV,hV.',1);
                h(1).FaceColor= [0.1,0.1,1];
                h(2).FaceColor= [1,0.1,0.1];
                set(gca,'XTick',xV,'XTickLabelRotation',90,'YScale','log')
                title(['Speed in Motile phase:(',num2str(mV,'%2.2f,'),') [',Prop.XUnit,'/',Prop.TUnit,']'])
                legend({'Forward','Backward'})
                colormap(gca,jet(4))
                ylabel('# measure')
                grid on
                
            set(gcf,'CurrentAxe',Graph.Tab{iTab}.hAxe(5))
                h=bar(xL,hL.',1);
                h(1).FaceColor= [0.1,0.1,1];
                h(2).FaceColor= [1,0.1,0.1];
                set(gca,'XTick',xL,'XTickLabelRotation',90,'YScale','log')
                title(['Length of Motile phase:(',num2str(mL,'%2.2f,'),') [',Prop.XUnit,']'])
                legend({'Forward','Backward'})
                colormap(gca,jet(4))
                ylabel('# event')
                grid on

            set(gcf,'CurrentAxe',Graph.Tab{iTab}.hAxe(6))
                h=bar(xT,hT.',1);
                h(1).FaceColor= [0.1,0.1,1];
                h(2).FaceColor= [1,0.1,0.1];
                set(gca,'XTick',xT,'XTickLabelRotation',90,'YScale','log')
                title(['Duration of Motile phase:(',num2str(mT,'%2.2f,'),') [',Prop.TUnit,']'])
                legend({'Forward','Backward'})
                ylabel('# event')
                colormap(gca,jet(4))
                grid on
                
            if exist('OutXLS','var')
                OutXLS=ExportXLS(OutXLS,'SpeedMotPh',[xV.',hV.'],mV,LVF,LVB,Prop);
                %OutXLS=ExportXLS(OutXLS,'LengthMotPh',[xL.',hL.'],mL,LLB,LLF,Prop);
                %OutXLS=ExportXLS(OutXLS,'DurationMotPh',[xT.',hT.'],mT,LTB,LTF,Prop);
            end
%==========================================================================
        case 'FullData'
            OutXLS=struct;
            Graph=varargin{3};
            Dendrite=varargin{4};        
            
            MDistrib=cat(1,Dendrite(:).Distrib);
            [NL,NC]=size(MDistrib);

            
            set(gcf,'CurrentAxe',Graph.Tab{1}.hAxe(1))
            M0=sum(MDistrib);
            M=M0./sum(M0)*100;
            for i1=1:4
                bar(i1,M(i1),'FaceColor',Prop.MotileMap(i1,:))
                hold on
            end
            hold off
            set(gca,'XTick',1:4,'XTickLabel',Prop.ListMode)
            ylabel('# vesicles');
            title(['Distribution of Modes % ',mat2str(M0)])

            set(gcf,'CurrentAxe',Graph.Tab{1}.hAxe(2))
            for i1=1:NL
                for i2=NC:-1:1
                   h=barh(i1,sum(MDistrib(i1,1:i2)),'stacked','FaceColor',Prop.MotileMap(i2,:));
                   h.ButtonDownFcn=['Graphs(''Dendrite'',',num2str(i1),');'];
                   hold on
                end
            end
            hold off
            set(gca,'YTick',1:length(Dendrite),'YTickLabel',{Dendrite(:).FileName});
            xlabel('# vesicles');
            title('Distribution of Modes')
            colormap(gca,Prop.MotileMap)

%            OutXLS=ExportXLS(OutXLS,'DistribMode',MDistrib,Prop,{Dendrite(:).FileName});

%==== GRAPH ===============================================================
% Courbes des temps de pause 
            set(gcf,'CurrentAxe',Graph.Tab{1}.hAxe(3))
            nd=length(Dendrite);
            LPause=[];
            for i1=1:nd
                LPause=cat(1,LPause*Prop.ut,cat(1,Dendrite(i1).Track(:).ListPause));
            end

            [~,edge]=histcounts(LPause,20);
            hx=0.5*(edge(1:end-1)+edge(2:end));
            hy=histcounts(LPause,edge);
            bar(hx,hy,'FaceColor',0.5*[0.5,1,0.5]);
            
            title(['Pause duration=',num2str(mean(LPause),3),Prop.TUnit])
%            OutXLS=ExportXLS(OutXLS,'PauseDuration',LPause); % ERROR
            
            [~,OutXLS]=Graphs('MotilePhase',Prop,Graph,Dendrite,OutXLS);            
        otherwise
            Graph=[];
            
    end
end

%==========================================================================
function Out=ExportXLS(varargin)
    Out=varargin{1};
    Mode=varargin{2};
    switch Mode
        case 'DistribMode'
% TO CHECK            
            Data=varargin{3};
            Prop=varargin{4};
            ListFileName=varargin{5}.';
            
            Tmp=cat(1,{'GLOBAL','','','',''},cat(2,{''},Prop.ListMode),cat(2,{'Sum'},num2cell(sum(Data))));
            Tmp=cat(1,Tmp,{'per Dendrite','','','',''},cat(2,{''},Prop.ListMode),cat(2,ListFileName,num2cell(Data)));
            Out.Name=Mode;
            Out.Data=Tmp;
        case 'PauseDuration'
% TO CHECK BUG ?           
            Data=varargin{3};
            n=max(1,round(sqrt(length(Data))/10))*10;
            [hy,hx]=hist(Data,n);
            Tmp=cat(1,{'img','Distrib Pause Duration'},cat(2,num2cell(hx.'),num2cell(hy.')),...
                {'MEAN',mean(Data)});
            Out(2).Name=Mode;
            Out(2).Data=Tmp;
            Out(3).Name='RawPauseData';
            Out(3).Data=Data;
            
        case 'SpeedMotPh'%checked
            Data=varargin{3};%[histox,histoVForW,histoVBackW]
            Prop=varargin{7};
            Tmp=cat(1,{['V [',Prop.XUnit,'/',Prop.TUnit,']'],'# Forward','# Backward'},num2cell(Data),...
                      cat(2,'MEAN',num2cell(varargin{4})));
            Out(4).Name=Mode;
            Out(4).Data=Tmp;
            
            if length(varargin{5})>length(varargin{6})
                Tmp2=num2cell(varargin{5});
                Tmp2(1:length(varargin{6}),2)=num2cell(varargin{6});
            else
                Tmp2=num2cell(varargin{6});
                Tmp2(1:length(varargin{5}),2)=num2cell(varargin{5});
                Tmp2=fliplr(Tmp2);
            end
            Tmp2=cat(1,{'Forward','Backward'},Tmp2);
            Out(5).Name=[Mode,'Raw'];
            Out(5).Data=Tmp2;
        case 'LengthMotPh'
% TO CHECK            
            Data=varargin{3};
            Tmp=cat(1,{'M[µm]','# Back','# Back'},num2cell(Data),...
                      cat(2,'MEAN',num2cell(varargin{4})));
            Out(6).Name=Mode;
            Out(6).Data=Tmp;
            
            if length(varargin{5})>length(varargin{6})
                Tmp2=num2cell(varargin{5});
                Tmp2(1:length(varargin{6}),2)=num2cell(varargin{6});
            else
                Tmp2=num2cell(varargin{6});
                Tmp2(1:length(varargin{5}),2)=num2cell(varargin{5});
                Tmp2=fliplr(Tmp2);
            end
            Tmp2=cat(1,{'Backward','Forward'},Tmp2);
            Out(7).Name=[Mode,'Raw'];
            Out(7).Data=Tmp2;
        
        case 'DurationMotPh'
% TO CHECK BUG ?            
            Data=varargin{3};
            Tmp=cat(1,{'T[s]','# Back','# Back'},num2cell(Data),...
                      cat(2,'MEAN',num2cell(varargin{4})));
            Out(8).Name=Mode;
            Out(8).Data=Tmp;
            
            if length(varargin{5})>length(varargin{6})
                Tmp2=num2cell(varargin{5});
                Tmp2(1:length(varargin{6}),2)=num2cell(varargin{6});
            else
                Tmp2=num2cell(varargin{6});
                Tmp2(1:length(varargin{5}),2)=num2cell(varargin{5});
                Tmp2=fliplr(Tmp2);
            end
            Tmp2=cat(1,{'Backward','Forward'},Tmp2);
            Out(9).Name=[Mode,'Raw'];
            Out(9).Data=Tmp2;
    end
end
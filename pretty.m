% TOOLS FOR GENERATING GOOD FIGURES

classdef pretty

    methods(Static)
 
function [MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=default()

 MarkerSize=20;   
 FontName='Helvetica';   
 TickFontSize=22;
 LabelFontSize=22;
 TitleFontSize=22;
 set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
 set(groot, 'defaultLegendInterpreter','latex'); 

end

%% Generates pretty figures
function plot(x,y,varargin)

[MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default();    
    
    if nargin==2
plot(x,y,'.k','MarkerSize',MarkerSize)
    elseif nargin>2
plot(x,y,varargin{:})      
    end
    
set(gca,'FontName',FontName,'FontSize',TickFontSize)
ax=gca;
ax.FontSize=TickFontSize;

% ,'FontWeight','bold')
%set(gcf,'OuterPosition',[0,0,1024,720]/1.2,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','A5')
%set(gcf,'OuterPosition',[0,0,720,620]/1.2,'PaperPositionMode','auto')
    
end

%% 
function loglog(x,y,varargin)

[MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default(); 

    if nargin==2
loglog(x,y,'.k','MarkerSize',MarkerSize)
    elseif nargin==3
        
loglog(x,y,varargin{:},'MarkerSize',MarkerSize)      
    end
set(gca,'FontName',FontName,'FontSize',TickFontSize)  % ,'FontWeight','bold')
%set(gcf,'OuterPosition',[0,0,720,620]/1.2,'PaperPositionMode','auto')

end



%%

function semilogx(x,y,varargin)
    
[MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default(); 

    if nargin==2
semilogx(x,y,'.k','MarkerSize',MarkerSize)
    elseif nargin==3
semilogx(x,y,varargin{:},'MarkerSize',MarkerSize)      
    end
set(gca,'FontName',FontName,'FontSize',TickFontSize)  % ,'FontWeight','bold')
% set(gcf,'OuterPosition',[0,0,1024,720]/1.2,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','A5')
%set(gcf,'OuterPosition',[0,0,720,620]/1.2,'PaperPositionMode','auto')
end

%%
function errorbar(x,y,e,varargin)
    
    [MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default();

    if nargin==2
errorbar(x,y,'.k','MarkerSize',MarkerSize)
    elseif nargin==3
errorbar(x,y,e,varargin{:},'MarkerSize',MarkerSize)      
    end
set(gca,'FontName',FontName,'FontSize',TickFontSize)  % ,'FontWeight','bold')
%set(gcf,'OuterPosition',[0,0,1024,720]/1.2,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','A5')

end
%%
% option=1 ---> plot
% option=2 ---> loglog
% option=3 ---> semilogx

function subplot(x,y, option, varargin)

    [MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default();
    
    switch option
        case 1 
    if nargin==3
plot(x,y,'.k','MarkerSize',20)
    elseif nargin==4
plot(x,y,varargin{:},'MarkerSize',MarkerSize)      
    end
        case 2
    if nargin==3
loglog(x,y,'.k','MarkerSize',MarkerSize)
    elseif nargin==4
loglog(x,y,varargin{:},'MarkerSize',MarkerSize)      
    end
    
        case 3
    if nargin==3
semilogx(x,y,'.k','MarkerSize',MarkerSize)
    elseif nargin==4
semilogx(x,y,varargin{:},'MarkerSize',MarkerSize)      
    end    
    end
    set(gca,'FontName',FontName,'FontSize',TickFontSize)  % ,'FontWeight','bold')
   
end
%%
function suberrorbar(x,y,e,varargin)
    
[MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default();
    
    if nargin==3
errorbar(x,y,'.k','MarkerSize',MarkerSize)
    elseif nargin==4
errorbar(x,y,e,varargin{:},'MarkerSize',MarkerSize)      
    end
set(gca,'FontName',FontName,'FontSize',TickFontSize)  % ,'FontWeight','bold')
end


%% 
function prepsubplot()

%set(gcf,'WindowStyle','Normal','OuterPosition',[0,0,1024,720]/1.2,'PaperOrientation','landscape','PaperPositionMode','auto','PaperType','A5');
end
%%
    function xylabel(str1,str2)
        
 [MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default();       
 xlabel(str1,'FontSize',LabelFontSize,'FontAngle','Italic','FontWeight','bold','FontName',FontName,'Interpreter','LaTex')
 ylabel(str2,'FontSize',LabelFontSize,'FontAngle','Italic','FontWeight','bold','FontName',FontName,'Interpreter','LaTex')
    end
    
%% 
function title(str)
    [MarkerSize,FontName,TickFontSize,LabelFontSize,TitleFontSize]=pretty_m.default();     
title(str,'FontSize',TitleFontSize,'FontWeight','bold','FontAngle','Italic','FontName',FontName,'Interpreter','LaTex')

end

end

end


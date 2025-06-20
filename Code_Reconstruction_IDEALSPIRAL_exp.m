
%%

% Code MATLAB pour la Reconstruction des Images Metaboliques
% à partir des acquisitions enregistrées avec la séquence IDEAL SPIRAL
% nommée XXX

% A utiliser avec une version
% MATLAB au-delà de 2017

close all
clearvars -except
clc

% Figures Ancrées ou Non
docked_fig = 0;
if (docked_fig==1)
    set(0,'DefaultFigureWindowStyle','docked')
else
    set(0,'DefaultFigureWindowStyle','normal')
end

%% I - INPUT UTILISATEUR


% Pour commencer, il faut choisir quelques conditions et paramétres:

% Fichier contenant tout le fichier code de MATLAB IDEAL SPIRAL, veuillez
% le saisir:
folder_code='/home/wyoussfi/Data/IDEALspiral/matlab/IDEALSpiral_v250619';
addpath(folder_code)
cd(folder_code)

% Nombre de métabolites étudiés:Lactate, Alanine, ac-pyruv, Urée
Nm=3;
% Spectre acquis durant la même experience que les images (nouveau
% protocol) (1) ou acquis à part (0)
% 0:non  1:oui
spectre_compris=1;

% Correction de phase O1 due au SliceOffset
% 1: correction durant l'acquisition (par default pour la séquence XXX)
% 2: pas de correction voulu et/ou nécessaire
% 3: correction post-traitement (à partir du signal)
CorrType_array={'No O1 correction','Acqu O1 correction',...
    'Post O1 correction, si nécessaire'};
CorrctPhase=2;
CorrType=CorrType_array{CorrctPhase};

% Spectre identifié une 1ére fois peut étre fixé (1) pour les exp qui
% suivent ou modifié (0)
% 0:non  1:oui
spectre_fixe=0;

% Enregistrement des figures
% 0:non  1:oui
sauvegarde_figure=0;

% Application d'un filtre numérique aux images
% 0:non  1:oui
Filtre=0;


% FIN INPUT UTILISATEUR





%% II - SECTIONS

%| Ce programme comporte les sections suivantes:
%|
%| I    - INPUT UTILISATEUR
%| II   - SECTIONS
%| III  - INTRODUCTION
%| IV   - DIALOG
%| V    - ATTRIBUTION DES VALEURS AUX PARAMETRES NECESSAIRES
%| VI   - SPIRALE
%| VII  - FANTOME
%| VIII - VECTEURS TEMPORELS
%| IX   - CALCUL DES TERMES ET EQUATIONS NECESSAIRES
%| X    - RECONSTRUCTION DES ESPACES DES K
%| XI   - RECONSTRUCTION DES ESPACES IMAGES PAR LA METHODE DES MOINDRE CARRES
%| XII  - RECONSTRUCTION DE L'ESPCAE DES K ET D'IMAGE PAR LE " GRIDDING "
%| XIII - POINT SPREAD FUNCTION
%| XIV  - ENREGISTREMENT DES FIGURES
%| XV   - END



%% III - INTRODUCTION

%| VERSION FRANCAISE:
%|
%| Ce programme est un code réalisant, é partir des données acquises
%| in vitro ou numériquement, des reconstructions d'image métabolique
%| avec un balayage spiral.
%|
%| Pour cela, il y a un nombre de paramétre é insérer, et des conditions
%| d'acquisition é choisir.
%|
%| Pour une acquisition réelle, nous avons besoin de :
%| - données du balayage utilisé par l'appareil (imageur),
%| - le signal recueilli,
%| - fixer les fréquences des métabolites (spectre ou manuellement).
%|
%|
%| ENGLISH VERSION:
%|
%| This program is a code realizing, from the data acquired in vitro or
%| numerically, reconstructions of metabolic image with a spiral spatial
%| encoding.
%|
%| For this, there is a number of parameters to insert, and acquisition
%| conditions to choose.
%|
%| For a digital acquisition, we need:
%| - a spiral for spatial codding of the Fourier space,
%| - a digital phantom,
%| - a signal to calculate.
%|
%| For a real acquisition, we need:
%| - the spatial codding samples used by the device (imager),
%| - the collected signal,
%| - set the frequencies of metabolites (from a spectrum or manually).
%|
%|
%| in
%|
%|  %%% SPIRAL
%|	Np  		   int        Nombre d'échantillons de la spirale          [points]
%|	Nt  	   	   int 	      Nombre de tours                              [tours]
%|  kmax           double     Taille de la spirale                         [mm-1]
%|  Tntot          double     Durée totale de l'échantillonnage spiral     [s]
%|  %%% FANTOME
%|  N			   int        Taille de l'image désirée                    [pas d'unité]
%|  Nm             int        Nombre des métabolites                       [pas d'unité]
%|  Amp(i)         double     Amplitude initiale du signal du métabolite(i)[pas d'unité]
%|  E(i)           [1  4]     Ellipse du fantome é reconstruire            [pas d'unité]
%|  P(i)           [N  N]     Fantome initial du métabolite(i)             [pas d'unité]
%|  freq(i)        double     Fréquence de résonance du métabolite(i)      [Hz]
%|  T2e(i)         double     Temps de relaxation transversale T2*         [s]
%|  %%% ACQUISITION
%|  NTE            int        Nombre de temps d'écho                       [pas d'unité]
%|  deltaTE        double     Décalage de temps ^TE aprés chaque impulsion [s]
%|  FOV            int        Champ de vue (field of view)                 [mm]
%|  BW             double     Bande Passante                               [Hz]
%|  tetta          double     Flip Angle                                   [é]
%|
%|		kspace and FOV must be in reciprocal units!
%|
%|
%| out
%|
%|	esp(i)		   [Np  1]	  Espace de Fourier du métabolite(i)
%|	x(i)		   [N*N 1]	  Espace image du métabolite(i) par la méthode des moindres carrés
%|	im(i)          [N   N]	  Espace image du métabolite(i) dans la taille désirée
%|	k(i)_gridding  [Np  1]	  Espace de Fourier du métabolite(i) par la méthode GRIDDING
%|	im(i)_gridding [N*N 1]	  Espace image du métabolite(i) par la méthode GRIDDING
%|	perf           double     Valeur de performance d'une méthode de reconstruction
%|  PSNR           double     Valeur de performance d'une méthode de reconstruction
%|  PSF            [N   N]    Point Spread Function
%|
%|
%|
%| NB !!! Pour avoir une idée sur les équations utilisées dans ce code, veuillez consulter la publication :
%|
%| F. Wiesinger et al., IDEAL Spiral CSI for Dynamic Metabolic MR Imaging of Hyperpolarized [1-13C]Pyruvate.
%| Wiley Periodicals, Inc, 2012.
%|
%| Les equations dans le code sont numérotées selon la publication.
%|
%|
%| Copyright 2018-10-01 to 2021-09-30, Nour EL SABBAGH, Plateforme AgroResonance, INRA Theix, Clermont-Ferrand.





%% IV - DIALOG

%%%%%%%%%%%%%%%% Dialog %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Metabolites name and T2
prompt={'Metabolites Names (spaced out):','T2* in ms (spaced out):'};
filename = [folder_code '\specter_vs_frq.mat'];
if exist(filename, 'file')==2
    load('specter_vs_frq.mat')
    def=specter_vs_frq';
    clear specter_vs_frq
else
    def={'Lac Ala Urée','15 8 15'};
end
dlgTitle='Metabolites Identification';
lineNo=1;
AddOpts.Resize='off';
AddOpts.WindowStyle='modal';
AddOpts.Interpreter='tex';
AddOpts.Default = 'OK';
specter_vs_frq=newid(prompt,dlgTitle,lineNo,def,AddOpts);                 % Fenêtre de dialogue + bouton 'ENTER' pour 'OK'
save('specter_vs_frq.mat','specter_vs_frq')

% T2*??????????,,,,

%IDEAL SPIRAL SPECTRUM file
if(isfile('folder_exp_path.mat')==1)
    load('folder_exp_path.mat')
    folder_exp_path = uigetdir([folder_exp_path,'/..'],...
        'Select the IDEAL Spiral acquisition''s file directory');
else
    folder_exp_path = uigetdir(folder_code,...
        'Select the IDEAL Spiral acquisition''s file directory');
end
save('folder_exp_path.mat','folder_exp_path')

specter_choice=0;

%% V - ATTRIBUTION DES VALEURS AUX PARAMETRES NECESSAIRES

% Barre de Progression
prog_bar = waitbar(0,'Setting Parameters...','Position',[690 424 270 56]);
frames = java.awt.Frame.getFrames();
frames(end).setAlwaysOnTop(1);
pause(.25)
drawnow


% Dans cette partie, il faut attribuer les valeurs choisies aux variables
% correspondantes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% |:|:|:|:|:|:|:|:|:|:|:| Experimentation |:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|
% Pour l'expérimentation, deux cas sont possibles pour indiquer les
% fréquences des métabolites:
% (le cas est précisé selon la valeur de "specter_choice" dans VI - DIALOG
% || Expérimentation)
% 1- indiquer les valeurs manuellement (specter_choice=0)
% 2- utiliser un spectre RMN déjà enregistré (specter_choice=1)
%
% ETAPE A - Tout d'abord, pour les deux cas, la collecte des fichiers des
% paramètres correpondants à l'expérience à étudier est faite.
%
% ETAPE B - Ensuite, pour le deuxième cas (specter_choice=1), la
% determination des valeurs des fréquences se fait à partir de l'étude du
% spectre RMN aquis (compris dans la séquence d'imagerie ou à part).
%
% ETAPE C - Pour le premier cas (specter_choice=0), où les valeurs des
% fréquences sont déjà indiquées manuellement avant.
%
% ETAPE D - Enfin, tous les autres paramètres sont attribués normallement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% ETAPE A : start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(folder_exp_path)                    % Path du fichier de l'exp IDEAL SPIRAL                 []
% Qacqp_acquisition=ReadParameterFile('acqp'); % Fichier ACQP de l'expérience                    []
% Qmethod_acquisition=ReadParameterFile('method'); % Fichier Method de l'expérience              []
out=regexp(folder_exp_path,'/','split');
folder_exp=out{end};                   % Numéro de l'expérience                                []
folder_date_colle=out{end-1};          % Date de l'exp et Nom du fantéme                       []
folder_date=[folder_date_colle(7:8),...% Date de l'expérience DD-MM-YYY                        []
    '-',folder_date_colle(5:6),'-',folder_date_colle(1:4)];

%%%%%%%%%%%%%%%% ETAPE A : end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.33,prog_bar,'Setting Parameters...');drawnow
%%%%%%%%%%%%%%%% ETAPE B : start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.5,prog_bar,'Setting Parameters...');drawnow
type3_names=split(specter_vs_frq{1});
type3_names=char(type3_names);

%%%%%%%%%%%%%%%% ETAPE C : start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% ETAPE C : end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.67,prog_bar,'Setting Parameters...');drawnow

%%%%%%%%%%%%%%%% ETAPE D : start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(folder_exp_path)                    % Path du fichier de l'exp IDEAL SPIRAL                 []

% ACQUISITION:
A = regexp(fileread('method'),'\n','split');                              % Fichier METHOD
B = regexp(fileread('acqp'),'\n','split');                                % Fichier ACQP

whichline = find(contains(A,'PVM_Matrix='));                              % Index dans le fichier Method       []
N=cell2mat(textscan(A{whichline+1},'%f'));  N=N(1);                       % Taille de l'image désirée          []

whichline = find(contains(A,'PVM_FrqRef='));                              % Index dans le fichier Method       []
fref=cell2mat(textscan(A{whichline+1},'%f'));fref=fref(1)*1000000;        % Fréquence de référence             [Hz]

whichline = contains(A,'PVM_NSpiralIncrPP=');                       % Index dans le fichier Method       []
NTE=cell2mat(textscan(A{whichline}(22:end),'%f'));                        % Nombre de temps d'écho             []

whichline = contains(A,'PVM_RepetitionTime=');                      % Index dans le fichier Method       []
TR=cell2mat(textscan(A{whichline}(23:end),'%f'))/1000;                    % Temps de Répétition                [s]                                         

whichline = contains(A,'PVM_NRepetitions=');                        % Index dans le fichier Method       []
NTR=cell2mat(textscan(A{whichline}(21:end),'%f'));                        % Nombre de temps de répétition      []


whichline = contains(A,'PVM_NAverages=');                           % Index dans le fichier Method       []
NAv=cell2mat(textscan(A{whichline}(18:end),'%f'));                        % Nombre d'accumulation              []


whichline = contains(A,'PVM_EchoTimeIncrement=');                   % Index dans le fichier Method       []
deltaTE=cell2mat(textscan(A{whichline}(26:end),'%f'))/1000;               % Décalage de temps d'écho           [s]


whichline = find(contains(A,'EffectiveTESpiral='));                       % Index dans le fichier Method       []
whichline2 = find(contains(A,'PVM_DummyScans='));                         % Index dans le fichier Method       []
TEi_char=char(splitlines(strjoin(A(whichline+1:whichline2-1),'\n')));
TEi=zeros(10,100);
for indi=1:size(TEi_char,1)
    TEi(indi,1:length(cell2mat(textscan(TEi_char(indi,:),'%f'))'))=...
        cell2mat(textscan(TEi_char(indi,:),'%f'))';
end
TEi=reshape(TEi',1,1000);
TEi=TEi(1:NTE)/1000;                                                             % Train de  temps d'écho             [s]

whichline = contains(A,'PVM_EchoTime=');                            % Index dans le fichier Method       []
TEmin=cell2mat(textscan(A{whichline}(17:end),'%f'))/1000;                 % Temps d'écho (minimale)            [s]

whichline = contains(A,'PVM_EffSWh=');                              % Index dans le fichier Method       []
BW=cell2mat(textscan(A{whichline}(15:end),'%f'));                         % Bande passante d'acquisition       [Hz]

whichline = find(contains(A,'PVM_Fov='));                                 % Index dans le fichier Method       []
FOV=cell2mat(textscan(A{whichline+1},'%f'));FOV=FOV(1);                   % Champ de vue (field of vue)        [mm]

% kmax=N/(2*FOV);                                                           % Taille de la spirale               [mm-1]

whichline = contains(A,'PVM_EncNReceivers=');                       % Index dans le fichier Method       []
NRc=cell2mat(textscan(A{whichline}(22:end),'%f'));                        % Nombre de coil dans la sonde       []

whichline = contains(A,'PVM_Nucleus1Enum=<');                       % Index dans le fichier Method       []
Noyau=A{whichline}(22:24);                                                % Noyau Nucléaire étudié             []

whichline = contains(A,'PVM_SpiralFrequencyLimit=');                % Index dans le fichier Method       []
MaxFreq=cell2mat(textscan(A{whichline}(29:end),'%f'));                    % Frequence Limite des gradients     [Hz]

whichline = find(contains(A,'PVM_FrqWork='));                             % Index dans le fichier Method       []
FrqWork=cell2mat(textscan(A{whichline+1},'%f'));FrqWork=FrqWork(1)*10^6;  % Frequence de Résonance de travail  [Hz]

whichline = contains(A,'PVM_GradCalConst=');                        % Index dans le fichier Method       []
Grad_Amp_max=cell2mat(textscan(A{whichline}(21:end),'%f'));               % Amplitude Maximale des gradients   [Hz/mm]

whichline = contains(A,'PVM_SpiralEffectiveGradient=');             % Index dans le fichier Method       []
Grad_Amp_Eff=cell2mat(textscan(A{whichline}(32:end),'%f'));               % Amplitude effectif des gradients   [%]

whichline = contains(A,'PVM_SpiralSlewRateLimit=');                 % Index dans le fichier Method       []
Slew_Rate=cell2mat(textscan(A{whichline}(28:end),'%f'));                  % Slew Rate des gradients            [Hz/mm/ms]

whichline = contains(A,'RiseT=');                                   % Index dans le fichier Method       []
RiseT=cell2mat(textscan(A{whichline}(10:end),'%f'))/1000;                 % Rise time des gradients            [s]

whichline = contains(A,'PVM_SpiralPars=');                          % Index dans le fichier Method       []
Spiral_Pars=A{whichline}(23:27);                                          % Design de la spirale               []

if (contains(Spiral_Pars,'Three'))
    Spiral_Design='3D-3Contraints';
else
    Spiral_Design='2D-2Contraints';
end

whichline = contains(A,'ExcPul=(');                                 % Index dans le fichier Method       []
out=regexp(A{whichline}(12:end),',','split');BW_ExcPul=str2double(out{2});% Bande spectral d'excitation        [Hz]

whichline = find(contains(B,'ACQ_O1_list='));                             % Index dans le fichier Acqp         []
O1_list=cell2mat(textscan(B{whichline+1},'%f'));                          % O1 list                            [Hz]

whichline = find(contains(B,'ACQ_obj_order='));                           % Index dans le fichier Acqp         []
obj_order=cell2mat(textscan(B{whichline+1},'%f'));                        % Ordre des Objets                   []

whichline = contains(A,'PVM_SliceThick=');                          % Index dans le fichier Method       []
SliceThick=cell2mat(textscan(A{whichline}(19:end),'%f'));                 % Epaisseur de coupe                 [mm]

whichline = find(contains(A,'PVM_SliceOffset='));                         % Index dans le fichier Method       []
SliceOffset=cell2mat(textscan(A{whichline+1},'%f'));                      % Décalage de coupe                  [mm]

if (isempty(specter_vs_frq{2})==0)
    T2eCheck = 1;
    T2e=cell2mat(textscan(specter_vs_frq{2},'%f'))/1000;                                      % Temps de relaxation trans. T2*     [ms]
    invT2e=zeros(Nm,1);
    for i=1:Nm
        invT2e(i,1)=1/(T2e(i,1));
    end
else
    invT2e=zeros(Nm,1);
end


% SPIRAL:
whichline = contains(A,'PVM_SpiralSize=');                          % Index dans le fichier Method       []
Np=cell2mat(textscan(A{whichline}(19:end),'%f'));                         % Nb de points dans la spirale       []

whichline = contains(A,'PVM_SpiralPostSize=');                      % Index dans le fichier Method       []
Np_postsize=cell2mat(textscan(A{whichline}(23:end),'%f'));                % PostSize points dans la spirale    []

whichline = contains(A,'PVM_TrajSamples=');                         % Index dans le fichier Method       []
TrajSamples=cell2mat(textscan(A{whichline}(20:end),'%f'));                % Nb Total de points dans la spirale []

whichline = contains(A,'PVM_DigRes=');                              % Index dans le fichier Method       []
Extra_Samples=cell2mat(textscan(A{whichline}(15:end),'%f'));              % Nb de points sup dans la spirale   []


whichline = contains(A,'PVM_AcquisitionTime=');                     % Index dans le fichier Method       []
Tntot=cell2mat(textscan(A{whichline}(24:end),'%f'))/1000;                 % Durée totale de la spirale         [s]

whichline = contains(A,'PVM_SpiralNbOfGradientPoints=');            % Index dans le fichier Method       []
SpiralNbOfGradientPoints=cell2mat(textscan(A{whichline}(33:end),'%f'));   % Nb de point dans la forme gradiant []

whichline = contains(A,'PVM_SpiralNbOfInterleaves=');               % Index dans le fichier Method       []
Nint=cell2mat(textscan(A{whichline}(30:end),'%f'));                       % Nb d'interleaves                   []

whichline = find(contains(A,'NoGradRepetitionNumberPP='));                % Index dans le fichier Method       []
if (isempty(whichline))
    Spectr_GO=[];
else
    Spectr_GO=cell2mat(textscan(A{whichline+1},'%f'))';                       % Index Spectres durant la séquence  []
end
if (isempty(Spectr_GO)==1)
    Spectr_GO=0;
end
if (Spectr_GO~=0)
    Spectr_acquired=1;                                                    % Spectre acquis durant la séquence
else
    Spectr_acquired=0;                                                    % Spectre non acquis
end

whichline = find(contains(A,'NoGradRepetitionNumber='));                  % Index dans le fichier Method       []
if (isempty(whichline))
    Spectrums_achieved=[];
else
    Spectrums_achieved=length(cell2mat(textscan(A{whichline+1},'%f')));       % Nb de spectre acquis               []
end

if (isempty(Spectrums_achieved)==1)
    Spectrums_achieved=0;
end
NTR=NTR-Spectrums_achieved;                                               % Nb d'images acquis                 []

if (NTR>NTE)
    NTR=NTR/NTE;                                                          % Nb de block d'images acquis        []
end

%%%%%%%%%%%%%%%% ETAPE D : end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

waitbar(1,prog_bar,'Setting Parameters...');drawnow


%% VI - SPIRALE

%%%%%%%%%%%%%%%% Telecharger la spirale de l'imageur %%%%%%%%%%%%%%%%%%%%%%
waitbar(0,prog_bar,'Loading Spiral Trajectory...');drawnow
pause(.25)

fig1=figure2;
if (docked_fig~=1)
    fig1.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
end


%  Trajectoire théorique du fichier spiralKt (SpiralSize)

cd(folder_exp_path)
nametraj=[folder_date_colle(3:8),' spiralKt E',folder_exp];
nameGrad=[folder_date_colle(3:8),' spiralGt E',folder_exp];
if(exist([nametraj,'.mat'],'file')~=2)
    sp_bruker ='spiralKt';
    kt = regexp(fileread('spiralKt'),'\n','split');
    kt = cellfun(@str2num,kt,'un',0);
    kt=reshape(cell2mat(kt'),[length(kt)-1,2]);
    save(nametraj,'kt')
end
if(exist([nameGrad,'.mat'],'file')~=2)
    gx = regexp(fileread('spiralGx'),'\n','split');
    gx = cellfun(@str2num,gx,'un',0);
    gx=reshape(cell2mat(gx'),[length(gx)-1,1]);
    gy = regexp(fileread('spiralGx'),'\n','split');
    gy = cellfun(@str2num,gy,'un',0);
    gy=reshape(cell2mat(gy'),[length(gy)-1,1]);
    
    gt=gx+1i*gy;
    save(nameGrad,'gt','gx','gy')
end
load(nametraj)
load(nameGrad)
k(:,1)=kt(:,1)*N/FOV;
k(:,2)=kt(:,2)*N/FOV;
Gx=gx*Grad_Amp_max;
Gy=gy*Grad_Amp_max;
Gt=gt*Grad_Amp_max;


waitbar(.33,prog_bar,'Loading Spiral Trajectory...');drawnow

fig1;
subplot(2,3,2)


%%%%%%%%%%%%%  Trajectoire en disque  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     plot(k(1:end,1),k(1:end,2))
%     axis square
%     grid on
%     grid minor
%     xlim([min(k(1:end,1)) max(k(1:end,1))])
%     ylim([min(k(1:end,2)) max(k(1:end,2))])
%     xlabel('mm^{-1}')
%     ylabel('mm^{-1}')
%%%%%%%%%%%%%  Trajectoire en disque  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%  Trajectoire par rapport au temps  %%%%%%%%%%%%%%%%%%%%%%%%%%
sp_time=linspace(0,Tntot*1000,Np);
plot(sp_time,k(1:end,1))
hold on
plot(sp_time,k(1:end,2))
hold on
grid on
grid minor
title('Trajectoire Spirale vs temps')
xlabel('ms')
ylabel('mm^{-1}')
drawnow
%%%%%%%%%%%%%  Trajectoire par rapport au temps  %%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.67,prog_bar,'Loading Spiral Trajectory...');drawnow

%%%%%%%%%%%%%  Amplitude des Grandients par rapport au temps  %%%%%%%%%%%%%
fig1;
subplot(2,3,1)

tint=Tntot*1000*linspace(0,1,length(Gx));
yyaxis left
plot(tint,Gx)
hold on
plot(tint,Gy)
ylim([-max(Gx) max(Gx)])
xlimm=xlim;
ylabel('Amplitude de Gradient (Hz/mm)')
yyaxis right
ylim([-max(gx*100) max(gx*100)])
ylabel('% du Maximum d''Amplitude (%)')
grid on,grid minor
title('Amp Grad vs temps')
drawnow
%%%%%%%%%%%%%  Amplitude des Grandients par rapport au temps  %%%%%%%%%%%%%

waitbar(1,prog_bar,'Loading Spiral Trajectory...');drawnow

cd(folder_code)
% return

%% VII - FANTOME

waitbar(0,prog_bar,'Loading Acquisitions Data...');drawnow
pause(.25)
clc

% le fantôme dans ce cas correspond à l'image reconstruite à partir de la
% 1ère spirale acquise.

% TELECHARGEMENT DU SIGNAL DE L'IMAGEUR

cd (folder_exp_path)
data_set='fid';
id = fopen(data_set, 'r', 'l');            % Données enregistrées dans une seule colonne et parties
[donnee,~] = fread(id,  'int32');      % réelles et imaginaires de chaque échantillons entrelacées.
cd (folder_code)
data=zeros(length(donnee)/2,1);            % Rangement des données en forme complexe
for i=1:length(data)
    data(i,1)=donnee(2*i-1,1)+1i*donnee(2*i,1);
end
NTRR=size(data,1)/(TrajSamples+Extra_Samples);
datatxw=reshape(data(1:(TrajSamples+Extra_Samples)*NRc*(NTRR)),[(TrajSamples+Extra_Samples),NRc*(NTRR)]);
j=1;ss=1;iii=0;
clear newdatata
clear Spectr_GO_signal
waitbar(.33,prog_bar,'Loading Acquisitions Data...');drawnow

% SEPARATION DES SIGNAUX DES IMAGES ET DES SPECTRES, SI ACQUIS
if(Spectr_acquired~=0)
    spc_ind=1;
    spectr_at=Spectr_GO(spc_ind);
    imblc=0;
    if (docked_fig~=1)
        fig1.WindowState='minimize';                                 % Maximiser la fenétre de la figure     []
    end
    Spectr_GO_signal=zeros(TrajSamples,Spectrums_achieved);
    newdatata=zeros(TrajSamples,NTR*NTE);
    for i=1:NTRR                    % data en forme : TE par colonne
        if(i==spectr_at)
            Spectr_GO_signal(:,spc_ind)=datatxw(1:TrajSamples,i);
            iii=1;
            spc_ind=spc_ind+1;
            if(spc_ind<=length(Spectr_GO))
                spectr_at=Spectr_GO(spc_ind);
            end
            cprintf('key',['#' num2str(i) '  Spectre n° ' num2str(ss) '\n'])
            pause(0.25)
            ss=ss+1;
        else
            newdatata(:,j)=datatxw(1:TrajSamples,i);
            cprintf('text',['#' num2str(i) '  Image n° ' num2str(j)])
            if(imblc==0)
                cprintf('[1,0.5,0]','  \\ \n')
                imblc=imblc+1;
            else
                if(imblc==(NTE-1))
                    cprintf('[1,0.5,0]','  / \n')
                    imblc=0;
                else
                    if(imblc==fix(NTE/2))
                        if(mod(j,NTE)==0)
                            immet=fix(j/NTE);
                        else
                            immet=fix(j/NTE)+1;
                        end
                        cprintf('[1,0.5,0]',['  | Image Metabolique n° ' num2str(immet) '\n'])
                        imblc=imblc+1;
                    else
                        cprintf('[1,0.5,0]','  | \n')
                        imblc=imblc+1;
                    end
                end
            end
            pause(0.125)
            j=j+1;
        end
        waitbar(0.33+0.00825*i,prog_bar,'Loading Acquisitions Data...');drawnow
    end
    
    cprintf('text','\nEn total pour chaque métabolite: \n')
    cprintf('*magenta',[num2str(immet) ' image(s) '])
    cprintf('text','reconstruite(s) à different temps à partir de ')
    cprintf('*magenta',[num2str(NTE) ' acquisitions '])
    cprintf('text','chacune.\n')
    cprintf('text','\nL''évolution temporelle du spectre est étudiée sur ')
    cprintf('*magenta',[num2str(ss-1) ' points.\n\n'])
    
    pause(0.5)
    if (docked_fig~=1)
        fig1.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
    end
    
else
    newdatata=zeros(TrajSamples,NTRR);
    for i=1:NTRR                    % data en forme : TE par colonne
        newdatata(:,j)=datatxw(1:TrajSamples,i);
        j=j+1;
    end
end
waitbar(1,prog_bar,'Loading Acquisitions Data...');drawnow


% SIGNAL A FITTER POUR LES IMAGES METABOLIQUES
spiral_signal=newdatata;

% ETUDE DU SPECTRE
waitbar(0,prog_bar,'Finding Peaks from 1^{st} spectrum acquired...');drawnow
if(Spectr_acquired~=0 && spectre_compris==1)
    Y1=fftshift(fft(Spectr_GO_signal(:,1),TrajSamples))/TrajSamples;
    Y_module1 = abs((Y1(1:TrajSamples)));                      % Module de la TF(FID)                      [Hz]
    f1 = BW*linspace(-1/2,1/2,TrajSamples);                   % vecteur des fréquences                    [Hz]
    Minamp1=max(Y_module1)-5;
    NPics=0;
    figure2
    hold on
    Peak_Detector=yline(Minamp1,':','Label',['Peak Detector : ',num2str(Nm),' in total'],'LabelHorizontalAlignment','left');
    MinPeakDistanceppm=0;
    MinPeakDistanceHz=fref*10^-6*MinPeakDistanceppm;   % Distance min entre les pics en Hz         [Hz]
    findpeaks(Y_module1,...% Pics propriétés                       []
        f1,'MinPeakProminence',Minamp1,'Annotate','extents',...
        'MinPeakDistance',MinPeakDistanceHz);
    [pks_amp,pks_locs,pks_w,pks_p] = findpeaks(Y_module1,...% Pics propriétés                       []
        f1,'MinPeakProminence',Minamp1,'Annotate','extents',...
        'MinPeakDistance',MinPeakDistanceHz);
    legend('Location','northeast')
    xlim([-2000 2000])
    xlabel('Hz')
    drawnow
    wt=0.1;
    while(NPics~=Nm)
        children = get(gca, 'children');
        delete(children(1:end-1));
        set(gca,'ColorOrderIndex',1)
        Minamp1=Minamp1-3;
        [pks_amp,pks_locs,pks_w,pks_p] = findpeaks(Y_module1,...% Pics propriétés                       []
            f1,'MinPeakProminence',Minamp1,'Annotate','extents',...
            'MinPeakDistance',MinPeakDistanceHz);
        findpeaks(Y_module1,...% Pics propriétés                       []
            f1,'MinPeakProminence',Minamp1,'Annotate','extents',...
            'MinPeakDistance',MinPeakDistanceHz);
        legend('Location','northeast')
        Peak_Detector.Value=Minamp1;
        xlim([-2000 2000])
        title({['Spectrum in Hz exp #',folder_exp,' - peak(s) at'],[num2str(fix(pks_locs)),'  Hz']})
        drawnow
        NPics=length(pks_locs);
        wt=wt+0.01;
        if(wt<1)
            waitbar(wt,prog_bar,'Finding Peaks from 1^{st} spectrum acquired...');drawnow
        end
    end
    waitbar(1,prog_bar,'Finding Peaks from 1^{st} spectrum acquired...');drawnow
    pause(0.5)
    close
    pks_locss=pks_locs;
    pks_locss=sort(pks_locss,'descend');
    for i=1:Nm
        frq(i,1)=pks_locss(i);                              % Fréquences des métabolites            [Hz]
    end
    
    % NOM DES METABOLITES:
    frqq=num2str(round(frq,2));
    type3_names=split(specter_vs_frq{1});
    
    type3_names=char(type3_names);
    clear metab_names
    metab_names=char();
    for i=1:Nm
        metab_names(i,:)=[type3_names(i,:),': ',(frqq(i,:))];
    end
    metab_names_hor=blanks(Nm*2+numel(metab_names));
    for i=1:Nm
        metab_names_hor(1+(i-1)*(size(metab_names,2)+2):...
            i*(size(metab_names,2)+2))=[metab_names(i,:),'  '];
    end
    metab_names_hor=[metab_names_hor,' Hz '];
    while(contains(metab_names_hor,': ')==1 || contains(metab_names_hor,' :')==1)
        metab_names_hor=strrep(metab_names_hor,' :',':');
        metab_names_hor=strrep(metab_names_hor,': ',':');
    end
    
    % PLOT DU SPECTRE:
    waitbar(0,prog_bar,'Plotting all acquired spectrum...'); drawnow
    figure(fig1)
    
    subplot(2,3,3)
    hold on
    findpeaks(Y_module1,f1,'MinPeakProminence',Minamp1,'Annotate','extents','MinPeakDistance',MinPeakDistanceHz);
    legend('off')
    xlim([-2000 2000])
    grid on,grid minor
    text(pks_locs+100,pks_amp,flipud(type3_names))
    title(['Spctr Hz #',folder_exp,': ',metab_names_hor])
    xlabel('Hz')
    drawnow
    
    %%%% plot of other acquired spectrum during the experiment %%%%
    for spc=2:Spectrums_achieved
        Y2=fftshift(fft(Spectr_GO_signal(:,spc),TrajSamples))/TrajSamples;
        Y_module2 = abs((Y2(1:TrajSamples)));                      % Module de la TF(FID)
        plot(f1,Y_module2,'.','DisplayName',['Spectr ' num2str(spc)])
        drawnow
        Minamp2=max(Y_module2);
        NPics=0;
        while(NPics~=Nm)
            Minamp2=Minamp2-2;
            MinPeakDistanceHz=fref*10^-6*MinPeakDistanceppm;   % Distance min entre les pics en Hz         [Hz]
            [asp,bsp,csp,dsp] = findpeaks(Y_module1,...% Pics propriétés                       []
                f1,'MinPeakProminence',Minamp2,'Annotate','extents',...
                'MinPeakDistance',MinPeakDistanceHz);
            NPics=length(bsp);
        end
        pks_amp(:,spc)=asp;
        pks_locs(spc,:)=bsp;
        pks_w(spc,:)=csp;
        pks_p(:,spc)=dsp;
        waitbar(spc/(Spectrums_achieved),prog_bar,'Plotting all acquired spectrum...');drawnow
    end
    %%%% plot of other acquired spectrum during the experiment %%%%
    
    frq0=sort(frq,'ascend');
    %%%% Moyenner la valeur de la fréquence % à la largeur du pic %%%%%%%%%
    %     for peak_selection=1:Nm
    %         peakloc=pks_locs(1,peak_selection);peakamp=pks_amp(peak_selection,1);
    %         peak_indice=find(abs(f1-peakloc)==min(abs(f1-peakloc)));
    %         peak_find=true;
    %         peak_ind_Incr=0;
    % %         list_frq=zeros()
    %         list_frq(1)=f1(peak_indice);
    %         list_frq_ind=1;
    %         while(peak_find)
    %             peak_ind_Incr=peak_ind_Incr+1;
    %             if(Y_module1(peak_indice+peak_ind_Incr)>peakamp/2)
    %                 list_frq(list_frq_ind+1,1)=f1(peak_indice+peak_ind_Incr);
    %                 list_frq_ind=list_frq_ind+1;
    %             else
    %                 peak_find=false;
    %             end
    %         end
    %         peak_find=true;
    %         peak_ind_Incr=0;
    %         list_frq_ind=length(list_frq);
    %         while(peak_find)
    %             peak_ind_Incr=peak_ind_Incr+1;
    %             if(Y_module1(peak_indice-peak_ind_Incr)>peakamp/2)
    %                 list_frq(list_frq_ind+1,1)=f1(peak_indice-peak_ind_Incr);
    %                 list_frq_ind=list_frq_ind+1;
    %             else
    %                 peak_find=false;
    %             end
    %         end
    %         if(length(list_frq)>2)
    %             frq0(peak_selection)=mean(list_frq);
    %         end
    %         clear list_frq
    %     end
    %%%% Moyenner la valeur de la fréquence % à la largeur du pic %%%%%%%%%
    
    frq=sort(frq0,'descend');
end

% GRIDDING 1er TE
waitbar(0,prog_bar,'1^{st} Image reconstruction...');drawnow

%%%%%%%%%%%%%%%%%%  Correction de phase O1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(SliceOffset~=0 && CorrctPhase==3)
    for tr=1:NTR
        for i=1:NTE
            spiral_signal(1:Np,(tr-1)*NTE+i)=spiral_signal(1:Np,i)*exp(+1i*2*pi*O1_list*(i-1)*deltaTE);
        end
    end
end
%%%%%%%%%%%%%%%%%%  Correction de phase O1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[im_data_grid, k_data_grid, im_g, k_g] = mri_grid_linear(k,spiral_signal(1:Np,1),N, FOV,3); % GRIDDING
imagegrid=abs(im_data_grid);               % Image du 1er TE
fig1;
subplot(2,3,4)

imagesc(im_g{1},im_g{2},transpose(imagegrid))
hold on
plot(-FOV/2:FOV/2-2,zeros(length(-FOV/2:FOV/2-2),1),'w')
plot(zeros(length(-FOV/2:FOV/2),1),-FOV/2:FOV/2,'w')
axis square,%grid on,grid minor,
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap hot
cbar
title([Spiral_Design(1:2),'#',folder_exp,' Reco N=',num2str(N),'  MaxFrq=',num2str(MaxFreq),'Hz'])
hold off
drawnow
waitbar(1,prog_bar,'1^{st} Image reconstruction...');drawnow

waitbar(0,prog_bar,'Inverse space of the 1^{st} Image...');drawnow
fig1;
subplot(2,3,5)

imagesc(1:N,1:N,angle(k_data_grid))
title('Espace de Fourier')
axis square,%grid on,grid minor
set(gca,'xtick',[])
set(gca,'ytick',[])
cbar
title(colorbar,'rad')
waitbar(1,prog_bar,'Inverse space of the 1^{st} Image...');drawnow

waitbar(0,prog_bar,'Acquired Acquisition...');drawnow
fig1;
subplot(2,3,6)

t_ytest=linspace(0,Tntot*1000,length(spiral_signal));% Plot du signal
plot(t_ytest,abs(spiral_signal)),grid on,grid minor
xlabel('Time (ms)')
title('IDEAL SPIRAL FIDs')
drawnow
waitbar(1,prog_bar,'Acquired Acquisition...');drawnow



%% VIII - VECTEURS TEMPORELS

%%%%%%%%%%%%%%%% Vecteurs Temporels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(Nm==1)
    NTE=1;deltaTE=0;
end

% Vecteur temporel des échantillons de la spirale
%     tn=TEmin:Tntot/Np:TEmin+Tntot-Tntot/Np;
tn=0:Tntot/Np:Tntot-Tntot/Np;

% Vecteur temporel des incréments des spirales
tm=TEmin:deltaTE:TEmin+deltaTE*(NTE-1);
if (deltaTE==0)
    tm=TEmin*ones(1,NTE);
end




%% IX - CALCUL DES TERMES ET EQUATIONS NECESSAIRES

waitbar(0,prog_bar,'Compute Required Parameters for Wiesinger algorithm...');drawnow

%%%%%%%%%%%%%%%% Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Em,q=exp(i*wq*TEm)  ET  Em,q'                 Equation [1]
Eq=exp(1i*2*pi*frq*tm);
% Edag=Eq';
Edag=pinv(Eq);
waitbar(.33,prog_bar,'Compute Required Parameters for Wiesinger algorithm...');drawnow


%%%%%% exp(i*kn*rp)  ET  Fn,p=exp(i*kn*rp) ET Fn,p'  Equation [1]
rp=zeros(2,N^2);
intervv=linspace(-FOV,FOV,N);
for ii=1:N
    rp(1,1+(ii-1)*N:ii*N)=linspace(-FOV,FOV,N);
    rp(2,1+(ii-1)*N:ii*N)=intervv(ii);
end

waitbar(.67,prog_bar,'Compute Required Parameters for Wiesinger algorithm...');drawnow

Fq=exp(1i*pi*(k*rp));
Fdag0=Fq';
%     Fdag=pinv(Fq);

waitbar(1,prog_bar,'Compute Required Parameters for Wiesinger algorithm...');drawnow


%% X - RECONSTRUCTION DES ESPACES DES K

waitbar(0,prog_bar,'Running Wiesinger Fitting algorithm...');drawnow

signal_final=spiral_signal;

%%%%%%%%%%%%%%%% Espace des k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Eyn=Emq'ymn                                   Equation [2]
Eyn=ones(Np,Nm,NTR);
for tr=1:NTR
    for n=1:Np
        Eyn(n,:,tr)=signal_final(n,(tr-1)*NTE+1:tr*NTE)*Edag;
    end
end
waitbar(.33,prog_bar,'Running Wiesinger Fitting algorithm...');drawnow

%%%% esp.q(kn)=exp(-i*wq*tn*)(Eyn)q                Equation [2]

%%%% es : espace des k (dim=Np*1)
espacekk=cell(1,Nm);
for i=1:Nm
    espacekk{1,i}='es';
end
espacektot=genvarname(espacekk, 'es');
for i=1:Nm
    eval([espacektot{i} '= ones(Np,NTR);']);
end
for i=1:Nm
    for tr=1:NTR
        for n=1:Np
            eval([espacektot{i} '(n,tr)=exp(-1i*2*pi*(frq(i))*tn(1,n))*Eyn(n,i,tr);']);
        end
    end
    waitbar(.33+i*0.1675,prog_bar,'Running Wiesinger Fitting algorithm...');drawnow
end




%% XI - RECONSTRUCTION DES ESPACES IMAGES PAR LA METHODE DES MOINDRE CARRES

waitbar(0,prog_bar,'Metabolic Image Rec. with Least Square Estimation...');drawnow

%%%%%%%%%%%%%%%% Image par Moindres Carrés %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Correction dimension spatiale
facteur_corr=1;%(N/(2*FOV))/3;
ma=weight_vor(k(:,1), k(:,2),1)*facteur_corr;        % Correction de voronoi
Fdag=Fq';
% for i=1:Np
%     Fdag(:,i)=Fdag(:,i)*ma(i);
% end

%%% Espace des k Final

espacek=cell(1,Nm);
for i=1:Nm
    espacek{1,i}='esp';
end
espacektot2=genvarname(espacek, 'esp');
for i=1:Nm
    eval([espacektot2{i} '= zeros(Np,4,NTR);']);
end

%%% Pondération (default : signaux non pondérés)
ma(:,1)=ones(length(k),1);                                                 % Sans pondération       h=1
ma(:,2)=weight_vor(k(:,1), k(:,2),1)/max(weight_vor(k(:,1), k(:,2),1));    % Correction de voronoi  h=2
for i=1:Nm                                                                 % Pondération T2*        h=3
    ma(:,2+i)=transpose((exp(tn*invT2e(i))/100-1/100));%/max(transpose((exp(tn/T2e(i))/100-1/100)));
end
ma_legend=char('Sans pondération(m-c)','Correction de Voronoi','Pondération T2*(m-c)','la méthode GRIDDING');

im=zeros(N,N,4,Nm+1,NTR);                                                 % matrice image
for h=1%:3                                                                % Choix de la pondération désirée: 1, 2 ou 3
    x=zeros(N*N,1,Nm+1,NTR);                                              % matrice pré-image
    
    % Espace des k
    for i=1:Nm
        for tr=1:NTR
            eval([espacektot2{i} '(:,h,tr)=es' num2str(i) '(:,tr);']);     % Espace des k
        end
    end
    
    % Espace Image: xq=F'*esp.q(k)          Equation[3]
    if(h~=3)
        for i=1:Nm
            for tr=1:NTR
                eval(['x(:,1,i,tr)=(Fdag)*(esp' num2str(i) '(:,h,tr).*ma(:,h));']);  % espace image du métabolite (i)
                eval('x(:,1,Nm+1,tr)=x(:,1,Nm+1,tr)+x(:,1,i,tr);');                     % espace image de l'ensemble des métabolites
            end
        end
    else
        for i=1:Nm
            for tr=1:NTR
                eval(['x(:,1,i,tr)=(Fdag)*(esp' num2str(i) '(:,h,tr).*ma(:,2+i));']);% espace image du métabolite (i)
                eval('x(:,1,Nm+1,tr)=x(:,1,Nm+1,tr)+x(:,1,i,tr);');                     % espace image de l'ensemble des métabolites
            end
        end
    end
    
    % Reshape des Espace Image
    for i=1:Nm
        for tr=1:NTR
            eval('im(:,:,h,i,tr)=embed(x(:,1,i,tr),true(N));');                      % Image Finale du métabolite (i)
            eval('im(:,:,h,Nm+1,tr)=embed(x(:,1,Nm+1,tr),true(N));');                    % Image Finale de l'ensemble des métabolites
        end
    end
end


waitbar(0.33,prog_bar,'Metabolic Image Rec. with Least Square Estimation...');drawnow

pondr=char(' ','avec correction de Voronoi','avec pondération T2*');

%%%%%%%%%%%%%%%% Filtre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Filtre==1) % filtre directement appliqué aux images
    imagegrid_filtre=zeros(N,N,Nm,NTR);
    max_imagegrid=zeros(1,NTR);
    for q=1:Nm
        for tr=1:NTR
            imagegrid_filtre(:,:,tr)=abs(im(:,:,1,q,tr));
            max_imagegrid(1,tr)=max(max(imagegrid_filtre(:,:,tr)));
            for i=1:N
                for j=1:N
                    if(abs(im(i,j,1,q,tr))<max_imagegrid(1,tr)*0.60)
                        im(i,j,1,q,tr)=im(i,j,1,q,tr)/max_imagegrid(1,tr);
                    end
                end
            end
            imagegrid_filtre(:,:,tr)=abs(im(:,:,1,q,tr));
        end
    end
else % images filtrées enregistrées à part
    imagegrid_filtre=zeros(N,N,Nm,NTR);
    max_imagegrid=zeros(Nm,NTR);
    prctg_filtre=0.3;
    imagegrid_filtre(:,:,:,:)=abs(im(:,:,1,1:Nm,:));
    for q=1:Nm
        for tr=1:NTR
            max_imagegrid(q,tr)=max(max(imagegrid_filtre(:,:,q,tr)));
            for i=1:N
                for j=1:N
                    if(abs(imagegrid_filtre(i,j,q,tr))<max_imagegrid(q,tr)*prctg_filtre)
                        imagegrid_filtre(i,j,q,tr)=imagegrid_filtre(i,j,q,tr)/max_imagegrid(q,tr);
                    end
                end
            end
        end
    end
    imagegrid_filtre = abs(imagegrid_filtre);
end

%%%%%%%%%%%%%%%% Reconstruction par Moindres Carrés %%%%%%%%%%%%%%%%%%%%%%%
if (NTR==NTE)%(NTR==1)
    % Un seul TR => Etude Stationnaire
    
    sbplt=subplotin(Nm*2);
    
    for h=1%:2
        eval(['fig2' num2str(h) '=figure2;'])
        if (docked_fig~=1)
            eval(['fig2' num2str(h) '.WindowState=''maximize'';'])               % Maximiser la fenétre de la figure     []
        end
        drawnow
        % Images Métaboliques
        for i=1:Nm
            subplot(sbplt(1),sbplt(2),i)
            img_abs_tr_flip=abs(transpose(rot90(im(:,:,h,i),2)));
            imagesc(im_g{1},im_g{2},img_abs_tr_flip)
            axis square,grid on,grid minor,colormap hot
            title([metab_names(i,:),' Hz. #',folder_exp,', max=',num2str(max(max(img_abs_tr_flip))/10^6),'*10^6'],'FontSize',8)
            hold on
            plot(-FOV/2:FOV/2-2,zeros(length(-FOV/2:FOV/2-2),1),'w')
            plot(zeros(length(-FOV/2:FOV/2),1),-FOV/2:FOV/2,'w')
        end
        if(NTE>1)
            annotation('textbox', [0 0.9 1 0.1],'String', ['Slice Offset = ',num2str(SliceOffset) ,' mm -  ',CorrType,' - Images Métaboliques par la méthode des Moindres Carrés '],'Color','b','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
        end
        drawnow
        % Superposition avec les Images de Référence en 1H
        for i=1:Nm
            Map2=gray;
            eval(['ax1 = axes(''Parent'',fig2' num2str(h) ');'])
            eval(['ax2 = axes(''Parent'',fig2' num2str(h) ');'])
            subplot(sbplt(1),sbplt(2),Nm+i,ax1)
            subplot(sbplt(1),sbplt(2),Nm+i,ax2)
            set(ax1,'Visible','off');
            set(ax2,'Visible','off');
            Map1=hot;
            Map1(1,:)=[1 1 1];
            resol=N;
            eval(['load(''image_1H_bruker_' num2str(resol) '.mat'')'])
            eval('plotmetb=imshow(abs(transpose(fliplr(flipud((im(:,:,h,i)))))),''Parent'',ax2,''Colormap'',Map1,''DisplayRange'',[-max(max(abs(im(:,:,h,i))))/63 max(max(abs(im(:,:,h,i))))]);')
            axis square,grid on,grid minor
            hold on
            eval(['imshow(((transpose(image_1H_bruker_' num2str(resol) '))),''Parent'',ax1,''Colormap'',Map2,''DisplayRange'',[0 max(max(abs(image_1H_bruker_' num2str(resol) ')))]);'])
            axis square,grid on,grid minor
            title(['Superposition #',folder_exp,' du ',folder_date])
            plotmetb.AlphaData = 0.4;
            axis square,grid on,grid minor
            hold off
        end
        if(NTE>1)
            annotation('textbox', [0 0.42 1 0.1],'String', 'Superposition avec l''image de référence en 1H','Color','b','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
        end
        drawnow
    end
else
    
    % Plusieurs TR => Etude Temporelle
    resol=N;
    eval(['load(''image_1H_bruker_' num2str(resol) '.mat'')']) % reference 1H image
    sbplt=subplotin(NTR*2);
    metab_evol_max=zeros(Nm,NTR);
    for q=1:Nm
        for h=1%:2
            timeTR=(0:NTRR-1)*TR/60;
            timeSpct=TR*(Spectr_GO-1)/60;
            for ok=1:NTRR
                timeSpiral=timeTR(ismember(timeTR,timeSpct)==false);
            end
            timeSpiral=timeSpiral(NTE*(1:NTR)-NTE+1);
            
            f8=figure2;
            if (docked_fig~=1)
                f8.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
            end
            
            % Images Métaboliques
            for i=1:NTR
                subplot(sbplt(1),sbplt(2),i)
                img_abs_tr_flip=abs(transpose(rot90(im(:,:,h,q,i),2)));
                imagesc(im_g{1},im_g{2},img_abs_tr_flip)
                axis square,grid on,grid minor,colormap hot
                title(['t = ' num2str(timeSpiral(i)+TR*NTE/60/2) ' s.'])
                hold on
                plot(-FOV/2:FOV/2-2,zeros(length(-FOV/2:FOV/2-2),1),'w')
                plot(zeros(length(-FOV/2:FOV/2),1),-FOV/2:FOV/2,'w')
                %                 metab_evol_max(q,i)=abs(max(max(im(:,:,h,q,i))));
                a=reshape(maxk(abs(im(:,:,h,q,i)),5),5*N,1);
                a=sort(a,'descend');
                metab_evol_max(q,i)=mean(a(1:5));
            end
            if(NTE>1)
                annotation('textbox',[0.0031 0.9188 0.1255 0.0593],'String',...
                    ['\bf',type3_names(q,:),' ',num2str(fix(frq(q))),'Hz'],...
                    'Color','b','EdgeColor', 'none', 'HorizontalAlignment',...
                    'center','FontSize',15,'EdgeColor','b')
                annotation('textbox',[0.0078 0.8266 0.1026 0.0878],'String',...
                    ['Exp #',folder_exp,' du ',folder_date],...
                    'Color','black','EdgeColor', 'none', 'HorizontalAlignment',...
                    'center','FontSize',12)
                annotation('textbox', [0 0.875 1 0.1],'String',...
                    ['Images Métaboliques par la méthode des Moindres Carrés  - Slice Offset = ',num2str(SliceOffset) ,' mm'],...
                    'Color','black','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
            end
            drawnow
            
            % Superposition avec les Images de Référence en 1H
            for i=1:NTR
                Map2=gray;
                ax1 = axes('Parent',f8);
                ax2 = axes('Parent',f8);
                subplot(sbplt(1),sbplt(2),NTR+i,ax1)
                subplot(sbplt(1),sbplt(2),NTR+i,ax2)
                set(ax1,'Visible','off');
                set(ax2,'Visible','off');
                Map1=hot;
                Map1(1,:)=[1 1 1];
                resol=N;
                eval(['load(''image_1H_bruker_' num2str(resol) '.mat'')'])
                eval('plotmetb=imshow(abs(transpose(fliplr(flipud((im(:,:,h,q,i)))))),''Parent'',ax2,''Colormap'',Map1,''DisplayRange'',[-max(max(abs(im(:,:,h,q,i))))/63 max(max(abs(im(:,:,h,q,i))))]);')
                axis square,grid on,grid minor
                hold on
                eval(['imshow(transpose(image_1H_bruker_' num2str(resol) '),''Parent'',ax1,''Colormap'',Map2,''DisplayRange'',[0 max(max(abs(image_1H_bruker_' num2str(resol) ')))]);'])
                axis square,grid on,grid minor
                title(['t = ' num2str(timeSpiral(i)+TR*NTE/60/2) ' s.'])
                plotmetb.AlphaData = 0.4;
                axis square,grid on,grid minor
                hold off
            end
            %%% Filtrage & Superposition avec les Images de Référence en 1H
            %             for i=1:NTR
            %                 subplot(sbplt(1),sbplt(2),NTR+i)
            %                 im1F=transpose(rot90(imagegrid_filtre(:,:,q,i),2));
            %                 im1F=adapthisteq(im1F,'NumTiles',[2 2],'ClipLimit',0);
            %                 eval(['im2=transpose(image_1H_bruker_' num2str(resol) ');'])
            %                 im_imfuse = imfuse(im1F,im2,'falsecolor');
            %                 imagesc(im_imfuse)
            %                 axis square,grid on,grid minor
            %                 title(['t = ' num2str(timeSpiral(i)+TR*NTE/60/2) ' s.'])
            %             end
            
            if(NTE>1)
                annotation('textbox', [0 0.4 1 0.1],'String', 'Filtrage & Superposition avec l''image de référence en ^1H','Color','black','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
            end
            drawnow
        end
        drawnow
        waitbar(.33+q*0.1675,prog_bar,'Metabolic Image Rec. with Least Square Estimation...');drawnow
    end
    metab_evol_max=flipud(metab_evol_max);
    
    waitbar(0,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
    if (spectre_compris==1)
        fig_evolution=figure2;
        if (docked_fig~=1)
            fig_evolution.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
        end
        subplot(1,2,1)
        hold on,
        timeTR=(0:NTRR-1)*TR/60;
        timeSpct=TR*(Spectr_GO-1)/60;
        for ok=1:NTRR
            timeSpiral=timeTR(ismember(timeTR,timeSpct)==false);
        end
        timeSpiral=timeSpiral(NTE*(1:NTR)-NTE+1);
        met_evol_name=flipud(type3_names);
        for i=1:length(timeSpiral)
            aclrr=[0 0.4470 0.7410];
            if(i==1)
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock'...
                    ,'FaceColor',aclrr)
            else
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock',...
                    'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        
        for i=1:length(timeSpct)
            aclrr=[0.6350 0.0780 0.1840];
            if(i==1)
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'FaceColor',aclrr);
            else
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        
        waitbar(0.25,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
        
        for q=1:Nm
            str=sprintf(['peak ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            ppl=plot([0 timeSpct+TR/2/60],[0 pks_amp(q,:)]/max(max(pks_amp)),...
                'LineStyle','--','Marker','.','MarkerSize',20,'DisplayName',str,'LineWidth',1);
            c = get(ppl,'Color');
            str=sprintf(['image ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            plot(timeSpiral+TR*NTE/60/2,metab_evol_max(q,:)/max(max(metab_evol_max)),...
                'LineStyle','-','Marker','.','MarkerSize',20,'DisplayName',str, 'Color',c,'LineWidth',1)%,...
            %                 'HandleVisibility','off')
        end
        xlim([0 NTRR*TR/60])
        ylim([0 1])
        xlabel(['Repetition time (min) with TR = ' num2str(TR) ' s'])
        ylabel('Magnitude')
        title({['Peak (noramlized by ' num2str(round(max(max(pks_amp)))) ...
            ') &'],['Image Maximum (noramlized by ' num2str(round(max(max(metab_evol_max)))) ') Magntitude Evolution']})
        legend('Location','best')
        drawnow
        waitbar(0.5,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
        
        subplot(1,2,2)
        hold on,
        timeTR=(0:NTRR-1)*TR/60;
        timeSpct=TR*(Spectr_GO-1)/60;
        for ok=1:NTRR
            timeSpiral=timeTR(ismember(timeTR,timeSpct)==false);
        end
        timeSpiral=timeSpiral(NTE*(1:NTR)-NTE+1);
        met_evol_name=flipud(type3_names);
        for i=1:length(timeSpiral)
            aclrr=[0 0.4470 0.7410];
            if(i==1)
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock'...
                    ,'FaceColor',aclrr)
            else
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock',...
                    'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        for i=1:length(timeSpct)
            aclrr=[0.6350 0.0780 0.1840];
            if(i==1)
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'FaceColor',aclrr);
            else
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        
        waitbar(0.75,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
        
        for q=1:Nm
            str=sprintf(['peak ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            ppl=plot([0 timeSpct+TR/2/60],[0 pks_amp(q,:)./pks_w(:,q)']/max(max(pks_amp./pks_w')),...
                'LineStyle','--','Marker','.','MarkerSize',20,'DisplayName',str,'LineWidth',1);
            c = get(ppl,'Color');
            str=sprintf(['image ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            plot(timeSpiral+TR*NTE/60/2,metab_evol_max(q,:)/max(max(metab_evol_max)),...
                'LineStyle','-','Marker','.','MarkerSize',20,'DisplayName',str, 'Color',c,'LineWidth',1)%,...
            %                 'HandleVisibility','off')
        end
        xlim([0 NTRR*TR/60])
        ylim([0 1])
        xlabel(['Repetition time (min) with TR = ' num2str(TR) ' s'])
        ylabel('Magnitude')
        title({['Peak/FWHM (noramlized by ' num2str(round(max(max(pks_amp./pks_w')))) ...
            ') &'],['Image Maximum (noramlized by ' num2str(round(max(max(metab_evol_max)))) ') Magntitude Evolution']})
        legend('Location','best')
        annotation('textbox', [0.0109 0.8819 0.0740 0.1109],'String',...
            'Moindres Carrés','FontWeight','bold','EdgeColor', 'none',...
            'HorizontalAlignment', 'center','FontSize',16)
        drawnow
    end
    waitbar(1,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
end

drawnow

%%%%%%%%%%%%%%%%%%% Save en pdf pour papier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orient(fig21,'landscape')
% saveas(gcf,['P:\0Nour\3- Matlab\Code Imagerie Metabolique\Data set\Figures Papier\#',...
%     num2str(folder_exp),' ',num2str(SliceOffset) ,'mm ',CorrType,'.pdf'])
%%%%%%%%%%%%%%%%%%% Save en pdf pour papier %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% XII - RECONSTRUCTION DE L'ESPCAE DES K ET D'IMAGE PAR LE " GRIDDING "

% waitbar(1,prog_bar,'Finished');drawnow
% pause(2)
% close(prog_bar)
% return

waitbar(0,prog_bar,'Metabolite Image Rec. with Gridding...');drawnow

%%%%%%%%%%%%%%%% Reconstruction par GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Espace des k Final
for tr=1:NTR
    for i=1:Nm
        eval([espacektot2{i} '(:,4,tr)=es' num2str(i) '(:,tr);']);   % Espace des k
    end
end

% Espace Image: GRIDDING
for i=1:Nm
    for tr=1:NTR
        eval(['[im' num2str(i) '_gridding(:,:,tr), k' num2str(i) '_gridding(:,:,tr), im_g, k_g]' ...
            ' = mri_grid_linear(k, ' espacektot2{i} '(:,4,tr),N, FOV,3);']);
    end
end
for i=1:Nm
    for tr=1:NTR
        eval(['im(:,:,4,i,tr)=transpose(abs(im' num2str(i) '_gridding(:,:,tr)));']);
        im(:,:,4,Nm+1,tr)=im(:,:,4,Nm+1,tr)+im(:,:,4,i,tr);
    end
end
waitbar(.33,prog_bar,'Metabolite Image Rec. with Gridding...');drawnow

% Filtre
if (Filtre==1)
    for q=1:Nm
        for tr=1:NTR
            imagegrid_filtre(:,:,tr)=abs(im(:,:,4,q,tr));
            max_imagegrid(1,tr)=max(max(imagegrid_filtre(:,:,tr)));
            for i=1:N
                for j=1:N
                    if(abs(im(i,j,4,q,tr))<max_imagegrid(1,tr)*0.50)
                        im(i,j,4,q,tr)=im(i,j,4,q,tr)/max_imagegrid(1,tr);
                    end
                end
            end
            imagegrid_filtre(:,:,tr)=abs(im(:,:,4,q,tr));
        end
    end
end

% IMAGES METABOLIQUES
if(NTR==NTE)%(NTR==1)
    % Un seul TR => Etude Strationnaire
    
    sbplt=subplotin(Nm*2);
    fig3=figure2;
    if (docked_fig~=1)
        fig3.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
    end
    
    for i=1:Nm
        subplot(sbplt(1),sbplt(2),i)
        img_abs_tr_flip=abs(im(:,:,4,i));
        imagesc(im_g{1}, im_g{2},img_abs_tr_flip),axis square,grid on,grid minor,colormap hot,
        title([metab_names(i,:),' Hz. #',folder_exp,', intensité max = ',num2str(max(max(img_abs_tr_flip))/10^3),'*10^3'],'FontSize',8)
        caxis([-max(max(abs(im(:,:,4,i))))/63 max(max(abs(im(:,:,4,i))))])
        hold on
        plot(-FOV/2:FOV/2-2,zeros(length(-FOV/2:FOV/2-2),1),'w')
        plot(zeros(length(-FOV/2:FOV/2),1),-FOV/2:FOV/2,'w')
    end
    if(NTE>1)
        annotation('textbox', [0 0.9 1 0.1],'String', ['Slice Offset = ',num2str(SliceOffset) ,' mm -  ',CorrType,' - Images Métaboliques par la méthode de GRIDDING '],'Color','b','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
    end
    drawnow
    
    % Superposition avec les Images de Référence en 1H
    for i=1:Nm
        Map2=gray;
        ax1 = axes('Parent',fig3);
        ax2 = axes('Parent',fig3);
        subplot(sbplt(1),sbplt(2),Nm+i,ax1)
        subplot(sbplt(1),sbplt(2),Nm+i,ax2)
        set(ax1,'Visible','off');
        set(ax2,'Visible','off');
        Map1=hot;
        Map1(1,:)=[1 1 1];
        resol=N;
        eval(['load(''image_1H_bruker_' num2str(resol) '.mat'')'])
        eval('plotmetb=imshow(abs(im(:,:,4,i)),''Parent'',ax2,''Colormap'',Map1,''DisplayRange'',[-max(max(abs(im(:,:,4,i))))/63 max(max(abs(im(:,:,4,i))))]);')
        axis square,grid on,grid minor
        hold on
        eval(['imshow(transpose(image_1H_bruker_' num2str(resol) '),''Parent'',ax1,''Colormap'',Map2,''DisplayRange'',[0 max(max(abs(image_1H_bruker_' num2str(resol) ')))]);'])
        axis square,grid on,grid minor
        title(['Superposition #',folder_exp,' du ',folder_date])
        plotmetb.AlphaData = 0.4;
        axis square,grid on,grid minor
        hold off
    end
    if(NTE>1)
        annotation('textbox', [0 0.42 1 0.1],'String', 'Superposition avec l''image de référence en 1H','Color','b','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
    end
    drawnow
else
    % Plusieurs TR => Etude Temporelle
    
    metab_evol_max_gridding=zeros(Nm,NTR);
    for q=1:Nm
        sbplt=subplotin(NTR*2);
        fig3=figure2;
        if (docked_fig~=1)
            fig3.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
        end
        for i=1:NTR
            subplot(sbplt(1),sbplt(2),i)
            img_abs_tr_flip=abs(im(:,:,4,q,i));
            imagesc(im_g{1}, im_g{2},img_abs_tr_flip),axis square,grid on,grid minor,colormap hot,
            title(['t = ' num2str(timeSpiral(i)+TR*NTE/60/2) ' s.'])
            caxis([-max(max(abs(im(:,:,4,q,i))))/63 max(max(abs(im(:,:,4,q,i))))])
            hold on
            plot(-FOV/2:FOV/2-2,zeros(length(-FOV/2:FOV/2-2),1),'w')
            plot(zeros(length(-FOV/2:FOV/2),1),-FOV/2:FOV/2,'w')
            a=reshape(maxk(abs(im(:,:,4,q,i)),5),5*N,1);
            a=sort(a,'descend');
            metab_evol_max_gridding(q,i)=mean(a(1:5));
        end
        if(NTE>1)
            annotation('textbox',[0.0031 0.9188 0.1255 0.0593],'String',...
                ['\bf',type3_names(q,:),' ',num2str(fix(frq(q))),'Hz'],...
                'Color','b','EdgeColor', 'none', 'HorizontalAlignment',...
                'center','FontSize',15,'EdgeColor','b')
            annotation('textbox',[0.0078 0.8266 0.1026 0.0878],'String',...
                ['Exp #',folder_exp,' du ',folder_date],...
                'Color','black','EdgeColor', 'none', 'HorizontalAlignment',...
                'center','FontSize',12)
            annotation('textbox', [0 0.875 1 0.1],'String',...
                ['Images Métaboliques par la méthode Gridding  - Slice Offset = ',num2str(SliceOffset) ,' mm'],...
                'Color','black','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
        end
        drawnow
        
        % Superposition avec les Images de Référence en 1H
        for i=1:NTR
            Map2=gray;
            ax1 = axes('Parent',fig3);
            ax2 = axes('Parent',fig3);
            subplot(sbplt(1),sbplt(2),NTR+i,ax1)
            subplot(sbplt(1),sbplt(2),NTR+i,ax2)
            set(ax1,'Visible','off');
            set(ax2,'Visible','off');
            Map1=hot;
            Map1(1,:)=[1 1 1];
            resol=N;
            eval(['load(''image_1H_bruker_' num2str(resol) '.mat'')'])
            eval('plotmetb=imshow(abs(im(:,:,4,q,i)),''Parent'',ax2,''Colormap'',Map1,''DisplayRange'',[-max(max(abs(im(:,:,4,q,i))))/63 max(max(abs(im(:,:,4,q,i))))]);')
            axis square,grid on,grid minor
            hold on
            eval(['imshow(transpose(image_1H_bruker_' num2str(resol) '),''Parent'',ax1,''Colormap'',Map2,''DisplayRange'',[0 max(max(abs(image_1H_bruker_' num2str(resol) ')))]);'])
            axis square,grid on,grid minor
            title(['t = ' num2str(timeSpiral(i)+TR*NTE/60/2) ' s.'])
            plotmetb.AlphaData = 0.4;
            axis square,grid on,grid minor
            hold off
        end
        if(NTE>1)
            annotation('textbox', [0 0.4 1 0.1],'String', 'Superposition avec l''image de référence en ^1H','Color','black','EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12)
        end
        drawnow
        waitbar(.33+q*0.1675,prog_bar,'Metabolite Image Rec. with Gridding...');drawnow
    end
    
    waitbar(0,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
    metab_evol_max_gridding=flipud(metab_evol_max_gridding);
    if (spectre_compris==1)
        fig_evolution=figure2;
        if (docked_fig~=1)
            fig_evolution.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
        end
        subplot(1,2,1)
        hold on,
        timeTR=(0:NTRR-1)*TR/60;
        timeSpct=TR*(Spectr_GO-1)/60;
        for ok=1:NTRR
            timeSpiral=timeTR(ismember(timeTR,timeSpct)==false);
        end
        timeSpiral=timeSpiral(NTE*(1:NTR)-NTE+1);
        met_evol_name=flipud(type3_names);
        for i=1:length(timeSpiral)
            aclrr=[0 0.4470 0.7410];
            if(i==1)
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock'...
                    ,'FaceColor',aclrr)
            else
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock',...
                    'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        for i=1:length(timeSpct)
            aclrr=[0.6350 0.0780 0.1840];
            if(i==1)
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'FaceColor',aclrr);
            else
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        
        waitbar(.25,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
        
        for q=1:Nm
            str=sprintf(['peak ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            ppl=plot([0 timeSpct+TR/2/60],[0 pks_amp(q,:)]/max(max(pks_amp)),...
                'LineStyle','--','Marker','.','MarkerSize',20,'DisplayName',str,'LineWidth',1);
            c = get(ppl,'Color');
            str=sprintf(['image ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            plot(timeSpiral+TR*NTE/60/2,metab_evol_max_gridding(q,:)/max(max(metab_evol_max_gridding)),...
                'LineStyle','-','Marker','.','MarkerSize',20,'DisplayName',str, 'Color',c,'LineWidth',1)%,...
            %                 'HandleVisibility','off')
        end
        xlim([0 NTRR*TR/60])
        ylim([0 1])
        xlabel(['Repetition time (min) with TR = ' num2str(TR) ' s'])
        ylabel('Magnitude')
        title({['Peak (noramlized by ' num2str(round(max(max(pks_amp)))) ...
            ') &'],['Image Maximum (noramlized by ' num2str(round(max(max(metab_evol_max_gridding)))) ') Magntitude Evolution']})
        legend('Location','best')
        drawnow
        waitbar(0.5,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
        
        subplot(1,2,2)
        hold on,
        timeTR=(0:NTRR-1)*TR/60;
        timeSpct=TR*(Spectr_GO-1)/60;
        for ok=1:NTRR
            timeSpiral=timeTR(ismember(timeTR,timeSpct)==false);
        end
        timeSpiral=timeSpiral(NTE*(1:NTR)-NTE+1);
        met_evol_name=flipud(type3_names);
        for i=1:length(timeSpiral)
            aclrr=[0 0.4470 0.7410];
            if(i==1)
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock'...
                    ,'FaceColor',aclrr)
            else
                area([timeSpiral(i) timeSpiral(i)+NTE*TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','ImageBlock',...
                    'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        for i=1:length(timeSpct)
            aclrr=[0.6350 0.0780 0.1840];
            if(i==1)
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'FaceColor',aclrr);
            else
                area([timeSpct(i) timeSpct(i)+TR/60],[1,1],...
                    'FaceAlpha',0.25,'EdgeAlpha',0,...
                    'DisplayName','Spectrum'...
                    ,'HandleVisibility','off','FaceColor',aclrr)
            end
        end
        drawnow
        
        waitbar(0.75,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
        
        for q=1:Nm
            str=sprintf(['peak ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            ppl=plot([0 timeSpct+TR/2/60],[0 pks_amp(q,:)./pks_w(:,q)']/max(max(pks_amp./pks_w')),...
                'LineStyle','--','Marker','.','MarkerSize',20,'DisplayName',str,'LineWidth',1);
            c = get(ppl,'Color');
            str=sprintf(['image ' met_evol_name(q,:), '     ' , num2str(round(pks_locs(1,q))) ,' Hz ']);
            plot(timeSpiral+TR*NTE/60/2,metab_evol_max_gridding(q,:)/max(max(metab_evol_max_gridding)),...
                'LineStyle','-','Marker','.','MarkerSize',20,'DisplayName',str, 'Color',c,'LineWidth',1)%,...
            %                 'HandleVisibility','off')
        end
        xlim([0 NTRR*TR/60])
        ylim([0 1])
        xlabel(['Repetition time (min) with TR = ' num2str(TR) ' s'])
        ylabel('Magnitude')
        title({['Peak/FWHM (noramlized by ' num2str(round(max(max(pks_amp./pks_w')))) ...
            ') &'],['Image Maximum (noramlized by ' num2str(round(max(max(metab_evol_max_gridding)))) ') Magntitude Evolution']})
        legend('Location','best')
        annotation('textbox', [0.0109 0.9363 0.0740 0.0560],'String',...
            'Gridding','FontWeight','bold','EdgeColor', 'none',...
            'HorizontalAlignment', 'center','FontSize',16)
        drawnow
    end
    waitbar(1,prog_bar,'Metabolite Signal Evolution vs Time...');drawnow
end

drawnow

%% XIII - POINT SPREAD FUNCTION

% waitbar(1,prog_bar,'Finished');drawnow
% pause(2)
% close(prog_bar)
% return

% Pour l’évaluation du choix de l’échantillonnage spiral, on calcul la
% fonction d'étalement du point (PSF – Point Spread Function).
% Celle-ci décrit la réponse d'un système d'imagerie à une source
% ponctuelle ou à un objet ponctuel, et appartient au domaine spatial,
% à travers un produit de convolution. La PSF dans notre travail est la
% transformée de Fourier de nos échantillons dans l’espace-k. Plus le
% tracé de la surface en trois dimensions de cette fonction se rapproche
% d’une fonction de Dirac, moins il existera d’étalement spatial dans
% l’image reconstruite.
PSF = false;
if (PSF==true)
    waitbar(0,prog_bar,'Spiral Point Spread Function...');drawnow
    
    kmax=max(abs(k(:,1)+1i*k(:,2)));
    cordn=[k(:,1),k(:,2)];                                                     % Les coordonnées des échantillons de l'espace des k
    [X1,Y1] = ndgrid(-kmax:kmax/(N/2):+kmax);                                  % axes de la grille cartésienne choisie
    count=zeros(N+1,N+1);                                                      % La matrice de répartition des échantillons sur la grille
    % cartésienne.
    
    for t=1:Np
        x0=cordn(t,1);
        y0=cordn(t,2);
        
        for m=-kmax:kmax/(N/2):+kmax                                               % Dans cette boucle, je cherche pour un échantillons
            if(x0>=m && x0<(m+kmax/(N/2)))                                         % donné dans quelle case de la grille il appartient.
                x1=m;                                                              %
                x2=x1+kmax/(N/2);                                                  % La case est définie par un carré dans les quatres sommets
            end                                                                    % sont les points de la grille cartésienne.
            if(y0>=m && y0<(m+kmax/(N/2)))                                         %
                y1=m;                                                              % On cherche alors les coordonnées des sommets de la case
                y2=y1+kmax/(N/2);                                                  % correspondante formant un carré:
            end                                                                    % (x1,y1)-(x1,y2)-(x2,y1)-(x2,y2)
        end                                                                        %
        
        for kl=1:N+1                                                               % Dans cette boucle, je cherche le numéro de la case
            if(X1(kl,1)<=(x1+0.005) && X1(kl,1)>=(x1-0.005))                       % (indices de la case dans la matrice de répartition).
                il=kl;                                                             %
            end                                                                    % J'aoute une toute petite variation négligeable(0.005),
            % car sinon, je risque de compter les points sur les bords
            if(Y1(1,kl)<=(y1+0.005) && Y1(1,kl)>=(y1-0.005))                       % dans la mauvaise case.
                j=kl;                                                              %
            end                                                                    % Toujours, les points à l'arret droite appartient à la case
        end                                                                        % directement à gauche, et les points à l'arret au-dessous appartient
        count(j,il)=count(j,il)+1;                                                 % à la case directement au-dessous
        
    end
    
    waitbar(0.33,prog_bar,'Spiral Point Spread Function...');drawnow
    
    fig_PSF=figure2;
    if (docked_fig~=1)
        fig_PSF.WindowState='maximize';                                 % Maximiser la fenétre de la figure     []
    end
    
    subplot(2,2,1)
    axis square,grid on,grid minor,
    xlim([min(-kmax)   max(+kmax)])
    ylim([min(-kmax)   max(+kmax)])
    for i=1:N                                                                   % Répartition de densité des échantillons 2D
        for j=1:N
            patch([Y1(1,j),Y1(1,j),Y1(1,j+1),Y1(1,j+1)],[X1(i,1),X1(i+1,1),X1(i+1,1),X1(i,1)],count(i,j)/max(max(count)));
            axis ij
        end
    end
    colorbar
    title('Carte de densité d''échantillonnage 2D normalisée')
    colormap cool
    drawnow
    
    subplot(2,2,3)
    surf(X1, Y1, count/max(max(count)))                                                        % Répartition de densité des échantillons 3D
    axis square,grid on,grid minor,
    xlim([min(-kmax)   max(+kmax)])
    ylim([min(-kmax)   max(+kmax)])
    title('Carte de densité d''échantillonnage 3D normalisée')
    drawnow
    
    % subplot(2,2,3)
    % PSF_initiale=ifft2(count);                                                 % Point Spread Function
    % surf(0:N, 0:N, abs(PSF_initiale)/max(max(abs(PSF_initiale))))
    % axis square,grid on,grid minor,
    % xlim([0   N])
    % ylim([0   N])
    % % Forme de la matrice
    % % [ d  c ]
    % % [ b  a ]
    % title('Point Spread Function')
    
    waitbar(0.67,prog_bar,'Spiral Point Spread Function...');drawnow
    
    subplot(1,2,2)
    PSF_initiale=ifft2(count);                                                 % Point Spread Function
    PSF_finale=zeros(N,N);                                                     % Point Spread Function, mieux réparti
    PSF_finale(1:N/2+1,1:N/2)=abs(PSF_initiale(N/2+1:N+1,N/2+2:N+1))+PSF_finale(1:N/2+1,1:N/2); %a
    PSF_finale(1:N/2+1,N/2:N)=abs(PSF_initiale(N/2+1:N+1,1:N/2+1))+PSF_finale(1:N/2+1,N/2:N); %b
    PSF_finale(N/2+1:N,1:N/2)=abs(PSF_initiale(1:N/2,N/2+2:N+1))+PSF_finale(N/2+1:N,1:N/2); %c
    PSF_finale(N/2+1:N,N/2:N)=abs(PSF_initiale(1:N/2,1:N/2+1))+PSF_finale(N/2+1:N,N/2:N); %d
    surf(1:N, 1:N, abs(PSF_finale)/max(max(abs(PSF_finale))))
    axis square,grid on,grid minor
    xlim([0   N])
    ylim([0   N])
    % Forme de la matrice
    % [ a  b ]
    % [ c  d ]
    title('Normalized Point Spread Function')
    drawnow
    
    waitbar(1,prog_bar,'Spiral Point Spread Function...');drawnow
end

%% XIV - ENREGISTREMENT DES FIGURES

%%%%%%%%%%%%%%%% Enregistrement des figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(sauvegarde_figure==1)
    waitbar(0,prog_bar,'Saving Figures...');drawnow
    formatOut = 'YYmmDD hh_MM_ss';
    currDate=datestr(now,formatOut);
    folder_save=[folder_code,'\Archive 20',currDate(1:2),'\',currDate];
    
    prompt='Save File Name';
    dlgTitle='Save File Name';
    lineNo=1;
    AddOpts.Resize='on';
    AddOpts.WindowStyle='normal';
    AddOpts.Interpreter='tex';
    AddOpts.Default = 'OK';
    SaveFileName = newid(prompt,dlgTitle,lineNo,{''},AddOpts);             % Fenêtre pour rentrer les paramétres:
    
    folder_save_reco=[folder_save,'\',SaveFileName{1}];
    mkdir(folder_save_reco)
    cd(folder_save_reco)
    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 2:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = get(FigHandle, 'Name');
        FigNum    = get(FigHandle, 'Number');
        savefig(FigHandle, fullfile(folder_save_reco, [num2str(FigNum) '  ' FigName '.fig']));
    end
    
end




%% XV - END

waitbar(0,prog_bar,'Finishing');drawnow

cd(folder_code)
save('frequencies_13C.mat','frq')
% clc

frq_metab=num2str(frq');
cnt=true;
while(cnt)
    frq_metab=replace(frq_metab,'  ',' ');
    if(contains(frq_metab,'  ')==1)
        cnt=true;
    else
        cnt=false;
    end
end

%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nouns=specter_vs_frq{1};

waitbar(0.5,prog_bar,'Finishing');drawnow
if(spectre_compris==1)
    sc = ' (spectre(s) compris)';
    scc = [' \n ' num2str(length(Spectr_GO)) ' spectre(s) acquis'];
else
    sc='';
    scc='';
end
fprintf(['\n Expérience #' num2str(folder_exp) sc ...
    ' du ' folder_date ...
    ' \n ' num2str(Nm) ' métabolites: ' Nouns...
    ' \n Image métabolique à partir de'...
    ' ' num2str(NTE) ' temps d''écho, incrément de ' num2str(deltaTE*1000) ' ms'...
    scc...
    ' \n ' num2str(NTR) ' répétition(s) d''image, avec un TR = ' num2str(TR) ' ms'...
    ' \n Slice Offset    = ',num2str(SliceOffset) ,' mm '...
    ' \n Slice Thickness = ',num2str(SliceThick) ,' mm '...
    ' \n O1              = ',num2str(O1_list) ,' Hz '...
    ' \n ------- ' CorrType ' ------- \n\n'])
if(Filtre==1)
    fprintf(' ------- Images Filtrées ---------\n')
end

waitbar(1,prog_bar,'Finished');drawnow
pause(2)
close(prog_bar)

return




%%




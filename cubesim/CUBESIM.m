%*******************************************************************************
%*  CUBESIM	CubeSat Orbit and Attitude Simulator
%*
%*    version 2.7, February 14, 2005.  
%*    
%*    Copyright (C) 2002-2005 Jean-Francois Levesque  
%*    
%*    This software is provided 'as-is', without any express or implied  
%*    warranty.  In no event will the authors be held liable for any damages  
%*    arising from the use of this software.  
%*  
%*    Permission is granted to anyone to use this software for any purpose,  
%*    including commercial applications, and to alter it and redistribute it  
%*    freely, subject to the following restrictions:  
%*  
%*    1. The origin of this software must not be misrepresented; you must not  
%*       claim that you wrote the original software. If you use this software  
%*       in a product, an acknowledgment in the product documentation would be  
%*       appreciated but is not required.  
%*    2. Altered source versions must be plainly marked as such, and must not be  
%*       misrepresented as being the original software.  
%*    3. This notice may not be removed or altered from any source distribution.  
%*
%*    Jean-Francois Levesque, MS.
%*    jflev@yahoo.ca
%*
%*******************************************************************************
%
% CUBESIM is a Matlab scripts and functions package for spacecraft mission
% analysis and design.  It is intented to help CubeSat developers to design
% and simulate their system for different orbits, spacecraft configurations
% and ground station parameters.
%
% The simulator computes orbits, attitudes, ground station visibility,
% radiation, solar power, and much more using given orbit elements and
% spacecraft configuration.
%
%*******************************************************************************
%*  REQUIRED PACKAGE FUNCTIONS AND SCRIPTS:
%*  ---------------------------------------
%*  constants         		(global variable declaration)           Yes
%*  Rx, Ry, Rz          	(rotation cosine matrix)                -
%*  LL2XYZ, XYZ2LL      	(coordinate transformation)             Yes
%*  AziEl2XYZ, XYZ2AziEl	(coordinate transformation)             
%*  ROT2EUL, EUL2ROT		(rotation transformation)
%*  ROT2QUAT, QUAT2ROT		(rotation transformation)
%*  QUL2QUAT, QUAT2EUL		(rotation transformation)
%*  LEN                 	(norm of a vector array)
%*  plotTrace	        	(ground trace plot function)
%*  image2vector               	(convert binary image to XY vectors)
%*  getEphemeris, getElements, getTrace, plotTrace            	
%*  getMagField, getIntersect, getPtError, getRange     (magnetic stabilization)	
%*  identSC, trackSC, getDoppler, getRange, getLink     (G/S tracking)
%*  getShadow, getEclipse, getAtmDensity	
%*  plotMagField, plotHysteresis (simulations)
%*  T_ode2, T_ode7           (Temperature ODE simulation files)
%*  magStab3.mdl             (magnetic stabilization simulation)
%*  current.mat, default.mat (input parameter fields files)
%*
%*******************************************************************************
% 
% FEATURES
% --------
% - Orbit propagation and ephemeris calculation
% - Passive magnetic attitude stabilization using hysteresis materials
% - 3 magnetic hysteresis models
% - S/C antenna pointing calculation and trace 
% - Ground trace plots
% - Ground station tracking information and S/C visibility reports
% - Communication link occurences
% - Eclipse computation and reports
% - Solar cells power generation
% - Thermal simulations for 2nd and 7th order models
% - 22 simulation graphics and plots
% - Simulation results and data log on a text window
% 
% MAJOR VERSION UPDATES
% ---------------------
% - v1.5 Magnetic stabilization 3-D Simulink model based on rigid body dynamics
% - v1.5 Added moment of inertia tensor components input fields
% - v1.5 Added passive magnetic components description input fields
% - v1.5 Added initial attitude and rate input fields
% - v1.6 Corrected atmospheric density computation with exponential interpolation
% - v1.6 Quaternions and rates computed as attitude outputs
% - v1.6 Added magnetic angle and attitude graphs
% - v1.6 Orbit propagation using more accurate algorythms for eccentric orbits
% - v1.6 Added internal heat generation input field
% - v1.6 Revised solar Albedo calculations
% - v1.6 Display text box scrollable
% - v1.6 Added info buttons for field inputs
% - v1.6 Added check boxes for graph plotting
% - v1.7 Added antenna pointing deviation graphs
% - v1.8 Corrected UT and GST calculation
% - v1.8 Detailed magnetic hysteresis stabilization explanations
% - v1.9 Added 2 new models for magnetic hysteresis
% - v1.9 Parameters corrected in Simulink Model
% - v2.0 Added plot of angular rate
% - v2.1 Accurate World map
% - v2.1 Graph titles updated
% - v2.2 Field input parameters are now saved and restored automatically
% - v2.2 Restore default with RESET button
% - v2.2 Added plots of hysteresis models
% - v2.2 Correction in sampling frequency and time text fields
% - v2.3 Ground trace Legend tag names display bugs
% - v2.4 Replaced variable i (used as Matlab complex symbol) by variable k
% - v2.4 Correction of Maximum eclipse calculation
% - v2.5 Added simulation progression bar
% - v2.5 Added control window menubar
% - v2.5 Added simulation popup messages
% - v2.5 Changed tag names for legends
% - v2.5 Added graphs saving function to file and button
% - v2.6 Correction of bugs for invalid object handles
% - v2.6 Updated info question marks
% - v2.7 Correction of conflicts with new function names of Matlab v7
%
% NOTES
% -----
% - The graphical interface is set for display 1280x1024.
% - The simulation takes a lot of processing power and it is recommended
%   to run the simulation with a smaller number of samples than default
%   when running on a slow machine.
% - Use the question mark buttons on the control interface to get
%   more details about the simulation input parameters.
% 
% REFERENCES
% ----------
% - J-F. LEVESQUE, Passive Magnetic Attitude Stabilization using Hysteresis
%     Materials, Universite de Sherbrooke, SIgMA-PU-006-UdeS, Sept 2003.
% - LARSON, WERTZ, Space Mission Analysis and Design, Microcosm Press,
%     3rd edition, 1999.
% - WERTZ, Spacecraft Attitude Determination and Control, Microcosm, 1978.
% - BRYSON, Control of Spacecraft and Aircraft, Princeton University Press,
%     1994.
% - INCROPERA, DEWITT, Fundamentals of Heat and Mass Transfer, Wiley,
%     4th edition, 1996.
% - MUNSON, YOUNG, OKIISHI, Fundamental of Fluid Mechanics, Wiley, 
%     3rd edition, 1998.
% - David R. LIDE, CRC Handbook of Chemestry and Physics, CRC Press, 
%     82e Edition, 2001.
% - KREYZSIG, Advanced Engineering Mathematics, Wiley, 7th Edition, 1993
% - SADIKU, Elements of Electromagnetics, Saunders College Publishing, 1989.
% - Edward Della Torre, Magnetic Hysteresis, IEEE Press, New-York, 1999.
% - CUBESAT, http://cubesat.calpoly.edu
% - NARCISSAT, http://ssdl.stanford.edu/narcissat
% - AMSAT, http://www.amsat.org
%
%*******************************************************************************


function CUBESIM(action, parameter)

constants
format compact


if nargin < 1
   action = 'initialize';
end

if nargin < 2
   parameter = 'current';
end

FID = fopen('constants.m','r');
if FID < 0
     fprintf('Cannot locate cubesim package files, make sure all package files\n')
     fprintf('are present and set the Matlab path definition accordingly\n\n')
     return
end

%*******************************************************************************
% Control Window Formatting
%*******************************************************************************
if strcmp(action, 'initialize')
   close all
   ScreenSize = [1280 1024];
% ScreenWH = get(0,'ScreenSize'); 
% ScreenSize = [ScreenWH(1,3) ScreenWH(1,4)];

   fprintf('Control Window initialized for display %dx%d \n',ScreenSize)

   FIG = figure('Color',[0.8 0.8 0.8], ...  % 'menubar','none',...%'FileName','C:\MATLAB\work\CUBESIM.m', ...
      'PaperUnits','inches','PaperOrientation','landscape', 'PaperPositionMode','auto',... %'PaperPosition',[0 50 960 720],  ... 
      'Position',[0 35 ScreenSize-[0 64]], 'Tag','ControlWindow', 'ToolBar','none',...
      'Name','CUBESIM - CubeSat Orbit & Attitude Simulator v2.7');% 'menubar','none',...
     
     
   %============================================================
   % Graph Axes
   %============================================================
   %--------------------------------------------
   % Ground Trace
   %--------------------------------------------
   pos = get(FIG,'position');
   FrameSize = [pos(3)-380, (pos(3)-380)/2];
   FrameStart = [(pos(1) + 50), (pos(2) + pos(4) - (pos(3)-380)/2 - 60)];
   FrameEnd = FrameStart + FrameSize;

   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','on', 'CameraUpVector',[0 1 0], ...
      'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotGroundTrace', ...
      'XColor',[0 0 0], 'XLim',[-180 180], 'XLimMode','manual', 'XTickMode','manual', ...
      'XTick', [-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
         15 30 45 60 75 90 105 120 135 150 165 180], ...
      'YColor',[0 0 0], 'YLim',[-90 90], 'YLimMode','manual', 'YTickMode','manual', ...
      'YTick',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90], ...
      'ZColor',[0 0 0]);
   grid on;
   xlabel('Longitude [deg]','Tag', 'PlotGroundTrace', 'fontweight', 'bold')
   ylabel('Latitude [deg]','Tag', 'PlotGroundTrace', 'fontweight', 'bold')
   title('Orbital Ground Trace','Tag', 'PlotGroundTrace', 'fontweight', 'bold', 'fontsize', [12])
   axis([-180 180 -90 90]);
   
   % Compute and plot continents layout, takes lot of time to process
   WMimage = imread('wm3.bmp');	
   WM = image2vector(WMimage);
   plot(WM(2,:)/length(WMimage(1,:))*360-180, 90-WM(1,:)/length(WMimage(:,1))*180,...
      'k.','markersize',2, 'Tag','PlotGroundTraceMap', 'color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
   clear WMimage WM

   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''hideLeg'')', ...        
      'Position',[FrameEnd(1)-100 FrameStart(2)-35 100 15],'BackgroundColor',[0.8 0.8 0.8],...
      'String', 'Hide Legend', 'Style','checkbox', 'Tag','CheckBoxPlotHide');

   
   %--------------------------------------------
   % Orbit Shape
   %--------------------------------------------
   FrameSize = [pos(3)-380, pos(4)-(pos(3)-380)/2-195];
   FrameStart = [pos(1) + 50, 100];
   FrameEnd = FrameStart + FrameSize;

   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotOrbShape', 'UserData', 'Graph', 'XColor',[0 0 0], 'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Orbit Shape','Tag', 'PlotOrbShape', 'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Distance from Earth''s Center [km]','Tag', 'PlotOrbShape', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Distance from Earth''s Center [km]', 'Tag','PlotOrbShape', 'UserData', 'Graph', 'fontweight', 'bold')
   axis equal;

      
   %--------------------------------------------
   % Ground Station Tracking
   %--------------------------------------------
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotTrkAngle', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'XTick', [-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
         15 30 45 60 75 90 105 120 135 150 165 180], 'YColor',[0 0 0],...
      'YTick', [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180], 'YLim', [-180 180], 'ZColor',[0 0 0]);
   grid;
   title('Ground Station - Spacecraft Tracking Parameters','Tag', 'PlotTrkAngle',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);
   xlabel('Longitude [deg]','Tag', 'PlotTrkAngle', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Angles [deg]','Tag', 'PlotTrkAngle', 'UserData', 'Graph', 'fontweight', 'bold');
   axis([-180 180 -180 180]);
   
      
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotTrkRange', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'XTick', [-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
         15 30 45 60 75 90 105 120 135 150 165 180], 'XLim', [-180 180], ...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Ground Station - Spacecraft Tracking Ranges','Tag', 'PlotTrkRange',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Longitude [deg]','Tag', 'PlotTrkRange', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Ranges [km]', 'Tag', 'PlotTrkRange', 'UserData', 'Graph', 'fontweight', 'bold')
   %axis([-180 180 0 1.2*S_maxa]);
   
         
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotTrkDoppler', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'XTick', [-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
         15 30 45 60 75 90 105 120 135 150 165 180], 'XLim', [-180 180], ...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Ground Station - Spacecraft Link Doppler Shift','Tag', 'PlotTrkDoppler',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Longitude [deg]','Tag', 'PlotTrkDoppler', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Doppler Shift [kHz]', 'Tag','PlotTrkDoppler', 'UserData', 'Graph', 'fontweight', 'bold');
   %axis([-180 180 f0*.99995/1000 f0*1.00005/1000]);
   
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotTrkAngleTime', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YTick', [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180], 'YLim', [-180 180], 'ZColor',[0 0 0]);
   grid;
   title('Ground Station - Spacecraft Tracking Parameters','Tag', 'PlotTrkAngleTime',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);
   xlabel('Time [min]','Tag', 'PlotTrkAngleTime', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Angles [deg]','Tag', 'PlotTrkAngleTime', 'UserData', 'Graph', 'fontweight', 'bold');
   %axis([-180 180 -180 180]);
   
      
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotTrkRangeTime', 'UserData', 'Graph', 'XColor',[0 0 0],...
       'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Ground Station - Spacecraft Tracking Ranges','Tag', 'PlotTrkRangeTime',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotTrkRangeTime', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Ranges [km]', 'Tag', 'PlotTrkRangeTime', 'UserData', 'Graph', 'fontweight', 'bold');
   %axis([-180 180 0 1.2*S_maxa]);

   
         
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotTrkDopplerTime', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;


   title('Ground Station - Spacecraft Link Doppler Shift','Tag', 'PlotTrkDopplerTime',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotTrkDopplerTime', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Doppler Shift [kHz]', 'Tag','PlotTrkDopplerTime', 'UserData', 'Graph', 'fontweight', 'bold');
   %axis([-180 180 f0*.99995/1000 f0*1.00005/1000]);
   
   
   %--------------------------------------------
   % Radiation, Solar Power & Temperature
   %--------------------------------------------
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotRadArea', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Radiation - Area Exposed to Different Sources','Tag', 'PlotRadArea',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotRadArea', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Area [m²]', 'Tag','PlotRadArea', 'UserData', 'Graph', 'fontweight', 'bold');   
   

   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotRadPower', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Radiation - Power Generated from Solar Cells','Tag', 'PlotRadPower',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotRadPower', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Power [W]', 'Tag','PlotRadPower', 'UserData', 'Graph', 'fontweight', 'bold');   
   

   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotRadInput2', 'YLim', [-30 30], 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Radiation - Energy Absobed from Different Sources, 2nd Order Model','Tag', 'PlotRadInput2',....
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotRadInput2', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Power [W]', 'Tag','PlotRadInput2', 'UserData', 'Graph', 'fontweight', 'bold');     
   
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotRadTemp2', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'YTick', [-50 -25 0 25 50 75 100 125 150], 'YLim', [-50 150], 'ZColor',[0 0 0]);
   grid;
   title('Radiation - S/C Temperature History, 2nd Order Model','Tag', 'PlotRadTemp2',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotRadTemp2', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Temperature [°C]', 'Tag','PlotRadTemp2', 'UserData', 'Graph', 'fontweight', 'bold');     
   
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotRadInput7', 'YLim', [-30 30], 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Radiation - Energy Absobed from Different Sources, 7th Order Model',...
       'Tag', 'PlotRadInput7', 'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotRadInput7', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Power [W]', 'Tag','PlotRadInput7', 'UserData', 'Graph', 'fontweight', 'bold');     
   
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag','PlotRadTemp7', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'YTick', [-50 -25 0 25 50 75 100 125 150], 'YLim', [-50 150],'ZColor',[0 0 0]);
   grid;
   title('Radiation - S/C Temperature History, 7th Order Model','Tag', 'PlotRadTemp7',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12]);      
   xlabel('Time [min]','Tag', 'PlotRadTemp7', 'UserData', 'Graph', 'fontweight', 'bold');
   ylabel('Temperature [°C]', 'Tag','PlotRadTemp7', 'UserData', 'Graph', 'fontweight', 'bold');     
 
   
   %--------------------------------------------
   % Magnetic Stabilization Simulation
   %--------------------------------------------
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotMagSim', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'XTick', [-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
         15 30 45 60 75 90 105 120 135 150 165 180], 'YColor',[0 0 0],...
      'YTick', [0 15 30 45 60 75 90 105 120 135 150 165 180], 'ZColor',[0 0 0]);
   grid;
   title('Passive Magnetic Attitude Stabilization - Magnet Angle to Earth''s Field', 'Tag', 'PlotMagSim',...
       'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12])
   xlabel('S/C Longitude [deg]', 'Tag', 'PlotMagSim', 'UserData', 'Graph', 'fontweight', 'bold')
   ylabel('Angle [deg]', 'Tag', 'PlotMagSim', 'UserData', 'Graph', 'fontweight', 'bold')
   axis([-180 180 0 180]);
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotMagSimTime', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'YTick', [0 15 30 45 60 75 90 105 120 135 150 165 180], 'ZColor',[0 0 0]);
   grid;
   title('Passive Magnetic Attitude Stabilization - Magnet Angle to Earth''s Field',...
      'Tag', 'PlotMagSimTime', 'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12])
   xlabel('Time [min]', 'Tag', 'PlotMagSimTime', 'UserData', 'Graph', 'fontweight', 'bold')
   ylabel('Angle [deg]', 'Tag', 'PlotMagSimTime', 'UserData', 'Graph', 'fontweight', 'bold')
   %axis([-180 180 0 180]);
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotRateTime', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'ZColor',[0 0 0]);
   grid;
   title('Passive Magnetic Attitude Stabilization - Angular Rate',...
      'Tag', 'PlotMagSimTime', 'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12])
   xlabel('Time [min]', 'Tag', 'PlotRateTime', 'UserData', 'Graph', 'fontweight', 'bold')
   ylabel('Angle Rate [deg/s]', 'Tag', 'PlotRateTime', 'UserData', 'Graph', 'fontweight', 'bold')
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotAntPoint', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'XTick', [-180 -165 -150 -135 -120 -105 -90 -75 -60 -45 -30 -15 0 ...
         15 30 45 60 75 90 105 120 135 150 165 180], 'YColor',[0 0 0],...
      'YTick', [0 15 30 45 60 75 90 105 120 135 150 165 180], 'ZColor',[0 0 0]);
   grid;
   title('Passive Magnetic Attitude Stabilization - Antenna Pointing Deviation to Ground Stations', ...
      'Tag', 'PlotAntPoint', 'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12])
   xlabel('S/C Longitude [deg]', 'Tag', 'PlotAntPoint', 'UserData', 'Graph', 'fontweight', 'bold')
   ylabel('Angle [deg]', 'Tag', 'PlotAntPoint', 'UserData', 'Graph', 'fontweight', 'bold')
   axis([-180 180 0 180]);
   
   h1 = axes('Parent',FIG, 'Units','pixels', 'Box','off',...
      'CameraUpVector',[0 1 0], 'CameraUpVectorMode','manual', 'Color',[1 1 1], 'NextPlot','add', ...
      'Position',[FrameStart FrameSize], 'Tag', 'PlotAntPointTime', 'UserData', 'Graph', 'XColor',[0 0 0],...
      'YColor',[0 0 0], 'YTick', [0 15 30 45 60 75 90 105 120 135 150 165 180], 'ZColor',[0 0 0]);
   grid;
   title('Passive Magnetic Attitude Stabilization - Antenna Pointing Deviation to Ground Stations',...
      'Tag', 'PlotAntPointTime', 'UserData', 'Graph', 'fontweight', 'bold', 'fontsize', [12])
   xlabel('Time [min]', 'Tag', 'PlotAntPointTime', 'UserData', 'Graph', 'fontweight', 'bold')
   ylabel('Angle [deg]', 'Tag', 'PlotAntPointTime', 'UserData', 'Graph', 'fontweight', 'bold')
   %axis([-180 180 0 180]);
   
   
   %--------------------------------------------
   % Graph Select Title Bar List
   %--------------------------------------------   
   labelList=[' Orbit Shape| GS Tracking Angles| GS Tracking Ranges| GS Doppler Shifts|',...
         ' GS Tracking Angles (time scale)| GS Tracking Ranges (time scale)| GS Doppler Shifts (time scale)|',...
         ' Attitude - Magnetic Sabilization| Attitude - Magnetic Sabilization (time scale)|', ...
         ' Attitude - S/C Angular Rate (time scale)|', ...
         ' Attitude - Antenna Pointing| Attitude - Antenna Pointing (time scale)|', ...
         ' Radiation - Exposed Area| Radiation - Solar Power Generation|',...
         ' Radiation - Energy Absorbed, 2-DOF Model| Radiation - Temperature, 2-DOF Model|',...
         ' Radiation - Energy Absorbed, 7-DOF Model| Radiation - Temperature, 7-DOF Model|',...
         ' Ground Trace 3D - Geo Referential| Ground Trace 3D - Sun Referential',...
         ];
   cmdList=str2mat( ...
      ' set(findobj(FIG, ''Tag'', ''PlotOrbShape''), ''Visible'',''on'');', ...
      ' set(findobj(FIG, ''Tag'', ''PlotTrkAngle''), ''Visible'',''on'');', ...
      ' set(findobj(FIG, ''Tag'', ''PlotTrkRange''), ''Visible'',''on'');', ...
      ' set(findobj(FIG, ''Tag'', ''PlotTrkDoppler''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotTrkAngleTime''), ''Visible'',''on'');', ...
      ' set(findobj(FIG, ''Tag'', ''PlotTrkRangeTime''), ''Visible'',''on'');', ...
      ' set(findobj(FIG, ''Tag'', ''PlotTrkDopplerTime''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotMagSim''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotMagSimTime''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotRateTime''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotAntPoint''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotAntPointTime''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotRadArea''), ''Visible'',''on'');', ...
      ' set(findobj(FIG, ''Tag'', ''PlotRadPower''), ''Visible'',''on'');', ...
      ' set(findobj(FIG, ''Tag'', ''PlotRadInput2''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotRadTemp2''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotRadInput7''), ''Visible'',''on'');',...
      ' set(findobj(FIG, ''Tag'', ''PlotRadTemp7''), ''Visible'',''on'');',...
      ' figure(findobj(''Tag'',''PlotTrace3Dgeo'',''Type'',''figure''))',...
      ' figure(findobj(''Tag'',''PlotTrace3Dsun'',''Type'',''figure''))');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.8 0.8 0.8], 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', ...
      'Position',[FrameEnd-[315 360] 90 15], 'String', 'Select Graph', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'BackgroundColor',[1 1 1], 'Callback','cubeSim(''plot'')', ...
      'Position',[FrameEnd-[220 360] 220 20], ...
      'String',labelList, 'UserData' ,cmdList,...
      'Style','popupmenu', 'Tag','PopupPlot', 'Value',1);  
   
   LocalUpdateGraphs;  
   
   %--------------------------------------------
   % Display TextBox and Ctrl Buttons
   %--------------------------------------------
   %FrameSize = [pos(3), 45];
   FrameStart = [pos(1) + 50, 10];
   FrameEnd = [pos(3), 55];
   %FrameEnd = FrameStart + FrameSize;
   
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameStart 575 45], 'String', {'CubeSat Orbit & Attitude Simulator'}, 'Style','list', 'Tag','TextDisplay');    

   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''info'')', ...        
      'Position',[8 15 35 35], 'String', '?', 'Style','pushbutton', 'Tag','ButtonInfo',...
      'fontweight', 'bold', 'fontsize', [20], 'BackgroundColor', [0.7 0 1]);
   
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''history'')', ...        
      'Position',[FrameStart(1)+580 FrameStart(2)+5 75 30], 'String', 'Data Log', 'Style','pushbutton', 'Tag','ButtonHistory',...
      'UserData', 'No data available yet.  Run a simulation first.');
%    h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback', 'cubeSim(''clearGraph'');',...
%       'Position',[FrameStart(1)+660 FrameStart(2)+5 80 30], 'String', 'Clear Graphs', 'Style','pushbutton', 'Tag','ButtonClr');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback', 'cubeSim(''saveGraphs'');',...
      'Position',[FrameStart(1)+660 FrameStart(2)+5 80 30], 'String', 'Save Graphs', 'Style','pushbutton', 'Tag','ButtonSave');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''plotMag'')', ...        
      'Position',[FrameStart(1)+745 FrameStart(2)+5 75 30], 'String', 'Plot L-Shell', 'Style','pushbutton', 'Tag','ButtonPlot');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''plotHyst'')', ...        
      'Position',[FrameStart(1)+825 FrameStart(2)+5 75 30], 'String', 'Plot Models', 'Style','pushbutton', 'Tag','ButtonPlotH');
%    h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback', 'set(findobj(FIG,''style'',''edit''),''string'','''');',...
%      'Position',[FrameStart(1)+825 FrameStart(2)+5 80 30], 'String', 'Clear Fields', 'Style','pushbutton', 'Tag','ButtonFld');
%    h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'cubesim(''setDefault'');',...
%      'Position',[FrameStart(1)+825 FrameStart(2)+5 80 30], 'String', 'Default', 'Style','pushbutton', 'Tag','ButtonFld');
   
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'BackgroundColor', [0.8 0.8 0.8], 'FontWeight', 'Bold', 'Callback',...
      'set(findobj(FIG,''tag'',''ButtonSW1''),''userdata'',1);set(findobj(FIG,''tag'',''ButtonSW2''),''value'',0);set(findobj(FIG,''tag'',''ButtonSW3''),''value'',0);cubeSim(''eval'');',...
      'Position',[FrameStart(1)+920 FrameStart(2)+24 130 12], 'String', 'Hysteresis model 1', 'Style','radiobutton', 'Tag','ButtonSW1',...
      'UserData',1, 'value', 1);
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'BackgroundColor', [0.8 0.8 0.8], 'FontWeight', 'Bold', 'Callback',...
      'set(findobj(FIG,''tag'',''ButtonSW1''),''userdata'',2);set(findobj(FIG,''tag'',''ButtonSW1''),''value'',0);set(findobj(FIG,''tag'',''ButtonSW3''),''value'',0);cubeSim(''eval'');',...
      'Position',[FrameStart(1)+920 FrameStart(2)+12 130 12], 'String', 'Hysteresis model 2', 'Style','radiobutton', 'Tag','ButtonSW2');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'BackgroundColor', [0.8 0.8 0.8], 'FontWeight', 'Bold', 'Callback',...
      'set(findobj(FIG,''tag'',''ButtonSW1''),''userdata'',3);set(findobj(FIG,''tag'',''ButtonSW1''),''value'',0);set(findobj(FIG,''tag'',''ButtonSW2''),''value'',0);cubeSim(''eval'');',...
      'Position',[FrameStart(1)+920 FrameStart(2)+0 130 12], 'String', 'Hysteresis model 3', 'Style','radiobutton', 'Tag','ButtonSW3');
   
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback', 'cubeSim(''close'');cubeSim(''initialize'',''default'');',...
      'Position',[FrameEnd(1)-170 FrameStart(2) 80 35], 'String', 'Reset', 'Style','pushbutton', 'Tag','ButtonDef');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''close'')', ...        
      'Position',[FrameEnd(1)-85 FrameStart(2) 80 35], 'BackgroundColor',[0.7 0 1],...
      'String', 'Close', 'Style','pushbutton', 'Tag','ButtonClose');
   


   
   
   %============================================================
   % Simulation Text Fields  
   %============================================================


   FrameSize = [105, 255];   
   FrameEnd = [pos(1), pos(2)] + [pos(3), pos(4)] - [5, 60];
   FrameStart = FrameEnd - FrameSize;   
%   FrameStart = [960 525];
%   FrameEnd = [1270 930];
   FrameColor = [1 1 0.3];
   TextColor = [0 0 0];
   TextBoxColor = FrameColor;
   
   FID1 = fopen('current.mat','r');
   FID2 = fopen('default.mat','r');
   
   if parameter == 'current' 
       if FID1 > 0
           load current.mat	%load parameter values for text fields
       elseif FID2 > 0
           load default.mat	%load parameter values for text fields
           fprintf('Unable to find field parameters file current.mat, using default.mat\n')
           warndlg('Unable to find field parameters file current.mat, using default.mat')
       else
           fprintf('Unable to find field parameters files current.mat or default.mat\n')
           warndlg('Unable to find field parameters files current.mat or default.mat')           
       end
   else
       if FID2 > 0
           load default.mat	%load parameter values for text fields
           warndlg('The interface has been reset with default parameters');
       elseif FID1 > 0
           load current.mat	%load parameter values for text fields
           fprintf('Unable to find field parameters file default.mat, using current.mat\n')
           warndlg('Unable to find field parameters file default.mat, using current.mat')
           return
       else
           fprintf('Unable to find field parameters files current.mat or default.mat\n')
           warndlg('Unable to find field parameters files current.mat or default.mat')
           return
       end
   
   end
    
      
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',FrameColor, 'Units','pixels', ...
      'Position',[FrameStart FrameSize], 'Style','frame', 'Tag','FrameSim');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.8 0.8 0.8], 'FontWeight', 'Bold', 'fontsize', [10],...
      'Position',[FrameEnd-[95 25] [90 25]-[10 10]], 'Style','text', 'String', 'Sim', 'Tag','StaticTextSim'); 
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''simInfo'')', ...        
      'Position',[FrameEnd-[25 25] 15 15], 'String', '?', 'Style','pushbutton', 'Tag','ButtonInfo');

    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setlast'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[55 65] 45 20], 'String',sprintf('%0.2g',sim_dur), ...
      'Style','edit', 'Tag','EditTextDur', 'UserData',sim_dur);  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setdur'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[55 85] 45 20], 'String',sprintf('%0.4g',sim_last/60), ...
      'Style','edit', 'Tag','EditTextLast', 'UserData',sim_last);  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setlast'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[55 125] 45 20], 'String',sprintf('%0.4g',sim_sample), ...
      'Style','edit', 'Tag','EditTextSample', 'UserData',sim_sample);  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setfs'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[55 145] 45 20], 'String',sprintf('%0.4g',fs), ...
      'Style','edit', 'Tag','EditTextfs', 'UserData',fs);      
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setTs'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[55 165] 45 20], 'String',sprintf('%0.4g',Ts), ...
      'Style','edit', 'Tag','EditTextTs', 'UserData',Ts);  
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[95 42] 85 15], 'String', 'Duration','BackgroundColor',FrameColor,'Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', ...
      'Position',[FrameEnd-[95 62] 40 15], 'String', 'T [rev] ', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', ...
      'Position',[FrameEnd-[95 82] 40 15], 'String', 'T [min] ', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[95 102] 85 15], 'String', 'Sampling','BackgroundColor',FrameColor,'Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[95 122] 40 15], 'String', 'N     []', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[95 142] 40 15], 'String', 'fs [Hs]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[95 162] 40 15], 'String', 'Ts  [s]', 'Style','text', 'Tag','StaticText1'); 
   
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','', ...        
      'Position',[FrameEnd-[95 185] 85 15], 'String', 'Plot Trace',...
      'Value',1,'Visible','off', ... %remove checkbox from interface
      'BackgroundColor',[0.8 0.8 0.8],'Style','checkbox', 'Tag','CheckBoxPlotTrace');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','', ...        
      'Position',[FrameEnd-[95 185] 85 15], 'String', 'Therm Sim', ...
      'BackgroundColor',[0.8 0.8 0.8],'Style','checkbox', 'Tag','CheckBoxTempSim');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','', ...        
      'Position',[FrameEnd-[95 200] 85 15], 'String', '2-D Plots',...
      'BackgroundColor',[0.8 0.8 0.8],'Style','checkbox', 'Tag','CheckBoxPlot2D');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','', ...        
      'Position',[FrameEnd-[95 215] 85 15], 'String', '3-D Plots',...
      'BackgroundColor',[0.8 0.8 0.8],'Style','checkbox', 'Tag','CheckBoxPlot3D');   
   
   %'global XYZfield_I umag uhyst1 uhyst2 Bmag Bs1 Bs2 Vmag Vhyst1 Vhyst2 uo Br Hc sw Jinv t fs H0 q0;'
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback',...
      'global XYZfield_I umag uhyst1 uhyst2 Bmag Bs1 Bs2 Vmag Vhyst1 Vhyst2 uo Br Hc sw Jinv t fs H0 q0; cubeSim(''clearGraph'');cubeSim(''run'');cubeSim(''plot'')',...
      'Position',[FrameEnd-[95 245] 85 25], 'String', 'Run', 'Style','pushbutton', 'Tag','ButtonRun');%, 'BackgroundColor',[0.3 1 0.3]);
   
 
   %============================================================
   % Orbit Elements Text Fields 
   %============================================================
   FrameSize = [205, 255];   
   FrameEnd = FrameEnd - [105 0];
   FrameStart = FrameEnd-FrameSize;   
%   FrameStart = [960 555];
%   FrameEnd = [1165 820];
   FrameColor = [0.4 0.5 1];
   TextColor = [0 0 0];
   TextBoxColor = FrameColor;
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.8 0.8 0.8],...
      'Position',[FrameStart(1) FrameEnd(2) 310 20], 'Style','text', 'FontSize', [6],...
      'String', 'CubeSim v2.7 - (C) 2005 Jeff Levesque ( jflev@yahoo.ca )', 'Tag','StaticTextSim'); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',FrameColor, 'Units','pixels', ...
      'Position',[FrameStart FrameSize], 'Style','frame', 'Tag','FrameOrbit');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.8 0.8 0.8], 'FontWeight', 'Bold', 'fontsize', [10],...
      'Position',[FrameEnd-[195 25] [195 25]-[10 10]], 'Style','text', 'String', 'Orbit Elements', 'Tag','StaticTextOrbit');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''orbitInfo'')', ...        
      'Position',[FrameEnd-[25 25] 15 15], 'String', '?', 'Style','pushbutton', 'Tag','ButtonInfo');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 55] 55 20], 'String',sprintf('%0.1g',EYear), 'Style','edit', 'Tag','EditTextYear'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 75] 55 20], 'String',sprintf('%0.6g',Epoch), 'Style','edit', 'Tag','EditTextEpoch');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 95] 55 20], 'String',sprintf('%0.4g',incl*180/pi), 'Style','edit', 'Tag','EditTextIncl');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 115] 55 20], 'String',sprintf('%0.4g',raan*180/pi), 'Style','edit', 'Tag','EditTextRAAN');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback',...
      'cubeSim(''setn'');cubeSim(''setlast'');cubeSim(''testfs'');cubeSim(''plotOrbit'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 135] 55 20], 'String',sprintf('%0.4g',e), 'Style','edit', 'Tag','EditTextEcc');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 155] 55 20], 'String',sprintf('%0.4g',argp*180/pi), 'Style','edit', 'Tag','EditTextARGP');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 175] 55 20], 'String',sprintf('%0.4g',Msc0*180/pi), 'Style','edit', 'Tag','EditTextMsc0');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback',...
   'cubeSim(''seth'');cubeSim(''setlast'');cubeSim(''testfs'');cubeSim(''plotOrbit'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 195] 55 20], 'String',sprintf('%0.6g',n), 'Style','edit', 'Tag','EditTextn', 'UserData',n);
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback',...
      'cubeSim(''setn'');cubeSim(''setlast'');cubeSim(''testfs'');cubeSim(''plotOrbit'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[60 215] 55 20], 'String',sprintf('%0.4g',h), 'Style','edit', 'Tag','EditTexth', 'UserData',h);
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels',  'Callback','cubeSim(''eval'')',...
      'Position',[FrameEnd-[60 235] 55 20], 'String',sprintf('%0.2g',0), 'Style','edit', 'Tag','EditTextDecay');
%    h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels',  'Callback','cubeSim(''eval'')',...
%       'Position',[FrameEnd-[60 235] 55 20], 'String',sprintf('%0.2g',abs(orb_decay)), 'Style','edit', 'Tag','EditTextDecay');
   %h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
   %   'Position',[FrameEnd-[60 255] 55 20], 'String','', 'Style','edit', 'Tag','EditTextLife'); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 52] 135 15], 'String', 'Epoch year        [year]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 72] 135 15], 'String', 'Epoch time         [day]', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 92] 135 15], 'String', 'Inclination          [deg]', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 112] 135 15], 'String', 'RA of node         [deg]', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 132] 135 15], 'String', 'Eccentricity             []',  'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 152] 135 15], 'String', 'Arg of perigee     [deg]', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 172] 135 15], 'String', 'Mean anomaly     [deg]', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 192] 135 15], 'String', 'Mean motion [rev/day]', 'Style','text', 'Tag','StaticText1');   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 212] 135 15], 'String', 'Alt at perigee       [km]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[195 232] 135 15], 'String', 'Decay rate   [rev/day²]', 'Style','text', 'Tag','StaticText1'); 
   %h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
   %   'Position',[FrameEnd-[195 252] 135 15], 'String', 'Orbit life               [rev]', 'Style','text', 'Tag','StaticText1'); 
   
   %============================================================
   % Orbit Parameters Text Fields  
   %============================================================     
   FrameSize = [310, 160];   
   FrameEnd = FrameStart + [310 0];
   FrameStart = FrameEnd - FrameSize;
   BackColor = [0.9 0.9 0.9];
%   FrameStart = [960 395];
%   FrameEnd = [1270 555];  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',FrameColor, 'Units','pixels', ...
      'Position',[FrameStart FrameEnd-FrameStart], 'Style','frame', 'Tag','FrameOrbit2');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[195 30] 40 20], 'String',sprintf('%0.5g',a), 'Style','text', 'Tag','EditTexta');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[195 50] 40 20], 'String',sprintf('%0.5g',rp), 'Style','text', 'Tag','EditTextrp');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[195 70] 40 20], 'String',sprintf('%0.5g',ra), 'Style','text', 'Tag','EditTextra');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[195 90] 40 20], 'String',sprintf('%0.5g',Vp), 'Style','text', 'Tag','EditTextVp');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[195 110] 40 20], 'String',sprintf('%0.5g',Va), 'Style','text', 'Tag','EditTextVa');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[195 130] 40 20], 'String',sprintf('%0.5g',Eclipse_max/60), 'Style','text', 'Tag','EditTextEclmax');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[195 150] 40 20], 'String',sprintf('%0.5g',Eclipse_min/60), 'Style','text', 'Tag','EditTextEclmin');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[45 30] 40 20], 'String',sprintf('%0.5g',Psc/60), 'Style','text', 'Tag','EditTextPsc');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[45 50] 40 20], 'String',sprintf('%0.5g',dLON*180/pi), 'Style','text', 'Tag','EditTextdLON');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[45 70] 40 20], 'String',sprintf('%0.5g',Draan*Te*180/pi), 'Style','text', 'Tag','EditTextDraan');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[45 90] 40 20], 'String',sprintf('%0.4g',Dargp*Te*180/pi), 'Style','text', 'Tag','EditTextDargp');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[45 110] 40 20], 'String',sprintf('%0.5g',-Drag*1e+9), 'Style','text', 'Tag','EditTextDrag');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[45 130] 40 20], 'String',sprintf('%0.1g',-Da_rev/Psc*Psun),'Style','text', 'Tag','EditTextDecay2');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',BackColor, 'Units','pixels', ...
      'Position',[FrameEnd-[45 150] 40 20], 'String',sprintf('%3.0f',OrbitLife), 'Style','text', 'Tag','EditTextLife');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 27] 105 15], 'String','Semi-major    [km]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 47] 105 15], 'String','Perigee        [km]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 67] 105 15], 'String','Apogee        [km]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 87] 105 15], 'String','Perigee vel[km/s]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 107] 105 15], 'String','Apogee vel[km/s]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 127] 105 15], 'String','Max eclipse [min]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 147] 105 15], 'String','Min eclipse  [min]', 'Style','text', 'Tag','StaticText1');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[150 27] 100 15], 'String','Period         [min]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[150 47] 105 15], 'String','Node spac [°/rev]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[150 67] 105 15], 'String','Node prec [°/day]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[150 87] 105 15], 'String','Rate ArgP [°/day]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[150 107] 105 15], 'String','Air Drag        [nN]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[150 127] 105 15], 'String','Decay rate[km/yr]', 'Style','text', 'Tag','StaticText1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[150 147] 105 15], 'String','Orbit Life     [day]', 'Style','text', 'Tag','StaticText1');
   
   %============================================================
   % Spacecraft Text Fields
   %============================================================
   FrameSize = [310, 165];   
   FrameEnd = FrameStart + [310 0];
   FrameStart = FrameEnd - FrameSize;
%   FrameStart = [960 250];
%   FrameEnd = [1270 395];
   FrameColor = [1 0.7 0];
   TextColor = [0 0 0];
   TextBoxColor = FrameColor;   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',FrameColor, 'Units','pixels', ...
      'Position',[FrameStart FrameEnd-FrameStart], 'Style','frame', 'Tag','FrameSC');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.8 0.8 0.8], 'FontWeight', 'Bold', 'fontsize', [10],...
      'Position',[FrameEnd-[300 25] [300 25]-[10 10]], 'Style','text', 'String', 'Spacecraft Description', 'Tag','StaticTextSC');     
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''scInfo'')', ...        
      'Position',[FrameEnd-[25 25] 15 15], 'String', '?', 'Style','pushbutton', 'Tag','ButtonInfo');
   

   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 52] 20 15], 'String', 'Ixx [kgm²]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 72] 20 15], 'String', 'Iyy', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 92] 20 15], 'String', 'Izz', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 112] 20 15], 'String', 'Ixy', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 132] 20 15], 'String', 'Ixz',  'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 152] 20 15], 'String', 'Iyz',  'Style','text', 'Tag','StaticText1');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[280 55] 45 20], 'String',sprintf('%0.4g',Ixx*1e6), 'Style','edit', 'Tag','EditTextIxx');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[280 75] 45 20], 'String',sprintf('%0.4g',Iyy*1e6), 'Style','edit', 'Tag','EditTextIyy');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[280 95] 45 20], 'String',sprintf('%0.4g',Izz*1e6), 'Style','edit', 'Tag','EditTextIzz');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[280 115] 45 20], 'String',sprintf('%0.4g',Ixy*1e6), 'Style','edit', 'Tag','EditTextIxy');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels',  'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[280 135] 45 20], 'String',sprintf('%0.4g',Ixz*1e6), 'Style','edit', 'Tag','EditTextIxz');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels',  'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[280 155] 45 20], 'String',sprintf('%0.4g',Iyz*1e6), 'Style','edit', 'Tag','EditTextIyz');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[225 52] 65 15], 'String', 'Mass  [kg]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[225 72] 65 15], 'String', 'Width  [m]', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[225 92] 65 15], 'String', 'Length [m]', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[225 112] 65 15], 'String', 'Face1 [m²]', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[225 132] 65 15], 'String', 'Face2 [m²]',  'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[225 152] 65 15], 'String', 'Wall  [mm]',  'Style','text', 'Tag','StaticText1');   
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubesim(''setCd''),cubeSim(''eval'')', ...
      'Position',[FrameEnd-[160 55] 40 20], 'String',sprintf('%0.3g',msc), 'Style','edit', 'Tag','EditTextmsc');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubesim(''setA1sc''),cubesim(''setA2sc''),cubesim(''setCd''),cubeSim(''eval'')', ...
      'Position',[FrameEnd-[160 75] 40 20], 'String',sprintf('%0.3g',side), 'Style','edit', 'Tag','EditTextcsc');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback',...
      'cubesim(''setA2sc''),cubesim(''setCd''),cubeSim(''eval'')', ...
      'Position',[FrameEnd-[160 95] 40 20], 'String',sprintf('%0.3g',lsc), 'Style','edit', 'Tag','EditTextlsc');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[160 115] 40 20], 'String',sprintf('%0.4g',A1sc), 'Style','edit', 'Tag','EditTextA1sc','UserData',0.01);
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels',  'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[160 135] 40 20], 'String',sprintf('%0.4g',A2sc), 'Style','edit', 'Tag','EditTextA2sc','UserData',0.01);
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels',  'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[160 155] 40 20], 'String',sprintf('%0.4g',wallsc*1000), 'Style','edit', 'Tag','EditTextwallsc');

   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[110 52] 70 15], 'String', 'Drag coef []',  'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[110 72] 70 15], 'String', 'Cb   [kg/m²]',  'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[110 92] 70 15], 'String', 'PV eff.     %', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[110 112] 70 15], 'String', 'PV pack1 %', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[110 132] 70 15], 'String', 'PV pack2 %', 'Style','text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[110 152] 70 15], 'String', 'Heat     [W]',  'Style','text', 'Tag','StaticText1');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubesim(''setCd''),cubeSim(''eval'')', ...
      'Position',[FrameEnd-[40 55] 35 20], 'String',sprintf('%0.3g',Cd), 'Style','edit', 'Tag','EditTextCd','UserData',2);
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubesim(''setCb''),cubeSim(''eval'')', ...
      'Position',[FrameEnd-[40 75] 35 20], 'String',sprintf('%0.3g',Cb), 'Style','edit', 'Tag','EditTextCb','UserData',50);
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[40 95] 35 20], 'String',sprintf('%0.3g',PVeff*100), 'Style','edit', 'Tag','EditTextPVeff');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[40 115] 35 20], 'String',sprintf('%0.3g',PVpack1*100), 'Style','edit', 'Tag','EditTextPVpack1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[40 135] 35 20], 'String',sprintf('%0.3g',PVpack2*100), 'Style','edit', 'Tag','EditTextPVpack2');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[40 155] 35 20], 'String',sprintf('%0.3g',Qin), 'Style','edit', 'Tag','EditTextQin');   
   
   
   %============================================================
   % Magnetic Attitude Stabilization Text Fields  
   %============================================================
   FrameSize = [310, 185];   
   FrameEnd = FrameStart + [310 0];
   FrameStart = FrameEnd - FrameSize;
%   FrameStart = [960 145];
%   FrameEnd = [1270 250];
   FrameColor = [1 0.7 0];
   TextColor = [0 0 0];
   TextBoxColor = FrameColor;
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',FrameColor, 'Units','pixels', ...
      'Position',[FrameStart FrameEnd-FrameStart], 'Style','frame', 'Tag','FrameStab');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.8 0.8 0.8], 'FontWeight', 'Bold', 'fontsize', [10],...
      'Position',[FrameEnd-[300 25] [300 25]-[10 10]], 'Style','text', 'String', ...
      'Magnetic Attitude Stabilization', 'Tag','StaticTextStab');
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''magInfo'')', ...        
      'Position',[FrameEnd-[25 25] 15 15], 'String', '?', 'Style','pushbutton', 'Tag','ButtonInfo');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[250 50] 70 20], 'String', 'Direction','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[180 50] 40 20], 'String', 'm [g]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[140 50] 50 20], 'String', 'V [cm³]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[90 50] 40 20], 'String', 'Bs [T]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[45 50] 45 20], 'String', 'M [Am²]','Style', 'text', 'Tag','StaticText1');   
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 62] 150 15], 'String', 'Alnico-5','Style', 'text', 'Tag','StaticText1'); 
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[250 62] 5 15], 'String', '[','Style', 'text', 'Tag','StaticText1');   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[245 65] 20 20], 'String',sprintf('%0.1g',umagx), 'Style','edit', 'Tag','EditTextumagx');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[225 65] 20 20], 'String',sprintf('%0.1g',umagy), 'Style','edit', 'Tag','EditTextumagy');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[205 65] 20 20], 'String',sprintf('%0.1g',umagz), 'Style','edit', 'Tag','EditTextumagz');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[185 62] 5 15], 'String', ']','Style', 'text', 'Tag','StaticText1'); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setmmag'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[180 65] 40 20], 'String',sprintf('%0.4g',mmag*1e3), 'Style','edit', 'Tag','EditTextmmag', 'UserData',mmag);     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setVmag'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[135 65] 40 20], 'String',sprintf('%0.4g',Vmag*1e6), 'Style','edit', 'Tag','EditTextVmag', 'UserData',Vmag);    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setmmag'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[90 65] 40 20], 'String',sprintf('%0.3g',Bmag), 'Style','edit', 'Tag','EditTextBmag');    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback', ...
      'cubeSim(''setMmag'');cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[45 65] 40 20], 'String',sprintf('%0.3g',Mmag), 'Style','edit', 'Tag','EditTextMmag', 'UserData',Mmag); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 82] 50 15], 'String', 'HyMu80','Style', 'text', 'Tag','StaticText1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[250 82] 5 15], 'String', '[','Style', 'text', 'Tag','StaticText1');   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[245 85] 20 20], 'String',sprintf('%0.1g',uhyst1x), 'Style','edit', 'Tag','EditTextuhyst1x');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[225 85] 20 20], 'String',sprintf('%0.1g',uhyst1y), 'Style','edit', 'Tag','EditTextuhyst1y');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[205 85] 20 20], 'String',sprintf('%0.1g',uhyst1z), 'Style','edit', 'Tag','EditTextuhyst1z');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[185 82] 5 15], 'String', ']','Style', 'text', 'Tag','StaticText1'); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setmhyst1'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[180 85] 40 20], 'String',sprintf('%0.4g',mhyst1*1e3), 'Style','edit', 'Tag','EditTextmhyst1', 'UserData',mhyst1);    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setVhyst1'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[135 85] 40 20], 'String',sprintf('%0.4g',Vhyst1*1e6), 'Style','edit', 'Tag','EditTextVhyst1', 'UserData',Vhyst1);    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setmhyst1'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[90 85] 40 20], 'String',sprintf('%0.3g',Bs1), 'Style','edit', 'Tag','EditTextBs1');    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setMhyst1'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[45 85] 40 20], 'String',sprintf('%0.3g',Mhyst1), 'Style','edit', 'Tag','EditTextMhyst1', 'UserData',Mhyst1); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 102] 50 15], 'String', 'HyMu80','Style', 'text', 'Tag','StaticText1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[250 102] 5 15], 'String', '[','Style', 'text', 'Tag','StaticText1');   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[245 105] 20 20], 'String',sprintf('%0.1g',uhyst2x), 'Style','edit', 'Tag','EditTextuhyst2x');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[225 105] 20 20], 'String',sprintf('%0.1g',uhyst2y), 'Style','edit', 'Tag','EditTextuhyst2y');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[205 105] 20 20], 'String',sprintf('%0.1g',uhyst2z), 'Style','edit', 'Tag','EditTextuhyst2z');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[185 102] 5 15], 'String', ']','Style', 'text', 'Tag','StaticText1'); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setmhyst2'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[180 105] 40 20], 'String',sprintf('%0.4g',mhyst2*1e3), 'Style','edit', 'Tag','EditTextmhyst2', 'UserData',mhyst2);    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setVhyst2'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[135 105] 40 20], 'String',sprintf('%0.4g',Vhyst2*1e6), 'Style','edit', 'Tag','EditTextVhyst2', 'UserData',Vhyst2);    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setmhyst2'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[90 105] 40 20], 'String',sprintf('%0.3g',Bs2), 'Style','edit', 'Tag','EditTextBs2');    
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''setMhyst2'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[45 105] 40 20], 'String',sprintf('%0.3g',Mhyst2), 'Style','edit', 'Tag','EditTextMhyst2', 'UserData',Mhyst2); 
   
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
%      'Position',[FrameEnd-[195 125] 40 20], 'String','0.5', 'Style','edit', 'Tag','EditTextwspin');     
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
%      'Position',[FrameEnd-[45 125] 40 20], 'String','15', 'Style','edit', 'Tag','EditTextgamma'); 
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
%      'Position',[FrameEnd-[45 145] 40 20], 'String','0.1', 'Style','edit', 'Tag','EditTextwnu');    
   
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
%      'Position',[FrameEnd-[45 175] 40 20], 'String','11.356', 'Style','text', 'Tag','EditTextStabRes'); 
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', ...
%      'Position',[FrameEnd-[60 152] 55 15], 'String', 'Stab [°]', 'Style','text', 'Tag','StaticText1'); 
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
%      'Position',[FrameEnd-[300 122] 100 15], 'String', 'Spin rate  [°/s]','Style', 'text', 'Tag','StaticText1');
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
%         'Position',[FrameEnd-[150 122] 100 15], 'String', 'Nutation angle [°]', 'Style','text', 'Tag','StaticText1');
%   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
%      'Position',[FrameEnd-[150 142] 100 15], 'String', 'Nutation freq [°/s]', 'Style','text', 'Tag','StaticText1'); 

   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 127] 120 15], 'String', 'Stab. resolution  [°]','Style', 'text', 'Tag','StaticText1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[180 127] 40 15], 'String','11.38', 'Style','text', 'Tag','EditTextStabRes');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[135 127] 95 15], 'String', 'Dyn freq [rad/s]','Style', 'text', 'Tag','StaticText1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.9 0.9 0.9], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[45 127] 40 15], 'String','0.2482', 'Style','text', 'Tag','EditTextwn');     


   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 152] 100 15], 'String', 'Init. rate  [rad/s]','Style', 'text', 'Tag','StaticText1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[200 152] 5 15], 'String', '[','Style', 'text', 'Tag','StaticText1');   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[195 155] 40 20], 'String',sprintf('%0.3g',wbx), 'Style','edit', 'Tag','EditTextwbx');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[155 155] 40 20], 'String',sprintf('%0.3g',wby), 'Style','edit', 'Tag','EditTextwby');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''testfs'');cubeSim(''eval'')', ...
      'Position',[FrameEnd-[115 155] 40 20], 'String',sprintf('%0.3g',wbz), 'Style','edit', 'Tag','EditTextwbz');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[75 152] 5 15], 'String', ']','Style', 'text', 'Tag','StaticText1'); 
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 172] 100 15], 'String', 'Init. quaternions','Style', 'text', 'Tag','StaticText1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[200 172] 5 15], 'String', '[','Style', 'text', 'Tag','StaticText1');   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[195 175] 40 20], 'String',sprintf('%0.3g',q1), 'Style','edit', 'Tag','EditTextq1');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[155 175] 40 20], 'String',sprintf('%0.3g',q2), 'Style','edit', 'Tag','EditTextq2');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[115 175] 40 20], 'String',sprintf('%0.3g',q3), 'Style','edit', 'Tag','EditTextq3');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[75 175] 40 20], 'String',sprintf('%0.3g',q4), 'Style','edit', 'Tag','EditTextq4');     
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[35 172] 5 15], 'String', ']','Style', 'text', 'Tag','StaticText1'); 
   
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback',...
      'cubeSim(''qnorm'');', 'Position',[FrameEnd-[25 175] 20 20], 'String', '|q|', 'Style','pushbutton', 'Tag','ButtonRun');
   


   %============================================================
   % Ground Station Text Fields  

   %============================================================
   FrameSize = [310, 115];   
   FrameEnd = FrameStart + [310 0];
   FrameStart = FrameEnd - FrameSize;
%   FrameStart = [960 10];
%   FrameEnd = [1270 145];
   FrameColor = [0.3 1 0.3];
   TextColor = [0 0 0];
   TextBoxColor = FrameColor;
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',FrameColor, 'Units','pixels', ...
      'Position',[FrameStart FrameEnd-FrameStart], 'Style','frame', 'Tag','FrameGS');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[0.8 0.8 0.8], 'FontWeight', 'Bold', 'fontsize', [10],...
      'Position',[FrameEnd-[300 25] [300 25]-[10 10]], 'Style','text', 'String', 'Ground Stations', 'Tag','StaticTextGS');     
   h1 = uicontrol('Parent',FIG, 'Units','pixels', 'FontWeight', 'Bold', 'Callback','cubeSim(''gsInfo'')', ...        
      'Position',[FrameEnd-[25 25] 15 15], 'String', '?', 'Style','pushbutton', 'Tag','ButtonInfo');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[245 50] 40 20], 'String', 'Name','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[200 50] 50 20], 'String', 'LON [°]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[145 50] 50 20], 'String', 'LAT [°]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[95 50] 50 20], 'String', 'Elev [°]','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', ...
      'Position',[FrameEnd-[45 50] 40 20], 'String', '[MHz]','Style', 'text', 'Tag','StaticText1');
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[245 65] 40 20], 'String',NAME_GS1(1:3), 'Style','edit', 'Tag','EditTextGS1');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[200 65] 50 20], 'String',sprintf('%0.4g',LON_GS1*180/pi), 'Style','edit', 'Tag','EditTextGS1lon');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[145 65] 50 20], 'String',sprintf('%0.4g',LAT_GS1*180/pi), 'Style','edit', 'Tag','EditTextGS1lat');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[90 65] 40 20], 'String',sprintf('%0.3g',MEL_GS1*180/pi), 'Style','edit', 'Tag','EditTextGS1el');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[45 65] 40 20], 'String',sprintf('%0.4g',f1*1e-6), 'Style','edit', 'Tag','EditTextGS1fq');  
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[245 85] 40 20], 'String',NAME_GS2(1:3), 'Style','edit', 'Tag','EditTextGS2');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[200 85] 50 20], 'String',sprintf('%0.4g',LON_GS2*180/pi), 'Style','edit', 'Tag','EditTextGS2lon');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[145 85] 50 20], 'String',sprintf('%0.4g',LAT_GS2*180/pi), 'Style','edit', 'Tag','EditTextGS2lat');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[90 85] 40 20], 'String',sprintf('%0.4g',MEL_GS2*180/pi), 'Style','edit', 'Tag','EditTextGS2el');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[45 85] 40 20], 'String',sprintf('%0.4g',f2*1e-6), 'Style','edit', 'Tag','EditTextGS2fq');  
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[245 105] 40 20], 'String',NAME_GS3(1:3), 'Style','edit', 'Tag','EditTextGS3');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[200 105] 50 20], 'String',sprintf('%0.4g',LON_GS3*180/pi), 'Style','edit', 'Tag','EditTextGS3lon');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[145 105] 50 20], 'String',sprintf('%0.4g',LAT_GS3*180/pi), 'Style','edit', 'Tag','EditTextGS3lat');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[90 105] 40 20], 'String',sprintf('%0.4g',MEL_GS3*180/pi), 'Style','edit', 'Tag','EditTextGS3el');  
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
      'Position',[FrameEnd-[45 105] 40 20], 'String',sprintf('%0.4g',f3*1e-6), 'Style','edit', 'Tag','EditTextGS3fq');  
   
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 62] 55 15], 'String', 'Station 1','Style', 'text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 82] 55 15], 'String', 'Station 2', 'Style','text', 'Tag','StaticText1');
   h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
      'Position',[FrameEnd-[300 102] 55 15], 'String', 'Station 3', 'Style','text', 'Tag','StaticText1'); 
   
   %h1 = uicontrol('Parent',FIG, 'BackgroundColor',[1 1 1], 'Units','pixels', 'Callback','cubeSim(''eval'')', ...
   %   'Position',[FrameEnd-[60 130] 55 20], 'String','435.0', 'Style','edit', 'Tag','EditTextf0');  

   %h1 = uicontrol('Parent',FIG, 'BackgroundColor',TextBoxColor, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', ...
   %   'Position',[FrameEnd-[300 127] 235 15], 'String', 'Communication Frequency [MHz]','Style', 'text', 'Tag','StaticText1');
   
   %============================================================
   LocalUpdateGraphs;
   LocalDisplay('CubeSat Orbit & Attitude Simulator v2.7 - (C) 2005, Jean-Francois Levesque (jflev@yahoo.ca)', 'new')
   %LocalDisplay('Control Window initialized for display ','nl');
   %LocalDisplay(sprintf('%dx%d',ScreenSize), 'add');
   LocalDisplay(sprintf('Control Window initialized for display %dx%d',ScreenSize), 'nl');
   LocalDisplayInfo;
   
   
%******************************************************************************
% Text Field Updates
%******************************************************************************
elseif strcmp(action, 'qnorm')
   FIG = findobj('Tag','ControlWindow');
   q1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq1'),'String'), '%f');	%[] quaternion
   q2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq2'),'String'), '%f');	%[] quaternion
   q3 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq3'),'String'), '%f');	%[] quaternion
   q4 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq4'),'String'), '%f');	%[] quaternion
   q0 = [q1; q2; q3; q4]./sqrt(q1^2+q2^2+q3^2+q4^2); 		% intitial quaternion
   set(findobj(FIG, 'Tag', 'EditTextq1'),'String',sprintf('%0.3g',q0(1)));	
   set(findobj(FIG, 'Tag', 'EditTextq2'),'String',sprintf('%0.3g',q0(2)));	
   set(findobj(FIG, 'Tag', 'EditTextq3'),'String',sprintf('%0.3g',q0(3)));	
   set(findobj(FIG, 'Tag', 'EditTextq4'),'String',sprintf('%0.3g',q0(4)));	

elseif strcmp(action, 'testfs')
   FIG = findobj('Tag','ControlWindow');
   fs = sscanf(get(findobj(FIG, 'Tag', 'EditTextfs'),'String'), '%f'); 	      	%[s]
   wbx = sscanf(get(findobj(FIG, 'Tag', 'EditTextwbx'),'String'), '%f');%[rad/s] 
   wby = sscanf(get(findobj(FIG, 'Tag', 'EditTextwby'),'String'), '%f');%[rad/s] 
   wbz = sscanf(get(findobj(FIG, 'Tag', 'EditTextwbz'),'String'), '%f');%[rad/s]
   wb0 = [wbx; wby; wbz];		% initial rotation resolved in body coordinates
   Ixx = sscanf(get(findobj(FIG, 'Tag', 'EditTextIxx'),'String'), '%f')*1e-6;           %[kg.m²] 
   Iyy = sscanf(get(findobj(FIG, 'Tag', 'EditTextIyy'),'String'), '%f')*1e-6;           %[kg.m²]
   Izz = sscanf(get(findobj(FIG, 'Tag', 'EditTextIzz'),'String'), '%f')*1e-6;           %[kg.m²]
   Ixy = sscanf(get(findobj(FIG, 'Tag', 'EditTextIxy'),'String'), '%f')*1e-6;           %[kg.m²] product of inertia
   Ixz = sscanf(get(findobj(FIG, 'Tag', 'EditTextIxz'),'String'), '%f')*1e-6;           %[kg.m²] product of inertia
   Iyz = sscanf(get(findobj(FIG, 'Tag', 'EditTextIyz'),'String'), '%f')*1e-6;           %[kg.m²] product of inerti
   Isc = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];                                 	%[kg.m²] Moment of inertia tensor
   Mmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextMmag'),'String'), '%f');	%[Am^2]
   wn = sqrt(Mmag*B_EARTH/min(max(Isc)));
   set(findobj(FIG, 'Tag', 'EditTextwn'),'String',sprintf('%0.4g',wn));		%[rad/s]
   if fs < norm(wb0)/pi | fs < wn/pi	 	%Nyquist sampling criterion
      LocalDisplay('WARNING! Sampling frequency below initial attitude dynamics','nl')
      set(findobj(FIG, 'Tag', 'EditTextfs'),'ForegroundColor','r');     
      set(findobj(FIG, 'Tag', 'EditTextTs'),'ForegroundColor','r');
   else
      set(findobj(FIG, 'Tag', 'EditTextfs'),'ForegroundColor','k');     
      set(findobj(FIG, 'Tag', 'EditTextTs'),'ForegroundColor','k');
   end
   
elseif strcmp(action, 'setn')
   FIG = findobj('Tag','ControlWindow');
   h = sscanf(get(findobj(FIG, 'Tag', 'EditTexth'),'String'), '%f');	%[km] altitude at perigee
   e = sscanf(get(findobj(FIG, 'Tag', 'EditTextEcc'),'String'), '%f');	%[] eccentricity
   a = (h+Re) / (1-e);							%[km] semi-major axis
   n = sqrt(GMe/a^3)/(2*pi)*Te;						%[rev/day]
   set(findobj(FIG, 'Tag', 'EditTextn'),'String',sprintf('%0.6g',n));	%[rev/day]
   set(findobj(FIG, 'Tag', 'EditTextn'),'UserData',n);			%[rev/day]
   

elseif strcmp(action, 'seth')
   FIG = findobj('Tag','ControlWindow');
   n = sscanf(get(findobj(FIG, 'Tag', 'EditTextn'),'String'), '%f');    %[rev/day]
   e = sscanf(get(findobj(FIG, 'Tag', 'EditTextEcc'),'String'), '%f');	%[] eccentricity
   a = (sqrt(GMe)/(2*pi)*Te/n)^(2/3);					%[km] semi-major axis
   h = a*(1-e) - Re;							%[km] altitude at perigee
   set(findobj(FIG, 'Tag', 'EditTexth'),'String',sprintf('%0.5g',h));	%[km]
   set(findobj(FIG, 'Tag', 'EditTexth'),'UserData',h);			%[km]
   
elseif strcmp(action, 'setlast')
   FIG = findobj('Tag','ControlWindow');
   %global fs
   sim_sample =  sscanf(get(findobj(FIG, 'Tag', 'EditTextSample'),'String'), '%f');	%[]
   n = sscanf(get(findobj(FIG, 'Tag', 'EditTextn'),'String'), '%f');        	%[rev/day]
   sim_dur = sscanf(get(findobj(FIG, 'Tag', 'EditTextDur'),'String'), '%f');	%[rev]
   sim_last = sim_dur / n * Te;							%[s]
   Ts = sim_last / sim_sample;							%[s]
   fs = 1/ Ts;									%[Hz]
   set(findobj(FIG, 'Tag', 'EditTextLast'),'String',sprintf('%0.4g',sim_last/60));%[min]
   set(findobj(FIG, 'Tag', 'EditTextfs'),'String',sprintf('%0.4g',fs));		%[Hz]     
   set(findobj(FIG, 'Tag', 'EditTextTs'),'String',sprintf('%0.4g',Ts));		%[s]     
   set(findobj(FIG, 'Tag', 'EditTextDur'),'UserData',sim_dur);               	%[rev]
   set(findobj(FIG, 'Tag', 'EditTextLast'),'UserData',sim_last);               	%[s]
   set(findobj(FIG, 'Tag', 'EditTextfs'),'UserData',fs);               		%[Hz]     
   set(findobj(FIG, 'Tag', 'EditTextTs'),'UserData',Ts);               		%[s] 
   
elseif strcmp(action, 'setdur')
   FIG = findobj('Tag','ControlWindow');
   %global fs   
   sim_sample =  sscanf(get(findobj(FIG, 'Tag', 'EditTextSample'),'String'), '%f');	%[]
   n = sscanf(get(findobj(FIG, 'Tag', 'EditTextn'),'String'), '%f');		%[rev/day]
   sim_last = sscanf(get(findobj(FIG, 'Tag', 'EditTextLast'),'String'), '%f')*60; %[sec]
   sim_dur = n / Te * sim_last;							%[rev]
   Ts = sim_last / sim_sample;							%[s]
   fs = 1/ Ts;									%[Hz]
   set(findobj(FIG, 'Tag', 'EditTextDur'),'String',sprintf('%0.4g',sim_dur));	%[rev]   
   set(findobj(FIG, 'Tag', 'EditTextfs'),'String',sprintf('%0.4g',fs));		%[Hz]     
   set(findobj(FIG, 'Tag', 'EditTextTs'),'String',sprintf('%0.4g',Ts));		%[s]     
   set(findobj(FIG, 'Tag', 'EditTextDur'),'UserData',sim_dur);               	%[rev]

   set(findobj(FIG, 'Tag', 'EditTextLast'),'UserData',sim_last);               	%[s]
   set(findobj(FIG, 'Tag', 'EditTextfs'),'UserData',fs);               		%[Hz]     
   set(findobj(FIG, 'Tag', 'EditTextTs'),'UserData',Ts);               		%[s]     
   
elseif strcmp(action, 'setfs')
   FIG = findobj('Tag','ControlWindow');
   %global fs   
   sim_sample =  sscanf(get(findobj(FIG, 'Tag', 'EditTextSample'),'String'), '%f');	%[]
   Ts = sscanf(get(findobj(FIG, 'Tag', 'EditTextTs'),'String'), '%f'); 	      	%[s]
   n = sscanf(get(findobj(FIG, 'Tag', 'EditTextn'),'String'), '%f');		%[rev/day]
   fs = 1/ Ts;									%[Hz]
   sim_last = Ts * sim_sample;							%[s]
   sim_dur = n / Te * sim_last;							%[rev]
   set(findobj(FIG, 'Tag', 'EditTextLast'),'String',sprintf('%0.4g',sim_last/60));%[min]
   set(findobj(FIG, 'Tag', 'EditTextDur'),'String',sprintf('%0.4g',sim_dur));	%[rev]     
   set(findobj(FIG, 'Tag', 'EditTextfs'),'String',sprintf('%0.4g',fs));		%[Hz]     
   set(findobj(FIG, 'Tag', 'EditTextDur'),'UserData',sim_dur);               	%[rev]
   set(findobj(FIG, 'Tag', 'EditTextLast'),'UserData',sim_last);               	%[s]
   set(findobj(FIG, 'Tag', 'EditTextfs'),'UserData',fs);               		%[Hz]     
   set(findobj(FIG, 'Tag', 'EditTextTs'),'UserData',Ts);               		%[s]     
   
elseif strcmp(action, 'setTs')
   FIG = findobj('Tag','ControlWindow');
   %global fs
   sim_sample =  sscanf(get(findobj(FIG, 'Tag', 'EditTextSample'),'String'), '%f');	%[]
   fs = sscanf(get(findobj(FIG, 'Tag', 'EditTextfs'),'String'), '%f'); 	      	%[s]
   n = sscanf(get(findobj(FIG, 'Tag', 'EditTextn'),'String'), '%f');		%[rev/day]
   Ts = 1/ fs;									%[Hz]
   sim_last = Ts * sim_sample;							%[s]
   sim_dur = n / Te * sim_last;							%[rev]
   set(findobj(FIG, 'Tag', 'EditTextLast'),'String',sprintf('%0.4g',sim_last/60));%[min]
   set(findobj(FIG, 'Tag', 'EditTextDur'),'String',sprintf('%0.4g',sim_dur));	%[rev]     
   set(findobj(FIG, 'Tag', 'EditTextTs'),'String',sprintf('%0.4g',Ts));		%[s]     
   set(findobj(FIG, 'Tag', 'EditTextDur'),'UserData',sim_dur);               	%[rev]
   set(findobj(FIG, 'Tag', 'EditTextLast'),'UserData',sim_last);               	%[s]
   set(findobj(FIG, 'Tag', 'EditTextfs'),'UserData',fs);               		%[Hz]     
   set(findobj(FIG, 'Tag', 'EditTextTs'),'UserData',Ts);               		%[s]     
   
   
elseif strcmp(action, 'setCb')
   FIG = findobj('Tag','ControlWindow');
   Cb = sscanf(get(findobj(FIG, 'Tag', 'EditTextCb'),'String'), '%f');		%[km/m^2] ballistic coef
   msc = sscanf(get(findobj(FIG, 'Tag', 'EditTextmsc'),'String'), '%f');	%[km] mass
   side = sscanf(get(findobj(FIG, 'Tag', 'EditTextcsc'),'String'), '%f');	%[m] width
   lsc = sscanf(get(findobj(FIG, 'Tag', 'EditTextlsc'),'String'), '%f');	%[m] length
   Cd = msc/side/lsc/Cb;								%[] drag coef
   set(findobj(FIG, 'Tag', 'EditTextCb'),'UserData',Cb);			%[km/m^2]   
   set(findobj(FIG, 'Tag', 'EditTextCd'),'String',sprintf('%0.3g',Cd));		%[]
   set(findobj(FIG, 'Tag', 'EditTextCd'),'UserData',Cd);			%[]
   
elseif strcmp(action, 'setCd')
   FIG = findobj('Tag','ControlWindow');
   Cd = sscanf(get(findobj(FIG, 'Tag', 'EditTextCd'),'String'), '%f');		%[] drag coef
   msc = sscanf(get(findobj(FIG, 'Tag', 'EditTextmsc'),'String'), '%f');	%[km] mass
   side = sscanf(get(findobj(FIG, 'Tag', 'EditTextcsc'),'String'), '%f');	%[m] width
   lsc = sscanf(get(findobj(FIG, 'Tag', 'EditTextlsc'),'String'), '%f');	%[m] length
   Cb = msc/side/lsc/Cd;								%[km/m^2] ballistic coef
   set(findobj(FIG, 'Tag', 'EditTextCd'),'UserData',Cd);			%[]
   set(findobj(FIG, 'Tag', 'EditTextCb'),'String',sprintf('%0.3g',Cb));		%[km/m^2]
   set(findobj(FIG, 'Tag', 'EditTextCb'),'UserData',Cb);			%[km/m^2]   
   
elseif strcmp(action, 'setA2sc')
   FIG = findobj('Tag','ControlWindow');
   side = sscanf(get(findobj(FIG, 'Tag', 'EditTextcsc'),'String'), '%f');	%[m] width
   lsc = sscanf(get(findobj(FIG, 'Tag', 'EditTextlsc'),'String'), '%f');	%[m] length
   A2sc = side*lsc;								%[m^2] long face area
   set(findobj(FIG, 'Tag', 'EditTextA2sc'),'String',sprintf('%0.4g',A2sc));	%[m^2]
   set(findobj(FIG, 'Tag', 'EditTextA2sc'),'UserData',A2sc);			%[m^2]     
   
elseif strcmp(action, 'setA1sc')
   FIG = findobj('Tag','ControlWindow');
   side = sscanf(get(findobj(FIG, 'Tag', 'EditTextcsc'),'String'), '%f');	%[m] width
   A1sc = side*side;								%[m^2] base face area
   set(findobj(FIG, 'Tag', 'EditTextA1sc'),'String',sprintf('%0.4g',A1sc));	%[m^2]
   set(findobj(FIG, 'Tag', 'EditTextA1sc'),'UserData',A1sc);			%[m^2]   
   
   
%---Mag
elseif strcmp(action, 'setVmag')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGalnico
   Vmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextVmag'),'String'), '%f')*1e-6;	%[m^3] Volume
   Bmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextBmag'),'String'), '%f');	%[Tesla] Magnetic induction (field strength)
   Mmag = Bmag * Vmag / uo;							%[A.m^2] Magnetic moment
   mmag = Vmag * (SGalnico * 1000);						%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextVmag'),'UserData',Vmag);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextmmag'),'UserData',mmag);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextmmag'),'String',sprintf('%0.4g',mmag*1000));%[g] mass 
   set(findobj(FIG, 'Tag', 'EditTextMmag'),'UserData',Mmag);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextMmag'),'String',sprintf('%0.4g',Mmag));	%[A.m^2] Magnetic moment
   
elseif strcmp(action, 'setMmag')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGalnico

   Mmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextMmag'),'String'), '%f');	%[A.m^2] Magnetic moment
   Bmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextBmag'),'String'), '%f');	%[Tesla] Magnetic induction (field strength)

   Vmag = Mmag * uo / Bmag;							%[m^3] Volume
   mmag = Vmag * (SGalnico * 1000);						%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextMmag'),'UserData',Mmag);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextmmag'),'UserData',mmag);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextmmag'),'String',sprintf('%0.4g',mmag*1000));%[g] mass 
   set(findobj(FIG, 'Tag', 'EditTextVmag'),'UserData',Vmag);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextVmag'),'String',sprintf('%0.4g',Vmag*1e+6));%[cm^3] Volume  
   
elseif strcmp(action, 'setmmag')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGalnico
   mmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextmmag'),'String'), '%f')/1000;	%[kg] mass
   Bmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextBmag'),'String'), '%f');	%[Tesla] Magnetic induction (field strength)
   Vmag = mmag / (SGalnico * 1000);						%[m^3] Volume
   Mmag = Bmag * Vmag / uo;							%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextmmag'),'UserData',mmag);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextVmag'),'UserData',Vmag);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextVmag'),'String',sprintf('%0.4g',Vmag*1e+6));%[cm^3] Volume  
   set(findobj(FIG, 'Tag', 'EditTextMmag'),'UserData',Mmag);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextMmag'),'String',sprintf('%0.4g',Mmag));	%[A.m^2] Magnetic moment
   
  
%---Hyst1  
elseif strcmp(action, 'setVhyst1')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGhymu
   Vhyst1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextVhyst1'),'String'), '%f')*1e-6;	%[m^3] Volume
   Bs1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs1'),'String'), '%f');		%[Tesla] Magnetic induction (field strength)
   Mhyst1 = Bs1 * Vhyst1 / uo;							%[A.m^2] Magnetic moment
   mhyst1 = Vhyst1 * (SGhymu * 1000);							%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextVhyst1'),'UserData',Vhyst1);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextmhyst1'),'UserData',mhyst1);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextmhyst1'),'String',sprintf('%0.4g',mhyst1*1000));	%[g] mass 
   set(findobj(FIG, 'Tag', 'EditTextMhyst1'),'UserData',Mhyst1);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextMhyst1'),'String',sprintf('%0.4g',Mhyst1));		%[A.m^2] Magnetic moment
   
elseif strcmp(action, 'setMhyst1')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGhymu
   Mhyst1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextMhyst1'),'String'), '%f');		%[A.m^2] Magnetic moment
   Bs1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs1'),'String'), '%f');		%[Tesla] Magnetic induction (field strength)
   Vhyst1 = Mhyst1 * uo / Bs1;							%[m^3] Volume
   mhyst1 = Vhyst1 * (SGhymu * 1000);							%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextMhyst1'),'UserData',Mhyst1);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextmhyst1'),'UserData',mhyst1);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextmhyst1'),'String',sprintf('%0.4g',mhyst1*1000));	%[g] mass 
   set(findobj(FIG, 'Tag', 'EditTextVhyst1'),'UserData',Vhyst1);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextVhyst1'),'String',sprintf('%0.4g',Vhyst1*1e+6));	%[cm^3] Volume  
   
elseif strcmp(action, 'setmhyst1')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGhymu
   mhyst1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextmhyst1'),'String'), '%f')/1000;	%[kg] mass
   Bs1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs1'),'String'), '%f');		%[Tesla] Magnetic induction (field strength)
   Vhyst1 = mhyst1 / (SGhymu * 1000);							%[m^3] Volume
   Mhyst1 = Bs1 * Vhyst1 / uo;							%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextmhyst1'),'UserData',mhyst1);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextVhyst1'),'UserData',Vhyst1);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextVhyst1'),'String',sprintf('%0.4g',Vhyst1*1e+6));	%[cm^3] Volume  
   set(findobj(FIG, 'Tag', 'EditTextMhyst1'),'UserData',Mhyst1);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextMhyst1'),'String',sprintf('%0.4g',Mhyst1));		%[A.m^2] Magnetic moment
   
%---Hyst2  
elseif strcmp(action, 'setVhyst2')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGhymu
   Vhyst2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextVhyst2'),'String'), '%f')*1e-6;	%[m^3] Volume
   Bs2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs2'),'String'), '%f');		%[Tesla] Magnetic induction (field strength)
   Mhyst2 = Bs2 * Vhyst2 / uo;							%[A.m^2] Magnetic moment
   mhyst2 = Vhyst2 * (SGhymu * 1000);							%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextVhyst2'),'UserData',Vhyst2);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextmhyst2'),'UserData',mhyst2);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextmhyst2'),'String',sprintf('%0.4g',mhyst2*1000));	%[g] mass 
   set(findobj(FIG, 'Tag', 'EditTextMhyst2'),'UserData',Mhyst2);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextMhyst2'),'String',sprintf('%0.4g',Mhyst2));		%[A.m^2] Magnetic moment
   
elseif strcmp(action, 'setMhyst2')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGhymu
   Mhyst2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextMhyst2'),'String'), '%f');		%[A.m^2] Magnetic moment
   Bs2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs2'),'String'), '%f');		%[Tesla] Magnetic induction (field strength)
   Vhyst2 = Mhyst2 * uo / Bs2;							%[m^3] Volume
   mhyst2 = Vhyst2 * (SGhymu * 1000);							%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextMhyst2'),'UserData',Mhyst2);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextmhyst2'),'UserData',mhyst2);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextmhyst2'),'String',sprintf('%0.4g',mhyst2*1000));	%[g] mass 
   set(findobj(FIG, 'Tag', 'EditTextVhyst2'),'UserData',Vhyst2);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextVhyst2'),'String',sprintf('%0.4g',Vhyst2*1e+6));	%[cm^3] Volume  
   
elseif strcmp(action, 'setmhyst2')
   FIG = findobj('Tag','ControlWindow');
   %global uo SGhymu
   mhyst2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextmhyst2'),'String'), '%f')/1000;	%[kg] mass
   Bs2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs2'),'String'), '%f');		%[Tesla] Magnetic induction (field strength)
   Vhyst2 = mhyst2 / (SGhymu * 1000)							%[m^3] Volume
   Mhyst2 = Bs2 * Vhyst2 / uo							%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextmhyst2'),'UserData',mhyst2);			%[kg] mass
   set(findobj(FIG, 'Tag', 'EditTextVhyst2'),'UserData',Vhyst2);			%[m^3] Volume
   set(findobj(FIG, 'Tag', 'EditTextVhyst2'),'String',sprintf('%0.4g',Vhyst2*1e+6));	%[cm^3] Volume  
   set(findobj(FIG, 'Tag', 'EditTextMhyst2'),'UserData',Mhyst2);			%[A.m^2] Magnetic moment
   set(findobj(FIG, 'Tag', 'EditTextMhyst2'),'String',sprintf('%0.4g',Mhyst2));		%[A.m^2] Magnetic moment
   
   
 
   
%************************************************************
% Text Fields Process & Calculations
%************************************************************
elseif strcmp(action, 'eval')|strcmp(action, 'run')
   FIG = findobj('Tag','ControlWindow');
   %============================================================
   % Text Fields Read
   %============================================================
   %external simulation variable globalization for RUN callback
   global XYZfield_I umag uhyst1 uhyst2 Bmag Bs1 Bs2 Vmag Vhyst1 Vhyst2 uo Br Hc sw Jinv t fs H0 q0
   % === Simulation Parameters ===    
   %global fs
   sim_last = get(findobj(FIG, 'Tag', 'EditTextLast'),'UserData');	        	%[sec]
   sim_dur = get(findobj(FIG, 'Tag', 'EditTextDur'),'UserData');                    	%[rev]
   sim_sample =  sscanf(get(findobj(FIG, 'Tag', 'EditTextSample'),'String'), '%f');	%[]
   fs = get(findobj(FIG, 'Tag', 'EditTextfs'),'UserData');        	        	%[s]
   Ts = get(findobj(FIG, 'Tag', 'EditTextTs'),'UserData'); 	                	%[s]
   
   % === Orbit ===
   EYear = mod(sscanf(get(findobj(FIG, 'Tag', 'EditTextYear'),'String'), '%f'),2000); 	%[year2000]
   Epoch = sscanf(get(findobj(FIG, 'Tag', 'EditTextEpoch'),'String'), '%f');     	%[day]
   incl = sscanf(get(findobj(FIG, 'Tag', 'EditTextIncl'),'String'), '%f')*pi/180;	%[rad]
   raan = sscanf(get(findobj(FIG, 'Tag', 'EditTextRAAN'),'String'), '%f')*pi/180;	%[rad]
   e = sscanf(get(findobj(FIG, 'Tag', 'EditTextEcc'),'String'), '%f');     	%[]
   argp = sscanf(get(findobj(FIG, 'Tag', 'EditTextARGP'),'String'), '%f')*pi/180;	%[rad]
   Msc0 = sscanf(get(findobj(FIG, 'Tag', 'EditTextMsc0'),'String'), '%f')*pi/180;	%[rad]
   n = get(findobj(FIG, 'Tag', 'EditTextn'),'UserData'); 		            	%[rev/day]
   h = get(findobj(FIG, 'Tag', 'EditTexth'),'UserData');		             	%[km]
   
   % === Spacecraft ===   
   %global msc Ixx
   Ixx = sscanf(get(findobj(FIG, 'Tag', 'EditTextIxx'),'String'), '%f')*1e-6;           %[kg.m²] 
   Iyy = sscanf(get(findobj(FIG, 'Tag', 'EditTextIyy'),'String'), '%f')*1e-6;           %[kg.m²]
   Izz = sscanf(get(findobj(FIG, 'Tag', 'EditTextIzz'),'String'), '%f')*1e-6;           %[kg.m²]
   Ixy = sscanf(get(findobj(FIG, 'Tag', 'EditTextIxy'),'String'), '%f')*1e-6;           %[kg.m²] product of inertia
   Ixz = sscanf(get(findobj(FIG, 'Tag', 'EditTextIxz'),'String'), '%f')*1e-6;           %[kg.m²] product of inertia
   Iyz = sscanf(get(findobj(FIG, 'Tag', 'EditTextIyz'),'String'), '%f')*1e-6;           %[kg.m²] product of inerti
   Isc = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];                                 	%[kg.m²] Moment of inertia tensor
   Jinv = Isc^(-1);									% Inverse of tensor   msc = sscanf(get(findobj(FIG, 'Tag', 'EditTextmsc'),'String'), '%f');           	%[kg] mass
   msc = sscanf(get(findobj(FIG, 'Tag', 'EditTextmsc'),'String'), '%f');           	%[kg] total S/C mass
   side = sscanf(get(findobj(FIG, 'Tag', 'EditTextcsc'),'String'), '%f');           	%[m] base width
   lsc = sscanf(get(findobj(FIG, 'Tag', 'EditTextlsc'),'String'), '%f');		%[m] length
   A1sc = side^2;									%[m^2] base face area
   A2sc = side*lsc;									%[m^2] long face area
   wallsc = sscanf(get(findobj(FIG, 'Tag', 'EditTextwallsc'),'String'), '%f')/1000;	%[m]
   Cd = get(findobj(FIG, 'Tag', 'EditTextCd'),'UserData');              		%[] Drag coef
   Cb = get(findobj(FIG, 'Tag', 'EditTextCb'),'UserData');              		%[] Ballistic coef
   PVeff = sscanf(get(findobj(FIG, 'Tag', 'EditTextPVeff'),'String'), '%f')/100;   	%[] Solar cells efficiency
   PVpack1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextPVpack1'),'String'), '%f')/100;	%[] packing on base
   PVpack2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextPVpack2'),'String'), '%f')/100;	%[] packing on sides
   
   Qin = sscanf(get(findobj(FIG, 'Tag', 'EditTextQin'),'String'), '%f');
   
   Xbody = [1 0 0]';
   Ybody = [0 1 0]';
   Zbody = [0 0 1]';
   Nbody = [Xbody; Ybody; -Xbody; -Ybody; Zbody; -Zbody];         	%[]
   PVpacking = [PVpack2 PVpack2 PVpack2 PVpack2 PVpack1 PVpack1]';	%[]

   
   % === Stabilization ===   
   %global wspin Bmag Bs1 Bs2 Br Hc umag uhyst1 uhyst2 Vmag Vhyst1 Vhyst2% for simulink model

   %wspin = sscanf(get(findobj(FIG, 'Tag', 'EditTextwspin'),'String'), '%f')*pi/180; 	%[]
   %inclbody = sscanf(get(findobj(FIG, 'Tag', 'EditTextibody'),'String'), '%f')*pi/180;	%[]
   %wbody = sscanf(get(findobj(FIG, 'Tag', 'EditTextwbody'),'String'), '%f')*pi/180;   	%[]
   
   Bmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextBmag'),'String'), '%f');	%[Tesla] Magnetic induction (saturation field strength)
   Vmag = get(findobj(FIG, 'Tag', 'EditTextVmag'),'UserData');			%[m^3] Volume 
   umagx = sscanf(get(findobj(FIG, 'Tag', 'EditTextumagx'),'String'), '%f');	%[] Orientation in Body-frame components
   umagy = sscanf(get(findobj(FIG, 'Tag', 'EditTextumagy'),'String'), '%f');	%[] Orientation in Body-frame components
   umagz = sscanf(get(findobj(FIG, 'Tag', 'EditTextumagz'),'String'), '%f');	%[] Orientation in Body-frame components
   umag = [umagx; umagy ; umagz]./ norm([umagx; umagy ; umagz]);		%normalized   
   mmag = sscanf(get(findobj(FIG, 'Tag', 'EditTextmmag'),'String'), '%f')/1000;	%[kg] mass
 
   
   Bs1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs1'),'String'), '%f');	%[Tesla] Magnetic induction (saturation field strength)
   Vhyst1 = get(findobj(FIG, 'Tag', 'EditTextVhyst1'),'UserData');		%[m^3] Volume
   uhyst1x = sscanf(get(findobj(FIG, 'Tag', 'EditTextuhyst1x'),'String'), '%f');	%[] Orientation in Body-frame components
   uhyst1y = sscanf(get(findobj(FIG, 'Tag', 'EditTextuhyst1y'),'String'), '%f');	%[] Orientation in Body-frame components
   uhyst1z = sscanf(get(findobj(FIG, 'Tag', 'EditTextuhyst1z'),'String'), '%f');	%[] Orientation in Body-frame components
   uhyst1 = [uhyst1x; uhyst1y; uhyst1z]./ norm([uhyst1x; uhyst1y; uhyst1z]);		%normalized   
   mhyst1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextmhyst1'),'String'), '%f')/1000;	%[kg] mass
   
   Bs2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs2'),'String'), '%f');	%[Tesla] Magnetic induction (saturation field strength)
   Vhyst2 = get(findobj(FIG, 'Tag', 'EditTextVhyst2'),'UserData');		%[m^3] Volume
   uhyst2x = sscanf(get(findobj(FIG, 'Tag', 'EditTextuhyst2x'),'String'), '%f');	%[] Orientation in Body-frame components
   uhyst2y = sscanf(get(findobj(FIG, 'Tag', 'EditTextuhyst2y'),'String'), '%f');	%[] Orientation in Body-frame components
   uhyst2z = sscanf(get(findobj(FIG, 'Tag', 'EditTextuhyst2z'),'String'), '%f');	%[] Orientation in Body-frame components
   uhyst2 = [uhyst2x; uhyst2y; uhyst2z]./ norm([uhyst2x; uhyst2y; uhyst2z]);		%normalized   
   mhyst2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextmhyst2'),'String'), '%f')/1000;	%[kg] mass
   
   Br = 0.35;	%Remanence induction
   Hc = 1.59; 	%Coercive field
   sw = get(findobj(FIG, 'Tag', 'ButtonSW1'),'UserData');	%Hysteresis model switch {1=delay+ramp Br/H,   2=delay+ramp Bs/2H,   3=delay only}
   
   q1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq1'),'String'), '%f');	%[] quaternion
   q2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq2'),'String'), '%f');	%[] quaternion
   q3 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq3'),'String'), '%f');	%[] quaternion
   q4 = sscanf(get(findobj(FIG, 'Tag', 'EditTextq4'),'String'), '%f');	%[] quaternion
   q0 = [q1; q2; q3; q4]./sqrt(q1^2+q2^2+q3^2+q4^2); % intitial quaternion
   
   wbx = sscanf(get(findobj(FIG, 'Tag', 'EditTextwbx'),'String'), '%f');%[rad/s] 
   wby = sscanf(get(findobj(FIG, 'Tag', 'EditTextwby'),'String'), '%f');%[rad/s] 
   wbz = sscanf(get(findobj(FIG, 'Tag', 'EditTextwbz'),'String'), '%f');%[rad/s]
   wb0 = [wbx; wby; wbz];		% initial rotation resolved in body coordinates

   H0 = Isc * wb0;			% initial momentum resolvein in body coordinates
        
   
   % === Ground Stations ===   
   NAME_GS1 = get(findobj(FIG, 'Tag', 'EditTextGS1'),'String');                       	%[] 
   LON_GS1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS1lon'),'String'), '%f')*pi/180;	%[rad]
   LAT_GS1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS1lat'),'String'), '%f')*pi/180;	%[rad]
   MEL_GS1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS1el'),'String'), '%f')*pi/180;	%[rad]
   f1 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS1fq'),'String'), '%f')*1000000;     	%[Hz]
   
   NAME_GS2 = get(findobj(FIG, 'Tag', 'EditTextGS2'),'String');                        	%[]
   LON_GS2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS2lon'),'String'), '%f')*pi/180;	%[rad]
   LAT_GS2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS2lat'),'String'), '%f')*pi/180;	%[rad]
   MEL_GS2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS2el'),'String'), '%f')*pi/180;	%[rad]
   f2 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS2fq'),'String'), '%f')*1000000;     	%[Hz]
   
   NAME_GS3 = get(findobj(FIG, 'Tag', 'EditTextGS3'),'String');                        	%[]
   LON_GS3 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS3lon'),'String'), '%f')*pi/180;	%[rad]
   LAT_GS3 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS3lat'),'String'), '%f')*pi/180;	%[rad]
   MEL_GS3 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS3el'),'String'), '%f')*pi/180;	%[rad]
   f3 = sscanf(get(findobj(FIG, 'Tag', 'EditTextGS3fq'),'String'), '%f')*1000000;     	%[Hz]
   
   if length(NAME_GS1)< 3
      NAME_GS1=[NAME_GS1,' ',' ','1'];
   end
   if length(NAME_GS2)< 3
      NAME_GS2=[NAME_GS2,' ',' ','2'];
   end
   if length(NAME_GS3)< 3
      NAME_GS3=[NAME_GS3,' ',' ','3'];
   end
   NAME_GS = [NAME_GS1(1:3); NAME_GS2(1:3); NAME_GS3(1:3)];		%[string] ground station tag name
   LON_GS = [LON_GS1; LON_GS2; LON_GS3];				%[rad] Longitude
   LAT_GS = [LAT_GS1; LAT_GS2; LAT_GS3];				%[rad] Latitude
   MEL_GS = [MEL_GS1; MEL_GS2; MEL_GS3];				%[rad] Minimum Elevation 
   N_GS = length(LON_GS);
   
   
   %save cubesim/current.mat %,...
      %sim_last, sim_dur, sim_sample, fs, Ts, ...
      %EYear, Epoch, incl, raan, e, argp, Msc0, n, h, ...
      %Ixx, Iyy, Izz, Ixy, Ixz, Iyz, msc, side, lsc, ...
      %A1sc, A2sc, wallsc, Cd, Cb, PVeff, PVpack1, PVpack2, ...
      %Bmag, Vmag, umagx, umagy, umagz, umag, ...
      %Bs1, Vhyst1, uhyst1x, uhyst1y, uhyst1z, uhyst1, ...
      %Bs2, Vhyst2, uhyst2x, uhyst2y, uhyst2z, uhyst2, ...
      %sw, q1, q2, q3, q4, wbx, wby, wbz,
      %NAME_GS1, LON_GS1, LAT_GS1, MEL_GS1, f1, ...
      %NAME_GS2, LON_GS2, LAT_GS2, MEL_GS2, f2, ...
      %NAME_GS3, LON_GS3, LAT_GS3, MEL_GS3, f3;

   
   %============================================================
   % S/C Orbit Properties (at Epoch)
   %============================================================
   a = (sqrt(GMe)/(2*pi)*Te/n)^(2/3);	%[km] semi-major axis
   rp = a*(1-e);                 	%[km] perigee
   ra = a*(1+e);                 	%[km] apogee
   c = a - rp;                      	%[km] orbit center distance from focii
   b = sqrt(a^2-c^2);               	%[km] semi-minor axis
   h = rp - Re;                     	%[km] altitude at perigee

   Psc = 2*pi*sqrt(a^3/GMe);         	%[s] orbital period
   wsc = sqrt(GMe/a^3);             	%[rad/s] mean angular velocity

   Va = sqrt(GMe*(2/rp-1/a));        	%[km/s] velocity at apogee
   Vp = sqrt(GMe*(2/ra-1/a));        	%[km/s] velocity at perigee
   Vmean = sqrt(GMe/a);             	%[km/s] mean circular velocity
   Vesc =  sqrt(2*GMe/a);           	%[km/s] mean escape velocity

   Eclipse_max = asin(Re/rp)/pi * Psc;				%[s]
   Eclipse_min = asin(sqrt( (Re/ra)^2 - (sin(incl))^2 ))/pi * Psc;	%[s]
   if ~isreal(Eclipse_min)
      Eclipse_min = 0;
   end
   t_VE0 = (Epoch-VERNAL_EQX(EYear)) * Te;	%[sec] time from vernal equinox

   %============================================================
   % S/C Orbital Perturbations (at Epoch)
   %============================================================
   %Other Bodies
   Draan_sun = -0.00338*cos(incl)/n*(pi/180)/Te;			%[rad/s]
   Draan_moon = -0.00154*cos(incl)/n*(pi/180)/Te;			%[rad/s]
   Draan_geo = -1.5*wsc*J2*(Re/a)^2*cos(incl)/(1-e^2)^2;		%[rad/s]
   Draan = Draan_sun + Draan_moon + Draan_geo;
   
   Dargp_sun = 0.00169*(4-5*(sin(incl))^2)/n*(pi/180)/Te;		%[rad/s]
   Dargp_moon = 0.00077*(4-5*(sin(incl))^2)/n*(pi/180)/Te;		%[rad/s]
   Dargp_geo = 0.75*wsc*J2*(Re/a)^2*(4-5*(sin(incl))^2)/(1-e^2)^2;	%[rad/s]
   Dargp = Dargp_sun + Dargp_moon + Dargp_geo;
   
   %Air Drag
   [rhoa, ScaleHt, rhoa_max, rhoa_min] = getAtmDensity(h); %[kg/m³], [km]  !warning, dimensions are not consistent

   adrag = -0.5*rhoa/Cb*(Vp*1000)^2; 		%[m/s^2]
   Drag = adrag*msc;                         	%[N]
   orb_decay = adrag/(2*pi*a*1000)*(Te^2);       	%[rev/day^2]
   Da_rev = -2*pi/Cb*rhoa*(a*1000)^2/1000;	%[km/rev] variation of semi-major axis
   OrbitLife = -ScaleHt/Da_rev *Psc/Te;              %[day] estimated orbit lifetime
 
   adrag_max = -0.5*rhoa_max/Cb*(Vp*1000)^2; 		%[m/s^2]
   Drag_max = adrag_max*msc;                         	%[N]
   orb_decay_max = adrag_max/(2*pi*a*1000)*(Te^2);       	%[rev/day^2]
   Da_rev_max = -2*pi/Cb*rhoa_max*(a*1000)^2/1000;	%[km/rev] variation of semi-major axis
   Life_min = -ScaleHt/Da_rev_max *Psc/Te;              %[day] estimated orbit lifetime

   adrag_min = -0.5*rhoa_min/Cb*(Vp*1000)^2; 		%[m/s^2]
   Drag_min = adrag_min*msc;                         	%[N]
   orb_decay_min = adrag_min/(2*pi*a*1000)*(Te^2);       	%[rev/day^2]
   Da_rev_min = -2*pi/Cb*rhoa_min*(a*1000)^2/1000;	%[km/rev] variation of semi-major axis
   Life_max = -ScaleHt/Da_rev_min *Psc/Te;              %[day] estimated orbit lifetime

   %============================================================
   % Communication Link (At Epoch, circular orbit approx)
   %============================================================
   % Downlink ranges 
   S_maxp = sqrt(rp^2-Re^2);		%[km] maximum slant range perigee
   S_minp = rp - Re; 			%[km] minimum range perigee
   S_maxa = sqrt(ra^2-Re^2);		%[km] maximum slant range apogee
   S_mina = ra - Re; 			%[km] minimum range apogee

   SR_GS = getRange(MEL_GS, rp);		%[km] Nominal Slant Range

   % Downlink times
   [S, BetaLink] = getRange(0*pi/180, a);	%[km],[rad]
   DLmax = BetaLink /(wsc-we);			%[s] max downlink duration (without orbital perturbation correction)
   DLday = 2*BetaLink /(2*pi) *n;		%[1/day] downlink probability per day
   dLON = (-Draan+we)*Psc;			%[rad/rev] node spacing on earth (per S/C revolution)
      
   %============================================================
   % Magnetic Stabilization
   %============================================================
   % Alnico-5 Permanent Magnet (in bars 1/8 x 1/8 x 4")
   mmag = Vmag * SGalnico * 1000;	%[kg] mass 
   Mmag = Bmag * Vmag / uo;		%[A.m^2] Magnetic moment       
     
   % Hymu-80 Hysteresis Bars (in bars 1/8 x 1/8 x 4")
   mhyst1 = Vhyst1 * SGhymu * 1000;	%[kg] mass 
   Mhyst1 = Bs1 * Vhyst1 / uo;		%[A.m^2] Magnetic moment  
   mhyst2 = Vhyst2 * SGhymu * 1000;	%[kg] mass 
   Mhyst2 = Bs2 * Vhyst2 / uo;		%[A.m^2] Magnetic moment 

   % simulation
   %wo = wspin;
   %thetao = 4*pi/3;		%attitude angle to magnetic field
  
   Mres = Mmag.*umag + Mhyst1.*uhyst1 + Mhyst2.*uhyst2;
   Stab_err = acos(dot(Mres, umag)/norm(Mres));	%steady-state pointing error
   Stab_res = Stab_err + asin(Hc/B_EARTH*uo);   %pointing error accounting stabilization resolution


   %============================================================
   % Text Fields Update
   %============================================================
   set(findobj(FIG, 'Tag', 'EditTextDecay'),'String',sprintf('%0.2g',-orb_decay));
   set(findobj(FIG, 'Tag', 'EditTexta'),'String',sprintf('%0.5g',a));
   %set(findobj(FIG, 'Tag', 'EditTextb'),'String',sprintf('%0.5g',b));
   %set(findobj(FIG, 'Tag', 'EditTextc'),'String',sprintf('%0.5g',c));
   set(findobj(FIG, 'Tag', 'EditTextra'),'String',sprintf('%0.5g',ra));
   set(findobj(FIG, 'Tag', 'EditTextrp'),'String',sprintf('%0.5g',rp));
   set(findobj(FIG, 'Tag', 'EditTextVa'),'String',sprintf('%0.5g',Va));
   set(findobj(FIG, 'Tag', 'EditTextVp'),'String',sprintf('%0.5g',Vp));
   set(findobj(FIG, 'Tag', 'EditTextEclmax'),'String',sprintf('%0.5g',Eclipse_max/60));
   set(findobj(FIG, 'Tag', 'EditTextEclmin'),'String',sprintf('%0.5g',Eclipse_min/60));
   set(findobj(FIG, 'Tag', 'EditTextPsc'),'String',sprintf('%0.5g',Psc/60));
   set(findobj(FIG, 'Tag', 'EditTextdLON'),'String',sprintf('%0.5g',dLON*180/pi));
   set(findobj(FIG, 'Tag', 'EditTextDraan'),'String',sprintf('%0.5g',Draan*Te*180/pi));
   set(findobj(FIG, 'Tag', 'EditTextDargp'),'String',sprintf('%0.4g',Dargp*Te*180/pi));
   set(findobj(FIG, 'Tag', 'EditTextDrag'),'String',sprintf('%0.5g',-Drag*1e+9));
   set(findobj(FIG, 'Tag', 'EditTextDecay2'),'String',sprintf('%0.1g',-Da_rev/Psc*Psun));
   set(findobj(FIG, 'Tag', 'EditTextLife'),'String',sprintf('%3.0f',OrbitLife));
   set(findobj(FIG, 'Tag', 'EditTextStabRes'),'String',sprintf('%0.4g',Stab_res*180/pi));
   
   save C:\Users\spaceconcordia\Desktop\cubesim\current.mat 
   
   fprintf('Parameters updated and saved in cubesim/current.mat\n')
   LocalDisplay('Parameters updated and saved in cubesim/current.mat', 'nl')
   
   %************************************************************
   % Elliptic Orbit Simulation (including perturbations)
   %************************************************************
   if strcmp(action, 'run')

      fprintf('\n Simulation Process |')
      LocalDisplay('Simulation started. May take a while.', 'nl');
      LocalDisplay('Simulation progress |', 'nl');
      WB = WAITBAR(0.0,'Simulation Progress');
   
      
      %============================================================
      % Time and Sampling
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.01,WB);
      
      %global t_VE
      circle = linspace(-pi,pi,91);			%[rad]
      t = linspace(0, sim_last, sim_sample);		%[sec] sample time vector
      t_VE = t_VE0 + t;					%[sec] time from vernal equinox
      
      %SolarYearPhase = mod((EYear+Epoch/365.24 - 2000.7), 11.1)/11.1 .*2*pi;
      SolarDayPhase = (t_VE./Te - APHELION(EYear))/365.24 .*2*pi;
      %SolarFlux = (Gsolar_min + (Gsolar_max-Gsolar_min)*cos(SolarYearPhase)) / (1.0004 + 0.0334*cos(SolarDayPhase))
      SolarFlux = Gsolar ./ (1.0004 + 0.0334.*cos(SolarDayPhase));

      %============================================================
      % S/C Orbital Calculations
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.05,WB);
      
      Msc = Msc0 + wsc.*t;					%[rad] Mean anomaly
      NUsc = Msc + 2*e.*sin(Msc) + 1.25*e^2.*sin(2*Msc);	%[rad] True anomaly (approximation for small eccentricity)
      %Esc = acos(e+cos(NUsc))./(1+e.*cos(NUsc));		%[rad] Eccentric anomaly (based on approx)
      %Esc = 2*atan(tan(NUsc./2)./sqrt((1+e)/(1-e)));		%[rad] Eccentric anomaly (based on approx)
      %Rsc = a*(1 - e.*cos(Esc));				%[km] Distance from focus (based on approx)
      %Rsc = (rp*(1+e))./(1 + e.*cos(NUsc));			%[km] Distance from focus (based on approx)
      %Vsc = sqrt(GMe.*(2./Rsc-1/a));				%[km/s] Velocity
      %usc = argp + Dargp.*t + NUsc;				%[rad] Angle to ascending node
      
      if Epoch < 10
         EpochStr = sprintf('%2.0f00%0.6f',EYear,Epoch);
      elseif Epoch < 100
         EpochStr = sprintf('%2.0f0%0.6f',EYear,Epoch);
      else
         EpochStr = sprintf('%2.0f%0.6f',EYear,Epoch);
      end
      ElementStruct = struct('Epoch_time', EpochStr,...
         'Inclination',sprintf('%0.8g',incl*180/pi),...
         'RA_of_node',sprintf('%0.8g',raan*180/pi),...
         'Eccentricity',sprintf('%0.8g',e),...
         'Arg_of_perigee',sprintf('%0.8g',argp*180/pi),...
         'Mean_anomaly',sprintf('%0.8g',Msc0*180/pi),...
         'Mean_motion',sprintf('%0.8g',n));
      
      LocalDisplay(' Ephemeris ', 'add')      
      XYZsc_geo = getEphemeris(ElementStruct,Cb,t_VE);		%Propagate orbit elements
      IJKsc_geo = XYZsc_geo./([1;1;1]*LEN(XYZsc_geo));
      
      %============================================================
      % Sun Orbital Calculations (Sun viewed as orbiting object around Earth)

      %============================================================
      
      

      
      
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.10,WB);
      
      %Msun0 = pi/2;
      %Msun = Msun0 + wsun.*t_VE;                                    	%[rad] Mean anomaly
      
      %NUsun = Msun + 2*ECCsun.*sin(Msun) + 1.25*ECCsun^2.*sin(2*Msun);	%[rad] True anomaly (approximation for small eccentricity)
      %Esun = acos(ECCsun+cos(NUsun))./(1+ECCsun.*cos(NUsun);        	%[rad] Eccentric anomaly (based on approx)
      %Esun = 2*atan(tan(NUsun./2)./sqrt((1+ECCsun)/(1-ECCsun)));      	%[rad] Eccentric anomaly (based on approx)
      %Rsun = a*(1 - ECCsun.*cos(Esun));                             	%[km] Distance from focus (based on approx)
      %Rsun = (Rpsun*(1+ECCsun))./(1 + ECCsun.*cos(NUsun));           	%[km] Distance from focus (based on approx)
      %usun = pi/2 + NUsun;                                           	%[rad] Angle to ascending node
      
      % Greenwch sideral time (WERTZ p.803)
      [D1900 C1900]= getD1900(EYear, Epoch + t./Te);	%[day] Julian days since Jan 0.5 1900
      UT = (t./Te + mod(Epoch,1))*360;			%[deg] Universal time
      GST = mod((99.6910 + 36000.7689.*C1900 + 0.0004.*C1900.^2 + UT),360) .*pi/180;	%[rad]
       
      % Compute longitude and anomaly measured in the ecliptic from mean equinox of date (WERTZ p.141)
      LONsun_ecl = (279.696678 + 0.9856473354.*D1900 + 2.267e-13.*D1900.^2) .*pi/180; %[rad] Mean longitude
      Msun_ecl = (358.475845 + 0.985600267.*D1900 - 1.12e-13.*D1900.^2 -7e-20.*D1900.^3) .*pi/180; %[rad]
      dLONsun_ecl = (1.918.*sin(Msun_ecl) + 0.02.*sin(2.*Msun_ecl)) .*pi/180; %[rad] True longitude correction

      % True anomaly (small eccentricity approx) and distance Earth-sun
      NUsun_ecl = Msun_ecl + 2*ECCsun.*sin(Msun_ecl) + 1.25*ECCsun^2.*sin(2.*Msun_ecl);	%[rad] True anomaly (approximation for small eccentricity)
      Rsun = AU*(1-ECCsun^2) ./ (1 + e.*cos(NUsun_ecl));	%[km]
      
      % Position of sun resolven in ecliptic and Earth referential
      XYZsun_ecl = LL2XYZ([LONsun_ecl;zeros(1,sim_sample)], Rsun);
      IJKsun_ecl = XYZsun_ecl./ ([1;1;1]*LEN(XYZsun_ecl));
      XYZsun_geo = zeros(3,sim_sample);  	% spacecraft vector
      for j=1:sim_sample
         ROTecl_geo = Rz(GST(j))*Rx(-incle);
         XYZsun_geo(:,j) = ROTecl_geo * XYZsun_ecl(:,j);
      end
      IJKsun_geo = XYZsun_geo./ ([1;1;1]*LEN(XYZsun_geo));         
      
      
      %============================================================
      % Ground Traces (2-D plot)
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.15,WB);
      % Spacecraft
      [LLsc Rsc] = XYZ2LL(XYZsc_geo);
      %incl, raan + Draan.*t, usc ,t_VE);	%[rad]

      LATsc = LLsc(2,:);                                  	%[rad]
      LONsc = LLsc(1,:);                                  	%[rad] including perturbations
      %LONsc = mod(LONsc+pi,2*pi*ones(size(LONsc)))-pi;    	%[rad] normalisation @ ±Pi
      %LLsc = [LONsc; LATsc];                             	%[rad]
      
      % Sun
      %LLsun = getTrace(-incle, 0, usun, t_VE);          	%[rad]
      LLsun = XYZ2LL(XYZsun_geo);				%[rad]
      LATsun = LLsun(2,:);                                	%[rad]
      LONsun = LLsun(1,:);                                	%[rad] 
      %LONsun = mod(LONsun+pi,2*pi*ones(size(LONsun)))-pi;	%[rad] normalisation @ ±Pi
      %LLsun = [LONsun; LATsun];                          	%[rad]
      
      % Eclipse Shadow (Sun Ray Vector)
      LATray = -LLsun(2,:);                              	%[rad]
      LONray = LLsun(1,:) +pi ;                          	%[rad] 
      LLray = [LONray; LATray];                          	%[rad]
      
      % Orbit Pole (Normal to Orbital Plane)
      LATorbpole = (pi/2-incl);                                      	%[rad]
      Lorbpole = -pi/2;                                              	%[rad]
      LONorbpole = raan + Lorbpole + Draan.*t - we.*(t_VE);        	%[rad]
      LONorbpole = mod(LONorbpole+pi,2*pi*ones(size(LONorbpole)))-pi;	%[rad] normalisation @ ±Pi
      LLorbpole = [LONorbpole; LATorbpole*ones(1,sim_sample)];       	%[rad]
      clear Lorbpole
     
      % Ground Station Circles (corrected for latitude map stretching deformation)
      S_maxp_V = S_maxp*ones(1,length(circle));  		%[km]
      SR_GS_V = SR_GS*ones(1,length(circle));  	 		%[km]
      LLgs_cir0 = zeros(2*N_GS,length(circle));  		%[rad] GS circle, Elev=0, Slant Range perigee
      LLgs_cirMEL = zeros(2*N_GS,length(circle));		%[rad] GS circle, MEL, Slant Range semi-major axis
     for k=1:N_GS
        LLgs_cir0(2*k-1:2*k,:) = identSC([LON_GS(k);LAT_GS(k)], [circle; 0*circle], S_maxp_V);			%[rad]
         LLgs_cirMEL(2*k-1:2*k,:) = identSC([LON_GS(k);LAT_GS(k)], [circle; MEL_GS(k)*ones(1,length(circle))], SR_GS_V(k,:));	%[rad]
     end
      
                  
      %============================================================
      % Attitude simulation (passive magnetic stabilization)
      %============================================================  
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.20,WB,'Simulation Progress (Simulink model)');
      %=== Compute Earth magnetic field vector at S/C location ===
      XYZfield_geo = getMagField(XYZsc_geo);
      
      %ADDED
      disp('LEN');
      disp(LEN(XYZfield_geo));
          
      disp('Argument 1');
      disp(-IJKsc_geo);
      
      disp('Argument 2');
      (XYZfield_geo./([1;1;1]*LEN(XYZfield_geo)))
      %END
      
      
      
      %alpha = acos(dot(-IJKsc_geo, (XYZfield_geo./([1;1;1]*LEN(XYZfield_geo))))); %field angle relative to Nadir
      
      alpha = acos(dot(-IJKsc_geo, (XYZfield_geo./([1;1;1]*LEN(XYZfield_geo))))); %field angle relative to Nadir
      
      [LLmagpt Rmagpt] = getIntersect(LLsc, Rsc, XYZfield_geo);
      XYZmagtr_geo = LL2XYZ(LLmagpt) * Re;	% 3-D pointing trace
       
      XYZfield_I = zeros(3,sim_sample);
      for k=1:sim_sample
         ROTinertial_geo = Rz(we*t_VE(k))*Rx(incle);
         XYZfield_I(:,k) = ROTinertial_geo^(-1) * XYZfield_geo(:,k);
      end
      XYZfield_I;
      
      %=== Simulate S/C attitude
      LocalDisplay(' Simulink ', 'add')      
      [t_sim X_sim OUT_sim] = sim('magSim3', [0 t(2:sim_sample-1) t(sim_sample)]);
      
      if length(t_sim) < sim_sample
         t_sim = [t_sim; NaN];
         X_sim = [X_sim; NaN*zeros(size(X_sim(1,:)))];
         OUT_sim = [OUT_sim; NaN*zeros(size(OUT_sim(1,:)))];
      end
            
      QUATsc_I = OUT_sim(:,1:4)';
      wb_body = OUT_sim(:,5:7)';	%[deg/s]
      wb_len = OUT_sim(:,8)';		%[deg/s]
      wspin = OUT_sim(:,9)';		%[deg/s]
      Hb_I = OUT_sim(:,10:12)';
      Xbody_I = OUT_sim(:,13:15)';
      Ybody_I = OUT_sim(:,16:18)';
      Zbody_I = OUT_sim(:,19:21)';
           
      %=== Compute magnet angle relative to Earth field
      umag_I = zeros(3,sim_sample);
      for k=1:sim_sample
         ROTbody_I = [Xbody_I(:,k), Ybody_I(:,k), Zbody_I(:,k)];
         umag_I(:,k) = ROTbody_I * umag;
      end
      phi = acos( dot(XYZfield_I./([1;1;1]*LEN(XYZfield_I)), umag_I));	%angle between Field and magnet
     
      %=== Rotate attitude and Field to Earth coordinates
      Zbody_geo = zeros(3,sim_sample);
      Ybody_geo = zeros(3,sim_sample);
      Xbody_geo = zeros(3,sim_sample);
      umag_geo = zeros(3,sim_sample);
      for k=1:sim_sample
         ROTinertial_geo = Rz(we*t_VE(k))*Rx(incle);
         Xbody_geo(:,k) = ROTinertial_geo * Xbody_I(:,k);
         Ybody_geo(:,k) = ROTinertial_geo * Ybody_I(:,k);
         Zbody_geo(:,k) = ROTinertial_geo * Zbody_I(:,k);
         umag_geo(:,k) = ROTinertial_geo * umag_I(:,k);
      end
      Nbody_geo = [Xbody_geo; Ybody_geo; -Xbody_geo; -Ybody_geo; Zbody_geo; -Zbody_geo];
      
      %=== Compute magnet angle relative to Nadir, and pointing trace
      disp('--------->XYZsc<----------');
      disp(size(XYZsc_geo));
      beta = acos(dot(-IJKsc_geo, umag_geo));
      
      [LLattpt Rattpt] = getIntersect(LLsc, Rsc, umag_geo); %Magnet pointing trace
      XYZatttr_geo = LL2XYZ(LLattpt) * Re;	%  3-D pointing trace
      
      %============================================================
      % Ground Station Traces (3-D plot)
      %============================================================      
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.40,WB,'Simulation Progress');
      % Ground Stations
      XYZgs_geo = zeros(3*N_GS,1);                    	% G/S location
      XYZgs_cir0_geo = zeros(3*N_GS,length(circle));   	% G/S circle, Elev=0, Slant Range perigee
      XYZgs_cirMEL_geo = zeros(3*N_GS,length(circle)); 	% G/S circle, MEL, Slant Range semi-major axis
      for k=1:N_GS
         XYZgs_geo(3*k-2:3*k,:) = LL2XYZ([LON_GS(k); LAT_GS(k)]) * Re;
         XYZgs_cir0_geo(3*k-2:3*k,:) = LL2XYZ(LLgs_cir0(2*k-1:2*k,:)) * Re;
         XYZgs_cirMEL_geo(3*k-2:3*k,:) = LL2XYZ(LLgs_cirMEL(2*k-1:2*k,:)) * Re;
      end
      
      % Earth References (X:= [0, 0], Y:= [90, 0], Z:= [0, 90] (North Pole)
      XYZequa_geo = [Re*cos(circle); Re*sin(circle); zeros(size(circle))]; 	% Equator (LAT = 0)
      XYZgw_geo = [abs(Re*cos(circle)); zeros(size(circle)); Re*sin(circle)];	% Greenwich (LON = 0)
      XYZzero_geo = Re*[1; 0; 0];						% [LON, LAT] = [0, 0]
            
      %============================================================
      % Position Vectors relative to Sun
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.45,WB);
      Res_graph = 3*a;
      XYZgeo_sun = -Res_graph*IJKsun_ecl;	% earth vectorfrom fictive close sun

      XYZzero_sun = zeros(3,sim_sample);	% [LON, LAT] = [0, 0]
      XYZsc_sun = zeros(3,sim_sample);  	% spacecraft vector
      for k=1:sim_sample
         ROTecl_geo = Rz(GST(k))*Rx(-incle);
         XYZzero_sun(:,k) = (ROTecl_geo' * XYZzero_geo) + XYZgeo_sun(:,k);
         XYZsc_sun(:,k) = (ROTecl_geo' * XYZsc_geo(:,k)) + XYZgeo_sun(:,k);
      end
      
      
      %============================================================
      % Ground Station Tracking Elements
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.50,WB);
      Azimut_GS = 	zeros(N_GS, sim_sample);
      Elevation_GS = 	zeros(N_GS, sim_sample);
      Range_GS = 	zeros(N_GS, sim_sample);
      AGLgs2ant = 	zeros(N_GS, sim_sample);
      Rgs2tg =    	zeros(N_GS, sim_sample);
      Doppler_GS = 	zeros(N_GS, sim_sample);
      DpRate_GS = 	zeros(N_GS, sim_sample);
      f0 = [f1 f2 f3];
      
      filter = 1;		% filter data when out-of view if flag = 1
      for k=1:N_GS
         %G/S tracking data
         [AziEl_GSi, Range_GSi, XYZgs2sc_i] = trackSC([LON_GS(k); LAT_GS(k)], [LONsc; LATsc], Rsc, filter);
         Azimut_GS(k,:) = AziEl_GSi(1,:);
         Elevation_GS(k,:) = AziEl_GSi(2,:);
         Range_GS(k,:) = Range_GSi;
         
         %S/C antenna magnetic pointing error to ground station location
         %[AGLgs2ant_i, Rgs2ant_i] = getPtError([LON_GS(k); LAT_GS(k)], LLmagpt, Rmagpt, Range_GSi);
         [AGLgs2ant_i, Rgs2ant_i] = getPtError([LON_GS(k); LAT_GS(k)], LLattpt, Rattpt, Range_GSi);
         AGLgs2ant(k,:) = AGLgs2ant_i;
         Rgs2ant(k,:) = Rgs2ant_i;
         
         %Link Doppler shift
         [Doppler_i, DpRate_i] = getDoppler(XYZgs2sc_i, fs, f0(k));
         Doppler_GS(k,:) = Doppler_i;
         DpRate_GS(k,:) = DpRate_i;
         % filter for LOS data only
         for j=1:sim_sample
            if Elevation_GS(k,j) >= 0
            else
               Doppler_GS(k,j) = NaN;
               DpRate_GS(k,j) = NaN;
            end
         end
  
      end
      clear AziEl_GSi Range_GSi XYZgs2sc_i AGLgs2ant_i Rgs2tg_i Doppler_i DpRate_i Rgs2ant_i
      
 
      %============================================================
      % Ground Station Communication Link Occurences / Visibility
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.55,WB);
      tLink_GS = zeros(N_GS,10);
      Pass_GS = zeros(N_GS,10);
      Link_GS = zeros(N_GS,10);
      Link_state = zeros(N_GS,sim_sample);
      LinkT_GS = zeros(N_GS,1);
      PassT_GS = zeros(N_GS,1);

      Beam = 120*pi/180;
      
      for k=1:N_GS
         [tLink_GSi, Pass_GSi, Link_state(k,:), Link_GSi] = ...
            getLink(Elevation_GS(k,:), MEL_GS(k), Ts, AGLgs2ant(k,:), Beam);
         
         tLink_GS(k,:) = tLink_GSi(1:10);
         Pass_GS(k,:) = Pass_GSi(1:10);
         Link_GS(k,:) = Link_GSi(1:10);         
         
         for j=1:length(Link_GSi)
            LinkT_GS(k) = LinkT_GS(k) + Link_GSi(j);
         end
         
         for j=1:length(Pass_GSi)
            PassT_GS(k) = PassT_GS(k) + Pass_GSi(j);
         end         
      end
      clear tLink_GSi Link_GSi Pass_GSi      
      
      %============================================================
      % Eclipses Occurences
      %============================================================
      [beta_ecl theta_ecl] = getShadow(LLray(:,1), LLorbpole(:,1), a);

      Eclipse_approx = theta_ecl/(2*pi) * Psc;				%[s]
      
      [Eclipse_state, VFsun, Eclipse_time, Eclipse] = getEclipse(XYZsun_geo, XYZsc_geo, Ts);
      EclipseT = 0;
      for k=1:length(Eclipse_approx)
         EclipseT = EclipseT + Eclipse(k);  %sum of all eclipse times
      end
      

      
      %============================================================
      % Faces Exposition to Radiations (Sun, Earth, Albedo)
      %============================================================
      %--------------------------------------------
      % --- Cube faces exposition to solar ray ---
      %--------------------------------------------
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.60,WB);
      Obody_sun = zeros(6,sim_sample);
      for j=1:6
         Obody_sun(j,:) = acos(dot(Nbody_geo(3*j-2:3*j,:), IJKsun_geo));
         for k=1:sim_sample
            Obody_sun(j,k) = min(Obody_sun(j,k), pi/2);		% Back side of face is not exposed (angle>90°)
         end
      end
      Abody_sun = ([A2sc;A2sc;A2sc;A2sc;A1sc;A1sc]*ones(1,sim_sample)).*cos(Obody_sun);
      Atot_sun = Abody_sun(1,:) + Abody_sun(2,:) + Abody_sun(3,:) + Abody_sun(4,:) + Abody_sun(5,:) + Abody_sun(6,:);
      
      % --- Max (sqrt(3)*Aface)---
      Obody_sun_max = zeros(6,1);
      for j=1:6
         Obody_sun_max(j,:) = acos(dot(Nbody(3*j-2:3*j,1), -([1; 1; 1]./sqrt(3)) ));
         Obody_sun_max(j,:) = min(Obody_sun_max(j,:), pi/2);		% Back side of face is not exposed (angle>90°)
      end
      Abody_sun_max = [A2sc;A2sc;A2sc;A2sc;A1sc;A1sc].*cos(Obody_sun_max);
      Atot_sun_max = Abody_sun_max(1,:)+ Abody_sun_max(2,:) + Abody_sun_max(3,:) + Abody_sun_max(4,:) + Abody_sun_max(5,:) + Abody_sun_max(6,:);
      
      % --- Min (1*Aface)---
      Obody_sun_min = zeros(6,1);
      for j=1:6
         Obody_sun_min(j,:) = acos( dot( Nbody(3*j-2:3*j,1), -([0; 0; 1]) ));
         Obody_sun_min(j,:) = min(Obody_sun_min(j,:), pi/2);		% Back side of face is not exposed (angle>90°)
      end
      Abody_sun_min = [A2sc;A2sc;A2sc;A2sc;A1sc;A1sc].*cos(Obody_sun_min);
      Atot_sun_min = Abody_sun_min(1,:)+ Abody_sun_min(2,:) + Abody_sun_min(3,:) + Abody_sun_min(4,:) + Abody_sun_min(5,:) + Abody_sun_min(6,:);
      
      
      %--------------------------------------------
      % --- Cube faces exposition to Earth ---
      %--------------------------------------------
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.65,WB);
      rho_earth = asin(Re./Rsc); 				%[rad] radius angle
      Obody_geo = zeros(6,sim_sample);
      Obody_geo_F = Obody_geo;
      for j=1:6
         Obody_geo(j,:) = acos(dot(Nbody_geo(3*j-2:3*j,:), -IJKsc_geo));
         for k=1:sim_sample
            Obody_geo_F(j,k) = Obody_geo(j,k);
            %Obody_geo_F(j,k) = min(Obody_geo_F(j,k), pi/2+rho_earth(k));
            Obody_geo(j,k) = max(Obody_geo(j,k)-rho_earth(k), 0);   	% Front side of face is fully exposed (angle<0°)
            Obody_geo(j,k) = min(Obody_geo(j,k)-rho_earth(k), pi/2);	% Back side of face is not exposed (angle>90°)
         end
      end
      Abody_geo = ([A2sc;A2sc;A2sc;A2sc;A1sc;A1sc]*ones(1,sim_sample)) .* cos(Obody_geo);
      Atot_geo = Abody_geo(1,:) + Abody_geo(2,:) + Abody_geo(3,:) + Abody_geo(4,:) + Abody_geo(5,:) + Abody_geo(6,:);
      %VFearth = ([1;1;1;1;1;1]*(sin(rho_earth)).^2) .* ...
      %   (cos((pi/2-([1;1;1;1;1;1]*rho_earth)-Obody_geo_F)*pi/4./([1;1;1;1;1;1]*rho_earth))).^2;%[] View Factor (plate)
      %view factor of a finite sphere from a perpendicular infinitesimal plane 
      VFeartFIG = (1-[1;1;1;1;1;1]*cos(rho_earth));
      VFearth = VFeartFIG .* 0.5.*(1+ cos(max(min(pi/2./([1;1;1;1;1;1]*rho_earth).*(Obody_geo_F+([1;1;1;1;1;1]*rho_earth)-pi/2),pi),0)));
      %correction factor for not complete sphere in field of view
 
      % --- Max ---
      Obody_geo_max = zeros(6,sim_sample);
      for j=1:6
         Obody_geo_max(j,:) = acos(dot(Nbody(3*j-2:3*j,1), -([1; 1; 1]./sqrt(3)) ));
         for k=1:sim_sample   
            Obody_geo_max(j,k) = acos(dot(Nbody(3*j-2:3*j,1), -([1; 1; 1]./sqrt(3)) ));
            Obody_geo_max(j,k) = max(Obody_geo_max(j,k), 0);   		% Front side of face is fully exposed (angle<0°)
            Obody_geo_max(j,k) = min(Obody_geo_max(j,k), pi/2);		% Back side of face is not exposed (angle>90°)
         end
      end
      Abody_geo_max = ([A2sc;A2sc;A2sc;A2sc;A1sc;A1sc]*ones(1,sim_sample)).*cos(Obody_geo_max);
      Atot_geo_max = Abody_geo_max(1,:)+ Abody_geo_max(2,:) + Abody_geo_max(3,:) + Abody_geo_max(4,:) + Abody_geo_max(5,:) + Abody_geo_max(6,:);
      
      % --- Min ---
      Obody_geo_min = zeros(6,sim_sample);
      for j=1:6
         for k=1:sim_sample
            Obody_geo_min(j,k) = acos( dot( Nbody(3*j-2:3*j,1), -([0; 0; 1]) ));
            Obody_geo_min(j,k) = max(Obody_geo_min(j,k), 0);   		% Front side of face is fully exposed (angle<0°)
            Obody_geo_min(j,k) = min(Obody_geo_min(j,k), pi/2);		% Back side of face is not exposed (angle>90°)
         end
      end
      Abody_geo_min = ([A2sc;A2sc;A2sc;A2sc;A1sc;A1sc]*ones(1,sim_sample)).*cos(Obody_geo_min);
      Atot_geo_min = Abody_geo_min(1,:)+ Abody_geo_min(2,:) + Abody_geo_min(3,:) + Abody_geo_min(4,:) + Abody_geo_min(5,:) + Abody_geo_min(6,:);
      
      Obody_geo_min2 = zeros(6,sim_sample);
      for j=1:6
         for k=1:sim_sample
            Obody_geo_min2(j,k) = acos( dot( Nbody(3*j-2:3*j,1), -([cos(rho_earth(k)); sin(rho_earth(k)); 0]) ));
            Obody_geo_min2(j,k) = max(Obody_geo_min2(j,k), 0);   		% Front side of face is fully exposed (angle<0°)
            Obody_geo_min2(j,k) = min(Obody_geo_min2(j,k), pi/2);		% Back side of face is not exposed (angle>90°)
         end
      end
      Abody_geo_min2 = ([A2sc;A2sc;A2sc;A2sc;A1sc;A1sc]*ones(1,sim_sample)).*cos(Obody_geo_min2);
      Atot_geo_min2 = Abody_geo_min2(1,:)+ Abody_geo_min2(2,:) + Abody_geo_min2(3,:) + Abody_geo_min2(4,:) + Abody_geo_min2(5,:) + Abody_geo_min2(6,:);
      
      
      
      %--------------------------------------------
      % --- Cube faces exposition to Albedo ---
      %--------------------------------------------
      fprintf('*')
      LocalDisplay('*', 'add')
            
      theta = acos(dot(XYZsc_geo./([1;1;1]*Rsc), XYZsun_geo./([1;1;1]*Rsun) ));
      VFalbedo = VFearth .* 0.5.*(1+ cos(max(min(pi/2./([1;1;1;1;1;1]*(pi/2-rho_earth)).*([1;1;1;1;1;1]*(theta-rho_earth)),pi),0)));
     		%correction factor for position of S/C and field of view versus albedo hemisphere  
      
      Kalbedo = 0.664 + 0.521.*rho_earth - 0.203.*rho_earth.^2; 
      
      
      %--------------------------------------------
      % --- Cube faces exposition to Deep Space ---
      %--------------------------------------------
      fprintf('*')
      LocalDisplay('*', 'add')
      Abody_space = (1-VFearth).*([A2sc;A2sc;A2sc;A2sc;A1sc;A1sc]*ones(1,sim_sample));
      Atot_space = Abody_space(1,:) + Abody_space(2,:) + Abody_space(3,:) + Abody_space(4,:) + Abody_space(5,:) + Abody_space(6,:);
      
      
      %============================================================
      % Energy Input on Faces (Sun, Earth, Albedo)
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      %global Gsolar

      %Gsolar = (Gsolar_max + Gsolar_min)/2 + (Gsolar_max - Gsolar_min)/2.*sin(usun);

      
      emi_wall = PVpacking * emi_cell + (1-PVpacking) * emi_alu;
      abs_wall = PVpacking * abs_cell + (1-PVpacking) * abs_alu;
      
      % --- Earth IR ---
      q_earth = GIR * (emi_wall*ones(1,sim_sample)) .* VFearth;
      Qearth = q_earth .* Abody_geo;				%[W]
      
      % --- Sun ---
      q_sun = (abs_wall * (SolarFlux .* VFsun)); 			%[W/m²]
      Qsun = q_sun .* Abody_sun;				%[W]		
      
      % --- Sun Albedo ---
      Ka = 0.664 + 0.521 .* rho_earth - 0.203 .* rho_earth.^2;	%[]
      q_albedo =  (abs_wall * (Albedo * SolarFlux .* Ka)) .* VFalbedo;	%[W/m²]
      Qalbedo = q_albedo .* Abody_geo;%Abody_albedo;			%[W]	
      
      % --- Deep Space ---
      q_radspace = (Boltzmann .* emi_wall ) * ( (300*ones(1,sim_sample)).^4 - 4^4 );
      Qspace = q_radspace .* Abody_space;				%[W]
      Qspacesum = Qspace(1,:) + Qspace(2,:) + Qspace(3,:) + Qspace(4,:) + Qspace(5,:) + Qspace(6,:);
      
      % --- Summations ---
      Qface = Qearth + Qsun + Qalbedo;
      Qearthsum = Qearth(1,:) + Qearth(2,:) + Qearth(3,:) + Qearth(4,:) + Qearth(5,:) + Qearth(6,:);
      Qsunsum = Qsun(1,:) + Qsun(2,:) + Qsun(3,:) + Qsun(4,:) + Qsun(5,:) + Qsun(6,:);
      Qalbedosum = Qalbedo(1,:) + Qalbedo(2,:) + Qalbedo(3,:) + Qalbedo(4,:) + Qalbedo(5,:) + Qalbedo(6,:);
      Qtotal = Qearthsum + Qsunsum + Qalbedosum - Qspacesum;
      
      
      %============================================================
      % Solar Power Generation (Sun, Earth, Albedo)
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.70,WB);
      
      Pgen_geo = (PVpacking' * Qearth) * PVeff;
      Pgen_sun = (PVpacking' * Qsun) * PVeff;
      Pgen_albedo = (PVpacking' * Qalbedo) * PVeff;
      Pgen = Pgen_geo + Pgen_sun + Pgen_albedo;
      
      
      %============================================================
      % Thermal Analysis, Worst Cases, Single Body (SMAD p447)
      %============================================================
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.70,WB);
      
      absorb = 0.805;  		%[] (solar cells)
      emi = 0.825;   		%[] (solar cells)
      
      % --- Heat generation ---
      Qgenmax = 1;  		%W
      Qgenmin = 0.1;   		%W
      
      % --- Earth IR ---
      rho_eartFIG = asin(Re/(Re+h));   	%[rad]
      FeartFIG = (1-cos(rho_eartFIG))/2;	%[]
      q_eartFIG = FeartFIG*GIR*emi;  	%[W/m²]
      
      % --- Sun ---
      q_sun0 = Gsolar_max*absorb;        	%[W/m²]
      
      % --- Sun Albedo ---
      Ka0 = 0.664+0.521*rho_eartFIG-0.203*rho_eartFIG^2;   	%[]
      q_albedo0 = FeartFIG*Albedo*Gsolar_max*Ka0*absorb;	%[W/m²]
      
      % --- Cube in space with ideal conductivity ---
      Qearth1 = 3*A1sc*q_eartFIG;     		%[W]
      Qsun1 = 1.8*A1sc*q_sun0;    		%[W]
      Qalbedo1 = 3*A1sc*q_albedo0; 		%[W]
      Tmax1 = ((Qgenmax/6 + Qearth1 + Qalbedo1 + Qsun1)/(6*A1sc*Boltzmann*emi))^(0.25) - 273.15;	%[degC]
      Tmin1 = ((Qgenmin/6 + Qearth1)/(6*A1sc*Boltzmann*emi))^(0.25) - 273.15;			%[degC]
      
      % --- Single plane surface with no conduction, adiabatic backplane ---
      Qearth2 = A1sc*q_eartFIG;  			%[W]
      Qsun2 = A1sc*q_sun0;      			%[W]
      Qalbedo2 = A1sc*q_albedo0;  		%[W]
      Tmax2 = ((Qgenmax + Qearth2 + Qalbedo2 + Qsun2)/(A1sc*Boltzmann*emi))^(0.25) - 273.15; 	%[degC]
      Tmin2 = ((Qgenmin + Qearth2)/(A1sc*Boltzmann*emi))^(0.25) - 273.15;				%[degC]
      
      % --- Single plane surface with conduction, adiabatic backplane ---
      Qearth3 = A1sc*q_eartFIG;  			%[W]
      Qsun3 = A1sc*q_sun0;      			%[W]
      Qalbedo3 = A1sc*q_albedo0;  		%[W]
      Qrad3 = - 3*A1sc * Boltzmann * emi * 173^4;	%[W]
      Tmax3 = ((Qgenmax + Qearth3 + Qalbedo3 + Qsun3 + Qrad3)/(A1sc*Boltzmann*emi))^(0.25) - 273.15; %[degC]
      Tmin3 = ((Qgenmin + Qearth3 + Qrad3)/(A1sc*Boltzmann*emi))^(0.25) - 273.15;			%[degC]
      
      
      
      fprintf('*')
      LocalDisplay('*', 'add')
      waitbar(0.75,WB);
      %============================================================
      % Thermal Analysis, ODE45, n order Cube Model (SMAD p447)
      %============================================================
      T2 = [NaN NaN];
      T7 = [NaN NaN NaN NaN NaN NaN NaN];
      if get(findobj(FIG,'Tag','CheckBoxTempSim'),'Value')
         waitbar(0.75,WB,'Simulation Progress (Thermal simulation)');      
         
         T0_wall = 40;
         T0_in = 40;
         T0 = [T0_wall T0_wall T0_wall T0_wall T0_wall T0_wall T0_in] + 273;		%[K]
         %Qin = 0.2;				%[W]
       
         
         
         %--------------------------------------------
         % --- 2nd Order Model ---
         %--------------------------------------------      
         fprintf('*')
         LocalDisplay(' Thermal 2-DOF', 'add')
         waitbar(0.75,WB);
         global Qtot2
         Qtot2 = zeros(2,sim_sample);
         
         [t2,T2] = ode45('T_ode2', t(sim_sample), [T0_wall T0_in]+273, odeset('RelTol',1e-6,'AbsTol',[1 1]*1e-6), t, Qin, Rsc, ...
            Qearth, Qsun, Qalbedo, VFearth, PVpacking, wsc, A1sc, side, wallsc, msc);
         
         for k=1:2
            for j=1:length(Qtot2)  
               if Qtot2(k,j) ==0
                  Qtot2(k,j) = NaN;
               end
            end
         end
         
         % --- Deep Space ---
         q_radspace2 = zeros(6,sim_sample);
         for k=1:sim_sample
            j=find(t2>=t(k));
            T2i = T2(j(1),1);
            q_radspace2(:,k) = (Boltzmann .* emi_wall ) * (T2i.^4 - 4^4);
         end
         Qspace2 = q_radspace2 .* Abody_space;	%[W]
         Qspacesum2 = Qspace2(1,:) + Qspace2(2,:) + Qspace2(3,:) + Qspace2(4,:) + Qspace2(5,:) + Qspace2(6,:);
         Qtotal2 = Qearthsum + Qsunsum + Qalbedosum - Qspacesum2;
         
         if min(min(T2)) < 0 | max(max(T2)) > 1000
            LocalDisplay('WARNING! 2-DOF temperature simulation didn''t converge. Refine sampling','nl')
            LocalDisplay(' *', 'add')
         end
         
         
         %--------------------------------------------
         % --- 7th Order Model ---
         %--------------------------------------------      
         fprintf('*')
         LocalDisplay(' 7-DOF', 'add')
         waitbar(0.85,WB);
         global Qtot7
         global Qtot7b
         Qtot7 = zeros(7,sim_sample);
         Qtot7b = zeros(7,sim_sample);
         
         [t7,T7] = ode45('T_ode7', t(sim_sample), T0, odeset('RelTol',1e-6,'AbsTol',[1 1 1 1 1 1 1]*1e-6), t, Qin, Rsc, ...%Eclipse_state,...
            Qearth, Qsun, Qalbedo, VFearth, PVpacking, wsc, A1sc, side, wallsc, msc);
         
         for k=1:7
            for j=1:length(Qtot7)  
               if Qtot7(k,j) ==0
                  Qtot7(k,j) = NaN;

               end
               if Qtot7b(k,j) ==0
                  Qtot7b(k,j) = NaN;
               end
            end
         end
         
         % --- Deep Space ---
         q_radspace7 = zeros(6,sim_sample);
         for k=1:sim_sample
            j=find(t7>=t(k));
            T7i = T7(j(1),1);
            q_radspace7(:,k) = (Boltzmann .* emi_wall ) * (T7i.^4 - 4^4);
         end
         Qspace7 = q_radspace7 .* Abody_space;	%[W]
         Qspacesum7 = Qspace7(1,:) + Qspace7(2,:) + Qspace7(3,:) + Qspace7(4,:) + Qspace7(5,:) + Qspace7(6,:);
         Qtotal7 = Qearthsum + Qsunsum + Qalbedosum - Qspacesum7;   
         
         if min(min(T7)) < 0 | max(max(T7)) > 1000
            LocalDisplay('WARNING! 7-DOF temperature simulation didn''t converge. Refine sampling','nl')
            LocalDisplay(' *', 'add')
         end
         
      end
      
      
      %************************************************************
      % Data Log Display Window      
      %************************************************************
      LocalDisplay('*', 'add') 
      waitbar(0.95,WB,'Simulation Progress');
      
      DataStr0={'CUBESAT ORBIT SIMULATION RESULTS';
         '';
         sprintf('=== Orbit Elements =================================================');
         sprintf(' Altitude at Perigee (h)     	%0.5g [km]',h);
         sprintf(' Epoch Year                  	%0.5g [yr]',EYear+2000);
         sprintf(' Epoch Time                  	%0.5g [day]',Epoch);
         sprintf(' Orbit Inclination (i)       	%0.5g [deg]',incl*180/pi);
         sprintf(' RA of Node (raan)           	%0.5g [deg]',raan*180/pi);
         sprintf(' Orbit Eccentricity (e)      	%0.5g []',e);
         sprintf(' Argument of Perigee         	%0.5g [deg]',argp*180/pi);
         sprintf(' Mean Anomaly (M)            	%0.5g [deg]',Msc0*180/pi);
         sprintf(' Orbital Mean Motion (n)     	%0.5g [rev/day]',n);
         sprintf('====================================================================');
         '';
         sprintf('=== Orbital Characteristics ========================================');
         sprintf(' Semi-major axis (a)         	%0.5g [km]',a);
         sprintf(' Semi-minor axis (b)         	%0.5g [km]',b);
         sprintf(' Focal distance to center (c)	%0.5g [km]',c);
         sprintf(' Apogee (ra)                 	%0.5g [km]',ra);
         sprintf(' Perigee (rp)                	%0.5g [km]',rp);
         sprintf(' Velocity at apogee (Va)     	%0.5g [km/s]',Va);
         sprintf(' Velocity at perigee (Vp)    	%0.5g [km/s]',Vp);
         sprintf(' Orbital period (Psc)        	%0.5g [min]',Psc/60);
         sprintf(' Orbital mean motion (n)     	%0.5g [rev/day]',n);
         sprintf(' Node spacing                	%0.5g [deg/rev]',dLON*180/pi);
         sprintf(' Node precession (secular)   	%0.5g [deg/day]',Draan*Te*180/pi);
         sprintf(' Arg of perigee rate (secular)%0.5g [deg/day]',Dargp*Te*180/pi);
         sprintf(' Eclipse duration din        	%0.5g [min]',Eclipse_min/60);
         sprintf(' Eclipse duration dax        	%0.5g [min]',Eclipse_max/60);
         sprintf('====================================================================');
         '';      
         sprintf('=== Orbital Perturbations ==========================================');
         sprintf('  Spacecraft rotation	     	%0.5g [rad/s]',wsc);
         sprintf('  Earth''s rotation  	      	%0.5g [rad/s]',we);
         '';      
         sprintf(' Secular variations of RA of node');       
         sprintf('  Sun''s gravity     	      	%0.5g [rad/s]',Draan_sun);
         sprintf('  Moon''s gravity    	      	%0.5g [rad/s]',Draan_moon);
         sprintf('  Geopotential function J2 	%0.5g [rad/s]',Draan_geo);
         '';      
         sprintf(' Secular variations of Arg of perigee');       
         sprintf('  Sun''s gravity     	      	%0.5g [rad/s]',Dargp_sun);
         sprintf('  Moon''s gravity    	      	%0.5g [rad/s]',Dargp_moon);
         sprintf('  Geopotential function J2 	%0.5g [rad/s]',Dargp_geo);
         '';      
         sprintf(' Atmospheric drag');       
         sprintf('  Atmospheric density min    	%0.5g [kg/m^3]', rhoa_min);
         sprintf('  Atmospheric density max    	%0.5g [kg/m^3]', rhoa_max);
         sprintf('  Atmospheric drag min       	%0.5g [N]', Drag_min);
         sprintf('  Atmospheric drag max       	%0.5g [N]', Drag_max);
         sprintf('  Orbit decay rate min       	%0.5g [rev/day²]',orb_decay_min);
         sprintf('  Orbit decay rate max       	%0.5g [rev/day²]',orb_decay_max);
         sprintf('  Altitude decay rate min    	%0.5g [km/year]',Da_rev_min/Psc*Psun);
         sprintf('  Altitude decay rate max    	%0.5g [km/year]',Da_rev_max/Psc*Psun);
         sprintf('  Orbit lifetime min         	%0.5g [day]',Life_min);
         sprintf('  Orbit lifetime max         	%0.5g [day]',Life_max);
         sprintf('====================================================================');
         '';
      };
      DataStr1={
         sprintf('=== Spacecraft Parameters ==========================================');
         sprintf(' Mass                        	%0.5g [kg]',msc);         
         sprintf(' Moment of inertia tensor    	[%10.5f %10.5f %10.5f ] [kg.mm²]',Ixx*1e+6, Ixy*1e+6, Ixz*1e+6);
         sprintf('                              [%10.5f %10.5f %10.5f ]',Ixy*1e+6, Iyy*1e+6, Iyz*1e+6);
         sprintf('                              [%10.5f %10.5f %10.5f ]',Ixz*1e+6, Iyz*1e+6, Izz*1e+6);
         sprintf(' Base length                 	%0.5g [m]',side);
         sprintf(' Side length                 	%0.5g [m]',lsc);
         sprintf(' Base face area              	%0.5g [m²]',A1sc);
         sprintf(' Long face area              	%0.5g [m²]',A2sc);
         sprintf(' Total area                  	%0.5g [m²]',2*A1sc+4*A2sc);
         sprintf(' Wall thickness              	%0.5g [mm]',wallsc*1000);
         sprintf(' Solar cells efficiency      	%0.5g []',PVeff);
         sprintf(' Cells packing factor (base) 	%0.5g []',PVpack1);
         sprintf(' Cells packing factor (sides)	%0.5g []',PVpack2);
         sprintf(' Atmospheric drag coefficient	%0.5g []',Cd);
         sprintf(' Ballistic coefficient       	%0.5g [kg/m²]',Cb);
         sprintf('====================================================================');
         '';
         sprintf('=== Magnetic Stabilization Parameters ==============================');
         %sprintf(' Spacecraft Mass             	±%0.4g [kg]',msc);
         %sprintf(' Spacecraft Moment of Inertia	±%0.6g [kg.m²]',Ixx);
         sprintf(' Mass magnet Alnico-5        	±%0.4g [g]',mmag*1000);
         sprintf(' Mass hysteresis Hymu-80     	±%0.4g [g]',(mhyst1+mhyst2)*1000);
         sprintf(' System pointing resolution  	±%0.4g [°]',Stab_res*180/pi);
         sprintf('====================================================================');
         '';
         sprintf('=== Downlink Characteristics =======================================');
         sprintf(' Maximum possible downlink duration 	%0.5g [min]',DLmax/60);
         sprintf(' Nb. downlinks per day (estimate)  	%0.5g [day-1]',DLday);
         sprintf(' Downlink Slant range perigee      	%0.6g [km]',S_maxp);
         sprintf(' Downlink Zenith range perigee     	%0.6g [km]',S_minp);
         sprintf(' Downlink Slant range apogee       	%0.6g [km]',S_maxa);
         sprintf(' Downlink Zenith range apogee      	%0.6g [km]',S_mina);
         sprintf('====================================================================');
      };
      DataStr2={        
         '';



         sprintf('=== Communication Link Simulation Results for %4.1f Orbits ==========',sim_dur);
         sprintf([' ',NAME_GS(1,:)]);

         sprintf('  Time [min]     	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',tLink_GS(1,1:10)/60);
         sprintf('  Pass [±%0.1gmin]	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Ts/60/2,Pass_GS(1,1:10)/60);
         sprintf('  Link [±%0.1gmin]	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Ts/60/2,Link_GS(1,1:10)/60);
         sprintf('  Compiled Passes [min]	%1.1f ±%1.1f over %0.4g  (%0.4g%%)',...
            PassT_GS(1)/60, Ts/60/2*length(tLink_GS(1,:)), sim_last/60, PassT_GS(1)/sim_last*100);
         sprintf('  Compiled Links  [min]	%1.1f ±%1.1f over %0.4g  (%0.4g%%)',...
            LinkT_GS(1)/60, Ts/60/2*length(tLink_GS(1,:)), sim_last/60, LinkT_GS(1)/sim_last*100);
         sprintf([' ',NAME_GS(2,:)]);
         sprintf('  Time [min]     	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',tLink_GS(2,1:10)/60);
         sprintf('  Pass [±%0.1gmin]	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Ts/60/2,Pass_GS(2,1:10)/60);
         sprintf('  Link [±%0.1gmin]	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Ts/60/2,Link_GS(2,1:10)/60);
         sprintf('  Compiled Passes [min]	%1.1f ±%1.1f over %0.4g  (%0.4g%%)',...
            PassT_GS(2)/60, Ts/60/2*length(tLink_GS(2,:)), sim_last/60, PassT_GS(2)/sim_last*100);
         sprintf('  Compiled Links  [min]	%1.1f ±%1.1f over %0.4g  (%0.4g%%)',...
            LinkT_GS(2)/60, Ts/60/2*length(tLink_GS(2,:)), sim_last/60, LinkT_GS(2)/sim_last*100);
         sprintf([' ',NAME_GS(3,:)]);
         sprintf('  Time [min]     	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',tLink_GS(3,1:10)/60);
         sprintf('  Pass [±%0.1gmin]	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Ts/60/2,Pass_GS(3,1:10)/60);
         sprintf('  Link [±%0.1gmin]	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Ts/60/2,Link_GS(3,1:10)/60);
         sprintf('  Compiled Passes [min]	%1.1f ±%1.1f over %0.4g  (%0.4g%%)',...
            PassT_GS(3)/60, Ts/60/2*length(tLink_GS(3,:)), sim_last/60, PassT_GS(3)/sim_last*100);
         sprintf('  Compiled Links  [min]	%1.1f ±%1.1f over %0.4g  (%0.4g%%)',...
            LinkT_GS(3)/60, Ts/60/2*length(tLink_GS(3,:)), sim_last/60, LinkT_GS(3)/sim_last*100);
         sprintf('====================================================================');
         '';
         sprintf('=== Absolute Maximum Doppler Shifts  ===============================');
         sprintf([' ',NAME_GS(1,:), '		f = %0.4g [MHz], Doppler = %0.4g [kHz], Rate = %0.4g [kHz/s]'],...
            f0(1)/1000000,...
            max(abs(min(Doppler_GS(1,:))),max(Doppler_GS(1,:)))/1000,...
            max(abs(min(DpRate_GS(1,:))),max(DpRate_GS(1,:)))/1000);
         sprintf([' ',NAME_GS(2,:), '		f = %0.4g [MHz], Doppler = %0.4g [kHz], Rate = %0.4g [kHz/s]'],...
            f0(2)/1000000,...
            max(abs(min(Doppler_GS(2,:))),max(Doppler_GS(2,:)))/1000,...
            max(abs(min(DpRate_GS(2,:))),max(DpRate_GS(2,:)))/1000) ;        
         sprintf([' ',NAME_GS(3,:), '		f = %0.4g [MHz], Doppler = %0.4g [kHz], Rate = %0.4g [kHz/s]'],...
            f0(3)/1000000,...
            max(abs(min(Doppler_GS(3,:))),max(Doppler_GS(3,:)))/1000,...
            max(abs(min(DpRate_GS(3,:))),max(DpRate_GS(3,:)))/1000);
         sprintf('====================================================================');
         '';
         sprintf('=== Eclipse Characteristics (circular orbit approx Rsc=a ) =========');
         sprintf(' Solar incidence on orbit plane     	%0.5g°',beta_ecl*180/pi);
         sprintf(' Orbital portion in eclipse         	%0.4g°  (%0.4g%%)',theta_ecl*180/pi, theta_ecl/2/pi*100);
         sprintf(' Eclipse duration                   	%0.4g over %0.4g [min]',Eclipse_approx/60, Psc/60);
         sprintf('====================================================================');
         '';
         sprintf('=== Eclipse Simulation Results for %4.1f Orbits Simulation ==========',sim_dur);
         sprintf(' Time [min]    	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Eclipse_time(1:10)/60);
         sprintf(' Dur [±%0.1gmin]	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f,	%1.1f, ...',Ts/60/2,Eclipse(1:10)/60);
         sprintf(' Compiled [min]	%1.1f ±%1.1f over %0.4g  (%0.4g%%)',...
            EclipseT/60, Ts/60/2*length(Eclipse_time), sim_last/60, EclipseT/sim_last*100);
         sprintf('====================================================================');
         '';
         sprintf('=== Thermal Analysis (simple analytic models) ======================');
         sprintf(' Cube, perfect conduction:   	Tmax = %0.4g [°C]   Tmin = %0.4g [°C]',Tmax1, Tmin1);
         sprintf(' Plane face, no conduction:  	Tmax = %0.4g [°C]   Tmin = %0.4g [°C]',Tmax2, Tmin2);
         sprintf(' Plane face, some conduction:	Tmax = %0.4g [°C]   Tmin = %0.4g [°C]',Tmax3, Tmin3);
         sprintf('====================================================================');
         '';         
         sprintf('=== Thermal Simulation Results for %4.1f Orbits  ====================',sim_dur);
         sprintf(' External temperature	2-DOF model	Tmax = %0.4g [°C] Tmin = %0.4g [°C]',max(T2(:,1))-273, min(T2(:,1))-273);
         sprintf(' Internal temperature	2-DOF model	Tmax = %0.4g [°C] Tmin = %0.4g [°C]',max(T2(:,2))-273, min(T2(:,2))-273);
         sprintf(' External temperature	7-DOF model	Tmax = %0.4g [°C] Tmin = %0.4g [°C]',max(max(T7(:,1:6)))-273, min(min(T7(:,1:6)))-273);
         sprintf(' Internal temperature	7-DOF model	Tmax = %0.4g [°C] Tmin = %0.4g [°C]',max(T7(:,7))-273, min(T7(:,7))-273);
         sprintf('====================================================================');
         '';         
         sprintf('=== Power Generation Simulation Results for %4.1f Orbits ============',sim_dur);
         sprintf(' Average PV power generation from Sun	       	%0.4g [W]',mean(Pgen_sun));
         sprintf(' Average PV power generation from all sources	%0.4g [W]',mean(Pgen));
         sprintf('====================================================================');
         '';         
      };
      DataStr = [DataStr0;DataStr1;DataStr2];
      set(findobj(FIG,'Tag','ButtonHistory'),'UserData',DataStr);
      
      %************************************************************
      % Graphical Results
      %************************************************************      
      fprintf('*')
      LocalDisplay('*|  completed.', 'add')
      waitbar(1,WB);
      
      LocalDisplay('Graph progress |', 'nl')
      LocalDisplay('*', 'add')
      close(WB);
      WB = WAITBAR(0.01,'Graphs Generation Progress');
      %============================================================
      % Text Strings Definition (Graphs & Display)
      %============================================================
      str_alt = sprintf('h = %0.4gkm',h);
      str_incl = sprintf('i = %0.4g°',incl*180/pi);
      str_raan = sprintf('raan = %0.4g°',raan*180/pi);
      str_e = sprintf('e = %0.4g',e);
      str_argp = sprintf('argp = %0.4g°',argp*180/pi);
      str_Msc0 = sprintf('Mo = %0.4g°',Msc0*180/pi);
      str_time = sprintf('Epoch = %2.0f%6.2f',EYear,Epoch);
      str_n = sprintf('n = %0.4g[rev/day]',n);
      str_a = sprintf('a = %0.5gkm',a);
      str_b = sprintf('b = %0.5gkm',b);
      str_c = sprintf('c = %0.5gkm',c);
      str_rp = sprintf('rp = %0.5gkm',rp);
      str_ra = sprintf('ra = %0.5gkm',ra);
      %str_w = sprintf('w=%0.4g°/s',wo*180/pi);
      %str_theta = sprintf('theta=%0.4g°',thetao*180/pi);
      str_mag = sprintf('%0.4gg Alnico-5',mmag*1000);
      str_hyst = sprintf('%0.4gg Hymu-80',(mhyst1+mhyst2)*1000);      
      
      if get(findobj(FIG,'Tag','CheckBoxPlot2D'),'Value')
         
         %============================================================
         % Orbit Shape Plot
         %============================================================
         fprintf('*')
         LocalDisplay('*', 'add')
         figure(FIG);
         delete(findobj(FIG,'UserData', 'Legend', 'Type', 'axes'))
%          waitbar(0.03,WB); figure(FIG);
         
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotOrbShape', 'Type', 'axes', 'UserData', 'Graph'))
%          subplot(findobj(FIG,'Tag','PlotOrbShape', 'Type', 'axes'))
         axis([-2*ra 2*rp -1.2*b 1.2*b]);
         
         x_ellipse = linspace(-a, a, 100);
         y_ellipse = b.*sqrt(1-x_ellipse.^2/a^2);
         x_ellipse2 = [x_ellipse, x_ellipse]-c;
         y_ellipse2 = [y_ellipse, -y_ellipse];
         
         axis([-1.5*ra 1.5*rp -1.2*b 1.2*b]);
         H = plot(Re.*cos(circle), Re.*sin(circle), 'g-', 'linewidth', 2, 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = plot(x_ellipse2, y_ellipse2,'r-', 'linewidth', 2, 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = plot([-c -c], [-b b], 'b--', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = plot([-ra rp], [0 0], 'b--', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = plot(0,0,'r*',-2*c,0,'g*', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = plot(-ra,0,'r*',rp,0,'r*', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = text(-2*c/3, b/12, str_a, 'color', [0 0 1], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = text(-c+300, 3*b/4, str_b, 'color', [0 0 1], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = text(-2*c/3, -b/12, str_c, 'color', [0 0 1], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = text(rp+300, b/12, str_rp, 'color', [1 0 0], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         H = text(-ra+300, b/12, str_ra, 'color', [1 0 0], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
         %delete(findobj(FIG,'UserData', 'Legend', 'Type', 'axes'))
         H = legend('Earth & Focii', ['Orbit (',str_e,')']);
         set(H, 'Tag', 'PlotOrbShape', 'UserData', 'Legend');                            
         set(get(H, 'Children'), 'Tag', 'PlotOrbShape', 'UserData', 'Legend');       
         
         fprintf('*')
         LocalDisplay('*', 'add')
%          waitbar(0.05,WB); figure(FIG);
         %============================================================
         % Ground Station & Spacecraft Tracking Plot
         %============================================================

         fprintf('*')
         LocalDisplay('*', 'add')
         % === Plot angular tracking data ===
         LocalDisplay('*', 'add')
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotTrkAngle', 'Type', 'axes'))
%          waitbar(0.05,WB); figure(FIG);
         for k=1:N_GS
            H = plotTrace(180/pi*LONsc, 180/pi*Azimut_GS(k,:),'--');
            set(H, 'Tag', 'PlotTrkAngle', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
            H = plotTrace(180/pi*LONsc, 180/pi*Elevation_GS(k,:),'.');
            set(H, 'Tag', 'PlotTrkAngle', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         %   H = plotTrace(180/pi*LONsc, 180/pi*AGLgs2ant(k,:),'-');
         %   set(H, 'Tag', 'PlotTrkAngle', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         end
         plot(circle*180/pi, zeros(size(circle)), 'k+','Tag', 'PlotTrkAngle', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotTrkAngle', 'Type', 'axes'));
         H = legend(...
            ['G/S Azimut ', NAME_GS(1,:)],['G/S Elevation ', NAME_GS(1,:)],...
            ['G/S Azimut ', NAME_GS(2,:)],['G/S Elevation ', NAME_GS(2,:)],...
            ['G/S Azimut ', NAME_GS(3,:)],['G/S Elevation ', NAME_GS(3,:)]);
         %   ['G/S Azimut ', NAME_GS(1,:)],['G/S Elevation ', NAME_GS(1,:)],['Mag pointing deviation ', NAME_GS(1,:)],...
         %   ['G/S Azimut ', NAME_GS(2,:)],['G/S Elevation ', NAME_GS(2,:)],['Mag pointing deviation ', NAME_GS(2,:)],...
         %   ['G/S Azimut ', NAME_GS(3,:)],['G/S Elevation ', NAME_GS(3,:)],['Mag pointing deviation ', NAME_GS(3,:)]);
         set(H, 'Tag', 'PlotTrkAngle', 'UserData', 'Legend');
         set(get(H, 'Children'), 'Tag', 'PlotTrkAngle', 'UserData', 'Legend');                            
         
         
         % === Plot tracking ranges ===          
         LocalDisplay('*', 'add')
%          waitbar(0.10,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotTrkRange', 'Type', 'axes'))
         plot(circle*180/pi, S_maxa*ones(size(circle)), 'k+', circle*180/pi, S_minp*ones(size(circle)), 'k--',...
            'Tag', 'PlotTrkRange', 'UserData', 'Graph');
         for k=1:N_GS
            H = plotTrace(180/pi*LONsc, Range_GS(k,:),'.-');
            set(H, 'Tag', 'PlotTrkRange', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         %   H = plotTrace(180/pi*LONsc, Rgs2ant(k,:),'--');
         %   set(H, 'Tag', 'PlotTrkRange', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         end
         str_range1 = sprintf('Slant range = %0.5gkm',S_maxa);
         str_range2 = sprintf('Azimut range = %0.5gkm',S_minp);
         %delete(findobj(FIG,'Tag','PlotTrkRange', 'Type', 'axes'));
         H = legend(str_range1, str_range2,...
             ['G/S link range ', NAME_GS(1,:)], ['G/S link range ', NAME_GS(2,:)],['G/S link range ', NAME_GS(3,:)]);
         %   ['G/S link range ', NAME_GS(1,:)], ['Mag pointing gnd deviation ', NAME_GS(1,:)],...
         %   ['G/S link range ', NAME_GS(2,:)], ['Mag pointing gnd deviation ', NAME_GS(2,:)],...
         %   ['G/S link range ', NAME_GS(3,:)], ['Mag pointing gnd deviation ', NAME_GS(3,:)]);
         set(H, 'Tag', 'PlotTrkRange', 'UserData', 'Legend');                            
         set(get(H, 'Children'), 'Tag', 'PlotTrkRange', 'UserData', 'Legend');                            
         
         
         % === Plot doppler shift data ===
         LocalDisplay('*', 'add')
%          waitbar(0.15,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotTrkDoppler', 'Type', 'axes'))
         %plot(circle*180/pi, S_maxp*ones(size(circle)), 'k+', circle*180/pi, S_minp*ones(size(circle)), 'k--')
         for k=1:N_GS
            H = plotTrace(180/pi*LONsc, Doppler_GS(k,:)./1000,'-');
            set(H, 'Tag','PlotTrkDoppler', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);

         end
         %delete(findobj(FIG,'Tag','PlotTrkDoppler', 'Type', 'axes'));
         H = legend(['Doppler shift ', NAME_GS(1,:)],['Doppler shift ', NAME_GS(2,:)],['Doppler shift ', NAME_GS(3,:)]);
         set(H, 'Tag', 'PlotTrkDoppler', 'UserData', 'Legend');    
         set(get(H, 'Children'), 'Tag', 'PlotTrkDoppler', 'UserData', 'Legend');                            
         
         
         % === Plot angular tracking data ===
         LocalDisplay('*', 'add')
%          waitbar(0.20,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotTrkAngleTime', 'Type', 'axes'))
         for k=1:N_GS
            H = plot(t./60, 180/pi*Azimut_GS(k,:),'--');
            set(H, 'Tag', 'PlotTrkAngleTime', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
            H = plot(t./60, 180/pi*Elevation_GS(k,:),'.');
            set(H, 'Tag', 'PlotTrkAngleTime', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         %   H = plot(t./60, 180/pi*AGLgs2ant(k,:),'-');
         %   set(H, 'Tag', 'PlotTrkAngleTime', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 1-(k-1)/max(1,(N_GS-1)) 1]);
         end
         plot([0 t(sim_sample)/60], [0 0], 'k+--','Tag', 'PlotTrkAngleTime', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotTrkAngleTime', 'Type', 'axes'));
         H = legend(...
            ['G/S Azimut ', NAME_GS(1,:)],['G/S Elevation ', NAME_GS(1,:)],...
            ['G/S Azimut ', NAME_GS(2,:)],['G/S Elevation ', NAME_GS(2,:)],...
            ['G/S Azimut ', NAME_GS(3,:)],['G/S Elevation ', NAME_GS(3,:)]);
         %   ['G/S Azimut ', NAME_GS(1,:)],['G/S Elevation ', NAME_GS(1,:)],['Mag pointing deviation ', NAME_GS(1,:)],...
         %   ['G/S Azimut ', NAME_GS(2,:)],['G/S Elevation ', NAME_GS(2,:)],['Mag pointing deviation ', NAME_GS(2,:)],...
         %   ['G/S Azimut ', NAME_GS(3,:)],['G/S Elevation ', NAME_GS(3,:)],['Mag pointing deviation ', NAME_GS(3,:)]);
         set(H, 'Tag', 'PlotTrkAngleTime', 'UserData', 'Legend');
         set(get(H, 'Children'), 'Tag', 'PlotTrkAngleTime', 'UserData', 'Legend');                            
         
         
         % === Plot tracking ranges ===          
         LocalDisplay('*', 'add')
%          waitbar(0.25,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotTrkRangeTime', 'Type', 'axes'))
         plot([t(1) t(sim_sample)]/60, S_maxa*[1 1], 'k+--', [t(1) t(sim_sample)]/60, S_minp*[1 1], 'k.--',...
            'Tag', 'PlotTrkRangeTime', 'UserData', 'Graph');
         for k=1:N_GS
            H = plot(t./60, Range_GS(k,:),'.-');
            set(H, 'Tag', 'PlotTrkRangeTime', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         %   H = plot(t./60, Rgs2ant(k,:),'--');
         %   set(H, 'Tag', 'PlotTrkRangeTime', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 1-(k-1)/max(1,(N_GS-1)) 0.5]);
         end
         str_range1 = sprintf('Slant range = %0.5gkm',S_maxa);
         str_range2 = sprintf('Azimut range = %0.5gkm',S_minp);
         %delete(findobj(FIG,'Tag','PlotTrkRangeTime', 'Type', 'axes'));
         H = legend(str_range1, str_range2,...
             ['G/S link range ', NAME_GS(1,:)], ['G/S link range ', NAME_GS(2,:)],['G/S link range ', NAME_GS(3,:)]);
         %   ['G/S link range ', NAME_GS(1,:)], ['Mag pointing gnd deviation ', NAME_GS(1,:)],...
         %   ['G/S link range ', NAME_GS(2,:)], ['Mag pointing gnd deviation ', NAME_GS(2,:)],...
         %   ['G/S link range ', NAME_GS(3,:)], ['Mag pointing gnd deviation ', NAME_GS(3,:)]);
         set(H, 'Tag', 'PlotTrkRangeTime', 'UserData', 'Legend');                            
         set(get(H, 'Children'), 'Tag', 'PlotTrkRangeTime', 'UserData', 'Legend');                            
         
         
         % === Plot doppler shift data ===
         LocalDisplay('*', 'add')
%          waitbar(0.30,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotTrkDopplerTime', 'Type', 'axes'))
         %plot(circle*180/pi, S_maxp*ones(size(circle)), 'k+', circle*180/pi, S_minp*ones(size(circle)), 'k--')
         for k=1:N_GS
            H = plot(t./60, Doppler_GS(k,:)./1000,'-');
            set(H, 'Tag','PlotTrkDopplerTime', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         end
         %delete(findobj(FIG,'Tag','PlotTrkDopplerTime', 'Type', 'axes'));
         H = legend(['Doppler shift ', NAME_GS(1,:)],['Doppler shift ', NAME_GS(2,:)],['Doppler shift ', NAME_GS(3,:)]);
         set(H, 'Tag', 'PlotTrkDopplerTime', 'UserData', 'Legend');    
         set(get(H, 'Children'), 'Tag', 'PlotTrkDopplerTime', 'UserData', 'Legend'); 
         
         
         
         %============================================================
         % Radiation Plot
         %============================================================
         fprintf('*')
         LocalDisplay('*', 'add')
         % === Radiation Area ===
         LocalDisplay('*', 'add')
%          waitbar(0.35,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotRadArea', 'Type', 'axes'))
         H = plot(t./60, Atot_geo, 'b-', 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         H = plot(t./60, Atot_sun.*VFsun, 'r-', 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         %H = plot(t./60, Atot_albedo, 'g-', 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         H = plot(t./60, Atot_space, 'k--', 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         H = plot([0 t(sim_sample)]/60, [A1sc A1sc]*6, 'm-', 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         H = plot(t./60, VFsun*0.05, 'k-', 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         H = plot(t./60, Atot_geo_min, 'b--', t./60, Atot_geo_max, 'b--', t./60, Atot_geo_min2, 'b--',...
            'Tag', 'PlotRadArea', 'UserData', 'Graph');
         H = plot([0 t(sim_sample)]/60, [Atot_sun_min Atot_sun_min], 'r--', [0 t(sim_sample)]/60, [Atot_sun_max Atot_sun_max], ...
            'r--', 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         H = plot(0,0, 'Tag', 'PlotRadArea', 'UserData', 'Graph');
         %plot(t./60, Atot_albedo_min, t./60, 'g--',  Atot_albedo_max, 'g--');
         %delete(findobj(FIG,'Tag','PlotRadArea', 'Type', 'axes'));
         H = legend('Earth (disk)', 'Sun (point)', 'Space',...
            ['Total S/C area ', sprintf('%0.2gm²', 2*A1sc+4*A2sc)], 'Eclispe state (1/20)');      
         set(H, 'Tag', 'PlotRadArea', 'UserData', 'Legend');    
         set(get(H, 'Children'), 'Tag', 'PlotRadArea', 'UserData', 'Legend');     
         
         
         % === Radiation Power ===
         LocalDisplay('*', 'add')
%          waitbar(0.40,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotRadPower', 'Type', 'axes'))
         H = plot(t./60, Pgen_geo, 'b-', 'Tag', 'PlotRadPower', 'UserData', 'Graph');
         H = plot(t./60, Pgen_sun, 'r-', 'Tag', 'PlotRadPower', 'UserData', 'Graph');
         H = plot(t./60, Pgen_albedo, 'g-', 'Tag', 'PlotRadPower', 'UserData', 'Graph');
         H = plot(t./60, Pgen, 'm-', 'Tag', 'PlotRadPower', 'UserData', 'Graph');
         H = plot(t./60, VFsun, 'k-', 'Tag', 'PlotRadPower', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotRadPower', 'Type', 'axes'));
         H = legend('Earth IR', 'Sun', 'Sun Albedo', 'Total', 'Eclispe state');
         set(H, 'Tag', 'PlotRadPower', 'UserData', 'Legend');    
         set(get(H, 'Children'), 'Tag', 'PlotRadPower', 'UserData', 'Legend');       
         
         if get(findobj(FIG,'Tag','CheckBoxTempSim'),'Value')
            % === Radiation Energy, 2-DOF Model ===
            LocalDisplay('*', 'add')
%             waitbar(0.45,WB); figure(FIG);
            set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotRadInput2', 'Type', 'axes'))
            %for k=1:1
            %   H = plot(t./60, Qface(k,:), 'r-', 'Tag', 'PlotRadInput2, 'UserData', 'Graph');
            %   set(H, 'color', [(k-1)/max(1,(6-1)) 1-(k-1)/max(1,(6-1)) 0.5]);
            %end
            H = plot(t./60, Qearthsum, 'b-', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            H = plot(t./60, Qsunsum, 'r-', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            H = plot(t./60, Qalbedosum, 'g-', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            H = plot(t./60, -Qspacesum2, 'k-', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            H = plot(t./60, Qtotal2, 'm.-', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            H = plot(t./60, Qtot2(1,:), 'y-', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            %H = plot([0 t(sim_sample)/60], [1 1]*GIR*Asc*emi_wall(1),'b.-.', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            %H = plot([0 t(sim_sample)/60], [1 1]*Gsolar(1)*Asc*abs_wall(1),'r.-.', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            %H = plot([0 t(sim_sample)/60], [1 1]*Gsolar(1)*Albedo*Asc*abs_wall(1),'g.-.', 'Tag', 'PlotRadInput2', 'UserData', 'Graph');
            %H = plot(t./60, VFsun, 'k-', 'Tag', 'PlotRadInput2, 'UserData', 'Graph');
            %delete(findobj(FIG,'Tag','PlotRadInput2', 'Type', 'axes'));
            H = legend('Earth IR', 'Sun', 'Sun Albedo', 'Space', 'Total (computed)', 'Total (ODE sim)');%,'Nominal IR', 'Nominal Sun', 'Nominal Albedo')%, 'Eclipse State');
            set(H, 'Tag', 'PlotRadInput2', 'UserData', 'Legend');    
            set(get(H, 'Children'), 'Tag', 'PlotRadInput2', 'UserData', 'Legend');       
            
            
            % === Radiation Temperature, 2-DOF Model ===
            LocalDisplay('*', 'add')
%             waitbar(0.50,WB); figure(FIG);
            set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotRadTemp2', 'Type', 'axes'))
            H = plot(t2./60, T2(:,1)-273, 'm.-', 'Tag', 'PlotRadTemp2', 'UserData', 'Graph');
            H = plot(t2./60, T2(:,2)-273, 'g.-', 'Tag', 'PlotRadTemp2', 'UserData', 'Graph');
            H = plot(t./60, VFsun*100, 'k-', 'Tag', 'PlotRadTemp2', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [125 125],'r.-', 'Tag', 'PlotRadTemp2', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [-40 -40],'b.-', 'Tag', 'PlotRadTemp2', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [70 70],'r.-.', 'Tag', 'PlotRadTemp2', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [-5 -5],'b.-.', 'Tag', 'PlotRadTemp2', 'UserData', 'Graph');
            %delete(findobj(FIG,'Tag','PlotRadTemp2', 'Type', 'axes'));
            H = legend('Faces temperature', 'Internal temperature', 'Eclispe state (100x)');
            set(H, 'Tag', 'PlotRadTemp2', 'UserData', 'Legend');    
            set(get(H, 'Children'), 'Tag', 'PlotRadTemp2', 'UserData', 'Legend');       
            
            
            % === Radiation Energy, 7-DOF Model ===      
            LocalDisplay('*', 'add')
%             waitbar(0.55,WB); figure(FIG);
            set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotRadInput7', 'Type', 'axes'))
            %for k=1:6
            %   H = plot(t./60, Qface(k,:), 'r-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            %   set(H, 'color', [(k-1)/max(1,(6-1)) 1-(k-1)/max(1,(6-1)) 0.5]);
            %end
            H = plot(t./60, Qearthsum, 'b-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            H = plot(t./60, Qsunsum, 'r-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            H = plot(t./60, Qalbedosum, 'g-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            H = plot(t./60, -Qspacesum7, 'k-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            H = plot(t./60, Qtotal7, 'm.-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            H = plot(t./60, Qtot7(1,:), 'y-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            %H = plot(t./60, Qtot7b(1,:), 'b-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            %H = plot([0 t(sim_sample)/60], [1 1]*GIR*Asc*emi_wall(1),'b.-.', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            %H = plot([0 t(sim_sample)/60], [1 1]*Gsolar(1)*Asc*abs_wall(1),'r.-.', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            %H = plot([0 t(sim_sample)/60], [1 1]*Gsolar(1)*Albedo*Asc*abs_wall(1),'g.-.', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');

            %H = plot(t./60, VFsun, 'k-', 'Tag', 'PlotRadInput7', 'UserData', 'Graph');
            %H = legend('Face 1', 'Face 2', 'Face 3', 'Face 4', 'Face 5', 'Face 6', 'Total IR', 'Total Sun',...
            %   'Total Albedo', '-Total Space', 'Total', 'Nominal IR', 'Nominal Sun', 'Nominal Albedo')%, 'Eclipse State');
            %axis([0 120 -20 20])
            %delete(findobj(FIG,'Tag','PlotRadInput7', 'Type', 'axes'));
            H = legend('Earth IR', 'Sun', 'Sun Albedo', 'Space', 'Total (computed)', 'Total (ODE sim)');%, 'Eclipse State');
            set(H, 'Tag', 'PlotRadInput7', 'UserData', 'Legend');    
            set(get(H, 'Children'), 'Tag', 'PlotRadInput7', 'UserData', 'Legend');  
            
            
            % === Radiation Temperature, 7-DOF Model ===      
            LocalDisplay('*', 'add')
%             waitbar(0.60,WB); figure(FIG);
            set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotRadTemp7', 'Type', 'axes'))
            for k=1:6
               H = plot(t7./60, T7(:,k)-273, 'r-', 'Tag', 'PlotRadTemp7', 'UserData', 'Graph');
               set(H, 'color', [(k-1)/max(1,(6-1)) 1-(k-1)/max(1,(6-1)) 0.5]);
            end
            H = plot(t7./60, T7(:,7)-273, 'g.-', 'Tag', 'PlotRadTemp7', 'UserData', 'Graph');
            H = plot(t./60, VFsun*100, 'k-', 'Tag', 'PlotRadTemp7', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [125 125],'r.-', 'Tag', 'PlotRadTemp7', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [-40 -40],'b.-', 'Tag', 'PlotRadTemp7', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [70 70],'r.-.', 'Tag', 'PlotRadTemp7', 'UserData', 'Graph');
            H = plot([0 t(sim_sample)/60], [-5 -5],'b.-.', 'Tag', 'PlotRadTemp7', 'UserData', 'Graph');
            %delete(findobj(FIG,'Tag','PlotRadTemp7', 'Type', 'axes'));
            H = legend('Face 1', 'Face 2', 'Face 3', 'Face 4', 'Face 5', 'Face 6',...
               'Internal', 'Eclispe state (100x)');
            set(H, 'Tag', 'PlotRadTemp7', 'UserData', 'Legend');    
            set(get(H, 'Children'), 'Tag', 'PlotRadTemp7', 'UserData', 'Legend');   
            
         end
         
         %============================================================
         % Magnetic Stabilization
         %============================================================
         % === Angle to field ===
         LocalDisplay('*', 'add')
%          waitbar(0.65,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotMagSim', 'Type', 'axes'))
         %hold on;
    
         H = plotTrace(180/pi*LONsc,180/pi*alpha,'m'); set(H, 'Tag', 'PlotMagSim', 'UserData', 'Graph');
         H = plotTrace(180/pi*LONsc,min(180/pi*(alpha+Stab_res),180),'m--'); set(H, 'Tag', 'PlotMagSim', 'UserData', 'Graph');
         H = plotTrace(180/pi*LONsc,180/pi*beta,'g'); set(H, 'Tag', 'PlotMagSim', 'UserData', 'Graph');
         H = plotTrace(180/pi*LONsc,180/pi*phi,'b'); set(H, 'Tag', 'PlotMagSim', 'UserData', 'Graph');
         H = plotTrace(180/pi*LONsc,max(180/pi*(alpha-Stab_res),0),'m--'); set(H, 'Tag', 'PlotMagSim', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotMagSim', 'Type', 'axes'));
         H = legend('Field-Nadir angle', 'System resolution', 'Magnet-Nadir angle','Magnet-Field angle');
         set(H, 'Tag', 'PlotMagSim', 'UserData', 'Legend');    
         set(get(H, 'Children'), 'Tag', 'PlotMagSim', 'UserData', 'Legend');       
                    
         LocalDisplay('*', 'add')
%          waitbar(0.70,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotMagSimTime', 'Type', 'axes'))
         %hold on;
    
         H = plotTrace(t./60,180/pi*alpha,'m'); set(H, 'Tag', 'PlotMagSimTime', 'UserData', 'Graph');
         H = plotTrace(t./60,min(180/pi*(alpha+Stab_res),180),'m--'); set(H, 'Tag', 'PlotMagSimTime', 'UserData', 'Graph');
         H = plotTrace(t./60,180/pi*beta,'g'); set(H, 'Tag', 'PlotMagSimTime', 'UserData', 'Graph');
         H = plotTrace(t./60,180/pi*phi,'b'); set(H, 'Tag', 'PlotMagSimTime', 'UserData', 'Graph');
         H = plotTrace(t./60,max(180/pi*(alpha-Stab_res),0),'m--'); set(H, 'Tag', 'PlotMagSimTime', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotMagSimTime', 'Type', 'axes'));
         H = legend('Field-Nadir angle', 'System resolution', 'Magnet-Nadir angle','Magnet-Field angle');
         set(H, 'Tag', 'PlotMagSimTime', 'UserData', 'Legend');    
         set(get(H, 'Children'), 'Tag', 'PlotMagSimTime', 'UserData', 'Legend');       
         
         % === Angular rate ===
         LocalDisplay('*', 'add')
%         waitbar(0.75,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotRateTime', 'Type', 'axes'))
         %hold on;
    
         H = plotTrace(t./60,wb_len,'g.-'); set(H, 'Tag', 'PlotRateTime', 'UserData', 'Graph');
         H = plotTrace(t./60,wspin,'m'); set(H, 'Tag', 'PlotRateTime', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotRateTime', 'Type', 'axes'));
         H = legend('Total angular rate', 'Angular rate about magnet axis (nutation)');
         set(H, 'Tag', 'PlotRateTime', 'UserData', 'Legend');    
         set(get(H, 'Children'), 'Tag', 'PlotRateTime', 'UserData', 'Legend');
         
         % === Angle to Ground Station ===
         LocalDisplay('*', 'add')
%          waitbar(0.80,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotAntPoint', 'Type', 'axes'))
         for k=1:N_GS
            H = plot(180/pi*LONsc, 180/pi*AGLgs2ant(k,:),'-');
            set(H, 'Tag', 'PlotAntPoint', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         end
         plot(circle*180/pi, Beam/2*90/pi*ones(size(circle)), 'k+--','Tag', 'PlotAntPoint', 'UserData', 'Graph');
         plot(circle*180/pi, Beam*90/pi*ones(size(circle)), 'r+--','Tag', 'PlotAntPoint', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotAntPoint', 'Type', 'axes'));
         H = legend(['Mag pointing angle to ', NAME_GS(1,:)],['Mag pointing angle to ', NAME_GS(2,:)],...
             ['Mag pointing angle to ', NAME_GS(3,:)], '-3dB Beam width', '-6dB Beam width');
         set(H, 'Tag', 'PlotAntPoint', 'UserData', 'Legend');
         set(get(H, 'Children'), 'Tag', 'PlotAntPoint', 'UserData', 'Legend');          
         
         LocalDisplay('*', 'add')
%          waitbar(0.85,WB); figure(FIG);
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotAntPointTime', 'Type', 'axes'))
         for k=1:N_GS
            H = plot(t./60, 180/pi*AGLgs2ant(k,:),'-');
            set(H, 'Tag', 'PlotAntPointTime', 'UserData', 'Graph', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         end
         plot([0 t(sim_sample)/60], Beam/2*90/pi*[1 1], 'k+--','Tag', 'PlotAntPointTime', 'UserData', 'Graph');
         plot([0 t(sim_sample)/60], Beam*90/pi*[1 1], 'r+--','Tag', 'PlotAntPointTime', 'UserData', 'Graph');
         %delete(findobj(FIG,'Tag','PlotAntPointTime', 'Type', 'axes'));
         H = legend(['Mag pointing deviation from ', NAME_GS(1,:)],['Mag pointing deviation from ', NAME_GS(2,:)],...
             ['Mag pointing deviation from ', NAME_GS(3,:)], '-3dB Beam width', '-6dB Beam width');
         set(H, 'Tag', 'PlotAntPointTime', 'UserData', 'Legend');
         set(get(H, 'Children'), 'Tag', 'PlotAntPointTime', 'UserData', 'Legend');           
      end
      waitbar(0.50,WB); figure(FIG);
      %============================================================
      % Ground Trace 2-D Plot
      %============================================================
      if get(findobj(FIG,'Tag','CheckBoxPlotTrace'),'Value')
         fprintf('*')
         LocalDisplay('*', 'add')
%          waitbar(0.90,WB); figure(FIG);
         LONsc_link = ones(N_GS,1)*LONsc;
         LATsc_link = ones(N_GS,1)*LATsc;
         LONsc_ecl = LONsc;
         LATsc_ecl = LATsc;
         for k=1:sim_sample

            if Eclipse_state(k) == 0
               LONsc_ecl(k) = NaN;
               LATsc_ecl(k) = NaN;
            end
            for j=1:N_GS
               if Link_state(j,k) == 0
                  LONsc_link(j,k) = NaN;
                  LATsc_link(j,k) = NaN;
               end
            end
         end
         
         LocalDisplay('*', 'add')
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotGroundTrace', 'Type', 'axes'))%,'XGrid','on'));%'Box','off'));%,'ButtonDownFcn',''));
         % Plot ground traces
         H = plotTrace(180/pi*LONsc, 180/pi*LATsc,'g.-'); set(H, 'Tag','PlotGroundTrace','markersize',4,'linewidth',4);
         if(isnan(LATsc_ecl))
             H = plot(180/pi*-1, 0, 'k-'); set(H, 'Tag','PlotGroundTrace','markersize',4,'linewidth',4);
         else
             H = plotTrace(LONsc_ecl*180/pi, LATsc_ecl*180/pi, 'k-'); set(H, 'Tag','PlotGroundTrace','markersize',4,'linewidth',4);
         end
         H = plotTrace(180/pi*LLorbpole(1,:), 180/pi*LLorbpole(2,:),'g--'); set(H, 'Tag','PlotGroundTrace','markersize',4,'linewidth',2);
         H = plotTrace(180/pi*LONsun, 180/pi*LATsun,'y+-'); set(H, 'Tag','PlotGroundTrace','markersize',4,'linewidth',2);
         H = plot(180/pi*LON_MAG_SOUTH, 180/pi*LAT_MAG_SOUTH,'b+'); set(H, 'Tag','PlotGroundTrace','markersize',10,'linewidth',2);
         H = plot(180/pi*LON_MAG_NORTH, 180/pi*LAT_MAG_NORTH,'r+'); set(H, 'Tag','PlotGroundTrace','markersize',10,'linewidth',2);
         %H = plot(180/pi*LONray(1), 180/pi*LATray(1),'ko','markersize',4);


         H = plotTrace(180/pi*LLmagpt(1,:), 180/pi*LLmagpt(2,:), '--'); set(H, 'Tag','PlotGroundTrace', 'markersize',4,'color', [1 0 1]);
         
         LLattpt_trans = LLattpt;
         LLattpt_stab = LLattpt;
         for k=1:sim_sample
            if phi(k) > Stab_res + pi/180
               LLattpt_stab(:,k) = NaN;
            else
               LLattpt_trans(:,k) = NaN;
            end
         end
             
         H = plot(180/pi*LLattpt_trans(1,:), 180/pi*LLattpt_trans(2,:), '*'); set(H, 'Tag','PlotGroundTrace', 'markersize',2,'color', [1 0 0]);
         if(isnan(LLattpt_stab))
             H = plot(180/pi*-1, 0, '*'); set(H, 'Tag','PlotGroundTrace', 'markersize',2,'color', [0.5 1 0]);
         else
             H = plot(180/pi*LLattpt_stab(1,:), 180/pi*LLattpt_stab(2,:), '*'); set(H, 'Tag','PlotGroundTrace', 'markersize',2,'color', [0.5 1 0]);
         end
        
         for k=1:N_GS
            H = plot(LON_GS(k)*180/pi, LAT_GS(k)*180/pi, '.-', 'markersize',8,...
               'Tag','PlotGroundTrace', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         end
         %delete(findobj(FIG,'Tag','PlotGroundTraceLeg','UserData','Legend', 'Type', 'axes'));
         H=legend('Subsatellite trace','Trace in Eclipse','Orbit pole trace','Subsolar trace', 'Magnetic South', 'Magnetic North',...
            'Earth field pointing','S/C pointing (transient)','S/C pointing (stabilized)',...
            ['G/S contact ', NAME_GS(1,:)],['G/S contact ', NAME_GS(2,:)],['G/S contact ', NAME_GS(3,:)],4);
         set(H, 'Tag', 'PlotGroundTraceLeg','UserData','Legend','unit','pixels', 'position', [600 486 178 192]); %normalized 'Position',[0.45 0.505 0.14 0.20]);    
         %set(get(H, 'Children'), 'Tag', 'PlotGroundTrace'); 
%          set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotGroundTraceLeg','UserData','Legend', 'Type', 'axes'));         
         set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotGroundTrace', 'Type', 'axes'))%,'XGrid','on'));%'Box','off'));%,'ButtonDownFcn',''));
         
         % Compute and plot Solar Shade circles on graph (corrected for latitude deformation)
         %Rc1 = Re/tan(acos(Re/Rsc(1)))*ones(1,length(circle));
         %Rc2 = Re/tan(acos(Re/Rsc(fix(sim_sample/2))))*ones(1,length(circle));
         %Rc3 = Re/tan(acos(Re/Rsc(sim_sample)))*ones(1,length(circle));
         %LLray_cir1 = identSC(LLray(:,1), [circle; 0*circle], Rc1);			%[rad]
         %LLray_cir2 = identSC(LLray(:,fix(sim_sample/2)), [circle; 0*circle], Rc2);	%[rad]
         %LLray_cir3 = identSC(LLray(:,sim_sample), [circle; 0*circle], Rc3);		%[rad]
         %H = plotTrace(LLray_cir1(1,:)*180/pi, LLray_cir1(2,:)*180/pi, 'k--'); set(H, 'Tag','PlotGroundTrace');
         %H = plotTrace(LLray_cir2(1,:)*180/pi, LLray_cir2(2,:)*180/pi, 'k--'); set(H, 'Tag','PlotGroundTrace');
         %H = plotTrace(LLray_cir3(1,:)*180/pi, LLray_cir3(2,:)*180/pi, 'k--'); set(H, 'Tag','PlotGroundTrace');      
         
         %plot arrow markers
         H = plot(180/pi*[LONsc(1) LONsc(sim_sample)], 180/pi*[LATsc(1) LATsc(sim_sample)],'g>'); set(H, 'Tag','PlotGroundTrace', 'markersize',12);
         H = plot(180/pi*[LONsc(1) LONsc(sim_sample)], 180/pi*[LATsc(1) LATsc(sim_sample)],'g>'); set(H, 'Tag','PlotGroundTrace', 'markersize',10);
         H = plot(180/pi*[LONsc(1) LONsc(sim_sample)], 180/pi*[LATsc(1) LATsc(sim_sample)],'r>'); set(H, 'Tag','PlotGroundTrace', 'markersize',8);
         H = plot(180/pi*[LLorbpole(1,1) LLorbpole(1,sim_sample)], 180/pi*[LLorbpole(2,1) LLorbpole(2,sim_sample)],'g<'); set(H, 'Tag','PlotGroundTrace', 'markersize',12);
         H = plot(180/pi*[LLorbpole(1,1) LLorbpole(1,sim_sample)], 180/pi*[LLorbpole(2,1) LLorbpole(2,sim_sample)],'g<'); set(H, 'Tag','PlotGroundTrace', 'markersize',10);
         H = plot(180/pi*[LLorbpole(1,1) LLorbpole(1,sim_sample)], 180/pi*[LLorbpole(2,1) LLorbpole(2,sim_sample)],'r<'); set(H, 'Tag','PlotGroundTrace', 'markersize',8);
         H = plot(180/pi*[LONsun(1) LONsun(sim_sample)], 180/pi*[LATsun(1) LATsun(sim_sample)],'y<'); set(H, 'Tag','PlotGroundTrace', 'markersize',12);
         H = plot(180/pi*[LONsun(1) LONsun(sim_sample)], 180/pi*[LATsun(1) LATsun(sim_sample)],'y<'); set(H, 'Tag','PlotGroundTrace', 'markersize',10);
         H = plot(180/pi*[LONsun(1) LONsun(sim_sample)], 180/pi*[LATsun(1) LATsun(sim_sample)],'r<'); set(H, 'Tag','PlotGroundTrace', 'markersize',8);
         H = plot(180/pi*LON_MAG_SOUTH, 180/pi*LAT_MAG_SOUTH,'bo'); set(H, 'Tag','PlotGroundTrace','markersize',12,'linewidth',2);
         H = plot(180/pi*LON_MAG_NORTH, 180/pi*LAT_MAG_NORTH,'ro'); set(H, 'Tag','PlotGroundTrace','markersize',12,'linewidth',2);
         H = text(LON_MAG_SOUTH*180/pi-3, LAT_MAG_SOUTH*180/pi+6, 'S', 'color', [0 0 1], 'fontweight', 'bold', 'Tag','PlotGroundTrace');
         H = text(LON_MAG_NORTH*180/pi-3, LAT_MAG_NORTH*180/pi+6, 'N', 'color', [1 0 0], 'fontweight', 'bold', 'Tag','PlotGroundTrace');
                 
         for k=1:N_GS
            H = plot(LON_GS(k)*180/pi, LAT_GS(k)*180/pi, '+');
            set(H, 'Tag','PlotGroundTrace', 'markersize',8, 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
            H = plotTrace(LLgs_cir0(2*k-1,:)*180/pi, LLgs_cir0(2*k,:)*180/pi, 'b--');
            set(H, 'Tag','PlotGroundTrace','markersize',6,'linewidth',2, 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
            H = plotTrace(LLgs_cirMEL(2*k-1,:)*180/pi, LLgs_cirMEL(2*k,:)*180/pi, 'b-');
            set(H, 'Tag','PlotGroundTrace','markersize',6,'linewidth',2, 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
            H = plotTrace(LONsc_link(k,:)*180/pi, LATsc_link(k,:)*180/pi, '-'); 
            set(H, 'Tag','PlotGroundTrace','markersize',6,'linewidth',2, 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
            H = text(LON_GS(k)*180/pi+3, LAT_GS(k)*180/pi-2, NAME_GS(k,:), 'fontweight', 'bold');
            set(H, 'Tag','PlotGroundTrace', 'color', [(k-1)/max(1,(N_GS-1)) 0 1-(k-1)/max(1,(N_GS-1))/2]);
         end
      end
      
      waitbar(0.75,WB); figure(FIG);
      if get(findobj(FIG,'Tag','CheckBoxPlot3D'),'Value')
         %============================================================
         % Ground Trace 3-D Plot - Earth Referential (Geo)
         %============================================================
         fprintf('*')
         LocalDisplay('*', 'add')
%          waitbar(0.95,WB); figure(FIG);
         F3Dgeo = figure('position', [50 50 900 800],'Tag', 'PlotTrace3Dgeo');
         A3Dgeo = axes('Parent', F3Dgeo,'Tag','PlotTrace3Dgeo');      
         view(150,30);
         hold on;
         zoom(1.3)
         grid;
         axis('equal');
         axis(3*a*[-1 1 -1 1 -1 1])
         title(['Orbit System as wiewed in Earth Referential @ ',...
               str_alt,',  ',  str_incl,',  ', str_raan,',  ',str_argp,',  ', str_Msc0,',  ', str_time],...
            'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo', 'fontweight', 'bold');
         xlabel('X (Greenwich)','Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         ylabel('Y','Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         zlabel('Z (North Pole)','Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         
         DS = a*3; %Sun distace scale
         
         caxis('manual')
         caxis([-1 1])
         % Earth Sphere
         [Xe, Ye, Ze] = sphere(20);
         Xe = Xe*Re;
         Ye = Ye*Re;
         Ze = Ze*Re;
         surf(Xe, Ye, Ze, Ze/Re/10 - 0.7,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         
         % Equator & Greenwich
         plot(Re*cos(circle), Re*sin(circle), 'k-', 'linewidth',2,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         plot3((Re*cos(circle)), zeros(size(circle)), Re*sin(circle), 'k--', 'linewidth',2,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         plot3([0 0], [0 0], [-1 1]*Re*1.2, 'k-', 'linewidth',4,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         
         % Trace vector of objects
         plot3(XYZsc_geo(1,:),XYZsc_geo(2,:),XYZsc_geo(3,:),'g.-','Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         plot3(Re*IJKsc_geo(1,:),Re*IJKsc_geo(2,:),Re*IJKsc_geo(3,:),'g-', 'linewidth',2,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         plot3(XYZmagtr_geo(1,:),XYZmagtr_geo(2,:),XYZmagtr_geo(3,:),'r-', 'linewidth',2,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         plot3(DS*IJKsun_geo(1,:),DS*IJKsun_geo(2,:),DS*IJKsun_geo(3,:),'y.-','Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         
         % Sun Sphere
         [Xsun, Ysun, Zsun] = sphere(20);
         Xsun = Xsun*Re/4 + DS*IJKsun_geo(1,sim_sample);
         Ysun = Ysun*Re/4 + DS*IJKsun_geo(2,sim_sample);
         Zsun = Zsun*Re/4 + DS*IJKsun_geo(3,sim_sample);
         %surf(Xsun, Ysun, Zsun, Zsun/20 + 0.3*Re ,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo')
         surf(Xsun, Ysun, Zsun, -Zsun/Re/3 + 0.3 ,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo')

         % Sun Shadow
         %[Xecl, Yecl, Zecl] = cylinder(Re*[1 1], 60);
         %XYZshadow = (ROTsun_geo^(1))*[Zecl*Re*2; Xecl; Yecl]
         %mesh(XYZshadow(1,:),XYZshadow(2,:),XYZshadow(3,:),'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');  
         
         % Plot Ground Stations
         for k=1:N_GS
            plot3(XYZgs_cir0_geo(3*k-2,:),XYZgs_cir0_geo(3*k-1,:),XYZgs_cir0_geo(3*k,:),'y--',...
               'linewidth',2,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
            plot3(XYZgs_cirMEL_geo(3*k-2,:),XYZgs_cirMEL_geo(3*k-1,:),XYZgs_cirMEL_geo(3*k,:),'y-',...
               'linewidth',2,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
            plot3(XYZgs_geo(3*k-2),XYZgs_geo(3*k-1),XYZgs_geo(3*k),'y+',...
               'linewidth',2,'Parent', A3Dgeo,'Tag','PlotTrace3Dgeo');
         end
         
         hidden off
         shading interp
         rotate3d
         
         
         %============================================================
         % Ground Trace 3-D Plot - Sun Referential
         %============================================================
         fprintf('*')
         LocalDisplay('*', 'add')
%          waitbar(0.98,WB); figure(FIG);
         F3Dsun = figure('position', [75 75 900 800],'Tag', 'PlotTrace3Dsun');
         A3Dsun = axes('Parent', F3Dsun,'Tag','PlotTrace3Dsun');
         hold on;
         view(150,30);
         zoom(1.3);
         grid;
         axis('equal');
         axis(1.2*a*[-4 4 -4 4 -1 1]);
         title(['Orbit System as wiewed in Sun Referential @ ',...
               str_alt,',  ',  str_incl,',  ', str_raan,',  ',str_argp,',  ', str_Msc0,',  ', str_time],...
            'Parent', A3Dsun,'Tag','PlotTrace3Dsun', 'fontweight', 'bold');
         xlabel('Vernal Equinox ','Parent', A3Dsun,'Tag','PlotTrace3Dsun');
         ylabel('','Parent', A3Dsun,'Tag','PlotTrace3Dsun');
         zlabel('Normal to Eccliptic Plane','Parent', A3Dsun,'Tag','PlotTrace3Dsun');
         
         caxis('manual')
         caxis([-1 1])
         % Sun Sphere
         [Xsun, Ysun, Zsun] = sphere(20);
         Xsun = Xsun*Re;
         Ysun = Ysun*Re;
         Zsun = Zsun*Re;
         surf(Xsun, Ysun, Zsun, -Zsun/10/Re + 0.3 ,'Parent', A3Dsun,'Tag','PlotTrace3Dsun');
         
         % Other Traces
         plot3(XYZsc_sun(1,:),XYZsc_sun(2,:),XYZsc_sun(3,:),'g-','Parent', A3Dsun,'Tag','PlotTrace3Dsun');
         plot3(XYZgeo_sun(1,:),XYZgeo_sun(2,:),XYZgeo_sun(3,:),'b+','Parent', A3Dsun,'Tag','PlotTrace3Dsun');
         plot3(XYZzero_sun(1,:),XYZzero_sun(2,:),XYZzero_sun(3,:),'r-','Parent', A3Dsun,'Tag','PlotTrace3Dsun');
         
         hidden off
         shading interp
         rotate3d
         %============================================================
      end
      
      clear global Qtot2 Qtot7 Qtot7b
      
      waitbar(1,WB); figure(FIG);
      LocalDisplay('*|  completed.', 'add')
      
      fprintf('|\n Simulation Completed.\n')
      close(WB);
      msgbox('Simulation Completed');
   end
   %============================================================

   
   
   
%************************************************************
% Command Buttons & Other Actions 
%************************************************************   
elseif strcmp(action,'info')
   LocalDisplayInfo;
   
elseif strcmp(action,'simInfo')
   HelpStr={'';
      ' SIMULATION';
      ' ----------';
      ' Enter the simulation parameters.  The fields will update automatically';
      ' Check the appropriate boxes if graphs are to be processed and displayed.';
      ' Press the ''Run'' button to start the simulation. The simulation progress';
      ' is displayed in the text window at the bottom of the figure.';
      ' ';
      ' The simulation cannot be cancelled from the control interface.';
      ' Return to the Matlab command window and press Ctrl+C to terminate';
      ' the simulation before its end.';
      ' ';
      ' FIELDS DESCRIPTION';
      ' ------------------';
      ' T [rev] : simulation duration in orbit revolution';
      ' T [min] : simulation duration in minutes';
      ' N []    : number of samples used for simulation (see notes below)';
      ' fs [Hz] : sampling frequency in Hertz';
      ' Ts [s]  : sampling time in seconds';
      ' ';
      ' CHECK BOXES';
      ' -----------';
      ' These boxes are used to enable/disable some time-consuming processes';
      ' '; 
      ' - Check "Thermal Sim" to process the radiation and temperature ODE simulations';
      ' - Check "2-D Plot" to process the multiple simulation result graphs';
      ' - Check "3-D Plot" to process the additional 3-D orbits and ground traces';
      ' ';
      ' NOTES';
      ' -----';
      ' For attitude and thermal simulation purposes,';
      ' high initial dynamics might require high sampling';
      ' rates.  For instance, a high initial rate wb will';
      ' provoque high spin rates which must be sampled and';
      ' processed at a high frequency (in the order of the seconds)';
      ' insuring the stability of the simulations and the ';
      ' accuracy of the results.  But, high samplig rates';
      ' will significantly augment the simulation processing time';
      ' '; };
   helpwin(HelpStr, 'Simulation Info');
   
elseif strcmp(action,'orbitInfo')
   HelpStr={'';
      ' ORBIT';
      ' -----';
      ' The orbit is propagated over the simulation duration. The ephemeris';
      ' (S/C location) is computed taking into account Sun''s and Moon''s ';
      ' gravitational perturbation on the orbit as well as the geopotential';
      ' J2 parameter(secular variations).  The air drag is also computed and ';
      ' applied to the orbit perturbation forces.  Mean air density is used at ';
      ' every S/C position on the orbit to compute the drag and the orbit decay. ';
      ' ';
      ' Enter the orbit parameters at Epoch. The fields will update automatically';
      ' ';

      ' FIELDS DESCRIPTION';
      ' ------------------';
      ' (self explanatory)';
      ' ';
      ' ';
      ' ';
      ' ';
      ' ';
      ' ';
      ' ';
      ' ';
      ' ';
      ' ';
      ' '; };
   helpwin(HelpStr, 'Orbit Info');
   
elseif strcmp(action,'scInfo')
   HelpStr={'';
      ' SPACECRAFT DESCRIPTION';
      ' ----------------------';
      ' Enter the spacecraft parameters. ';
      ' The simulation assumes a rectangular-shaped spacecraft bus';
      ' with a square base.  The body frame z-axis is defined normal';
      ' to the square base of the bus.  The x and y-axes are defined';
      ' along the edges of the square base';
      ' ';
      ' The fields will update automatically';
      ' ';
      ' FIELDS DESCRIPTION';
      ' ------------------';
      ' Ixx, Iyy, Izz [kg.mm^2] : principal moments of inertia ';
      ' Ixy, Ixz, Iyz [kg.mm^2] : products of inertia (enter negative values)';
      ' Mass [kg]    : total spacecraft mass ';
      ' Width [m]    : spacecraft square base dimension';
      ' Length [m]   : length along the P-POD loading axis (also body z-axis)';
      ' Face 1 [m²]  : square base area';
      ' Face 2 [m²]  : long sides area';
      ' Wall [mm]    : external wall thickness for all 6 faces';
      ' Drag coef [] : atmospheric drag coefficient Cd (usually 1 to 4)';
      ' Cb [kg/m²]   : ballistic coefficient (mass/Cd/Area)';
      ' PV eff %     : solar cells (photovoltaics) efficiency';
      ' PV pack1 %    : solar cells (PV) packing factor on the square base';
      ' PV pack2 %    : solar cells (PV) packing factor on the long sides';
      ' Heat [W]     : constant internal heat generation in Watt';
      ' '; };
   helpwin(HelpStr, 'Spacecraft Info');
  
elseif strcmp(action,'magInfo')
   HelpStr={'';
      ' PASSIVE MAGNETIC ATTITUDE STABILIZATION';
      ' ---------------------------------------';
      ' About the magnetic stabilization, a system of permanent magnets';
      ' and hysteresis bars is simulated.  The magnets are placed to rotate';
      ' the spacecraft along the Earth''s magnetic field and the bars are';
      ' placed perpendicularily to the magnets.  The behavior of the libration';
      ' dynamics about the Earth''s field line is determined by the magnetic';
      ' moment of the 2 elements.   Alico-5 supermagnets and HyMu-80 bars are';
      ' the materials used for the simulation.  For given material induction';
      ' properties, the magnetic moment is given by the volume of material';
      ' and the permittivitiy of vacuum.  With proper density, the magnetic moment';
      ' computed and used for the simulation';
      ' ';
      ' The computation of the system magnetic moment assumes perfect dipole';
      ' distributions of the field around the magnets.  The magnetic moment is';
      ' based on the saturation field of the material and the relation';
      ' ';
      '                m = B.V/uo';
      ' where';
      '  m = magnetic moment [A.m^2]';
      '  B = rated field [Tesla]';
      '  V = material volume [m^3]';
      '  uo = permeability of vacuum = 4.Pi.1e-7 [1/A^2]';
      ' ';
      ' The hysteresis bars are made of high-permeability (soft) magnetic';
      ' materials in contrast to permanent (hard) magnetic materials. ';
      ' They work like permanent magnets but can have their magnetic';
      ' domain boundairies easily changed under induction from an external';
      ' magnetic field.  The displacement of these boundaries occurs with';
      ' molecular friction.  The hysteresis curve shows this phenomenon where the';
      ' area inside the curve is a measure of energy loss per cycle of unit volume';
      ' of material.  The energy is dissipated and transformed into heat.';
      ' ';
      ' A displacement of boundaries change the material magnetic'; 
      ' moment and may also reverse the polarity of the effective dipole.';
      ' In order to overcome the remanence (Br), and thus get a flip of the';
      ' magnetic polarization inside the hysteresis material, an external field';
      ' strength threshold equal to the material coercive field (Hc) must be met.';
      ' This field is assumed to come only from the Earth magnetic field. It has';
      ' an effective strength on the soft material depending on the location and';
      ' attitude of the spacecraft on its orbit.  Thus, the relative angle';
      ' between the hysteresis bars and the field lines (sinus of angle) is used';
      ' to compute the effective field value.  The Earth''s field is not strong';
      ' enough to act the same way on permanent magnets (their coercivity is much';
      ' higher than soft materials).  For optimal authority and performance,';
      ' the hyseresis bars must be placed perpendicularly to the permanent';
      ' magnets along both cross-axis';
      ' ';
      ' For simulation a L-Shell model of the Earth''s magnetic field is used';
      ' (SMAD p.213) The Earth''s magnetic field diagram is displayed when and';
      ' displayed when the "Mag Field" button is pressed.  The field vector';
      ' (strength and orientation) is computed at any location of the spacecraft';
      ' and for the whole simulation.  The data are used as input to a Simulink';
      ' attitude model. The model is used to simulate the non-linear behavior';
      ' of the hysteresis bars and magnets system dynamics. Three models of';
      ' hysteresis have been implemented in Simulink (press ''Plot Models''';
      ' to see):  ';
      ' ';
      '  Model 1: RAMP of (Br/Hc) + DELAY with Hc as hysteresis coercivity';
      '  Model 2: RAMP of (Bs/2Hc) + DELAY with Hc as hysteresis coercivity';
      '  Model 3: DELAY only, Hc as hysteresis coercivity, Bs as saturation level';
      '  (Bs=saturation induction, Br=remanence (H=0), Hc= coercive field)';
      ' ';
      ' The magnetic torque for each component of the stabilisation system uses';
      ' the usual vector product law for a material inside a magnetic field:';
      ' ';
      '                T = m x Be';
      ' where';
      '  Be = Earth magnetic field strength vector at S/C location [Tesla]';
      '  m =  magnetic moment [A.m^2]';
      '  T = magnetic torque vector [N.m]';
      ' ';
      ' The simulation is continuous in 3-D using single and rigid body dynamic';
      ' model based on momentums (The Simulink model is found under "magSim3.mdl")';
      ' Therefore, nutation and precession motions are completely described and ';
      ' observed by the model';
      ' ';
      ' The results usually display a damped oscillation about the Earth''s';
      ' magnetic field lines but at one critical point, the systems response';
      ' is edged to the magnetic field lines where no more damping occurs.  Since';
      ' the hysteresis bars are acting also as magnets, the dual-component system';
      ' has its resulting magnetic moment tilted off the desired axis which was';
      ' along the permanent magnets.  This new resulting axis is not fixed and';
      ' flips with the change of polarity inside the bars.  When the oscillations';
      ' gets bellow the point where the Earth''s effective field strength can no';
      ' longer reverse the polarity inside the hysteresis bars, the magnetic axis';
      ' become fixed and the system will continue its libration without energy';
      ' dissipation through friction from magnetic boundary displacement. The';
      ' stabilization resolution is easily computed for any vector combination';
      ' of permanent magnets (Alnico-5) and soft magnetic materials (HyMu-80)';
      ' ';
      ' Results for different volume (or mass) combinations of materials show that';
      ' the residual libration (resolution) of the system is a trade-off between';
      ' the volume (or mass) of HyMu and Alnico.  Therefore, the resolution, or';
      ' the stabilisation axis offset error cannot be brought to a near-zero value';
      ' unless the mass of magnets is greatly higher than of the hysteresis rods.';
      ' In such case, the response will be very undamped and the transient';
      ' oscillation will last for many orbits.  Also, in this case, the hysteresis';
      ' has a minor contribution on the system and other perturbations, such as';
      ' Eddy currents and internal cablings magnetic moment, might then have a';
      ' relevant impact on the systems dynamics.';
      ' ';
      ' ';
      ' see paper for more details';
      ' J-F. LEVESQUE, Passive Magnetic Attitude Stabilization using Hysteresis';
      '     Materials, Universite de Sherbrooke, SIgMA-PU-006-UdeS, Sept 2003.';
      ' ';
      ' ';
      ' FIELDS DESCRIPTION';
      ' ------------------';
      ' Alnico-5  : permanent magnet material used (Spec. Gravity = 7.01)';
      ' HyMu80    : hysteresis material used (Br=0.35, Hc=1.0, SG=8.74)';
      ' ';
      ' Direction : S/C body components of the magnetic moments direction,';
      '             (alignment of magnetic material North with respect to';
      '             S/C body frame).  The direction vector will be normalized';
      '             (made unitary) during process';
      ' m [g]     : mass of material used';
      ' V [cm^3]  : volume of material (computed with specific gravity)';
      ' Bs [T]    : saturation field strength in Tesla ([T] = 10^4[Gauss])';
      ' M [A.m^2] : magnetic moment computed using ';
      '             simple undeformed dipole model (M = V.Bs/uo)';
      ' ';
      ' Stab resolution [°] : computed stabilization final libration amplitude';
      ' Dyn freq [rad/s]    : natural frequency of the system about smallest';
      '                       principal moment of inertia assuming small angles';
      ' Init rate [rad/s]   : initial angular rate S/C body components [x; y; z]';
      ' Init quaternions    : initial S/C attitude [q1; q2; q3; q4]';
      '                       (use the button "|q|" to normalize the unit vector)';
      ' ';
      ' --------------------------------';
      ' (C) 2003, Jean-Francois Levesque'; 
      ' ';  };
   helpwin(HelpStr, 'Magnetic Attitude Stabilization Info');
  
elseif strcmp(action,'gsInfo')
   HelpStr={'';
      ' GROUND STATIONS';
      ' ---------------';
      ' Enter the location parameters for up to three ground stations.';
      ' ';
      ' FIELDS DESCRIPTION';
      ' ------------------';
      ' Name     : ground station acronyme (use 3 letters)';
      ' LON [°]  : East longitude of the ground stations location';
      ' LAT [°]  : North latitude of the ground stations location';
      ' Elev [°] : minimum horizon elevation of the ground station';
      ' [MHz]    : communication link frequency used to compute the Doppler shift';
      ' ';
      };
   helpwin(HelpStr, 'Ground Stations Info');
  
      
elseif strcmp(action,'history')
   FIG = findobj('Tag','ControlWindow');
   DataStr = get(findobj(FIG,'Tag','ButtonHistory'),'UserData');
   helpwin(DataStr, 'Sim History');            
      
elseif strcmp(action,'plot')
   LocalUpdateGraphs;
   
elseif strcmp(action,'saveGraphs')
   LocalSaveGraphs;
   
elseif strcmp(action,'close')
   clear global XYZfield_I umag uhyst1 uhyst2 Bmag Bs1 Bs2 Vmag Vhyst1 Vhyst2 uo Br Hc sw Jinv H0 t fs q0 Qtot2 Qtot7 Qtot7b angle tout;
   close(findobj('Tag','PlotTrace3Dgeo','Type','figure'));
   close(findobj('Tag','PlotTrace3Dsun','Type','figure'));   
   close(findobj('Tag','PlotMagField','Type','figure'));
   close(findobj('Tag','PlotHyst','Type','figure'));
   close(findobj('Tag','ControlWindow'));

elseif strcmp(action,'clearGraph')
   FIG = findobj('Tag','ControlWindow');
   A = findobj(FIG, 'Type', 'axes');
   for k=1:length(A)
      set(FIG,'CurrentAxes',A(k));
      cla;
   end
   close(findobj('Tag','PlotTrace3Dgeo','Type','figure'));
   close(findobj('Tag','PlotTrace3Dsun','Type','figure'));   
   close(findobj('Tag','PlotMagField','Type','figure'));
   close(findobj('Tag','PlotHist','Type','figure'));
   LocalDisplay('Graphs Cleared', 'nl');
   
%elseif strcmp(action,'setDefault')
 
   
elseif strcmp(action,'hideLeg')
   FIG = findobj('Tag','ControlWindow');    
   if get(findobj(FIG,'Tag','CheckBoxPlotHide'),'Value')
      set(findobj(FIG,'Tag','PlotGroundTraceLeg'),'visible','off')
   else
      set(findobj(FIG,'Tag','PlotGroundTraceLeg'),'visible','on')
   end
   
elseif strcmp(action,'plotMag')
   FIG = findobj('Tag','ControlWindow');
   incl = sscanf(get(findobj(FIG, 'Tag', 'EditTextIncl'),'String'), '%f')*pi/180;	%[rad]
   e = sscanf(get(findobj(FIG, 'Tag', 'EditTextEcc'),'String'), '%f');     	%[]
   n = get(findobj(FIG, 'Tag', 'EditTextn'),'UserData'); 		            	%[rev/day]
   a = (sqrt(GMe)/(2*pi)*Te/n)^(2/3);							%[km] semi-major axis
   
   close(findobj('Tag','PlotMagField','Type','figure'));   
   plotMagField(incl, a, e);   
   
elseif strcmp(action,'plotHyst')
   FIG = findobj('Tag','ControlWindow'); 
   Bs = sscanf(get(findobj(FIG, 'Tag', 'EditTextBs1'),'String'), '%f');	%[A/m^2]
   Br = 0.35;
   Hc = 1.59;
      
   close(findobj('Tag','PlotHyst','Type','figure'));   
   plotHysteresis(Br, Bs, Hc);   
   
elseif strcmp(action,'plotOrbit')
   FIG = findobj('Tag','ControlWindow');
   e = sscanf(get(findobj(FIG, 'Tag', 'EditTextEcc'),'String'), '%f');     	%[]
   n = get(findobj(FIG, 'Tag', 'EditTextn'),'UserData'); 		         %[rev/day]
   
   a = (sqrt(GMe)/(2*pi)*Te/n)^(2/3);	%[km] semi-major axis
   rp = a*(1-e);                 	%[km] perigee
   ra = a*(1+e);                 	%[km] apogee
   c = a - rp;                      	%[km] orbit center distance from focii
   b = sqrt(a^2-c^2);               	%[km] semi-minor axis
   
   circle = linspace(-pi,pi,91);			%[rad]
   x_ellipse = linspace(-a, a, 100);
   y_ellipse = b.*sqrt(1-x_ellipse.^2/a^2);
   x_ellipse2 = [x_ellipse, x_ellipse]-c;
   y_ellipse2 = [y_ellipse, -y_ellipse];
   
   str_e = sprintf('e = %0.4g',e);
   str_a = sprintf('a = %0.5gkm',a);
   str_b = sprintf('b = %0.5gkm',b);
   str_c = sprintf('c = %0.5gkm',c);
   str_rp = sprintf('rp = %0.5gkm',rp);
   str_ra = sprintf('ra = %0.5gkm',ra); 
   
   set(FIG,'CurrentAxes',findobj(FIG,'Tag','PlotOrbShape', 'UserData', 'Graph', 'Type', 'axes'));%,'ButtonDownFcn',''));
   cla;
   delete(findobj(FIG,'Tag','PlotOrbShape','UserData','Legend', 'Type', 'axes'))
   
   axis([-1.5*ra 1.5*rp -1.2*b 1.2*b]);
   H = plot(Re.*cos(circle), Re.*sin(circle), 'g-', 'linewidth', 2, 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = plot(x_ellipse2, y_ellipse2,'r-', 'linewidth', 2, 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = plot([-c -c], [-b b], 'b--', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = plot([-ra rp], [0 0], 'b--', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = plot(0,0,'r*',-2*c,0,'g*', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = plot(-ra,0,'r*',rp,0,'r*', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = text(-2*c/3, b/12, str_a, 'color', [0 0 1], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = text(-c+300, 3*b/4, str_b, 'color', [0 0 1], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = text(-2*c/3, -b/12, str_c, 'color', [0 0 1], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = text(rp+300, b/12, str_rp, 'color', [1 0 0], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = text(-ra+300, b/12, str_ra, 'color', [1 0 0], 'fontweight', 'bold', 'Tag', 'PlotOrbShape', 'UserData', 'Graph');
   H = legend('Earth & Focii', ['Orbit (',str_e,')']);
   set(H, 'Tag', 'PlotOrbShape', 'UserData', 'Legend');                            
   set(get(H, 'Children'), 'Tag', 'PlotOrbShape', 'UserData', 'Legend');      
   
   
   
   
%*******************************************************************************  
end
  
%if nargout > 0, fig = FIG; end



%*******************************************************************************
% Built-in Functions
%*******************************************************************************  
function LocalUpdateGraphs();
   FIG = findobj('Tag','ControlWindow');
   G = findobj(FIG, 'UserData', 'Graph');
   L = findobj(FIG, 'UserData', 'Legend');
   for k=1:length(G)
       set(G(k),'Visible','off'); %hide secondary graphs
   end
   for k=1:length(L)
       set(L(k),'Visible','off'); %hide secondary legends
   end
   Popup = findobj(FIG, 'Tag', 'PopupPlot');
   cmdNum = get(Popup,'Value');
   cmdList = get(Popup,'UserData');
   eval(cmdList(cmdNum,:));	%set visible current popuplist graph

   
%================================================================  
function LocalSaveGraphs();
   FIG = findobj('Tag','ControlWindow');
   Popup = findobj(FIG, 'Tag', 'PopupPlot');
   cmd = get(Popup,'UserData');
   
   G = findobj(FIG, 'UserData', 'Graph');
   L = findobj(FIG, 'UserData', 'Legend');
   N = 18;              %number of popup graphs
   
   % save directory and file format
   DIRR = uigetdir('','Select destination directory for graphic files');
   if DIRR == 0
       warndlg('Saving operation aborted.');
       return;
   end
   type = menu('Select an image format','BMP','JPEG','TIFF')
   if type == 1
     FORM = 'bmp';
     EXT = 'bmp';
   elseif type == 2
     FORM = 'jpeg'
     EXT = 'jpg';
   elseif type == 3
     FORM = 'tiff'
     EXT = 'tif';
   end
   
   % move primary graph (ground trace) to temporaty figure for exportation
   % to file
   AXE = findobj(FIG, 'Tag','PlotGroundTrace','type','axe');
   LEG = findobj(FIG, 'Tag','PlotGroundTraceLeg','type','axe');
   POSA = get(AXE,'Position');
   POSL = get(LEG,'Position');
   FIG1 = figure('Position',[100 300 POSA(3)+90 POSA(4)+90], 'ToolBar','none','PaperOrientation','landscape','PaperPositionMode','auto');   

   set(AXE, 'Parent', FIG1, 'Position', [55,45,POSA(3),POSA(4)]);
   if LEG
       set(LEG, 'Parent', FIG1, 'Position', [605,46,POSL(3),POSL(4)]);
   end
   fname = 'GroundTrace';
   print(FIG1, '-djpeg', [DIRR,'\',fname,'.',EXT])
   set(AXE, 'Parent', FIG, 'Position', POSA);
   if LEG
       set(LEG, 'Parent', FIG, 'Position', POSL);
   end
   close(FIG1)
   
   % move secondary axes to temporary figure for exportation to file
   FIG2 = figure('Position',[100 400 960 480], 'ToolBar','none','PaperOrientation','landscape','PaperPositionMode','auto');
   for g = 1:N
       for k=1:length(G)
          set(G(k),'Visible','off'); %hide secondary graphs
       end
       for k=1:length(L)
          set(L(k),'Visible','off'); %hide secondary legends
       end

       eval(cmd(g,:));	    %set visible current popuplist graph
       
       AXE = findobj(FIG, 'UserData', 'Graph','type','axe','Visible','on');
       LEG = findobj(FIG, 'UserData', 'Legend','type','axe','Visible','on');
       if AXE
           set(AXE, 'Parent', FIG2);
       end
       if LEG
           set(LEG,'unit','pixels')
           set(LEG, 'Parent', FIG2);
       end
              
       fname = get(AXE,'Tag');  %use axis name tag as filename
       print(FIG2, '-djpeg', [DIRR,'\',fname,'.',EXT])
       
       if AXE
           set(AXE, 'Parent', FIG);
       end
       if LEG
           set(LEG, 'Parent', FIG);
           set(LEG,'unit','normalized','UserData','Legend')
       end
   end
   close(FIG2)
   msgbox(['All graphs have been saved to ',DIRR]);
   
%================================================================   
function LocalDisplay(inStr, loc);
   FIG = findobj('Tag','ControlWindow');
   H = findobj(FIG,'Tag','TextDisplay');
   if strcmp(loc, 'add')
      dispStr = get(H, 'String');
      n = length(dispStr);
      if n==0
         dispStr = {inStr};
      elseif n==1
         lastStr = dispStr(1);
         dispStr = {[lastStr{:},inStr]};
      else
         lastStr = dispStr(1);
         dispStr = {[lastStr{:},inStr],dispStr{2:n}};
      end
   elseif strcmp(loc, 'nl')
      dispStr = get(H, 'String');
      %dispStr = str2mat(get(H, 'String'),inStr);
      dispStr = {inStr,dispStr{:}};
   else
      dispStr = {inStr};
   end
   set(H,'String', dispStr);
   
function  LocalDisplayInfo();
   HelpStr={
      '****************************************************************************';
      ' CUBESIM - CubeSat Orbit & Attitude Simulation';
      '    ';
      '    version 2.5, June 4, 2004.';
      '    ';
      '    Copyright (C) 2002-2004 Jean-Francois Levesque';
      '    ';
      ' This software is provided ''as-is'', without any express or implied';
      ' warranty.  In no event will the authors be held liable for any damages';
      ' arising from the use of this software.';
      ' ';
      ' Permission is granted to anyone to use this software for any purpose,';
      ' including commercial applications, and to alter it and redistribute it';
      ' freely, subject to the following restrictions:';
      '  ';
      ' 1. The origin of this software must not be misrepresented; you must not';
      '    claim that you wrote the original software. If you use this software';
      '    in a product, an acknowledgment in the product documentation would be';
      '    appreciated but is not required.';
      ' 2. Altered source versions must be plainly marked as such, and must not be';
      '    misrepresented as being the original software. ';
      ' 3. This notice may not be removed or altered from any source distribution.';
      ' ';
      ' Jean-Francois Levesque, MS.';
      ' jflev@yahoo.ca';
      ' ';
      '****************************************************************************';
      ' ';
      ' - Orbit propagation and ephemeris calculation';
      ' - Attitude determination using passive magnetic stabilization and hysteresis';
      ' - Antenna pointing computation ';
      ' - Ground trace plots';
      ' - Ground station tracking information and S/C visibility reports';
      ' - Communication link occurences';
      ' - Eclipse computation and reports';
      ' - Solar cells power generation';
      ' - Thermal simulations for 2nd and 7th order models';
      ' - 22 simulation graphics and plots';
      ' - Simulation results and data log on a text window';
      ' ';
      ' MAJOR VERSION UPDATES';
      ' ---------------------';
      ' - v1.5 Magnetic stabilization 3-D Simulink model based on rigid body dynamics';
      ' - v1.5 Added moment of inertia tensor components input fields';
      ' - v1.5 Added passive magnetic components description input fields';
      ' - v1.5 Added initial attitude and rate input fields';
      ' - v1.6 Corrected atmospheric density computation with exponential interpolation';
      ' - v1.6 Quaternions and rates computed as attitude outputs';
      ' - v1.6 Added magnetic angle and attitude graphs';
      ' - v1.6 Orbit propagation using more accurate algorythms for eccentric orbits';
      ' - v1.6 Added internal heat generation input field';
      ' - v1.6 Revised solar Albedo calculations';
      ' - v1.6 Display text box scrollable';
      ' - v1.6 Added info buttons for field inputs';
      ' - v1.6 Added check boxes for graph plotting';
      ' - v1.7 Added antenna pointing deviation graphs';
      ' - v1.8 Corrected UT and GST calculation';
      ' - v1.8 Detailed magnetic hysteresis stabilization explanations';
      ' - v1.9 Added 2 new models for magnetic hysteresis';
      ' - v1.9 Parameters corrected in Simulink Model';
      ' - v2.0 Added plot of angular rate';
      ' - v2.1 Accurate World map';
      ' - v2.1 Graph titles updated';
      ' - v2.2 Field input parameters are now saved and restored automatically'
      ' - v2.2 Restore default with RESET button'
      ' - v2.2 Added plots of hysteresis models'
      ' - v2.2 Correction in sampling frequency and time text fields'
      ' - v2.3 Ground trace Legend tag names display bugs'
      ' - v2.4 Replaced variable i (used as Matlab complex symbol) by variable k'
      ' - v2.4 Correction of Maximum eclipse calculation'
      ' - v2.5 Added simulation progression bar'
      ' - v2.5 Added control window menubar'
      ' - v2.5 Added simulation popup messages'
      ' - v2.5 Changed tag names for legends'
      ' - v2.5 Added graphs saving function to file and button'
      ' - v2.6 Correction of bugs for invalid object handles';
      ' - v2.6 Updated info question marks';
      ' - v2.7 Correction of conflicts with new function names of Matlab v7'
      ' ';
      ' NOTES';
      ' -----';
      ' - The graphical interface is set for display 1280x1024.';
      ' ';
      ' - The simulation takes a lot of processing power and it is recommended';
      '   to run the simulation with a smaller number of samples than default';
      '   when running on a slow machine.';
      ' ';
      ' - Use the question mark buttons on the control interface to get';
      '   more details about the simulation input parameters.';
      ' ';    
      '****************************************************************************';
      ' ';

      ' ORBIT PROPAGATION AND EPHEMERIS';
      ' -------------------------------';
      ' The orbit is propagated over the simulation duration. The ephemeris';
      ' (S/C location) is computed taking into account Sun''s and Moon''s ';
      ' gravitational perturbation on the orbit as well as the geopotential';
      ' J2 parameter(secular variations).  The air drag is also computed and ';
      ' applied to the orbit perturbation forces.  Mean air density is used at ';
      ' every S/C position on the orbit to compute the drag and the orbit decay. ';
      ' ';
      ' ';

      ' ATTITUDE SIMULATION';
      ' -------------------';
      ' A simulink model is used to compute the S/C attitude in 3-D at every';
      ' position. The single and rigid body model is used.  External moments';
      ' (torques) are assumed to come only from the magnetic moment.  Gravity';
      ' gradients and Eddy current disturbances are neglected. The magnetic ';
      ' stabilization inputs and the magnetic field vector are used as inputs';
      ' to the simulations.';
      ' ';
      ' About the magnetic stabilization, a system of permanent magnets';
      ' and hysteresis bars is simulated.  The magnets are placed to rotate';
      ' the spacecraft along the Earth''s magnetic field and the bars are';
      ' placed perpendicularly to the magnets.  The behavior of the libration';
      ' dynamics about the Earth''s field line is determined by the magnetic';
      ' moment of the 2 elements.   Alico-5 supermagnets and HyMu80 bars are';
      ' the materials used for the simulation.  For given material induction';
      ' properties, the magnetic moment is given by the volume of material';
      ' and the permeability of vacuum.  With proper density, the magnetic';
      ' moment computed and used for the simulation';
      ' ';
      ' The hysteresis bars works like magnets but can have their magnetic';
      ' domain boundairies changed under induction from an external magnetic';
      ' field.  The displacement of these boundaries is not without molecular'; 
      ' friction.  The hysteresis curve shows this phenomenon where the area';
      ' inside the curve represents the energy dissipated and transformed into';
      ' heat.  A displacement of boundaries change the materials magnetic'; 
      ' moment and may also reverse the polarity of the effective dipole.';
      ' In order to get a flip of the magnetic polarization inside the';
      ' hysteresis material, an external field strength threshold (Hc) must'
      ' be met.  This field is assumed to be from the Earth and has a strength';
      ' and orientation depending on the location of the spacecraft around';
      ' the Earth, but also from the relative angle between the hysteresis';
      ' bars and the field lines (sinus of angle).  ';

      ' ';
      ' For simulation a L-Shell model of the Earth''s magnetic field is used';
      ' (SMAD p.213) The Earth''s magnetic field diagram is displayed when and';
      ' displayed when the "Mag Field" button is pressed.  The field vector';
      ' (strength and orientation) is computed at any location of the spacecraft';
      ' and for the whole simulation.  The data are used as input to a Simulink';
      ' attitude model. The model is used to simulate the non-linear behavior';
      ' of the hysteresis bars and magnets system dynamics. Simple models';
      ' appriximating the hysteresis behavior has been made of DELAY, RAMP';
      ' & SATURATION.  The simulation is continuous in 3-D using single and ';
      ' rigid body dynamic model based on momentums (The Simulink model is found';
      ' under "magSim3.mdl")';
      ' ';
      ' (See Stabilization help for more details)';
      ' ';
      ' ';
      ' RADIATION AND SOLAR INPUT POWER';
      ' -------------------------------';
      ' The radiation on the spacecraft is computed for every faces.  Radiation';
      ' agents taken into acount are the sun, its albedo as direct reflective';
      ' from Earth''s horizon and the Earth infrared radiation.  The computation';

      ' is vectorial using transformation matrixes to determine directions,';
      ' orientations, and view factors.  The spacecrafts attitude for the';
      ' calculations is given from the previous magnetic stabilization data.';
      ' two packing factors are used, one for the top and bottom faces, the other';
      ' for the side walls of the cube.';
      ' ';
      ' ode45 is used with a 2nd and a 7th order model to simulate the thermal ';
      ' behavior of the spacecraft assuming a given generated power and ';
      ' materials.  (The models are programmed in the scripts "T_ode2.m" and';
      ' T_ode7.m")';
      ' ';
      '****************************************************************************';
      ' ';
      ' REFERENCES';

      ' ----------';
      ' - J-F. LEVESQUE, Passive Magnetic Attitude Stabilization using Hysteresis';
      '     Materials, Universite de Sherbrooke, SIgMA-PU-006-UdeS, Sept 2003.';
      ' - LARSON, WERTZ, Space Mission Analysis and Design, Microcosm Press,';
      '     3rd edition, 1999.';
      ' - WERTZ, Spacecraft Attitude Determination and Control, Microcosm, 1978.';
      ' - BRYSON, Control of Spacecraft and Aircraft, Princeton University Press,';
      '     1994.';
      ' - INCROPERA, DEWITT, Fundamentals of Heat and Mass Transfer, Wiley,';
      '     4th edition, 1996.';
      ' - MUNSON, YOUNG, OKIISHI, Fundamental of Fluid Mechanics, Wiley,';
      '     3rd edition, 1998.';
      ' - KREYZSIG, Advanced Engineering Mathematics, Wiley, 7th Edition, 1993.';
      ' - SADIKU, Elements of Electromagnetics, Saunders College Publishing, 1989.';
      ' - Edward Della Torre, Magnetic Hysteresis, IEEE Press, New-York, 1999.';
      ' - AMSAT, https://www.amsat.org';
      ' - CUBESAT, http://cubesat.calpoly.edu';
      ' - NARCISSAT, http://ssdl.stanford.edu/narcissat';
      ' ';
      '****************************************************************************';
      ' (C) 2002-2005, Jean-Francois Levesque';       
   };
   helpwin(HelpStr, 'CUBESIM help');
   

%*******************************************************************************   
%*  (C) 2002-2005, Jean-Francois Levesque

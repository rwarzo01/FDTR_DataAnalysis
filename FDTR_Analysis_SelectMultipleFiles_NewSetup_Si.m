% FDTR Analysis with ability to select multiple files 
% 23JUN19, Warzoha
clc; clear; close all

%% User Inputs
% Settings
delete=[0 1]; %Delete number of values at beginning and end. e.g. [2 3] deletes two in beginning, three at end
print_confidence=1; %Set to 1 to print confidence intervals
    alpha=0.1; %set alpha (.1 for 90% confidence)
print_to_excel=0; %Set to 1 to print results to excel
    outputname='Output.xlsx'; %Set output name
same_guess=1; %Set to 1 if using same guess values

% Guess Values
p1=[120];%Parameter 1
p2=[80];%Parameter 2
p3=[1e8];%Parameter 3
lb=[0 0 1e5]; %Lower Bound of Parameters
ub=[500 500 1e10]; %Upper Bound of Parameters

% Beam Sizes
W.w0 = (6.9e-6)/2; %Pump spot size
W.w1 = (6.4e-6)/2; %Probe spot size

%__________________________________________________________________________
%% Behind the Scenes
% User selects files
[f,pathname]=uigetfile('*.*','Select Data File','MultiSelect','off');
%[fref,pathnameref]=uigetfile('*.*','Select Data File','MultiSelect','off');
% f = fullfile(pathname,filename);
% if ~iscell(filename) %Changes to cell if it's just a string (if 1 file)
%      filename = cellstr(f);
% end

cd(pathname)
figure(1)
    in=importdata(f,';'); %Imports data from txt
   %in2=importdata(fref,';'); %Imports data from txt
    freq=in.data(delete(1)+1:end-delete(2),1); %Frequency
    phi_data=in.data(delete(1)+1:end-delete(2),2); %Data Phase Lag
    phi_ref=in.data(delete(1)+1:end-delete(2),3); %Reference Phase Lag
    phi_exp=phi_data-phi_ref; % Corrected phase lag between ref and data
    
    % Correct if off by 180 degrees
    if phi_exp<-180
        phi_exp=phi_exp+180;
    elseif phi_exp>0
        phi_exp=phi_exp-180;
    end
    

    params=[p1,p2,p3];
    
    % Plot Data
    hold on
    scatter(freq,phi_exp);
    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel('\Phi (\circ)')
    xlim([min(freq) max(freq)])
    
    phi_fit = @(X,freq,W)FDTR(X,freq,W);
    SSEphifit = @(X)sum((phi_exp - phi_fit(X,freq,W)).^2);
    options = optimset('MaxFunEvals',1200,'MaxIter',1200);
    % Non Linear regression
    [extparams,SSE]=fminsearch(SSEphifit,params,options);
    fprintf('Parameters\nk=%f\nC=%e\nG=%e\n',extparams);

    % Calculate function using the regressed parameters
    [phi_model]=FDTR(extparams,freq,W);
    
    % Plot model
    hold on
    plot(freq,phi_model,'-r');
    legend('Measured', 'Curve fit')
    
    % Export formatted outputs to excel
    if print_to_excel==1
        output1=[freq phi_exp phi_model];
        header={'Freq','Exp','Model','','Guess'};
        xlswrite(outputname,filename,'A1');
        xlswrite(outputname,header,'A2')
        xlswrite(outputname,output1,'A3')
        xlswrite(outputname,params,'E3')
        xlswrite(outputname,{'Extracted Params'},'E4')
        xlswrite(outputname,extparams,'E5')
        xlswrite(outputname,ci','E6')
        xlswrite(outputname,{'Pump','Probe'},'E10')
        xlswrite(outputname,[w0 w1],'E11')
    end


function [Phi]=FDTR(X,freq,W)
%This function returns the inverse tangent of the ratio between the real
% and imaginary parts of the integrated Hankel transform.

% Layer Properties
S.kr = [220 0 X(1)]; %Thermal conductivity in the radial (r-)direction
S.kz = [220 X(3) X(1)]; %Thermal conductivity in the z-direction; G is included here
S.Cp = [2.48 0 1.8].*1e6;%Specific heat capacity of heach layer
S.t = [143e-9 0 1];%Thickness of each layer                                                                     
k = logspace(2,9,200);

omega = freq.*pi.*2;
Phi = ones([1 length(omega)]);
PhiExt=0;

for i = 1:length(omega)
......................................................................................................
    H = Hankel(k,W,S,omega(i));
    Hint = (1/(2.*pi)).*trapz(k,H);%

    Phi(i) = atand(imag(Hint)/real(Hint))-PhiExt;

end
Phi=Phi';
end

function [Hint] = Hankel(k,W,S,omega)
%This function returns the Hankel transform to the integrating function
%FDTR, above. This is a matrix solution that permits the use of an infinite
%number of layers in the multilayer material system.

Hint = ones([1 length(k)]);

for i = 1:length(k)

    MM = [1 0;0 1];
    
    for j = length(S.t):-1:1
       
        
        q(j) = ((S.kr(j).*k(i).^2+S.Cp(j).*1i.*omega)./(S.kz(j))).^(1/2);
        
        if S.t(j) ~= 0           
            N(:,:,j) = [ 1  -tanh(q(j)*S.t(j))/S.kz(j)/q(j);  -S.kz(j)*q(j)*tanh(q(j)*S.t(j))  1 ] ;
        else            
            N(:,:,j) = [ 1  -S.kz(j)^-1; 0 1] ;
        end
        
        MM = MM*N(:,:,j);

    end
    
    NDC = -MM(2,2)/MM(2,1);
    
    Hint(i) = k(i).*(NDC).*exp((-(k(i).^2).*(W.w0^2+W.w1^2))./8);

end

end


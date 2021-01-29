% eegplugin_runimatl() - EEGLAB plugin for interfacing the NSG portal
%
% Usage:
%   >> eegplugin_imat(fig);
%   >> eegplugin_imat(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] EEGLAB figure
%   trystrs    - [struct]  "try" strings for menu callbacks. See notes. 
%   catchstrs  - [struct]  "catch" strings for menu callbacks. See notes. 
%
% Author: Johanna Wagner, Ramon Martinez-Cancino SCCN, INC, UCSD, 2020
%
% See also: eeglab()

% Copyright (C) 2020 Ramon Martinez-Cancino
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_imat(fig, trystrs, catchstrs)
    
    vers = 'imat0.1';
    if nargin < 3
        error('eegplugin_imat requires 3 arguments');
    end;
    
  
    % -----------------------
    if ~exist('pop_runimat')
        p = fileparts(which('eegplugin_imat'));
        addpath(p);
    end
    
    % find tools menu (Singe subject)
    % ---------------
    menu = findobj(fig, 'tag', 'tools'); 

    
    % menu callback commands (Singe subject)
    % ----------------------
    catchstrs.store_and_hist = [ ...
                        '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); ' ...
                        'catch, errordlg2(lasterr, ''EEGLAB error''); LASTCOM= ''''; clear EEGTMP; end;' ...
                        'h(LASTCOM); disp(''Done.''); eeglab(''redraw'');' ];
    comrunima = [ trystrs.check_data '[EEG, IMA] = pop_runIMA(EEG);' catchstrs.add_to_hist ]; 
    
    cb_plot1 = 'pop_plotspecdecomp(EEG);';
    cb_plot2 = 'pop_plotspecenv(EEG);';
    cb_plot3 = 'pop_plotIMtimecourse(EEG);';
    
    % create menus (Singe subject)
    % ------------
    submenu = uimenu( menu, 'Label', 'Decompose IC spectrograms by IMA');
    set(menu, 'enable', 'on', 'userdata', 'startup:on;study:on');
    uimenu( submenu, 'Label', 'Run IMA', 'CallBack', comrunima);
    
    subsubmenu = uimenu( submenu, 'Label', 'Plot IMA results');
    uimenu( subsubmenu, 'Label', 'Superimposed components', 'CallBack', cb_plot1); % Joa check name of menu
    uimenu( subsubmenu, 'Label', 'Spectral envelope', 'CallBack', cb_plot2);       % Joa check name of menu
    uimenu( subsubmenu, 'Label', 'Time courses', 'CallBack', cb_plot3);            % Joa check name of menu
    
    % find tools menu (STUDY)
    % ---------------
    studymenu = findobj(fig, 'tag', 'study');
    
    cb_std_menu1 = '[STUDY] = pop_runIMA_study(STUDY, ALLEEG);';
    
    cb_std_menu21 = 'pop_plotspecdecomp_study(STUDY)';
    cb_std_menu22 = 'pop_plotspecenv_study(STUDY)';
    cb_std_menu23 = 'pop_plotIMtimecourse_study(STUDY)';
    
    cb_std_menu31 = 'pop_collecttemplates(STUDY)';
    cb_std_menu32 = '[STUDY] = pop_clusterIMAtemplates(STUDY, ALLEEG)';
    cb_std_menu33 = 'pop_plotIMAcluster(STUDY)';
    
      % menu callback commands (STUDY)
    % ----------------------
    submenustudy = uimenu( studymenu, 'Label', 'STUDY IMA', 'separator', 'on','userdata', 'startup:off;study:on');
    uimenu( submenustudy, 'Label', 'Run STUDY IMA', 'CallBack', cb_std_menu1, 'userdata', 'startup:off;study:on');
    
    std_subsubsubmenu1 = uimenu( submenustudy, 'Label', 'Plot IMA results','separator', 'on','userdata', 'startup:off;study:on');
    uimenu( std_subsubsubmenu1, 'Label', 'IM decomposition', 'CallBack', cb_std_menu21, 'userdata', 'startup:off;study:on');
    uimenu( std_subsubsubmenu1, 'Label', 'Spectral envelope', 'CallBack', cb_std_menu22, 'userdata', 'startup:off;study:on');
    uimenu( std_subsubsubmenu1, 'Label', 'IM timecourse',     'CallBack', cb_std_menu23, 'userdata', 'startup:off;study:on');
    
    std_subsubsubmenu2 = uimenu( submenustudy, 'Label', 'Cluster IMs','separator', 'on','userdata', 'startup:off;study:on');
    uimenu( std_subsubsubmenu2, 'Label', 'Collect templates', 'CallBack', cb_std_menu31, 'userdata', 'startup:off;study:on');
    uimenu( std_subsubsubmenu2, 'Label', 'Cluster IMs', 'CallBack', cb_std_menu32, 'userdata', 'startup:off;study:on');
    uimenu( std_subsubsubmenu2, 'Label', 'Plot clusters', 'CallBack', cb_std_menu33, 'userdata', 'startup:off;study:on');
    
    




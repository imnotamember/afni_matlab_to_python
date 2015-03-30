function [image_matrix_corrected,PHASES,REGRESSORS,OTHER] = ...
    retroicor_rvhr(image_matrix,slice_order,TR,cardiac_trig_times,resp,OPTIONS);
%function [image_matrix_corrected,PHASES,REGRESSORS,OTHER] = ...
%    retroicor_rvhr(image_matrix,slice_order,TR,cardiac_trig_times,resp,OPTIONS);
% 
% * catie chang,   catie.chang@nih.gov
%  
% * original: 12/13/11
% * modified: 02/14/12 from retroicor_main.m. This version
%   optionally includes RVHRcor regressors too! (RV*RRF, HR*CRF,
%   + time derivatives). 
%
%
% See the following for background:
% Glover et al., 2000: MRM 44, 162–167.
% Birn et al., 2006: Neuroimage 31, 1536-1548.
% Chang et al., 2009: Neuroimage 47, 1448-1459 (appendix A)
% Chang et al., 2009: Neuroimage 44, 857-869
%
% ---------------------------
% INPUTS:
% ---------------------------
% * image_matrix: 4D matlab matrix of fmri data, ** in native space and
%       before slice-timing correction ***
% * slice order:  vector indicating order of slice acquisition
%      (e.g. [30 28 26, .... 29 27 ... 1] for 30 slices with
%      "interleaved down" order)
% * TR: in seconds
% * cardiac_trig_times: vector of cardiac (R-wave peak) times, in
% seconds, relative to the start of the fMRI scan.
% * resp: structure with field:
%           * resp.trig: vector of respiratory peak times, in sec.
%           ---- OR instead with both of the following 2 fields ----
%           * resp.wave: respiration belt signal (aligned with the
%           fMRI scan)
%                 *AND*
%           * resp.dt:  sampling interval between the points in respiration
%             signal (in seconds, e.g. resp.dt=0.02 for 50 Hz sampling)
% * OPTIONS: optional structure with optional fields
%           * OPTIONS.retro_mode =1 for cardiac only, =2 for respiration
%           only, =3 for both, =0 for neither (i.e. no retroicor) [default = 3]
%           * OPTIONS.rvhr_mode: =1 for HR*CRF only, =2 for RV*RRF
%           only, =3 for both, =0 for neither (i.e. no rvhrcor) [default = 3]
%           * OPTIONS.slice_delta: constant temporal offset (in sec) to add
%           to all slice acquisition times. Usually useful only for
%           testing... [default = 0] 
%           * OPTIONS.doCorr:  =1 to correct the data, =0 if you
%           only want the regressors and phase outputs & not alter
%           the data. [default = 1]
%           * OPTIONS.verbose: =1 to display progress messages
%           [default = 1]
%
%  (** setting cardiac_trig_times = [] will ignore cardiac in both corrections)
%  (** setting resp = [] will ignore respiration in both corrections)
%
%
% ---------------------------
% OUTPUTS:
% ---------------------------
% * image_matrix_corrected: 4D matrix of corrected fmri data
% * PHASES: cell array of cardiac & respiration phases for each
%      slice. PHASES{i}(:,1) contains the cardiac phase for
%      slice "i", and PHASES{i}(:,2) contains the resp phases for
%      slice "i".
% * REGRESSORS: the retroicor & rvhrcor regressors for each
%      slice. dimensions = #timepoints x #regressors x #
%      slices. I.e., the regressors for slice "i" are the columns
%      of REGRESSORS(:,:,i). 
% * OTHER: possibly useful stuff, see code
%
%


% defaults
delta = 0;
doCorr = 1;
verbose = 1;
rvhr_mode = 3; %1=card, 2=resp, 3=both
retro_mode = 3;  %1=card, 2=resp, 3=both
Twin = 6; % default 6-sec window for computing RV & HR

CARD = 1;
RESP = 1;

CARD_RVHR = 1;
RESP_RVHR = 1;
CARD_RET = 1;
RESP_RET = 1;

% --------------------------------------------------------------
% load & check options
% --------------------------------------------------------------
if ~isempty(OPTIONS)
  if isfield(OPTIONS,'slice_delta')
    delta = OPTIONS.slice_delta
  end
  if isfield(OPTIONS,'doCorr')
    doCorr = OPTIONS.doCorr;
  end
  if isfield(OPTIONS,'verbose')
    verbose = OPTIONS.verbose;
  end
  if isfield(OPTIONS,'retro_mode')
    retro_mode = OPTIONS.retro_mode;
  end
  if isfield(OPTIONS,'rvhr_mode')
    rvhr_mode = OPTIONS.rvhr_mode;
  end
end

if isempty(cardiac_trig_times)
  CARD = 0;
end
if isempty(resp)
  RESP = 0;
end
if (CARD+RESP==0)
  error('need resp and/or cardiac input')
end
if (retro_mode+rvhr_mode==0)
  error('seriously?! :-)');
end
if (retro_mode == 0)
  CARD_RET = 0; RESP_RET = 0;
elseif (retro_mode == 1);
  CARD_RET = 1; RESP_RET = 0;
elseif (retro_mode == 2);
  CARD_RET = 0; RESP_RET = 1;
elseif (retro_mode == 3);
  CARD_RET = 1; RESP_RET = 1;
else
  error('unknown retro_mode');
end
if (rvhr_mode == 0)
  CARD_RVHR = 0; RESP_RVHR = 0;
elseif (rvhr_mode == 1);
  CARD_RVHR = 1; RESP_RVHR = 0;
elseif (rvhr_mode == 2);
  CARD_RVHR = 0; RESP_RVHR = 1;
elseif (rvhr_mode == 3);
  CARD_RVHR = 1; RESP_RVHR = 1;
else
  error('unknown rvhr_mode');
end
if (CARD==0)
  CARD_RET = 0; CARD_RVHR = 0;
end
if (RESP==0)
  RESP_RET = 0; RESP_RVHR = 0;
end

% --------------------------------------------------------------
% process inputs
% --------------------------------------------------------------
% get image parameters
fmri_dims = size(image_matrix);
nslc = fmri_dims(3);
nframes = fmri_dims(4);
npix_x = fmri_dims(1);
npix_y = fmri_dims(2);

% cardiac input
if (CARD)
  etrig = cardiac_trig_times;
end

% respiration input
if (RESP)
  if isfield(resp,'trig')
    RESP_TRIGGERS = 1;
    rtrig = resp.trig;
  elseif isfield(resp,'wave')
    if ~isfield(resp,'dt')
      error('please specify resp.dt');
    end
    RESP_TRIGGERS = 0;
  else
    error('incorrect resp format');
  end
  
  if ~RESP_TRIGGERS
    % shift 
    respwave = resp.wave-min(resp.wave);
    % bin respiration signal into 100 values
    [Hb,bins] = hist(respwave,100);
    % calculate derivative
    % first, filter respiratory signal - just in case
    f_cutoff = 1; % max allowable freq
    fs = 1/resp.dt;
    wn = f_cutoff/(fs/2);
    ntaps = 20;
    b = fir1(ntaps,wn);
    respfilt = filtfilt(b,1,respwave);
    drdt = diff(respfilt);
  end
end

% --------------------------------------------------------------
% find cardiac and respiratory phase vectors
% (not yet accounting for slice acquisition order - that's next)
% --------------------------------------------------------------
PHASES_0 = {};
for jj=1:nslc
  % times at which ith slice was acquired (midpoint):
  slice_times = [(TR/nslc)*(jj-0.5):TR:TR*nframes];
  slice_times = slice_times(1:nframes);
  % incorporate potential shift
  slice_times = slice_times + delta;
  
  phases_thisSlice = [];
  for ii=1:length(slice_times)
    
    % cardiac
    if (CARD)
      prev_trigs = find(etrig<slice_times(ii));
      if isempty(prev_trigs)
        t1 = 0;
      else
        t1 = etrig(prev_trigs(end));
      end
      next_trigs = find(etrig>slice_times(ii));
      if isempty(next_trigs)
        t2 = nframes*TR;
      else
        t2 = etrig(next_trigs(1));
      end
      phi_cardiac = (slice_times(ii) - t1)*2*pi/(t2-t1);
    else
      phi_cardiac = [];
    end
    
    % respiration: method based on triggers
    if (RESP)
      if (RESP_TRIGGERS)
        prev_trigs = find(rtrig<slice_times(ii));
        if isempty(prev_trigs)
          t1 = 0;
        else
          t1 = rtrig(prev_trigs(end));
        end
        next_trigs = find(rtrig>slice_times(ii));
        if isempty(next_trigs)
          t2 = nframes*TR;
        else
          t2 = rtrig(next_trigs(1));
        end
        phi_resp = (slice_times(ii) - t1)*2*pi/(t2-t1);
      else
        % respiration: method based on amplitude histogram
        tslice = slice_times(ii);
        iphys = max(1,round(tslice/resp.dt)); % closest idx in resp waveform
        iphys = min(iphys,length(drdt));
        amp = respwave(iphys);
        dbins = abs(amp-bins);
        [blah,thisBin] = min(dbins);  %closest resp histo bin
        numer = sum(Hb(1:thisBin));
        phi_resp = pi*sign(drdt(iphys))*(numer/length(respfilt));
      end
    else
      phi_resp = [];
    end
    
    % store
    phases_thisSlice(ii,:) = [phi_cardiac, phi_resp];
  end
    
  PHASES_0{jj} = phases_thisSlice;
end
% --------------------------------------------------------------
% permute according to slice acquisition order
% --------------------------------------------------------------
PHASES = PHASES_0(slice_order);
% --------------------------------------------------------------
% generate slice-specific retroicor regressors
% --------------------------------------------------------------
REGRESSORS_RET = [];
if (retro_mode>0)
  
  for jj=1:nslc    
    if (CARD_RET & RESP_RET)
      phi_c = PHASES{jj}(:,1);
      phi_r = PHASES{jj}(:,2);
    elseif (CARD_RET)
      phi_c = PHASES{jj}(:,1);
    elseif (RESP_RET)
      phi_r = PHASES{jj}(:,1);
    end
    
    covs = [];
    % Fourier expansion of cardiac phase
    if (CARD_RET)
      c1_c = cos(phi_c);
      s1_c = sin(phi_c);
      c2_c = cos(2*phi_c);
      s2_c = sin(2*phi_c);
      covs = [c1_c, s1_c, c2_c, s2_c];
    end
    % Fourier expansion of respiratory phase
    if (RESP_RET)
      c1_r = cos(phi_r);
      s1_r = sin(phi_r);
      c2_r = cos(2*phi_r);
      s2_r = sin(2*phi_r);
      covs = [covs, c1_r, s1_r, c2_r, s2_r];
    end
    REGRESSORS_RET(:,:,jj) = covs;  
  end % slices
  
end % retro_mode
% --------------------------------------------------------------
% generate slice-specific rvhrcor regressors
% --------------------------------------------------------------
REGRESSORS_RVHR0 = {};   
REGRESSORS_RVHR = cell(1,nslc);
if (rvhr_mode>0)
  
  for jj=1:nslc
    % times at which ith slice was acquired (midpoint):
    slice_times = [(TR/nslc)*(jj-0.5):TR:TR*nframes];
    slice_times = slice_times(1:nframes);
    % incorporate potential shift
    slice_times = slice_times + delta;
    X = [];
    % make slice RV*RRF regressor
    if (RESP_RVHR)
      if isfield(resp,'resp.trigs')
        display(['You supplied respiratory trigger times, but RV ', ...
                 'needs the continuous waveform ----> **Ignoring the RV ', ...
                 'terms**']);
      else
        nresp = length(respwave);
        for kk = 1:nframes
          t = slice_times(kk);
          i1 = max(0,floor((t - Twin*.5)/resp.dt)); 
          i2 = min(nresp, floor((t + Twin*.5)/resp.dt));
          if (i2<i1)
            error('respiration data is shorter than length of scan');
          end
          if i1==0;  i1 = i1+1; end
          rv(kk) = std(respwave(i1:i2));
        end
        rv = rv(:);
        % conv(rv, rrf)
        rv = rv-mean(rv);
        t = [0:TR:40-TR]; % 40-sec impulse response
        R = 0.6*(t.^2.1).*exp(-t/1.6) - 0.0023*(t.^3.54).*exp(-t/4.25); 
        R = R/max(R);
        rv_rrf = conv(rv,R);
        rv_rrf = rv_rrf(1:length(rv));
        % time derivative
        rv_rrf_d = diff(rv_rrf);
        rv_rrf_d = [rv_rrf_d(1); rv_rrf_d];
        X = [X, rv_rrf, rv_rrf_d];
        if (jj==1)
          rv_save = rv;
        end
      end
    end % RESP
    
    % make slice HR*CRF regressor
    if (CARD_RVHR)
      ntrig = length(etrig);
      for kk = 1:nframes
        t = slice_times(kk); 
        inds = intersect(find(etrig>=(t-Twin*.5)), ...
                         find(etrig<=(t+Twin*.5)));
        i1 = inds(1); i2 = inds(end);
        hr(kk) = (i2-i1)*60/(etrig(i2) - etrig(i1));  % bpm 
      end
      hr = hr(:);
      % conv(hr, crf)
      hr = hr - mean(hr);
      t = [0:TR:40-TR];  % 40-sec impulse response
      H = 0.6*(t.^2.7).*exp(-t/1.6) - 16*normpdf(t,12,3);
      H = H/max(H);
      hr_crf = conv(hr,H);
      hr_crf = hr_crf(1:length(hr));
      % time derivative
      hr_crf_d = diff(hr_crf);
      hr_crf_d = [hr_crf_d(1); hr_crf_d];
      X = [X, hr_crf, hr_crf_d];
      if (jj==1)
        hr_save = hr;
      end
    end
  
    REGRESSORS_RVHR0{jj} = X;
  end % slices 
      % --------------------------------------------------------------
      % permute according to slice acquisition order
      % --------------------------------------------------------------
  REGRESSORS_RVHR = REGRESSORS_RVHR0(slice_order);
end  % rvhr_mode

% --------------------------------------------------------------
% final set of physio regressors
% --------------------------------------------------------------
REGRESSORS = [];
for jj=1:nslc
  if ~isempty(REGRESSORS_RET)
    REGRESSORS(:,:,jj) = quaddetrend_cols([REGRESSORS_RET(:,:,jj), REGRESSORS_RVHR{jj}]);
  else
    REGRESSORS(:,:,jj) = quaddetrend_cols(REGRESSORS_RVHR{jj});
  end
end


% --------------------------------------------------------------
% correct the image data: slice-wise
% --------------------------------------------------------------
PCT_VAR_REDUCED = zeros(npix_x,npix_y,nslc);
if (doCorr)
  image_matrix_corrected = zeros(size(image_matrix));
  for jj=1:nslc
    if verbose; fprintf('%d ... ',jj); end
    slice_data = squeeze(image_matrix(:,:,jj,:));
    Y_slice = (reshape(slice_data,npix_x*npix_y,nframes))'; %ntime x nvox
    t = [1:nframes]';
    % design matrix
    XX = [t, t.^2, REGRESSORS(:,:,jj)]; 
    XX = [ones(size(XX,1),1), zscore(XX)];
    if verbose; fprintf('regressing ... '); end;
    Betas = pinv(XX)*Y_slice;
    Y_slice_corr = Y_slice - XX(:,4:end)*Betas(4:end,:);
    %%% -above line uses "4:end" to  keep mean and trends. Use 2:end to also remove
    % low-order trends / drift terms %%%
    
    % calculate percent variance reduction
    var_reduced = (var(Y_slice,0,1) - var(Y_slice_corr,0,1))./var(Y_slice,0,1);
    PCT_VAR_REDUCED(:,:,jj) = reshape(var_reduced',npix_x,npix_y);
    % fill corrected volume
    V_slice_corr = Y_slice_corr';
    if verbose; fprintf('storing ... \n'); end;
    for ii=1:nframes
      image_matrix_corrected(:,:,jj,ii) = reshape(V_slice_corr(:,ii),npix_x,npix_y); 
    end
  end
  fprintf('\n');
else
  % if you don't want to run the correction
  display('NOT applying correction to data ... (returning regressors only).')
  image_matrix_corrected = image_matrix;
end



% return other possibly useful stuff
OTHER.PCT_VAR_REDUCED = PCT_VAR_REDUCED;
OTHER.REGRESSORS_RET = REGRESSORS_RET; % final retroicor
                                       % regressors, each slice
OTHER.REGRESSORS_RVHR = REGRESSORS_RVHR; % final rvhrcor
                                       % regressors, each slice
if (RESP)
  OTHER.drdt = drdt;
end
if(rvhr_mode>0) 
  OTHER.rv_slice1 = rv_save;
  OTHER.hr_slice1 = hr_save;
end


function Yq = quaddetrend_cols(Y)
x = [1:size(Y,1)]';
Yq = Y;
for j=1:size(Y,2)
  y = Y(:,j);
  p = polyfit(x,y,2);
  ytrend = polyval(p, x);
  y = y - ytrend;
  Yq(:,j) = y;
end
'''
function [Opt, R, E] = RetroTS(SN)
'''
import sys
import numpy as np
'''
if (nargin < 1),
   fprintf(2,'Need some input.\n');
   return;
end

R = struct([]);
E = struct([]);
'''
if len(sys.argv) < 2:
    print 'Need some input.\n'
    quit()
#########################
'''
if (~isstruct(SN)), %mode 1, toy mode
   iscan = 12;
   lll = zglobb({ sprintf('Resp*%d*',iscan),...
                  sprintf('ECG*%d*',iscan),...
                  sprintf('scan*%d*', iscan)});
'''
if (type(sys.argv[1]) != 'numpy.ndarray':
    iscan = 12
    lll = zglobb(sprintf('Resp*%d*',iscan), sprintf('ECG*%d*',iscan), sprintf('scan*%d*', iscan)})#what is zglobb object? is this just loading defaults?
    ####lll = {'Resp*%d*'%iscan:zglobb[1], 'ECG*%d*'%iscan:zglobb[2], 'scan*%d*'%iscan:zglobb[3]}######Potentially this is a good workaround
'''
   %Get some info from header file and set params
   f = fopen(lll(3).name, 'r');
   s = fscanf(f,'%c');
   fclose(f);
   ns = length(s);
   pat = 'RT Physio:\W*sampling\W*';
   Opt.PhysFS = 1000/str2num(strtok(s(regexp(s,pat,'end'):ns)));
   Opt.Nslices = 20;
   Opt.VolTR = 2;
   Opt.SliceMajor = 1;
   Opt.ResampFS = Opt.PhysFS; %upsampling frequency in Hz
   Opt.Quiet = 1;
   Opt.ResamKernel = 'linear';   %resampling filter for envelopes and phase
   Opt.FIROrder = 40;  %order of fir filter
   Opt.RVTshifts = [0:5:20];  %shifts, in seconds, applied to RVT curve
   Opt.Demo = 0;
   Opt.zerophaseoffset = 0;
   Opt.fcutoff = 3; %cut off frequency for filter
   Opt.RespCutoffFreq = 3;
   Opt.CardCutoffFreq = 3;
   Opt.Respfile = lll(1).name;
   Opt.Cardfile = lll(1).name;
   Opt.SliceOffset = ...
      [0:Opt.VolTR./Opt.Nslices:Opt.VolTR-Opt.VolTR./Opt.Nslices];
   Opt.Prefix = sprintf('%d',iscan);
   Opt.SepDups = 0;
   clear ('s');
   clear ('SN');
'''
'''
    ####This is just dummy data for testing purposes, don't worry about this yet
    #Get some info from header file and set params
    f = open(lll[files for files in lll], 'r')
    ####Need a way to measure matrix size of f, then pass this along to
    ####np.chararray to create array before filling with f data
    s = np.chararray(f)
    close(f);####Not sure if close actually closes file f (look it up)
    ns = length(s);
    pat = 'RT Physio:\W*sampling\W*';
    Opt.PhysFS = 1000/str2num(strtok(s(regexp(s,pat,'end'):ns)));
    Opt.Nslices = 20;
    Opt.VolTR = 2;
    Opt.SliceMajor = 1;
    Opt.ResampFS = Opt.PhysFS; %upsampling frequency in Hz
    Opt.Quiet = 1;
    Opt.ResamKernel = 'linear';   %resampling filter for envelopes and phase
    Opt.FIROrder = 40;  %order of fir filter
    Opt.RVTshifts = [0:5:20];  %shifts, in seconds, applied to RVT curve
    Opt.Demo = 0;
    Opt.zerophaseoffset = 0;
    Opt.fcutoff = 3; %cut off frequency for filter
    Opt.RespCutoffFreq = 3;
    Opt.CardCutoffFreq = 3;
    Opt.Respfile = lll(1).name;
    Opt.Cardfile = lll(1).name;
    Opt.SliceOffset = ...
      [0:Opt.VolTR./Opt.Nslices:Opt.VolTR-Opt.VolTR./Opt.Nslices];
    Opt.Prefix = sprintf('%d',iscan);
    Opt.SepDups = 0;
    clear ('s');
    clear ('SN');
'''
'''
else,
   Opt = SN; clear ('SN');
'''
else:
    Opt = SN
    SN = ''
    Opt.err = 1
    Opt.zerophaseoffset = 0
    if ((!Opt['Respfile']) or Opt['Respfile'] is none):
        Opt['Respfile'] = ''
        Opt['Resp_out'] = 0
        Opt['RVT_out'] = 0
    if
'''
   Opt.err = 1; Opt.zerophaseoffset = 0;
   if ( (~isfield(Opt,'Respfile') | isempty(Opt.Respfile))),
      Opt.Respfile = '';
      Opt.Resp_out = 0;
      Opt.RVT_out = 0;
   end
   if ( (~isfield(Opt,'Cardfile') | isempty(Opt.Cardfile))),
      Opt.Cardfile = '';
      Opt.Card_out = 0;
   end
   if ( (~isfield(Opt,'Respfile') | isempty(Opt.Respfile)) & (~isfield(Opt,'Cardfile') | isempty(Opt.Cardfile))),
      fprintf(2,'No Respfile or Cardfile\n');
      return;
   end
   if ( ~isfield(Opt,'PhysFS') | isempty(Opt.PhysFS)),
      fprintf(2,'Missing field PhysFS\n');
      return;
   end
   if ( ~isfield(Opt,'SliceMajor') | isempty(Opt.SliceMajor)),
      Opt.SliceMajor = 1;
   end
   if ( ~isfield(Opt,'Nslices') | isempty(Opt.Nslices)),
      fprintf(2,'Missing field Nslices\n');
      return;
   end
   if ( ~isfield(Opt,'VolTR') | isempty(Opt.VolTR)),
      fprintf(2,'Missing field VolTR\n');
      return;
   end
   if ( ~isfield(Opt,'RVTshifts') | isempty(Opt.RVTshifts)),
      Opt.RVTshifts=[0:5:20];
   end
   if ( ~isfield(Opt,'ResampFS') | isempty(Opt.ResampFS)),
      Opt.ResampFS=Opt.PhysFS;
   end
   if ( ~isfield(Opt,'RespCutoffFreq') | isempty(Opt.RespCutoffFreq)),
      Opt.RespCutoffFreq=3;
   end
   if ( ~isfield(Opt,'CardCutoffFreq') | isempty(Opt.CardCutoffFreq)),
      Opt.CardCutoffFreq=3;
   end
   if ( ~isfield(Opt,'ResamKernel') | isempty(Opt.ResamKernel)),
      Opt.ResamKernel='linear';
   end

   if ( ~isfield(Opt,'FIROrder') | isempty(Opt.FIROrder)),
      Opt.FIROrder=40;
   end
   if ( ~isfield(Opt,'Quiet') | isempty(Opt.Quiet)),
      Opt.Quiet=1;
   end
   if ( ~isfield(Opt,'Demo') | isempty(Opt.Demo)),
      Opt.Demo=0;
   end
   if ( ~isfield(Opt,'Prefix') | isempty(Opt.Prefix)),
      Opt.Prefix = 'oba';
   end
   if ( ~isfield(Opt,'Resp_out') | isempty(Opt.Resp_out)),
      Opt.Resp_out = 1;
   end
   if ( ~isfield(Opt,'Card_out') | isempty(Opt.Card_out)),
      Opt.Card_out = 1;
   end
   if ( ~isfield(Opt,'RVT_out') | isempty(Opt.RVT_out)),
      Opt.RVT_out = 1;
   end
   if ( ~isfield(Opt,'SepDups') | isempty(Opt.SepDups)),
      Opt.SepDups = 0;
   end

   dtt = Opt.VolTR/Opt.Nslices; tt = 0.0;

   % & ~isfield(Opt, 'SliceOffset')
   % & (Opt.SliceOrder ~= 'alt+z')

      % default slice offset times are for alt+z (alternating slice timing)
   if ( ~isfield(Opt,'SliceOffset') | isempty(Opt.SliceOffset))
      Opt.SliceOffset=zeros(Opt.Nslices,1);
   end
   if(~isfield(Opt,'SliceOrder'))
      Opt.SliceOrder = 'alt+z'
   end

   if (isfield(Opt,'SliceOrder'))
      Opt.SliceOffset=zeros(Opt.Nslices,1);
      if(strcmpi(Opt.SliceOrder,'alt+z'))
         for (i=1:2:Opt.Nslices),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
         for (i=2:2:Opt.Nslices),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
      elseif(strcmpi(Opt.SliceOrder, 'alt+z2'))
         for (i=2:2:Opt.Nslices),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
         for (i=1:2:Opt.Nslices),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
      elseif(strcmpi(Opt.SliceOrder, 'seq+z'))
         for (i=1:1:Opt.Nslices),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
      elseif(strcmpi(Opt.SliceOrder,'seq-z'))
         for (i=Opt.Nslices:-1:1),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
      elseif(strcmpi(Opt.SliceOrder,'alt-z'))
         for (i=Opt.Nslices:-2:1),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
         for (i=Opt.Nslices-1:-2:1),
            Opt.SliceOffset(i) = tt; tt = tt+dtt;
         end
      elseif(strcmpi(Opt.SliceOrder,'Custom'))
          % timing already set in Opt.SliceOffset, do nothing
      else
         % read in time offsets from a file (SliceOrder is actually a
         % filename)
         readopt.verb = 0;
         [err, Opt.SliceOffset] = Read_1D(Opt.SliceOrder,readopt);
         if(length(Opt.SliceOffset)~=Opt.Nslices)
            fprintf('Could not read enough slice offsets from file');
            exit(1);
         end
      end
   end
   if(~Opt.Quiet)
      fprintf('Slice timing:'); Opt.SliceOffset
   end
   if ( ~isfield(Opt,'ShowGraphs') | isempty(Opt.ShowGraphs)),
      Opt.ShowGraphs = 1; % show graphs by default
   end
end
'''

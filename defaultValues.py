__author__ = 'Joshua Zosky'

import sys

if len(sys.argv) < 2:
    print 'Need some input.\n'
    quit()
R = {}
E = {}

#########################
'''
if (~isstruct(SN)), %mode 1, toy mode
   iscan = 12;
   lll = zglobb({ sprintf('Resp*%d*',iscan),...
                  sprintf('ECG*%d*',iscan),...
                  sprintf('scan*%d*', iscan)});
'''
if (type(sys.argv[1]) != 'dict'):
    iscan = 12
    #lll = zglobb(sprintf('Resp*%d*',iscan), sprintf('ECG*%d*',iscan), sprintf('scan*%d*', iscan)})#what is zglobb object? is this just loading defaults?
    ####lll = {'Resp*%d*'%iscan:zglobb[1], 'ECG*%d*'%iscan:zglobb[2], 'scan*%d*'%iscan:zglobb[3]}######Potentially this is a good workaround
'''
   %Get some info from header file and set params
   f = fopen(lll(3).name, 'r');
   s = fscanf(f,'%c');
   fclose(f);
'''
   ####Get some info from header file and set params
fileStuff = ''
with open('README.md', 'r') as f:
    s = f.read(1)
    while s != '':
        fileStuff = fileStuff + s
        s = f.read(1)
    print fileStuff
ns = len(fileStuff)
pat = 'RT Physio:\W*sampling\W*'#####What does this mean:" \W* "?
print ns
print pat
'''
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

'''
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
    if ((!Opt['Cardfile']) or Opt['Cardfile'] is none):
        Opt['Cardfile'] = ''
        Opt['Card_out'] = 0
    if ((!Opt['Respfile']) or Opt['Respfile'] is none) and ((!Opt['Cardfile']) or Opt['Cardfile'] is none):
        print 'No Respfile or Cardfile\n'
'''
'''
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
'''
__author__ = 'Joshua Zosky'
########################################################################################################################
# Respfile: Respiration data file
# Cardfile: Cardiac data file
# PhysFS: Physioliogical signal sampling frequency in Hz.
# Nslices: Number of slices
# VolTR: Volume TR in seconds
########################################################################################################################
import sys

retroTS_Dictionary = {}

########################################################################################################################
########################################################################################################################
def retroTS(Respfile,Cardfile,PhysFS, Nslices, VolTR,
            Prefix = 'Output_File_Name',
            SliceOffset = [],#figure this default out - needs to be a matrix of 0's based on Nslices
            RVTshifts = range(0, 20, 5),
            RespCutoffFreq = 3,
            CardCutoffFreq = 3,
            ResamKernel = 'linear',
            FIROrder = 40,
            Quiet = 1,
            Demo = 0,
            RVT_out = 0,
            Card_out = 0,
            Resp_out = 0,
            SliceOrder = 'alt+z',
            ShowGraphs = 1
            ):
    #########Determining SliceOffset based upon SliceOrder, VolTR, and Nslices.
    tt = 0.0 #Default float value to start iterations
    dtt = float(VolTR)/float(Nslices) #Increments for iteration
    SliceOffset = [0] * Nslices #Initial value for SliceOffset
    SliceFileList = [] #List for using external file for SliceOffset values/Indicates if using external file in last loop
    if SliceOrder[0:3] == 'alt': #Alternating?
        for i in range(0, Nslices, 2):
            if i < (Nslices-1):
                for ii in range(i,i+2): SliceOffset[ii] = tt
            else:
                SliceOffset[i] = tt
            tt += dtt
    elif SliceOrder[0:3] == 'seq': #Sequential?
        for i in range(0, Nslices):
            SliceOffset[i] = tt
            tt += dtt
    elif SliceOrder == 'Custom': #Does nothing, unsure of it's purpose
        pass
    else: #Open external file specified in argument line, fill SliceFileList with values, then load into SliceOffset
        with open(SliceOrder, 'r') as f:
            for i in f.readlines():
                SliceFileList.append(int(i))
                if len(SliceFileList) != Nslices:
                    print 'Could not read enough slice offsets from file'
                    print 'File should have as many offsets as Nslices'
                    quit()
            SliceOffset = SliceFileList
    if SliceOrder[3] == '-' and SliceFileList == []: #Check for a minus to indicate a reversed offset list
        SliceOffset.reverse()
    if Quiet != 1: #Show the slice timing (P.S. Printing is very time consuming in python)
        print 'Slice timing:', Opt.SliceOffset

    #########Create option copy for each type of signal.
    OptR = {'fcutoff': RespCutoffFreq,
            'AmpPhase': 1,      #amplitude based phase for respiration
            'as_percover': 50,  #percent overlap of windows for fft
            'as_windwidth': 0,  #window width in seconds for fft, 0 for full window
            'as_fftwin': 0      #1 == hamming window. 0 == no windowing
            }
    OptE = {'fcutoff': CardCutoffFreq,
            'AmpPhase': 0       #time based phase for cardiac signal
            }
    print OptR
    print OptE

########################################################################################################################
def peakFinder():
# function [R, e] = PeakFinder(vvec, Opt)
# %Example: PeakFinder('Resp*.1D');
# % or PeakFinder(v) where v is a column vector
# % if v is a matrix, each column is processed separately.
# %
# %clear all but vvec (useful if I run this function as as script)
    # keep('vvec', 'Opt');
    # e = 0;
    # R = struct([]);
    #

# if (nargin < 2) Opt = struct(); end
# if (~isfield(Opt,'PhysFS')  | isempty(Opt.PhysFS)),
#    Opt.PhysFS= 1/0.025; %sampling frequency
# end
# if (~isfield(Opt,'zerophaseoffset') | isempty(Opt.zerophaseoffset) ),
#    Opt.zerophaseoffset = 0.5;  %Fraction of the period that corresponds
#                            %to a phase of 0
#                        %0.5 means the middle of the period, 0 means the 1st peak
# end
# if (~isfield(Opt,'Quiet') | isempty(Opt.Quiet)),
#    Opt.Quiet = 0;
# end
# if (~isfield(Opt,'ResampFS') | isempty(Opt.ResampFS)),
#    Opt.ResampFS = Opt.PhysFS;
# end
# if (~isfield(Opt,'fcutoff') | isempty(Opt.fcutoff)),
#    Opt.fcutoff = 10;
# end
# if (~isfield(Opt,'FIROrder') | isempty(Opt.fcutoff)),%BC ???
#    Opt.FIROrder = 80;
# end
# if (~isfield(Opt,'ResamKernel') | isempty(Opt.ResamKernel)),
#    Opt.ResamKernel = 'linear';
# end
# if (~isfield(Opt,'Demo') | isempty(Opt.Demo)),
#    Opt.Demo = 0;
# end
# if (~isfield(Opt,'as_windwidth') | isempty(Opt.as_windwidth)),
#    Opt.as_windwidth = 0;
# end
# if (~isfield(Opt,'as_percover') | isempty(Opt.as_percover)),
#    Opt.as_percover = 0;
# end
# if (~isfield(Opt,'as_fftwin') | isempty(Opt.as_fftwin)),
#    Opt.as_fftwin = 0;
# end
# if (~isfield(Opt,'SepDups') | isempty(Opt.SepDups)),
#    Opt.SepDups = 0;
# end
#
#
# if (Opt.Demo),
#    Opt.Quiet = 0;
# else
#    pause off
# end
#
# %some filtering
# fnyq = Opt.PhysFS./2;
# %w(1) = 0.1/fnyq;  %cut frequencies below 0.1Hz
# %w(2) = Opt.fcutoff/fnyq;    % upper cut off frequency normalized
# %b = fir1(Opt.FIROrder, w, 'bandpass');     %FIR filter of order 40
# w = Opt.fcutoff/fnyq;    % upper cut off frequency normalized
# b = fir1(Opt.FIROrder, w, 'low');     %FIR filter of order 40
#
# NoDups = 1; % remove duplicates that might come up when improving peak location
#
# if (ischar(vvec)),
#    l = zglobb(vvec);
#    nl = length(l);
#    if (isnumeric(l)),
#       fprintf(2,'File (%s) not found\n', vvec);
#       e = 1;
#    return;
#    end
# else
#    l = [];
#    nl = size(vvec,2);
#    if (nl < 1),
#       fprintf(2,'No vectors\n', nl);
#       e = 1;
#       return;
#    end
# end
#
# clear R; %must clear it. Or next line fails
# R(nl) = struct( 'vname', '',...
#             't', [], ...
#             'X', [],...
#             'iz', [],...   %zero crossing (peak) locations
#             'ptrace', [], 'tptrace', [],...
#             'ntrace', [], 'tntrace', [],...
#             'prd', [], 'tmidprd', [], 'ptracemidprd', [],...
#             'phz', [],...
#             'RV', [], 'RVT', [] ...
#              );
#
# for (icol = 1:1:nl),
#
#    if (~isempty(l) && ~l(icol).isdir),
#       R(icol).vname = sprintf('%s%s', l(icol).path, l(icol).name);
#       v = Read_1D(R(icol).vname);
#    else,
#       R(icol).vname = sprintf('vector input col %d', icol);
#       v = vvec(:,icol);
#    end
#
#    windwidth = 0.2; %window for adjusting peak location in seconds
#
#
#    %remove the mean
#    v = (v - mean(v));
#    R(icol).v = v;      %store it for debugging
#
#    %filter both ways to cancel phase shift
#    v = filter(b,1,v); v = flipud(v); v = filter(b,1,v); v = flipud(v);
#
#    %get the analytic signal
#    R(icol).X = analytic_signal(v, Opt.as_windwidth.*Opt.PhysFS,...
#                                Opt.as_percover, Opt.as_fftwin);
#          %using local version to illustrate, can use hilbert
#          %Doing ffts over smaller windows can improve peak detection
#          %in the few instances that go undetected but what value to use
#          %is not clear and there seems to be at times more errors introduced
#          %in the lower envelope .
#
#    nt = length(R(icol).X);
#    R(icol).t = [0:1/Opt.PhysFS:(nt-1)/Opt.PhysFS]; % FIX FIX FIX
#    iz = find( imag(R(icol).X(1:nt-1)).*imag(R(icol).X(2:nt)) <= 0);
#    polall = -sign(imag(R(icol).X(1:nt-1)) - imag(R(icol).X(2:nt)));
#
#    pk = real(R(icol).X(iz));
#    pol = polall(iz);
#    tiz = R(icol).t(iz);
#
#
#    ppp = find(pol>0);
#    ptrace = pk(ppp);
#    tptrace = tiz(ppp);
#    ppp = find(pol<0);
#    ntrace = pk(ppp);
#    tntrace = tiz(ppp);
#    if (~Opt.Quiet),
#       fprintf(2,[ '--> Load signal\n',...
#                   '--> Smooth signal\n',...
#                   '--> Calculate analytic signal Z\n',...
#                   '--> Find zero crossing of imag(Z)\n',...
#                   '\n']);
#
#       figure(1); clf
#       subplot(211);
#       plot (R(icol).t, real(R(icol).X),'g'); hold on
#       %plot (R(icol).t, imag(R(icol).X),'g');
#       plot (tptrace, ptrace, 'ro');
#       plot (tntrace, ntrace, 'bo');
#       %plot (R(icol).t, abs(R(icol).X),'k');
#
#       subplot (413);
#       vn = real(R(icol).X)./(abs(R(icol).X)+eps);
#       plot (R(icol).t, vn, 'g'); hold on
#       ppp = find(pol>0);
#       plot (tiz(ppp), vn(iz(ppp)), 'ro');
#       ppp = find(pol<0);
#       plot (tiz(ppp), vn(iz(ppp)), 'bo');
#
#          drawnow ;
#          if (Opt.Demo),
#             uiwait(msgbox('Press button to resume', 'Pausing', 'modal'));
#          end
#    end
#
#
#
#    %Some polishing
#    if (1),
#       nww = ceil(windwidth/2 * Opt.PhysFS);
#       pkp = pk;
#       R(icol).iz = iz;
#       for (i=1:1:length(iz)),
#          n0 = max(2,iz(i)-nww);
#          n1 = min(nt,iz(i)+nww);
#          if (pol(i) > 0),
#             [xx, ixx] = max((real(R(icol).X(n0:n1))));
#          else,
#             [xx, ixx] = min((real(R(icol).X(n0:n1))));
#          end
#          R(icol).iz(i) = n0+ixx-2;
#          pkp(i) = xx;
#       end
#       tizp = R(icol).t(R(icol).iz);
#
#       ppp = find(pol>0);
#       R(icol).ptrace = pkp(ppp);
#       R(icol).tptrace = tizp(ppp);
#       ppp = find(pol<0);
#       R(icol).ntrace = pkp(ppp);
#       R(icol).tntrace = tizp(ppp);
#
#       if (NoDups),
#       %remove duplicates
#          if (Opt.SepDups),
#             fprintf(2,'YOU SHOULD NOT BE USING THIS.\n');
#             fprintf(2,' left here for the record\n');
#             [R(icol).tptrace, R(icol).ptrace] = ...
#                         remove_duplicates(R(icol).tptrace, R(icol).ptrace, Opt);
#             [R(icol).tntrace, R(icol).ntrace] = ...
#                         remove_duplicates(R(icol).tntrace, R(icol).ntrace, Opt);
#          else,
#             [R(icol).tptrace, R(icol).ptrace,...
#              R(icol).tntrace, R(icol).ntrace] = ...
#                         remove_PNduplicates(R(icol).tptrace, R(icol).ptrace,...
#                         R(icol).tntrace, R(icol).ntrace, Opt);
#          end
#          if (length(R(icol).ptrace) ~= length(R(icol).ntrace)),
#             fprintf(1,'Bad news in tennis shoes. I''m outa here.\n');
#             e = 1;
#             return;
#          end
#       end
#
#       if (~Opt.Quiet),
#          fprintf(2,[ '--> Improved peak location\n',...
#                      '--> Removed duplicates \n',...
#                      '\n']);
#          subplot(211);
#          plot( R(icol).tptrace, R(icol).ptrace,'r+',...
#                R(icol).tptrace, R(icol).ptrace,'r');
#          plot( R(icol).tntrace, R(icol).ntrace,'b+',...
#                R(icol).tntrace, R(icol).ntrace,'b');
#          drawnow ;
#          if (Opt.Demo),
#             uiwait(msgbox('Press button to resume', 'Pausing', 'modal'));
#          end
#       end
#    else
#       tizp = tiz;
#       R(icol).iz = iz;
#       pkp = pk;
#       R(icol).ptrace = ptrace;
#       nR(icol).ptrace = nptrace;
#    end
#
#
#    %Calculate the period
#    nptrc = length(R(icol).tptrace);
#    R(icol).prd = (R(icol).tptrace(2:nptrc) - R(icol).tptrace(1:nptrc-1) );
#    R(icol).ptracemidprd = (   R(icol).ptrace(2:nptrc) ...
#                             + R(icol).ptrace(1:nptrc-1) ) ./2.0;
#    R(icol).tmidprd = (  R(icol).tptrace(2:nptrc) ...
#                       + R(icol).tptrace(1:nptrc-1)) ./2.0;
#    if (~Opt.Quiet),
#          fprintf(2,[ '--> Calculated the period (from beat to beat)\n',...
#                      '\n']);
#       plot (R(icol).tmidprd, R(icol).ptracemidprd,'kx');
#       for (i=1:1:length(R(icol).prd)),
#        text( R(icol).tmidprd(i), R(icol).ptracemidprd(i),...
#              sprintf('%.2f', R(icol).prd(i)));
#       end
#          drawnow ;
#          if (Opt.Demo),
#             uiwait(msgbox('Press button to resume', 'Pausing', 'modal'));
#          end
#    end
#
#    if (~isempty(Opt.ResamKernel)),
#       %interpolate to slice sampling time grid:
#       R(icol).tR = [0:1./Opt.ResampFS:max(R(icol).t)];
#       R(icol).ptraceR = interp1( R(icol).tptrace', R(icol).ptrace, ...
#                                  R(icol).tR,Opt.ResamKernel);
#       R(icol).ntraceR = interp1( R(icol).tntrace', R(icol).ntrace, ...
#                                  R(icol).tR,Opt.ResamKernel);
#       R(icol).prdR = interp1(R(icol).tmidprd, R(icol).prd, ...
#                              R(icol).tR,Opt.ResamKernel);
#       %you get NaN when tR exceeds original signal time, so set those
#       %to the last interpolated value
#       R(icol).ptraceR = clean_resamp(R(icol).ptraceR);
#       R(icol).ntraceR = clean_resamp(R(icol).ntraceR);
#       R(icol).prdR = clean_resamp(R(icol).prdR);
#    end
#
#  if (icol ~= nl), input ('Hit enter to proceed...','s'); end
#
# end
#    if (~Opt.Quiet),   plotsign2(1); end
#
# return;
########################################################################################################################

for item in sys.argv[1:7]:
    try:
        retroTS_Dictionary[(item.split('='))[0]] = int(item.split('=')[1])
    except:
        retroTS_Dictionary[(item.split('='))[0]] = item.split('=')[1]
for item in retroTS_Dictionary:
    print "%s:%s" % (item,retroTS_Dictionary[item])

'''Prove there's data in the data file
with open(retroTS_Dictionary['Respfile'], 'r') as f:
    for i in f.readlines():
        print i
with open(retroTS_Dictionary['Cardfile'], 'r') as f:
    for i in f.readlines():
        print i
'''
########################################################################################################################
retroTS(**retroTS_Dictionary)
########################################################################################################################
# Prefix: Prefix of output file
# SliceOffset: Vector of slice acquisition time offsets in seconds.
#       (default is equivalent of alt+z)
# RVTshifts: Vector of shifts in seconds of RVT signal.
#       (default is [0:5:20])
# RespCutoffFreq: Cut off frequency in Hz for respiratory lowpass filter
#       (default 3 Hz)
# CardCutoffFreq: Cut off frequency in Hz for cardiac lowpass filter
#       (default 3 Hz)
# ResamKernel: Resampling kernel.
#       (default is 'linear', see help interp1 for more options)
# FIROrder: Order of FIR filter.
#       (default is 40)
# Quiet: [1]/0  flag.
#       (defaut is 1) Show talkative progress as the program runs
# Demo: [1]/0 flag.
#       (default is 0)
# RVT_out: [1]/0 flag for writing RVT regressors
# Card_out: [1]/0 flag for writing Card regressors
# Resp_out: [1]/0 flag for writing Resp regressors
# SliceOrder:['alt+z']/'alt-z'/'seq+z'/'seq-z'/'Custom'/filename.1D
#              Slice timing information in seconds. The default is
#              alt+z. See 3dTshift help for more info. 'Custom' allows
#              the program to use the values stored in the
#              Opt.SliceOffset array. If a value is placed into the
#              SliceOrder field other than these, it is assumed to be
#              the name of a 1D / text file containing the times for
#              each slice (also in seconds).



quit()
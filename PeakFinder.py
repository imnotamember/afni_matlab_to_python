__author__ = 'Joshua Zosky'

"""
import scipy.flipud as flipud
import scipy.mean as mean
import scipy.zeros as zeros
"""
from scipy import *
from scipy.fftpack import fft
from scipy.signal import lfilter
from scipy.signal import firwin
import numpy
"""
function [r, e] = PeakFinder(var_vector, Opt)
In python, the above command looks like:
import PeakFinder as PF
(r, e) = PF.peak_finder(var_vector, **opt)
where opt is a dictionary with the parameters for the function
"""

def peak_finder(var_vector,
                phys_fs=(1 / 0.025),
                zero_phase_offset=0.5,
                quiet=0,
                resample_fs=(1 / 0.025),
                f_cutoff=10,
                fir_order=80,
                resample_kernel='linear',
                demo=0,
                as_window_width=0,
                as_percover=0,
                as_fftwin=0,
                sep_dups=0,
                ):
    """
    Example: PeakFinder('Resp*.1D')
    or PeakFinder(v) where v is a column vector
    if v is a matrix, each column is processed separately.
    :param var_vector: column vector--list of list(s)
    :param phys_fs: Sampling frequency
    :param zero_phase_offset: Fraction of the period that corresponds to a phase of 0
                                0.5 means the middle of the period, 0 means the 1st peak
    :param quiet:
    :param resample_fs:
    :param f_cutoff:
    :param fir_order: BC ???
    :param resample_kernel:
    :param demo:
    :param as_window_width:
    :param as_percover:
    :param fftwin:
    :param sep_dups:
    :return: [r, e] r = Peak of var_vector; e = error value
    """
    # #clear all but var_vector (useful if I run this function as as script)
    # keep('var_vector', 'opt')
    default_div = 1 / 0.025
    if (phys_fs != default_div) and (resample_fs == default_div):
        resample_fs = phys_fs
    if demo:
        quiet = 0
    else:
        pause = False  # pause off
    e = False  # default value for e
    r = {}
    # Some filtering
    filter_nyquist = phys_fs / 2
    # w[1] = 0.1/filter_nyquist                     # Cut frequencies below 0.1Hz
    # w[2] = f_cutoff/filter_nyquist                # Upper cut off frequency normalized
    # b = signal.firwin(fir_order, w, 'bandpass')   # FIR filter of order 40
    w = f_cutoff / filter_nyquist                   # Upper cut off frequency normalized
    b = firwin(fir_order, w, 'low')          # FIR filter of order 40
    no_dups = 1                                     # Remove duplicates that might come up when improving peak location
    if isinstance(var_vector, str):
        L = zglobb(var_vector)  # NEED TO CONVERT ZGLOBB INTO LIST MAKER OF FILE OBJECTS; I.E. type(L) == list
        nl = len(L)
        if isinstance(L, (int, long, float, complex)):
            print 'Error: File (%s) not found\n', var_vector
            e = True
            return e
    else:
        L = []
        nl = len(var_vector)
        if nl < 1:
            print 'Error: No vectors\n', nl
            e = True
            return e

    # del(r) # "Must clear it. Or next line fails" -- Probably unnecessary in Python
    r_list = []
    for i in range(nl):
        r = {'v_name': '',
             't': [],
             'x': [],
             'iz': [],   # zero crossing (peak) locations
             'p_trace': [],
             'tp_trace': [],
             'n_trace': [],
             'tn_trace': [],
             'prd': [],
             't_mid_prd': [],
             'p_trace_mid_prd': [],
             'phz': [],
             'RV': [],
             'RVT': []
             }
        r_list.append(r)
    for i_column in range(nl):
        if L and not os.path.isdir(L):
            r_list[i_column]['v_name'] = '%s%s' % (sys.path, L[i_column]['name'])
            with open(nl[i_column]['v_name'], "rb") as f:
                v = f.read()
                print v
        else:
            r_list[i_column]['v_name'] = 'vector input col %d' % i_column
            v = var_vector[i_column]

        window_width = 0.2  # Window for adjusting peak location in seconds
        # Remove the mean
        v = [i - mean(v) for i in v]
        r_list[i_column]['v'] = v  # Store it for debugging
        # Filter both ways to cancel phase shift
        v = lfilter(b, 1, v)
        v = flipud(v)
        v = lfilter(b, 1, v)
        v = flipud(v)
        # Get the analytic signal
        r_list[i_column]['x'] = analytic_signal(v, as_window_width*phys_fs, as_percover, as_fftwin)  # # RESOLVE THIS: as_window_width.*phys_fs
        # LEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERELEFT OFF HERE
        # Using local version to illustrate, can use hilbert
        # Doing ffts over smaller windows can improve peak detection in the few instances that go undetected but
        # what value to use is not clear and there seems to be at times more errors introduced in the lower envelope.
        nt = len(r_list[i_column]['x'])
        r_list[i_column]['t'] = range(0,(nt - 1) / phys_fs, 1 / phys_fs)  # # % FIX FIX FIX
        iz = find(imag(r_list[i_column]['x'](1:nt - 1)).*imag(r_list[i_column]['x'](2:nt)) <= 0)
        polall = -sign(imag(r_list[i_column]['x'](1:nt-1)) - imag(r_list[i_column]['x'](2:nt)))

        pk = real(r_list[i_column]['x'](iz))
        pol = polall(iz)
        tiz = r_list[i_column].t(iz)

        '''
        ppp = find(pol>0)
        ptrace = pk(ppp)
        tptrace = tiz(ppp)
        ppp = find(pol<0)
        ntrace = pk(ppp)
        tntrace = tiz(ppp)
        if (~Opt.Quiet),
          fprintf(2,[ '--> Load signal\n',...
                      '--> Smooth signal\n',...
                      '--> Calculate analytic signal Z\n',...
                      '--> Find zero crossing of imag(Z)\n',...
                      '\n'])

          figure(1) clf
          subplot(211)
          plot (r_list[i_column].t, real(r_list[i_column]['x']),'g') hold on
          %plot (r_list[i_column].t, imag(r_list[i_column]['x']),'g')
          plot (tptrace, ptrace, 'ro')
          plot (tntrace, ntrace, 'bo')
          %plot (r_list[i_column].t, abs(r_list[i_column]['x']),'k')

          subplot (413)
          vn = real(r_list[i_column]['x'])./(abs(r_list[i_column]['x'])+eps)
          plot (r_list[i_column].t, vn, 'g') hold on
          ppp = find(pol>0)
          plot (tiz(ppp), vn(iz(ppp)), 'ro')
          ppp = find(pol<0)
          plot (tiz(ppp), vn(iz(ppp)), 'bo')

             drawnow 
             if (Opt.Demo),
                uiwait(msgbox('Press button to resume', 'Pausing', 'modal'))
             end
       end



       %Some polishing
       if (1),
          nww = ceil(window_width/2 * Opt.PhysFS)
          pkp = pk
          r_list[i_column].iz = iz
          for (i=1:1:length(iz)),
             n0 = max(2,iz(i)-nww)
             n1 = min(nt,iz(i)+nww)
             if (pol(i) > 0),
                [xx, ixx] = max((real(r_list[i_column]['x'](n0:n1))))
             else,
                [xx, ixx] = min((real(r_list[i_column]['x'](n0:n1))))
             end
             r_list[i_column].iz(i) = n0+ixx-2
             pkp(i) = xx
          end
          tizp = r_list[i_column].t(r_list[i_column].iz)

          ppp = find(pol>0)
          r_list[i_column].ptrace = pkp(ppp)
          r_list[i_column].tptrace = tizp(ppp)
          ppp = find(pol<0)
          r_list[i_column].ntrace = pkp(ppp)
          r_list[i_column].tntrace = tizp(ppp)

          if (NoDups),
          %remove duplicates
             if (Opt.SepDups),
                fprintf(2,'YOU SHOULD NOT BE USING THIS.\n')
                fprintf(2,' left here for the record\n')
                [r_list[i_column].tptrace, r_list[i_column].ptrace] = ...
                            remove_duplicates(r_list[i_column].tptrace, r_list[i_column].ptrace, Opt)
                [r_list[i_column].tntrace, r_list[i_column].ntrace] = ...
                            remove_duplicates(r_list[i_column].tntrace, r_list[i_column].ntrace, Opt)
             else,
                [r_list[i_column].tptrace, r_list[i_column].ptrace,...
                 r_list[i_column].tntrace, r_list[i_column].ntrace] = ...
                            remove_PNduplicates(r_list[i_column].tptrace, r_list[i_column].ptrace,...
                            r_list[i_column].tntrace, r_list[i_column].ntrace, Opt)
             end
             if (length(r_list[i_column].ptrace) ~= length(r_list[i_column].ntrace)),
                fprintf(1,'Bad news in tennis shoes. I''m outa here.\n')
                e = True
                return
             end
          end

          if (~Opt.Quiet),
             fprintf(2,[ '--> Improved peak location\n',...
                         '--> Removed duplicates \n',...
                         '\n'])
             subplot(211)
             plot( r_list[i_column].tptrace, r_list[i_column].ptrace,'r+',...
                   r_list[i_column].tptrace, r_list[i_column].ptrace,'r')
             plot( r_list[i_column].tntrace, r_list[i_column].ntrace,'b+',...
                   r_list[i_column].tntrace, r_list[i_column].ntrace,'b')
             drawnow 
             if (Opt.Demo),
                uiwait(msgbox('Press button to resume', 'Pausing', 'modal'))
             end
          end
       else
          tizp = tiz
          r_list[i_column].iz = iz
          pkp = pk
          r_list[i_column].ptrace = ptrace
          nR(i_column).ptrace = nptrace
       end


       %Calculate the period
       nptrc = length(r_list[i_column].tptrace)
       r_list[i_column].prd = (r_list[i_column].tptrace(2:nptrc) - r_list[i_column].tptrace(1:nptrc-1) )
       r_list[i_column].ptracemidprd = (   r_list[i_column].ptrace(2:nptrc) ...
                                + r_list[i_column].ptrace(1:nptrc-1) ) ./2.0
       r_list[i_column].tmidprd = (  r_list[i_column].tptrace(2:nptrc) ...
                          + r_list[i_column].tptrace(1:nptrc-1)) ./2.0
       if (~Opt.Quiet),
             fprintf(2,[ '--> Calculated the period (from beat to beat)\n',...
                         '\n'])
          plot (r_list[i_column].tmidprd, r_list[i_column].ptracemidprd,'kx')
          for (i=1:1:length(r_list[i_column].prd)),
           text( r_list[i_column].tmidprd(i), r_list[i_column].ptracemidprd(i),...
                 sprintf('%.2f', r_list[i_column].prd(i)))
          end
             drawnow 
             if (Opt.Demo),
                uiwait(msgbox('Press button to resume', 'Pausing', 'modal'))
             end
       end

       if (~isempty(Opt.ResamKernel)),
          %interpolate to slice sampling time grid:
          r_list[i_column].tR = [0:1./Opt.ResampFS:max(r_list[i_column].t)]
          r_list[i_column].ptraceR = interp1( r_list[i_column].tptrace', r_list[i_column].ptrace, ...
                                     r_list[i_column].tR,Opt.ResamKernel)
          r_list[i_column].ntraceR = interp1( r_list[i_column].tntrace', r_list[i_column].ntrace, ...
                                     r_list[i_column].tR,Opt.ResamKernel)
          r_list[i_column].prdR = interp1(r_list[i_column].tmidprd, r_list[i_column].prd, ...
                                 r_list[i_column].tR,Opt.ResamKernel)
          %you get NaN when tR exceeds original signal time, so set those
          %to the last interpolated value
          r_list[i_column].ptraceR = clean_resamp(r_list[i_column].ptraceR)
          r_list[i_column].ntraceR = clean_resamp(r_list[i_column].ntraceR)
          r_list[i_column].prdR = clean_resamp(r_list[i_column].prdR)
       end

     if (i_column ~= nl), input ('Hit enter to proceed...','s') end

    end
       if (~Opt.Quiet),   plotsign2(1) end

    return
    '''
    
# function h = analytic_signal(vi, window_width, percover, win),


def analytic_signal(vi, window_width, percover, win):
    n_vi = len(vi)
    vi_size = (n_vi, vi[0].len)
    h = zeros(vi_size)
    (bli, ble, num) = fft_segments(window_width, percover, n_vi)
    '''
    for ii in range(bli.len):
        v = vi(bli(ii): ble(ii))
        nv = length(v)
        if win == 1:
            fv = fft(v.*hamming(nv))  # # Guesswork
        else:
            fv = fft(v)
        window = zeros(size(v))
        # Zero negative frequencies, double positive frequencies
        if not nv & 1:
            window([1 nv/2+1]) = 1  # keep DC
            window([2:nv/2]) = 2    # double pos. freq
        else:
            window([1]) = 1
            window([2:(nv+1)/2]) = 2
        h(bli(ii):ble(ii)) = h(bli(ii):ble(ii)) + ifft(fv.*window)
    h = h./num
    '''
    return h

# function [bli, ble, num] = fftsegs (ww, po, nv)


def fft_segments(ww, po, nv):
    """
    Returns the segements that are to be used for fft calculations.
    Example: (bli, ble, num) = fftsegs (100, 70, 1000);
    :param ww: Segment width (in number of samples)
    :param po: Percent segment overlap
    :param nv: Total number of samples in original symbol
    :return: (bli, ble, num): bli, ble: Two Nblck x 1 vectors defining the segments' starting and ending indices;
                                num: An nv x 1 vector containing the number of segments each sample belongs to
    """
    return_tuple = (ww, po, nv)
    return return_tuple
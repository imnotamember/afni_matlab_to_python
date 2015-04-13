__author__ = 'Joshua Zosky'


# function [r] = PhaseEstimator(r, Opt)

def phase_estimator(r, opt):
    if opt['amp_phase'] != None:
        for i_column in len(r):
            # Calculate the phase of the trace, with the peak 
            # to be the start of the phase
            np_trace = len(r(i_column).tptrace)
            r(i_column).phz=-2.*ones(size(r(i_column).t))
            i = 1
            j = 1
            while i <= (np_trace - 1):
                while r(i_column).t(j) < r(i_column).tptrace(i+1):
                    if (r(i_column).t(j) >= r(i_column).tptrace(i)):
                        # Note: Using a constant244 period for each interval causes slope discontinuity within a period.
                        # One should resample prd(i) so that it is estimated at each time in r(i_column).t(j)
                        # dunno if that makes much of a difference in the end however.
                        r(i_column).phz(j) = (r(i_column).t(j) - (r(i_column).tptrace(i))) ./ r(i_column).prd(i) + Opt.zerophaseoffset
                        # ./  Right array divide. A./B denotes element-by-element division.  A and B
                        # must have the same dimensions unless one is a scalar.
                        # A scalar can be divided with anything.
                        if r(i_column).phz(j) < 0:
                            r(i_column).phz(j) = -r(i_column).phz(j)
                        if r(i_column).phz(j) > 1:
                            r(i_column).phz(j) = r(i_column).phz(j)-1
                    j += 1
                i += 1
        # Remove the points flagged as unset
        r(i_column).phz(find(r(i_column).phz<-1)) = 0.0
        # Change phase to radians
        r(i_column).phz = r(i_column).phz.*2.*pi
    else:  # Phase based on amplitude
        for i_column in r:
            max_amplitude = max(r(i_column).ptrace) #  Scale to the max
            glover_r = zscale(r(i_column).v, max_amplitude, 0) #  scale, per Glover 2000's paper
            bins = range(1, 100)./100.*max_amplitude #  .* = Element-wise multiplication
            [hb,bbins] = hist(glover_r, bins)
    
        
    #  #################################################################################################################
    # else, %phase based on amplitude 
    #    for (i_column=1:1:length(r)),
    #       % at first scale to the max
    #       max_amplitude = max(r(i_column).ptrace);
    #       glover_r = zscale(r(i_column).v, max_amplitude, 0); %scale, per Glover 2000's paper
    #       bins = [1:1:100]./100.*max_amplitude;
    #       [hb,bbins] = hist(glover_r, bins);
    #       if(Opt.ShowGraphs)
    #           bar (bins, hb);
    #       end
    #       %find the polarity of each time point in v
    #       i = 1; itp = 1; inp = 1;
    #       while (  i <= length(r(i_column).v) & ...
    #                r(i_column).t(i) < r(i_column).tptrace(1) & ...
    #                r(i_column).t(i) < r(i_column).tntrace(1) ),
    #          r(i_column).phz_pol(i) = 0;
    #          i = i + 1;
    #       end
    #       if (r(i_column).tptrace(1) < r(i_column).tntrace(1)), 
    #          cpol=-1;    %expiring phase, peak behind us
    #          itp = 2;
    #       else 
    #          cpol = 1; %inspiring phase (bottom behind us)
    #          inp = 2;
    #       end
    #       r(i_column).phz_pol = zeros(size(r(i_column).v));
    #       %add a fake point to tptrace and tntrace to avoid ugly if statements
    #       r(i_column).tptrace = [r(i_column).tptrace r(i_column).t(end)];
    #       r(i_column).tntrace = [r(i_column).tntrace r(i_column).t(end)];
    #       while(i <= length(r(i_column).v)),
    #          r(i_column).phz_pol(i) = cpol;
    #          if (r(i_column).t(i) == r(i_column).tptrace(itp)),
    #             cpol = -1; itp = min([itp+1, length(r(i_column).tptrace)]);
    #          elseif (r(i_column).t(i) == r(i_column).tntrace(inp)),
    #             cpol = +1; inp = min([inp+1, length(r(i_column).tntrace)]);
    #          end
    #          %cpol, inp, itp, i, r
    #          i = i + 1;
    #       end
    #       r(i_column).tptrace = [r(i_column).tptrace(1:end-1)];
    #       r(i_column).tntrace = [r(i_column).tntrace(1:end-1)];
    #       if(Opt.ShowGraphs),
    #           clf;
    #           plot (r(i_column).t, glover_r,'b'); hold on
    #           ip = find(r(i_column).phz_pol>0);
    #           plot (r(i_column).t(ip), 0.55.*max_amplitude,'r.');
    #           in = find(r(i_column).phz_pol<0);
    #           plot (r(i_column).t(in),0.45.*max_amplitude,'g.');
    #       end
    #       %Now that we have the polarity, without computing sign(dR/dt) 
    #       % as in Glover et al 2000, calculate the phase per eq. 3 of that paper
    #       %first the sum in the numerator
    #       glover_r = round(glover_r/max_amplitude.*100)+1; glover_r(find(glover_r>100))=100; 
    #       shb = sum(hb);
    #       hbsum = zeros(1,100);
    #       hbsum(1)=hb(1)./shb;
    #       for (i=2:1:100),
    #          hbsum(i) = hbsum(i-1)+hb(i)./shb;
    #       end
    #       for(i=1:1:length(r(i_column).t)),
    #          r(i_column).phz(i) = pi.*hbsum(round(glover_r(i))).*r(i_column).phz_pol(i); 
    #       end
    #    end
    # end
    # ##################################################################################################################
    
    for i_column in len(r):
        r(i_column).tst = range(0, max(r(i_column) -  0.5*opt['vol_tr'], opt['vol_tr']))  # Time-Series-Time vector
        r(i_column).phz_slc = zeros(len(r(i_column).tst), opt.n_slices)
        r(i_column).phz_slc_reg = zeros(len(r(i_column).tst), 4, opt.n_slices)
        for isl in range(1, opt.n_slices + 1):
            tslc = r(i_column).tst + opt.slice_offset(isl)
            for i in range(1, len(r(i_column).tst)):
                [mi,imin] = min(abs(tslc(i)-r(i_column).t))
                r(i_column).phz_slc(i,isl) = r(i_column).phz(imin)
            # And make regressors for each slice
            r(i_column).phz_slc_reg(:,1, isl) = sin(r(i_column).phz_slc(:,isl))
            r(i_column).phz_slc_reg(:,2, isl) = cos(r(i_column).phz_slc(:,isl))
            r(i_column).phz_slc_reg(:,3, isl) = sin(2.*r(i_column).phz_slc(:,isl))
            r(i_column).phz_slc_reg(:,4, isl) = cos(2.*r(i_column).phz_slc(:,isl))
        if (opt['quiet'] == 0) and (opt['show_graphs'] == 1):
            print '--> Calculated phase\n\n'
            subplot (413)
            plot (r(i_column).t, r(i_column).phz./2./pi, 'm')
            if hasattr(r(i_column),'phzR'):
                plot (r(i_column).tR, r(i_column).phzR./2./pi, 'm-.')
            subplot (414)
            plot (r(i_column).tst, r(i_column).phz_slc(:,1), 'ro',
                r(i_column).tst, r(i_column).phz_slc(:,2), 'bo',
                r(i_column).tst, r(i_column).phz_slc(:,2), 'b-')
            hold on
            plot (r(i_column).t, r(i_column).phz, 'k')
            grid on
            # title it
            title (r(i_column).vname, 'Interpreter', 'None')
                   drawnow
            if opt['demo'] == 1:
                uiwait(msgbox('Press button to resume', 'Pausing', 'modal'))

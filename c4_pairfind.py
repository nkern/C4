# pair_find()

import numpy as np

def pair_find(rdata,vdata):
        # Set Binning Constants
        rlim = 3                # Mpc
        vlim = 4000             # km/s
        vstep = 1000            # km/s
        nbin = vlim/vstep*2

        # Do binning and find number of points in each bin
        bins_temp = np.zeros(nbin)
        for i in range(nbin):
                bins_temp[i] = np.where((rdata<rlim)&(vdata<(vlim-i*vstep))&(vdata>(vlim-(i+1)*vstep)))[0].size

        # Make bin gridding denser by factor of 2, interpolate for higher resolution
        bins = np.zeros(nbin*2+1)
        j = 0
        for i in range(len(bins)):
                if i == 0 or i == range(len(bins))[-1]:
                        bins[i] = 0
                elif i % 2 != 0:
                        bins[i] = bins_temp[j]
                else:
                        bins[i] = np.mean([bins_temp[j],bins_temp[j+1]])
                        j += 1

        # Create 3 double peaked models and one single peaked model
        v_range = np.linspace(vlim,-vlim,nbin*2+1)
        peak = np.max(bins)
        double1 = np.zeros(nbin*2+1)            # clusters at +/- 1500 km/s
        double1[nbin/2:nbin/2+3] = peak
        double1[-nbin/2-3:-nbin/2] = peak
        double2 = np.zeros(nbin*2+1)            # clusters at 0km/s and -3000 km/s
        double2[nbin-1:nbin+2] = peak
        double2[-4:-1] = peak
        double3 = np.zeros(nbin*2+1)            # clusters at 0 km/s and 3000 km/s
        double3[nbin-1:nbin+2] = peak
        double3[1:4] = peak
        single = np.zeros(nbin*2+1)
        single[nbin/2+3:nbin/2+6] = peak
        single[[nbin/2+2,nbin/2+6]] = peak/2

        # Do Simple Least Squares Fit
        d1_chi = sum((bins-double1)**2)
        d2_chi = sum((bins-double2)**2)
        d3_chi = sum((bins-double3)**2)
        s_chi = sum((bins-single)**2)

        # If any double model chi_sq is < s_chi*thresh, it is likely a pair
        thresh = 2
        if d1_chi < s_chi*thresh or d2_chi < s_chi*thresh or d3_chi < s_chi*thresh:
                pair = True
        else:
                pair = False

        return pair, d1_chi, d2_chi, d3_chi, s_chi, double1, double2, double3, single, v_range, bins









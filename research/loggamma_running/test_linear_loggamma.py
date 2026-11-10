"""Deterministic tests for the isolated log-Gamma folding prototype."""
import math
import random
import unittest

from linear_loggamma import Geometry, fit_pooled, log_peak_height


class LogGammaFoldingTests(unittest.TestCase):
    def test_shifted_peaks_have_shared_shape_and_closed_area(self):
        lam,shape=0.8,14.0625
        shifts=(-0.35,-0.18,0.04,0.27,0.41)
        for sign in (1,-1):
            entries=[]
            for i,shift in enumerate(shifts):
                g=Geometry.uniform(41,0.04,lam,sign)
                height=math.log(800+100*i)
                y=[log_peak_height(z,shift,height,shape,lam,sign) for z in g.z]
                entries.append((g,g.moments(y),3.5*i))
            fit=fit_pooled(entries)
            self.assertAlmostEqual(fit.shared_slope,sign*shape*lam,places=7)
            self.assertLess(fit.residual_sum_squares,2e-7)
            for i,pk in enumerate(fit.peaks):
                self.assertAlmostEqual(pk.apex_local,shifts[i],places=7)
                true_area=(800+100*i)*math.exp(shape)*math.gamma(shape)/(lam*shape**shape)
                self.assertAlmostEqual(math.exp(pk.log_area_local)/true_area,1.,places=8)

    def test_cached_moments_rolling_both_orientations(self):
        rng=random.Random(42)
        full=[rng.uniform(-1,1) for _ in range(350)]
        n=33
        for sign in (1,-1):
            g=Geometry.uniform(n,0.035,0.8,sign)
            if sign==1:
                starts=range(0,len(full)-n+1)
                first=0
            else:
                starts=range(len(full)-n,-1,-1)
                first=len(full)-n
            running=g.moments(full[first:first+n])
            for j in starts:
                if j!=first:
                    if sign==1:
                        running=g.slide_right(running,full[j-1],full[j+n-1])
                    else:
                        running=g.slide_left(running,full[j+n],full[j])
                direct=g.moments(full[j:j+n])
                for attr in ('y0','yz','ye','yy'):
                    self.assertAlmostEqual(getattr(running,attr),getattr(direct,attr),delta=1e-9)

    def test_variable_length_profiles(self):
        entries=[]
        for n,phase in ((15,0.11),(21,-0.2),(37,0.07)):
            g=Geometry.uniform(n,0.07,0.6,1)
            entries.append((g,g.moments([log_peak_height(z,phase,math.log(1000),10,0.6)
                                         for z in g.z]),0.0))
        fit=fit_pooled(entries)
        self.assertAlmostEqual(fit.shared_slope,6.0,places=8)
        self.assertLess(fit.residual_sum_squares,1e-7)

    def test_nonuniform_geometry(self):
        lam=.5
        offsets=[-.75,-.62,-.41,-.22,-.08,0.07,0.25,0.39,0.57,0.8]
        g=Geometry.from_offsets(offsets,lam,1)
        y=[log_peak_height(z,0.12,math.log(1000),18,lam) for z in g.z]
        f=fit_pooled([(g,g.moments(y),1.0)])
        self.assertAlmostEqual(f.peaks[0].apex_global,1.12,places=7)
        self.assertLess(f.residual_sum_squares,1e-7)
        with self.assertRaises(ValueError):
            g.slide_right(g.moments(y),1.,2.)

    def test_reject_degenerate_or_nonintegrable(self):
        with self.assertRaises(ValueError):
            Geometry.from_offsets([0.,0.,0.,0.],.8,1)
        g=Geometry.uniform(11,.2,.8,1)
        ys=[-z for z in g.z]
        with self.assertRaises(ValueError):
            fit_pooled([(g,g.moments(ys),0.)])


if __name__=='__main__':
    unittest.main()

import math

class lognormal:
  def __init__(self, n_tot, meanr, gstdv):
    assert len(n_tot) == len(meanr) == len(gstdv)
    self.meanr = meanr
    self.stdev = gstdv
    self.n_tot = n_tot

  def __call__(self, lnr):
    total = 0.
    for m in range(0, len(self.n_tot)):
      total += self.n_tot[m] * math.exp(
          -pow((lnr - math.log(self.meanr[m])), 2) / 2 / pow(math.log(self.stdev[m]),2)
      ) / math.log(self.stdev[m]) / math.sqrt(2*math.pi);
    return total


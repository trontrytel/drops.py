from libcloudphxx.common import R_d, c_pd, g, p_1000

# perfect gas for for dry air                    
def T_fun(p_d, th_d):
  return th_d * pow(p_d / p_1000, R_d / c_pd)         

def rhod_fun(p_d, th_d):
  return p_d / R_d / T_fun(p_d, th_d)


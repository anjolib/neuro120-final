def g_h(t):
  ghrelin_h = 4.11635458e+02 * np.cos(2.65248261e-01 * t - 1.94183153e+00) + 1.31190578e+03
  return ghrelin_h 
def g_p(t):
  w_real = 7.0077939 / 13.829261660647743 #x.std()
  phi_real = 0.62327368 - (7.0077939 * 23.018136335209505 / 13.829261660647743)

  ghrelin_p = 0.06468781 * np.cos(w_real * t + phi_real) + 1.03423138
  return ghrelin_p
def l(t):
  leptin = 1.7132457 * np.cos(0.25891851 * t + 1.90885731) + 3.68477418
  return leptin 
def s(t):
  c_fos = 4.97669924e+02 * np.cos(2.44009422e-01 * t - 9.08555848e-01) + 8.69592999e+02
  return c_fos

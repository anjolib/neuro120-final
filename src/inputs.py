import numpy as np

def ghrelin(t: float, source: str = "h", normalize: bool = True) -> float:
    """
    Returns ghrelin levels at time t,
    using either hypothalamus levels (source='h') or plasma levels (source='p').

    Arguments:
        t (float): current time in hours
        source (str): source of ghrelin levels
        normalize (bool): if True, normalize function

    Returns:
        float: current ghrelin level
    """
    match source:
        case "h":
            return ghrelin_hypothalamus(t, normalize)
        case "p":
            return ghrelin_plasma(t, normalize)
        case _:
            raise Exception("Source must be 'h', 'hnorm', or 'p'.")
            return None

def leptin(t: float, normalize: bool = True) -> float:
    """
    Returns leptin levels at time t.

    Arguments:
        t (float): current time in hours
        normalize (bool): if True, normalize function

    Returns:
        float: leptin level at time t
    """
    if normalize:
        return np.cos(0.25891851 * t + 1.90885731)
    else:
        return (1.7132457 * np.cos(0.25891851 * t + 1.90885731) + 3.68477418)

def scn(t: float, normalize: bool = True) -> float:
    """
    Returns SCN output level at time t.

    Arguments:
        t (float): current time in hours
        normalize (bool): if True, normalize function

    Returns:
        float: leptin level at time t
    """
    if normalize:
        return np.cos(2.44009422e-01 * t - 9.08555848e-01)
    else:
        return (4.97669924e+02 * np.cos(2.44009422e-01 * t - 9.08555848e-01) + 8.69592999e+02)
        
def ghrelin_hypothalamus(t: float, normalize: bool = True):
    if normalize:
        return np.cos(2.65248261e-01 * t - 1.94183153e+00)
    else:
        return (4.11635458e+02 * np.cos(2.65248261e-01 * t - 1.94183153e+00) + 1.31190578e+03)
    
def ghrelin_plasma(t: float, normalize: bool = True):
    w_real = 7.0077939 / 13.829261660647743 #x.std()
    phi_real = 0.62327368 - (7.0077939 * 23.018136335209505 / 13.829261660647743)
    if normalize:
        return np.cos(w_real * t + phi_real)
    else:
        return (0.06468781 * np.cos(w_real * t + phi_real) + 1.03423138)

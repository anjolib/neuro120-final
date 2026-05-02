from . import params as p

_GHRELIN_ON = [
    (p.gen(0),   p.gen(3)),
    (p.gen(6),   p.gen(27)),
    (p.gen(30),  p.gen(51)),
    (p.gen(54),  p.gen(75)),
    (p.gen(78),  p.gen(99)),
    (p.gen(102), p.gen(123)),
    (p.gen(126), p.gen(147)),
    (p.gen(150), p.gen(171)),
]

def _ghrelin_production_rate(t):
    """
    Switches between maximum ghrelin production rate depending on time of day.

    Arguments:
        t (float): time in seconds

    Returns:
        float: production rate p.beta or p.beta2 corresponding to time
    """
    for t0, t1 in _GHRELIN_ON:
        if t0 <= t < t1:
            return p.beta
    return p.beta2

def dGI1_dt(GI1, GI2):
    return -GI1 * p.k12 / (1.0 + p.sigma * GI2)


def dGI2_dt(GI1, GI2):
    return p.k12 / (1.0 + p.sigma * GI2) * GI1 - p.k2x * GI2


def dGH_dt(GH, GI2, t):
    return (1.0 / (1.0 + p.gamma * GI2)) * _ghrelin_production_rate(t) - p.khx * GH

def dIN_dt(GI2, IN):
    return (p.rho
        + p.theta * GI2**p.s / (p.delta**p.s + GI2**p.s)
        - p.kix * IN)

def dLP_dt(LP, lam):
    return (
        p.FM * (p.R * p.eta) / (1.0 + (lam / p.EC50LEP)**p.l)
        - p.epsilon * LP
    ) / p.BV


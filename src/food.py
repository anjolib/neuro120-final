from .params import gen, break_, lun, din


def TRE(t: float, r: float) -> float:
    windows = [
        (gen(8),  break_(8)),   (gen(13), lun(13)),   (gen(17), din(17)),
        (gen(32), break_(32)),  (gen(37), lun(37)),   (gen(41), din(41)),
        (gen(56), break_(56)),  (gen(61), lun(61)),   (gen(65), din(65)),
    ]
    for t0, t1 in windows:
        if t0 < t < t1:
            return r
    return 0.0

def CD(t: float, r: float, s: float) -> float:
    windows = [
        (gen(8),    break_(8),  r),  (gen(13),   lun(13),   r),
        (gen(16),   din(16),    s),  (gen(19),   din(19),   r),
        (gen(20.5), din(20.5),  s),
        (gen(32),   break_(32), r),  (gen(37),   lun(37),   r),
        (gen(40),   din(40),    s),  (gen(43),   din(43),   r),
        (gen(44.5), din(44.5),  s),
        (gen(56),   break_(56), r),  (gen(61),   lun(61),   r),
        (gen(64),   din(64),    s),  (gen(67),   din(67),   r),
        (gen(68.5), din(68.5),  s),
    ]
    for t0, t1, rate in windows:
        if t0 < t < t1:
            return rate
    return 0.0

def IFI(t: float, r: float) -> float:
    windows = [
        (gen(13), lun(13)),  (gen(17), din(17)),
        (gen(37), lun(37)),  (gen(41), din(41)),
        (gen(61), lun(61)),  (gen(65), din(65)),
    ]
    for t0, t1 in windows:
        if t0 < t < t1:
            return r
    return 0.0

SCHEDULES = {
    "TRE": TRE,
    "CD":  CD,
    "IFI": IFI,
}


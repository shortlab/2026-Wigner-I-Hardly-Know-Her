from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import medfilt, savgol_filter


def _safe_std(v: np.ndarray, default: float = 1e-3) -> float:
    s = float(np.nanstd(v))
    if not np.isfinite(s) or s < 1e-12:
        return default
    return s


def find_lefts_rights_mc(
    t: np.ndarray,
    q: np.ndarray,
    n_samples: int = 1000,
    seed: int | None = None,
    peak_window: tuple[float, float] | None = None,
) -> dict | None:
    t = np.asarray(t, dtype=float).ravel()
    q = np.asarray(q, dtype=float).ravel()
    valid = np.isfinite(t) & np.isfinite(q)
    t, q = t[valid], q[valid]
    if len(t) < 50:
        return None

    order = np.argsort(t)
    t, q = t[order], q[order]
    keep = np.concatenate(([True], np.diff(t) != 0))
    t, q = t[keep], q[keep]

    if peak_window is not None:
        w = (t >= peak_window[0]) & (t <= peak_window[1])
        if np.count_nonzero(w) >= 50:
            t, q = t[w], q[w]

    frame = max(5, int(np.ceil(len(t) * 0.10)))
    if frame % 2 == 0:
        frame += 1
    if frame >= len(t):
        frame = max(5, len(t) - (1 - len(t) % 2))

    q_smooth = savgol_filter(q, window_length=frame, polyorder=3)
    dqdt = np.diff(q_smooth) / np.diff(t)
    dqdt = medfilt(dqdt, kernel_size=5)
    t_mid = (t[:-1] + t[1:]) / 2.0

    i_max = int(np.argmax(dqdt))
    i_min = int(np.argmin(dqdt))
    left_end = min(i_max, i_min)
    right_start = max(i_max, i_min)

    d_left = dqdt[: max(left_end, 1)]
    d_right = dqdt[right_start:] if right_start < len(dqdt) else dqdt[-1:]
    mu_left, sd_left = float(np.mean(d_left)), _safe_std(d_left)
    mu_right, sd_right = float(np.mean(d_right)), _safe_std(d_right, sd_left)

    lrloc = None
    for m in range(1, 21):
        idx = np.where(np.abs(dqdt[: max(left_end, 1)] - mu_left) <= sd_left * m)[0]
        if len(idx):
            lrloc = int(idx[-1])
            break

    rlloc = None
    for m in range(1, 21):
        idx = np.where(np.abs(dqdt[right_start:] - mu_right) <= sd_right * m)[0]
        if len(idx):
            rlloc = int(idx[0] + right_start)
            break

    if lrloc is None:
        return None

    lr_center = float(t_mid[lrloc])
    rl_center = float(t_mid[rlloc]) if rlloc is not None else lr_center + 20.0
    if lr_center > rl_center:
        lr_center, rl_center = rl_center, lr_center

    rng = np.random.default_rng(seed)
    ll = rng.uniform(lr_center - 17, lr_center - 13, n_samples)
    lr = rng.uniform(lr_center - 2, lr_center + 2, n_samples)
    rl = rng.uniform(rl_center - 2, rl_center + 2, n_samples)
    rr = rng.uniform(rl_center + 13, rl_center + 17, n_samples)
    return {
        "ll": ll,
        "lr": lr,
        "rl": rl,
        "rr": rr,
        "dTmax": float(t_mid[i_max]),
        "dTmin": float(t_mid[i_min]),
        "lrcenter": lr_center,
        "rlcenter": rl_center,
    }


def spline_integral_orig(
    t_data: np.ndarray,
    q_data: np.ndarray,
    t_ll: np.ndarray | float,
    t_lr: np.ndarray | float,
    t_rl: np.ndarray | float,
    t_rr: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = np.asarray(t_data, dtype=float).ravel()
    q = np.asarray(q_data, dtype=float).ravel()
    valid = np.isfinite(t) & np.isfinite(q)
    t, q = t[valid], q[valid]

    ll = np.asarray(t_ll, dtype=float).ravel()
    lr = np.asarray(t_lr, dtype=float).ravel()
    rl = np.asarray(t_rl, dtype=float).ravel()
    rr = np.asarray(t_rr, dtype=float).ravel()
    n = max(len(ll), len(lr), len(rl), len(rr))
    if len(ll) == 1 and n > 1:
        ll = np.repeat(ll, n)
    if len(lr) == 1 and n > 1:
        lr = np.repeat(lr, n)
    if len(rl) == 1 and n > 1:
        rl = np.repeat(rl, n)
    if len(rr) == 1 and n > 1:
        rr = np.repeat(rr, n)

    hcrys = np.full(n, np.nan, dtype=float)
    tout = np.full((len(t), n), np.nan, dtype=float)
    qout = np.full((len(t), n), np.nan, dtype=float)

    for i in range(n):
        left = (t > ll[i]) & (t < lr[i])
        right = (t > rl[i]) & (t < rr[i])
        corr = (t > ll[i]) & (t < rr[i])
        if np.count_nonzero(corr) < 3:
            continue

        tc = t[corr]
        qc = q[corr]

        if np.count_nonzero(left) < 2 or np.count_nonzero(right) < 2:
            baseline = np.interp(tc, [tc[0], tc[-1]], [qc[0], qc[-1]])
        else:
            t_left, q_left = t[left], q[left]
            t_right, q_right = t[right], q[right]
            p_left = np.polyfit(t_left, q_left - 1.0, 1)
            p_right = np.polyfit(t_right, q_right - 1.0, 1)
            t_fit = np.concatenate([t_right, t_left])
            q_fit = np.concatenate([
                np.polyval(p_right, t_right),
                np.polyval(p_left, t_left),
            ])
            uniq_t, uniq_idx = np.unique(t_fit, return_index=True)
            uniq_q = q_fit[uniq_idx]
            if len(uniq_t) < 3:
                baseline = np.interp(tc, [tc[0], tc[-1]], [qc[0], qc[-1]])
            else:
                f = interp1d(uniq_t, uniq_q, kind="cubic", fill_value="extrapolate")
                baseline = f(tc)

        q_corr = np.abs(qc - baseline)
        hcrys[i] = np.trapezoid(q_corr, tc) if len(tc) > 1 else np.nan
        tout[: len(tc), i] = tc
        qout[: len(tc), i] = q_corr

    return hcrys, tout, qout


def hook_mc_qad(
    t: np.ndarray,
    q: np.ndarray,
    fig_on: bool = False,
    n: int = 1000,
    mode: str | None = None,
    seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray, list, np.ndarray, np.ndarray]:
    del fig_on
    t = np.asarray(t, dtype=float).ravel()
    q = np.asarray(q, dtype=float).ravel()
    valid = np.isfinite(t) & np.isfinite(q)
    t, q = t[valid], q[valid]
    order = np.argsort(t)
    t, q = t[order], q[order]
    if len(t) < 10:
        return np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), ["heat", np.nan, np.array([np.nan, np.nan])], np.array([]), np.array([])

    if mode is None:
        d_avg = np.mean(q)
        mode = "cool" if abs(q[0] - d_avg) < abs(q[-1] - d_avg) else "heat"

    if mode == "cool":
        t = np.flipud(t)
        q = np.flipud(q)

    rng = np.random.default_rng(seed)
    n_sp = rng.random(n)
    n_sw = rng.random(n)
    n_ip = rng.random(n)
    n_iw = rng.random(n)

    dt = (t[1:] + t[:-1]) / 2.0
    dq = (q[1:] - q[:-1]) / (t[1:] - t[:-1])

    if mode == "heat":
        noise_mask = dt > t[0] + 10
    else:
        noise_mask = dt < t[-1] - 10

    if np.count_nonzero(noise_mask) < 5:
        noise_mask = np.ones_like(dt, dtype=bool)

    p_noise = np.polyfit(dt[noise_mask], dq[noise_mask], 2)
    dq_noise = np.polyval(p_noise, dt[noise_mask])
    rmse = float(np.sqrt(np.mean((dq[noise_mask] - dq_noise) ** 2)))

    infl = int(np.argmin(dq))
    t_infl = float(dt[infl])
    diff = dq - np.polyval(p_noise, dt)
    if mode == "heat":
        settles = np.where((diff > 0) & (dt > t_infl))[0]
        if len(settles) == 0:
            settles = np.array([min(infl + 1, len(t) - 2)])
        t_sp_min = t[settles[0]]
        t_sp_max = t_sp_min + 5.0
    else:
        settles = np.where((diff > 0) & (dt < t_infl))[0]
        if len(settles) == 0:
            settles = np.array([max(infl - 1, 0)])
        t_sp_max = t[settles[-1]]
        t_sp_min = t_sp_max - 5.0

    v_sp = n_sp * (t_sp_max - t_sp_min) + t_sp_min
    ind_sp = np.array([int(np.argmin(np.abs(t - v))) for v in v_sp], dtype=int)
    t_sp = t[ind_sp]

    w_sw = n_sw * (10.0 - 5.0) + 5.0
    if mode == "heat":
        ind_sw = np.array([int(np.argmin(np.abs(t - (t_sp[i] + w_sw[i])))) for i in range(n)], dtype=int)
    else:
        ind_sw = np.array([int(np.argmin(np.abs(t - (t_sp[i] - w_sw[i])))) for i in range(n)], dtype=int)

    infls = dq < (np.min(dq) + 2 * rmse)
    if np.count_nonzero(infls) == 0:
        s_ip = 1.0
    else:
        s_ip = float(np.max(np.abs(dt[infls] - t_infl)) / 4.0)
        if s_ip <= 0:
            s_ip = 1.0

    t_ip_min, t_ip_max = t_infl - 2 * s_ip, t_infl + 2 * s_ip
    v_ip = n_ip * (t_ip_max - t_ip_min) + t_ip_min
    ind_ip = np.array([int(np.argmin(np.abs(t - v))) for v in v_ip], dtype=int)
    t_ip = t[ind_ip]

    y_base = np.polyval(p_noise, dt)
    hm = dq[infl] + (y_base[infl] - dq[infl]) * 3.0 / 4.0
    hm_idx = np.where(dq < hm)[0]
    if len(hm_idx) == 0:
        w_iw_max = 3.0
    else:
        ws = [abs(np.min(dt[hm_idx]) - t_infl), abs(np.max(dt[hm_idx]) - t_infl)]
        w_iw_max = float(min(ws))
    w_iw_min = max(2 * s_ip, 0.5)
    if w_iw_max <= w_iw_min:
        w_iw_max = w_iw_min + 0.5

    w_iw = n_iw * (w_iw_max - w_iw_min) + w_iw_min
    ind_il = np.array([int(np.argmin(np.abs(t - (t_ip[i] - w_iw[i])))) for i in range(n)], dtype=int)
    ind_ir = np.array([int(np.argmin(np.abs(t - (t_ip[i] + w_iw[i])))) for i in range(n)], dtype=int)

    y_intersects = np.full(n, np.nan, dtype=float)
    x_intersects = np.full(n, np.nan, dtype=float)
    for i in range(n):
        a, b = sorted((ind_sp[i], ind_sw[i]))
        if b - a < 1:
            continue
        t_sr = t[a : b + 1]
        q_sr = q[a : b + 1]
        p_sr = np.polyfit(t_sr, q_sr, 1)

        c, d = sorted((ind_il[i], ind_ir[i]))
        if d - c < 1:
            d = min(c + 1, len(t) - 1)
            c = max(0, d - 1)
        t_ir = t[c : d + 1]
        q_ir = q[c : d + 1]
        p_ir = np.polyfit(t_ir, q_ir, 1)

        if abs(p_sr[0] - p_ir[0]) < 1e-12:
            continue
        x_int = (p_ir[1] - p_sr[1]) / (p_sr[0] - p_ir[0])
        y_int = p_sr[0] * x_int + p_sr[1]
        x_intersects[i] = x_int
        y_intersects[i] = y_int

    clean_y = y_intersects[np.isfinite(y_intersects)]
    clean_x = x_intersects[np.isfinite(x_intersects)]
    if len(clean_y) == 0:
        h_settle = np.array([np.nan, np.nan], dtype=float)
    else:
        h_settle = np.array([float(np.mean(clean_y)), float(np.std(clean_y))], dtype=float)

    q_start = float(q[0]) if mode == "heat" else float(q[-1])
    h_start = np.array([q_start, rmse], dtype=float)

    x_out = np.array([float(np.mean(clean_x)) if len(clean_x) else np.nan, float(np.std(clean_x)) if len(clean_x) else np.nan], dtype=float)
    detail = [mode, np.nan, x_out]
    return h_start, h_settle, detail, ind_il, ind_ir

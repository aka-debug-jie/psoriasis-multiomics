"""Unique-variant generalized IVW and explicit delta-method uncertainty."""
import numpy as np
from scipy.stats import norm, chi2


def deduplicate_sources(frame):
    """Selection uses exposure evidence only; input order and duplicates do not matter."""
    return frame.sort_values(['F_statistic', 'cohort', 'assay', 'rsid'],
                             ascending=[False, True, True, True]).drop_duplicates(['uniprot', 'variant_key'])


def allele_sign(ea, oa, target_ea, target_oa):
    a = [str(s).upper() for s in (ea, oa, target_ea, target_oa)]
    if any(s not in 'ACGT' or len(s) != 1 for s in a) or a[0] == a[1] or a[2] == a[3]:
        return np.nan
    if set(a[:2]) in ({'A', 'T'}, {'C', 'G'}):
        return np.nan
    if a[:2] == a[2:]: return 1.
    if a[:2] == a[2:][::-1]: return -1.
    comp = dict(A='T', T='A', C='G', G='C')
    b = [comp[s] for s in a[2:]]
    if a[:2] == b: return 1.
    if a[:2] == b[::-1]: return -1.
    return np.nan


def freeze_ld_subset(r, max_abs_r, max_condition=1e4):
    """Rows are already ordered by exposure F. No outcome influences pruning."""
    r = np.asarray(r, float)
    if not np.isfinite(r).all(): return [0], 'single_IV_missing_LD'
    keep = []
    for j in range(len(r)):
        if not keep or np.all(np.abs(r[j, keep]) <= max_abs_r): keep.append(j)
    a = r[np.ix_(keep, keep)]
    if np.linalg.eigvalsh(a).min() <= 0 or np.linalg.cond(a) > max_condition:
        return [0], 'single_IV_unstable_LD'
    return keep, 'LD_pruned' if len(keep) < len(r) else 'all_unique_variants'


def givw(x, y, sx, sy, r, rho_xx=0., rho_xy=0.):
    x, y, sx, sy, r = [np.asarray(a, float) for a in (x, y, sx, sy, r)]
    assert len(x) == len(y) == len(sx) == len(sy) == len(r)
    assert np.isfinite(np.concatenate([x, y, sx, sy])).all() and np.all(sx > 0) and np.all(sy > 0)
    np.linalg.cholesky(r)  # No pseudoinverse, ridge, or nearest-PD substitution.
    s_y = np.outer(sy, sy) * r
    wx = np.linalg.solve(s_y, x)
    wy = np.linalg.solve(s_y, y)
    denominator = x @ wx
    assert denominator > 0
    estimate = x @ wy / denominator
    gy = wx / denominator
    gx = (wy - 2 * estimate * wx) / denominator
    h = np.full_like(r, rho_xx); np.fill_diagonal(h, 1.)
    s_x = np.outer(sx, sx) * r * h
    cross = rho_xy * np.outer(sx, sy) * r
    # Validate the joint correlation, not only separate covariance blocks.
    joint_cor = np.block([[r*h, rho_xy*r], [rho_xy*r, r]])
    if np.linalg.eigvalsh(joint_cor).min() < -1e-10:
        raise ValueError('joint_error_correlation_not_positive_semidefinite')
    variance_delta = gy @ s_y @ gy + gx @ s_x @ gx + 2 * gx @ cross @ gy
    assert variance_delta > 0
    se = np.sqrt(1 / denominator); delta = np.sqrt(variance_delta)
    residual = y - estimate*x
    q = residual @ np.linalg.solve(s_y, residual)
    return dict(beta=float(estimate), se=float(se), se_delta=float(delta),
                p=float(2*norm.sf(abs(estimate/se))), neg_log10_p=float(-(np.log(2)+norm.logsf(abs(estimate/se)))/np.log(10)),
                p_delta=float(2*norm.sf(abs(estimate/delta))),
                CI_log_low=float(estimate-norm.ppf(.975)*se), CI_log_high=float(estimate+norm.ppf(.975)*se),
                Q=float(q), Q_df=len(x)-1, Q_p=float(chi2.sf(q, len(x)-1)) if len(x)>1 else np.nan,
                n_unique_variants=len(x), LD_condition=float(np.linalg.cond(r)))


def p_equivalent_se(beta, p):
    if not np.isfinite(beta) or not np.isfinite(p) or beta == 0 or not 0 < p < 1:
        return np.nan
    return abs(beta)/norm.isf(p/2)


def bh_fixed(p, family_size):
    p = np.asarray(p, float); valid = np.isfinite(p); result = np.full(len(p), np.nan)
    assert family_size >= valid.sum()
    v = p[valid]; order = np.argsort(v)
    adjusted = np.minimum.accumulate((v[order]*family_size/np.arange(1, len(v)+1))[::-1])[::-1].clip(max=1)
    restored = np.empty(len(v)); restored[order] = adjusted; result[valid] = restored
    return result

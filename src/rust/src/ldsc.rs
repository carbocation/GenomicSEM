//! Native core of `ldsc()`: every heritability and genetic-covariance LD-score
//! regression with its block jackknife, run across worker threads.
//!
//! The arithmetic follows R/ldsc.R step by step (same weights, same block
//! boundaries, sequential f64 block sums as in reference BLAS, compensated
//! sums where R accumulates in long double, and 2x2 solves ordered like
//! LAPACK dgesv), so results agree with the R implementation to rounding.
//!
//! Inputs are per-trait columns concatenated across traits (`offsets[t]` ..
//! `offsets[t + 1]`), each trait's rows already merged with the LD scores,
//! filtered on chi-square and sorted by CHR, BP. `id` is each row's index into
//! the LD-score SNP list, so a pair (j, k) is trait j's rows whose SNP is in
//! trait k, in trait j's order, which is what `merge(y1, y2, sort = FALSE)`
//! returns for unique SNP keys.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

pub const LDSC_OK: i32 = 0;
pub const LDSC_SINGULAR: i32 = 1;
pub const LDSC_TOO_FEW_SNPS: i32 = 2;

const ABSENT: u32 = u32::MAX;

/// Neumaier-compensated sum, standing in for R's long double accumulators
/// (`sum`, `mean`, `colSums`).
#[derive(Clone, Copy, Default)]
struct Sum {
    s: f64,
    c: f64,
}

impl Sum {
    #[inline]
    fn add(&mut self, x: f64) {
        let t = self.s + x;
        if self.s.abs() >= x.abs() {
            self.c += (self.s - t) + x;
        } else {
            self.c += (x - t) + self.s;
        }
        self.s = t;
    }
    #[inline]
    fn value(&self) -> f64 {
        self.s + self.c
    }
}

/// R's mean(): sum / n, then one correction pass.
fn r_mean<I: Iterator<Item = f64> + Clone>(values: I, n: usize) -> f64 {
    let mut s = Sum::default();
    for v in values.clone() {
        s.add(v);
    }
    let m = s.value() / n as f64;
    let mut t = Sum::default();
    for v in values {
        t.add(v - m);
    }
    m + t.value() / n as f64
}

/// Solve the 2x2 system as LAPACK dgesv does (dgetrf2 with partial pivoting
/// and reciprocal scaling, then dgetrs). Column-major a = [a00, a10, a01, a11].
/// Returns None when R's solve() would stop (reciprocal condition < eps).
fn solve2(a: [f64; 4], b: [f64; 2]) -> Option<[f64; 2]> {
    let (mut a00, mut a10, mut a01, mut a11) = (a[0], a[1], a[2], a[3]);
    let (mut b0, mut b1) = (b[0], b[1]);
    if a10.abs() > a00.abs() {
        std::mem::swap(&mut a00, &mut a10);
        std::mem::swap(&mut a01, &mut a11);
        std::mem::swap(&mut b0, &mut b1);
    }
    if a00 == 0.0 {
        return None;
    }
    let l = a10 * (1.0 / a00);
    let u11 = a11 + (-a01) * l;
    if u11 == 0.0 {
        return None;
    }
    // Reciprocal condition number in the 1-norm (exact for 2x2), as R's check.
    let det = a[0] * a[3] - a[2] * a[1];
    let norm = (a[0].abs() + a[1].abs()).max(a[2].abs() + a[3].abs());
    let inv_norm = ((a[3].abs() + a[1].abs()).max(a[2].abs() + a[0].abs())) / det.abs();
    if !(1.0 / (norm * inv_norm) >= f64::EPSILON) {
        return None;
    }
    b1 -= b0 * l;
    let x1 = b1 / u11;
    let x0 = (b0 - x1 * a01) / a00;
    Some([x0, x1])
}

pub struct TraitColumns<'a> {
    pub l2: &'a [f64],
    pub wld: &'a [f64],
    pub n: &'a [f64],
    pub z: &'a [f64],
    pub a1: &'a [i32],
    pub id: &'a [i32],
    pub complete: &'a [i32],
}

#[derive(Clone, Default)]
pub struct Regression {
    pub reg_tot: f64,
    pub intercept: f64,
    pub n_bar: f64,
    pub n_snps: f64,
    pub mean_y: f64,
    pub intercept_se: f64,
    pub tot_se: f64,
    pub pseudo: Vec<f64>,
    pub status: i32,
}

/// Block start rows (0-based, length n_blocks + 1, last = n) as
/// floor(seq(1, n, length.out = n_blocks + 1)) - 1.
fn block_starts(n: usize, n_blocks: usize) -> Option<Vec<usize>> {
    let by = (n as f64 - 1.0) / n_blocks as f64;
    let mut starts = Vec::with_capacity(n_blocks + 1);
    starts.push(0usize);
    for i in 1..n_blocks {
        starts.push((1.0 + i as f64 * by).floor() as usize - 1);
    }
    starts.push(n);
    // R's from:to runs backwards for an empty block; require every block to
    // have at least one row so the partition is the one R uses.
    if starts.windows(2).any(|w| w[1] <= w[0]) {
        return None;
    }
    Some(starts)
}

/// Per-row weight before normalisation (R: initial.w).
#[inline]
fn initial_weight(tot_agg: f64, n: f64, l2: f64, wld: f64, m_tot: f64) -> f64 {
    let ld = l2.max(1.0);
    let w_ld = wld.max(1.0);
    let c = tot_agg * n / m_tot;
    let one_c = 1.0 + c * ld;
    let het_w = 1.0 / (2.0 * (one_c * one_c));
    let oc_w = 1.0 / w_ld;
    (het_w * oc_w).sqrt()
}

#[inline]
fn tot_agg(m_tot: f64, mean_chi: f64, mean_l2n: f64) -> f64 {
    let t = (m_tot * (mean_chi - 1.0)) / mean_l2n;
    t.max(0.0).min(1.0)
}

/// One regression's rows: (l2, wld, n1, z1, n2, z2, same_allele).
#[derive(Clone, Copy)]
struct Row {
    l2: f64,
    wld: f64,
    n1: f64,
    z1: f64,
    n2: f64,
    z2: f64,
    same: bool,
}

fn regress(
    rows: &dyn Fn(&mut dyn FnMut(Row)),
    n: usize,
    pair: bool,
    m_tot: f64,
    n_blocks: usize,
) -> Regression {
    let mut out = Regression {
        n_snps: n as f64,
        pseudo: vec![f64::NAN; n_blocks],
        ..Default::default()
    };
    let starts = match block_starts(n, n_blocks) {
        Some(s) => s,
        None => {
            out.status = LDSC_TOO_FEW_SNPS;
            return out;
        }
    };

    // Pass 1: the means the weights and N.bar use (R: mean()).
    let mut collected: Vec<Row> = Vec::with_capacity(n);
    rows(&mut |r| collected.push(r));
    let chi1 = |r: &Row| r.z1 * r.z1;
    let chi2 = |r: &Row| r.z2 * r.z2;
    let mean_chi1 = r_mean(collected.iter().map(chi1), n);
    let mean_l2n1 = r_mean(collected.iter().map(|r| r.l2 * r.n1), n);
    let agg1 = tot_agg(m_tot, mean_chi1, mean_l2n1);
    let mean_n1 = r_mean(collected.iter().map(|r| r.n1), n);
    let (agg2, n_bar) = if pair {
        let mean_chi2 = r_mean(collected.iter().map(chi2), n);
        let mean_l2n2 = r_mean(collected.iter().map(|r| r.l2 * r.n2), n);
        let mean_n2 = r_mean(collected.iter().map(|r| r.n2), n);
        (tot_agg(m_tot, mean_chi2, mean_l2n2), (mean_n1 * mean_n2).sqrt())
    } else {
        (0.0, mean_n1)
    };
    // The regressand: chi1 for heritability, Z1 (aligned to trait 2's A1) * Z2
    // for a covariance.
    let y = |r: &Row| {
        if pair {
            let z1 = if r.same { r.z1 } else { -r.z1 };
            r.z2 * z1
        } else {
            r.z1 * r.z1
        }
    };
    out.mean_y = if pair { r_mean(collected.iter().map(y), n) } else { mean_chi1 };

    // Pass 2: normalising sum of the weights (R: sum(), long double).
    let raw_w = |r: &Row| {
        let w1 = initial_weight(agg1, r.n1, r.l2, r.wld, m_tot);
        if pair {
            w1 + initial_weight(agg2, r.n2, r.l2, r.wld, m_tot)
        } else {
            w1
        }
    };
    let mut wsum = Sum::default();
    for r in &collected {
        wsum.add(raw_w(r));
    }
    let wsum = wsum.value();

    // Pass 3: per-block X'X and X'y (sequential f64, as BLAS does per block).
    // X = cbind(L2, 1) * w, y = y * w.
    let mut xtx_blocks = vec![[0.0f64; 4]; n_blocks];
    let mut xty_blocks = vec![[0.0f64; 2]; n_blocks];
    for b in 0..n_blocks {
        let (mut s00, mut s01, mut s11, mut t0, mut t1) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for r in &collected[starts[b]..starts[b + 1]] {
            let w = raw_w(r) / wsum;
            let xa = r.l2 * w;
            let xb = w;
            let yw = y(r) * w;
            s00 += xa * xa;
            s01 += xa * xb;
            s11 += xb * xb;
            t0 += xa * yw;
            t1 += xb * yw;
        }
        xtx_blocks[b] = [s00, s01, s01, s11];
        xty_blocks[b] = [t0, t1];
    }
    // Totals (R: colSums, long double).
    let mut xtx_sum = [Sum::default(); 4];
    let mut xty_sum = [Sum::default(); 2];
    for b in 0..n_blocks {
        for e in 0..4 {
            xtx_sum[e].add(xtx_blocks[b][e]);
        }
        for e in 0..2 {
            xty_sum[e].add(xty_blocks[b][e]);
        }
    }
    let xtx: [f64; 4] = [xtx_sum[0].value(), xtx_sum[1].value(), xtx_sum[2].value(), xtx_sum[3].value()];
    let xty: [f64; 2] = [xty_sum[0].value(), xty_sum[1].value()];
    let reg = match solve2(xtx, xty) {
        Some(r) => r,
        None => {
            out.status = LDSC_SINGULAR;
            return out;
        }
    };
    let nb = n_blocks as f64;
    let mut pseudo_int = vec![0.0f64; n_blocks];
    for b in 0..n_blocks {
        let a = [
            xtx[0] - xtx_blocks[b][0],
            xtx[1] - xtx_blocks[b][1],
            xtx[2] - xtx_blocks[b][2],
            xtx[3] - xtx_blocks[b][3],
        ];
        let d = match solve2(a, [xty[0] - xty_blocks[b][0], xty[1] - xty_blocks[b][1]]) {
            Some(d) => d,
            None => {
                out.status = LDSC_SINGULAR;
                return out;
            }
        };
        out.pseudo[b] = nb * reg[0] - (nb - 1.0) * d[0];
        pseudo_int[b] = nb * reg[1] - (nb - 1.0) * d[1];
    }
    out.intercept = reg[1];
    out.n_bar = n_bar;
    out.reg_tot = reg[0] / n_bar * m_tot;
    // Jackknife SEs, used only for the log.
    let var = |v: &[f64]| {
        let m = v.iter().sum::<f64>() / nb;
        v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (nb - 1.0)
    };
    out.intercept_se = (var(&pseudo_int) / nb).sqrt();
    let coef_cov = var(&out.pseudo) / nb / (n_bar * n_bar);
    out.tot_se = (coef_cov * (m_tot * m_tot)).sqrt();
    out.status = LDSC_OK;
    out
}

/// Run the regressions for (pairs_j[p], pairs_k[p]); j == k is a heritability.
pub fn run(
    traits: &[TraitColumns<'_>],
    n_universe: usize,
    pairs: &[(usize, usize)],
    m_tot: f64,
    n_blocks: usize,
    requested_threads: usize,
) -> Vec<Regression> {
    // Row position of each LD-score SNP in each trait.
    let positions: Vec<Vec<u32>> = traits
        .iter()
        .map(|t| {
            let mut pos = vec![ABSENT; n_universe];
            for (row, &id) in t.id.iter().enumerate() {
                pos[id as usize] = row as u32;
            }
            pos
        })
        .collect();

    let one = |&(j, k): &(usize, usize)| -> Regression {
        let tj = &traits[j];
        if j == k {
            let n = tj.l2.len();
            let rows = |f: &mut dyn FnMut(Row)| {
                for i in 0..n {
                    f(Row { l2: tj.l2[i], wld: tj.wld[i], n1: tj.n[i], z1: tj.z[i], n2: 0.0, z2: 0.0, same: true });
                }
            };
            regress(&rows, n, false, m_tot, n_blocks)
        } else {
            let tk = &traits[k];
            let pk = &positions[k];
            let keep = |i: usize| tj.complete[i] != 0 && pk[tj.id[i] as usize] != ABSENT;
            let n = (0..tj.l2.len()).filter(|&i| keep(i)).count();
            let rows = |f: &mut dyn FnMut(Row)| {
                for i in 0..tj.l2.len() {
                    if !keep(i) {
                        continue;
                    }
                    let q = pk[tj.id[i] as usize] as usize;
                    f(Row {
                        l2: tj.l2[i],
                        wld: tj.wld[i],
                        n1: tj.n[i],
                        z1: tj.z[i],
                        n2: tk.n[q],
                        z2: tk.z[q],
                        same: tj.a1[i] == tk.a1[q],
                    });
                }
            };
            regress(&rows, n, true, m_tot, n_blocks)
        }
    };

    let workers = requested_threads
        .max(1)
        .min(thread::available_parallelism().map(|c| c.get()).unwrap_or(1))
        .min(pairs.len().max(1));
    if workers <= 1 {
        return pairs.iter().map(one).collect();
    }
    let next = AtomicUsize::new(0);
    let mut results: Vec<Option<Regression>> = vec![None; pairs.len()];
    let collected: Vec<Vec<(usize, Regression)>> = thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut local = Vec::new();
                    loop {
                        let p = next.fetch_add(1, Ordering::Relaxed);
                        if p >= pairs.len() {
                            break;
                        }
                        local.push((p, one(&pairs[p])));
                    }
                    local
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("ldsc worker panicked")).collect()
    });
    for (p, r) in collected.into_iter().flatten() {
        results[p] = Some(r);
    }
    results.into_iter().map(|r| r.expect("missing regression")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_starts_match_r_seq() {
        // floor(seq(1, 1000, length.out = 201)) - 1 for the first entries.
        let s = block_starts(1000, 200).unwrap();
        assert_eq!(s[0], 0);
        assert_eq!(s[1], 4); // floor(1 + 999/200) - 1 = floor(5.995) - 1
        assert_eq!(*s.last().unwrap(), 1000);
        assert!(block_starts(150, 200).is_none());
    }

    #[test]
    fn solve2_matches_direct_solution() {
        let a = [4.0, 2.0, 1.0, 3.0]; // [[4, 1], [2, 3]] column-major
        let x = solve2(a, [6.0, 8.0]).unwrap();
        assert!((4.0 * x[0] + 1.0 * x[1] - 6.0).abs() < 1e-12);
        assert!((2.0 * x[0] + 3.0 * x[1] - 8.0).abs() < 1e-12);
        assert!(solve2([1.0, 2.0, 2.0, 4.0], [1.0, 2.0]).is_none());
    }
}

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::thread;

const STATUS_OK: i32 = 0;
const STATUS_NON_FINITE: i32 = 1;
const STATUS_SINGULAR: i32 = 2;
const STATUS_NUMERICAL: i32 = 3;
const MIN_ROWS_PER_WORKER: usize = 256;

#[derive(Clone, Copy)]
struct Inputs<'a> {
    betas: &'a [f64],
    ses: &'a [f64],
    loadings: &'a [f64],
    corr: &'a [f64],
    q_lu: &'a [f64],
    q_pivots: &'a [usize],
    n: usize,
    traits: usize,
    factors: usize,
}

struct Scratch {
    normal: Vec<f64>,
    rhs: Vec<f64>,
    weighted_loadings: Vec<f64>,
    corr_weighted_loadings: Vec<f64>,
    meat: Vec<f64>,
    inverse: Vec<f64>,
    factor_work: Vec<f64>,
    residual: Vec<f64>,
    trait_work: Vec<f64>,
}

impl Scratch {
    fn new(traits: usize, factors: usize) -> Self {
        Self {
            normal: vec![0.0; factors * factors],
            rhs: vec![0.0; factors],
            weighted_loadings: vec![0.0; traits * factors],
            corr_weighted_loadings: vec![0.0; traits * factors],
            meat: vec![0.0; factors * factors],
            inverse: vec![0.0; factors * factors],
            factor_work: vec![0.0; factors],
            residual: vec![0.0; traits],
            trait_work: vec![0.0; traits],
        }
    }
}

fn cholesky_in_place(matrix: &mut [f64], size: usize) -> bool {
    for row in 0..size {
        for col in 0..=row {
            let mut value = matrix[row * size + col];
            for inner in 0..col {
                value -= matrix[row * size + inner] * matrix[col * size + inner];
            }

            if row == col {
                if !value.is_finite() || value <= 0.0 {
                    return false;
                }
                matrix[row * size + col] = value.sqrt();
            } else {
                matrix[row * size + col] = value / matrix[col * size + col];
            }
        }

        for col in (row + 1)..size {
            matrix[row * size + col] = 0.0;
        }
    }
    true
}

fn cholesky_solve(chol: &[f64], values: &mut [f64], size: usize) {
    for row in 0..size {
        let mut value = values[row];
        for col in 0..row {
            value -= chol[row * size + col] * values[col];
        }
        values[row] = value / chol[row * size + row];
    }

    for row in (0..size).rev() {
        let mut value = values[row];
        for col in (row + 1)..size {
            value -= chol[col * size + row] * values[col];
        }
        values[row] = value / chol[row * size + row];
    }
}

fn lu_in_place(matrix: &mut [f64], pivots: &mut [usize], size: usize) -> bool {
    for col in 0..size {
        let mut pivot = col;
        let mut pivot_value = matrix[col * size + col].abs();
        for row in (col + 1)..size {
            let candidate = matrix[row * size + col].abs();
            if candidate > pivot_value {
                pivot = row;
                pivot_value = candidate;
            }
        }
        if !pivot_value.is_finite() || pivot_value == 0.0 {
            return false;
        }
        pivots[col] = pivot;
        if pivot != col {
            for inner_col in 0..size {
                matrix.swap(col * size + inner_col, pivot * size + inner_col);
            }
        }

        let diagonal = matrix[col * size + col];
        for row in (col + 1)..size {
            matrix[row * size + col] /= diagonal;
            let multiplier = matrix[row * size + col];
            for inner_col in (col + 1)..size {
                matrix[row * size + inner_col] -= multiplier * matrix[col * size + inner_col];
            }
        }
    }
    true
}

fn lu_solve(lu: &[f64], pivots: &[usize], values: &mut [f64], size: usize) {
    for col in 0..size {
        if pivots[col] != col {
            values.swap(col, pivots[col]);
        }
    }
    for row in 0..size {
        for col in 0..row {
            values[row] -= lu[row * size + col] * values[col];
        }
    }
    for row in (0..size).rev() {
        for col in (row + 1)..size {
            values[row] -= lu[row * size + col] * values[col];
        }
        values[row] /= lu[row * size + row];
    }
}

fn invert_from_cholesky(chol: &[f64], inverse: &mut [f64], work: &mut [f64], size: usize) {
    inverse.fill(0.0);
    for col in 0..size {
        work.fill(0.0);
        work[col] = 1.0;
        cholesky_solve(chol, work, size);
        for row in 0..size {
            inverse[row * size + col] = work[row];
        }
    }
}

fn compute_row(
    row: usize,
    inputs: Inputs<'_>,
    beta_out: &mut [f64],
    se_out: &mut [f64],
    q_out: &mut f64,
    scratch: &mut Scratch,
) -> i32 {
    let k = inputs.traits;
    let f = inputs.factors;

    scratch.normal.fill(0.0);
    scratch.rhs.fill(0.0);
    scratch.weighted_loadings.fill(0.0);

    for trait_index in 0..k {
        let beta = inputs.betas[row + inputs.n * trait_index];
        let se = inputs.ses[row + inputs.n * trait_index];
        let corr_diag = inputs.corr[trait_index * k + trait_index];
        if !beta.is_finite() || !se.is_finite() || se <= 0.0 || corr_diag <= 0.0 {
            return STATUS_NON_FINITE;
        }

        let weight = 1.0 / (corr_diag * se * se);
        for factor_a in 0..f {
            let loading_a = inputs.loadings[trait_index + k * factor_a];
            if !loading_a.is_finite() {
                return STATUS_NON_FINITE;
            }
            scratch.rhs[factor_a] += loading_a * weight * beta;
            scratch.weighted_loadings[trait_index * f + factor_a] = loading_a / (corr_diag * se);

            for factor_b in 0..=factor_a {
                let loading_b = inputs.loadings[trait_index + k * factor_b];
                scratch.normal[factor_a * f + factor_b] += loading_a * weight * loading_b;
            }
        }
    }

    for factor_a in 0..f {
        for factor_b in 0..factor_a {
            scratch.normal[factor_b * f + factor_a] = scratch.normal[factor_a * f + factor_b];
        }
    }

    if !cholesky_in_place(&mut scratch.normal, f) {
        return STATUS_SINGULAR;
    }

    cholesky_solve(&scratch.normal, &mut scratch.rhs, f);
    beta_out.copy_from_slice(&scratch.rhs);
    invert_from_cholesky(
        &scratch.normal,
        &mut scratch.inverse,
        &mut scratch.factor_work,
        f,
    );

    scratch.corr_weighted_loadings.fill(0.0);
    for trait_a in 0..k {
        for trait_b in 0..k {
            let corr_value = inputs.corr[trait_a * k + trait_b];
            for factor_index in 0..f {
                scratch.corr_weighted_loadings[trait_a * f + factor_index] +=
                    corr_value * scratch.weighted_loadings[trait_b * f + factor_index];
            }
        }
    }

    scratch.meat.fill(0.0);
    for factor_a in 0..f {
        for factor_b in 0..f {
            let mut value = 0.0;
            for trait_index in 0..k {
                value += scratch.weighted_loadings[trait_index * f + factor_a]
                    * scratch.corr_weighted_loadings[trait_index * f + factor_b];
            }
            scratch.meat[factor_a * f + factor_b] = value;
        }
    }

    for factor_index in 0..f {
        let mut variance = 0.0;
        for factor_a in 0..f {
            for factor_b in 0..f {
                variance += scratch.inverse[factor_index * f + factor_a]
                    * scratch.meat[factor_a * f + factor_b]
                    * scratch.inverse[factor_b * f + factor_index];
            }
        }
        if !variance.is_finite() || variance < -1e-12 {
            return STATUS_NUMERICAL;
        }
        se_out[factor_index] = variance.max(0.0).sqrt();
    }

    for trait_index in 0..k {
        let mut fitted = 0.0;
        for factor_index in 0..f {
            fitted += inputs.loadings[trait_index + k * factor_index] * beta_out[factor_index];
        }
        let beta = inputs.betas[row + inputs.n * trait_index];
        let se = inputs.ses[row + inputs.n * trait_index];
        scratch.residual[trait_index] = (beta - fitted) / se;
    }

    scratch.trait_work.copy_from_slice(&scratch.residual);
    lu_solve(inputs.q_lu, inputs.q_pivots, &mut scratch.trait_work, k);
    let q = scratch
        .residual
        .iter()
        .zip(scratch.trait_work.iter())
        .map(|(left, right)| left * right)
        .sum::<f64>();
    if !q.is_finite() {
        return STATUS_NUMERICAL;
    }
    *q_out = q;

    STATUS_OK
}

fn compute_block(
    start_row: usize,
    inputs: Inputs<'_>,
    beta_out: &mut [f64],
    se_out: &mut [f64],
    q_out: &mut [f64],
    status_out: &mut [i32],
) {
    let mut scratch = Scratch::new(inputs.traits, inputs.factors);
    for local_row in 0..q_out.len() {
        let beta_row = &mut beta_out[local_row * inputs.factors..(local_row + 1) * inputs.factors];
        let se_row = &mut se_out[local_row * inputs.factors..(local_row + 1) * inputs.factors];
        status_out[local_row] = compute_row(
            start_row + local_row,
            inputs,
            beta_row,
            se_row,
            &mut q_out[local_row],
            &mut scratch,
        );
    }
}

fn run_kernel(
    betas: &[f64],
    ses: &[f64],
    loadings: &[f64],
    corr_column_major: &[f64],
    q_corr_column_major: &[f64],
    n: usize,
    traits: usize,
    factors: usize,
    requested_threads: usize,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<i32>), i32> {
    let mut corr = vec![0.0; traits * traits];
    let mut q_corr = vec![0.0; traits * traits];
    for row in 0..traits {
        for col in 0..traits {
            let left = corr_column_major[row + traits * col];
            let right = corr_column_major[col + traits * row];
            if !left.is_finite() || !right.is_finite() {
                return Err(1);
            }
            let tolerance = 64.0 * f64::EPSILON * left.abs().max(right.abs()).max(1.0);
            if (left - right).abs() > tolerance {
                return Err(2);
            }
            corr[row * traits + col] = 0.5 * (left + right);

            let q_left = q_corr_column_major[row + traits * col];
            let q_right = q_corr_column_major[col + traits * row];
            if !q_left.is_finite() || !q_right.is_finite() {
                return Err(1);
            }
            let q_tolerance = 64.0 * f64::EPSILON * q_left.abs().max(q_right.abs()).max(1.0);
            if (q_left - q_right).abs() > q_tolerance {
                return Err(2);
            }
            q_corr[row * traits + col] = 0.5 * (q_left + q_right);
        }
    }

    let mut q_lu = q_corr;
    let mut q_pivots = vec![0; traits];
    if !lu_in_place(&mut q_lu, &mut q_pivots, traits) {
        return Err(3);
    }

    let inputs = Inputs {
        betas,
        ses,
        loadings,
        corr: &corr,
        q_lu: &q_lu,
        q_pivots: &q_pivots,
        n,
        traits,
        factors,
    };

    let mut beta_out = vec![f64::NAN; n * factors];
    let mut se_out = vec![f64::NAN; n * factors];
    let mut q_out = vec![f64::NAN; n];
    let mut status_out = vec![STATUS_OK; n];

    let available_threads = thread::available_parallelism()
        .map(|count| count.get())
        .unwrap_or(1);
    let workers_for_rows = n.div_ceil(MIN_ROWS_PER_WORKER).max(1);
    let thread_count = requested_threads
        .max(1)
        .min(available_threads)
        .min(workers_for_rows)
        .min(n.max(1));
    let rows_per_thread = n.div_ceil(thread_count);

    if thread_count == 1 || n == 0 {
        compute_block(
            0,
            inputs,
            &mut beta_out,
            &mut se_out,
            &mut q_out,
            &mut status_out,
        );
    } else {
        thread::scope(|scope| {
            for (chunk_index, (((beta_chunk, se_chunk), q_chunk), status_chunk)) in beta_out
                .chunks_mut(rows_per_thread * factors)
                .zip(se_out.chunks_mut(rows_per_thread * factors))
                .zip(q_out.chunks_mut(rows_per_thread))
                .zip(status_out.chunks_mut(rows_per_thread))
                .enumerate()
            {
                let start_row = chunk_index * rows_per_thread;
                scope.spawn(move || {
                    compute_block(
                        start_row,
                        inputs,
                        beta_chunk,
                        se_chunk,
                        q_chunk,
                        status_chunk,
                    );
                });
            }
        });
    }

    Ok((beta_out, se_out, q_out, status_out))
}

/// Compute analytic GenomicSEM estimates for a batch of SNPs.
///
/// All input matrices use R's column-major layout. Factor estimates and standard
/// errors are written back in the same layout. The function returns zero on
/// success, a positive global validation code, or -1 if a Rust panic was caught.
#[no_mangle]
pub unsafe extern "C" fn genomicsem_gls_batch(
    betas: *const f64,
    ses: *const f64,
    n: usize,
    traits: usize,
    loadings: *const f64,
    factors: usize,
    corr: *const f64,
    q_corr: *const f64,
    requested_threads: usize,
    beta_out: *mut f64,
    se_out: *mut f64,
    q_out: *mut f64,
    status_out: *mut i32,
) -> i32 {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let betas = slice::from_raw_parts(betas, n * traits);
        let ses = slice::from_raw_parts(ses, n * traits);
        let loadings = slice::from_raw_parts(loadings, traits * factors);
        let corr = slice::from_raw_parts(corr, traits * traits);
        let q_corr = slice::from_raw_parts(q_corr, traits * traits);

        run_kernel(
            betas,
            ses,
            loadings,
            corr,
            q_corr,
            n,
            traits,
            factors,
            requested_threads,
        )
    }));

    match result {
        Ok(Ok((beta_values, se_values, q_values, status_values))) => {
            for factor_index in 0..factors {
                for row in 0..n {
                    *beta_out.add(row + n * factor_index) =
                        beta_values[row * factors + factor_index];
                    *se_out.add(row + n * factor_index) = se_values[row * factors + factor_index];
                }
            }
            std::ptr::copy_nonoverlapping(q_values.as_ptr(), q_out, n);
            std::ptr::copy_nonoverlapping(status_values.as_ptr(), status_out, n);
            0
        }
        Ok(Err(code)) => code,
        Err(_) => -1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_one_factor_has_expected_estimate() {
        let betas = vec![1.0, 2.0, 3.0];
        let ses = vec![1.0, 1.0, 1.0];
        let loadings = vec![1.0, 1.0, 1.0];
        let corr = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let (beta, se, q, status) =
            run_kernel(&betas, &ses, &loadings, &corr, &corr, 1, 3, 1, 1).unwrap();
        assert_eq!(status, vec![STATUS_OK]);
        assert!((beta[0] - 2.0).abs() < 1e-12);
        assert!((se[0] - (1.0_f64 / 3.0).sqrt()).abs() < 1e-12);
        assert!((q[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn multiple_threads_are_deterministic() {
        let betas = vec![1.0, 2.0, 1.5, 2.5, 2.0, 3.0];
        let ses = vec![1.0; 6];
        let loadings = vec![1.0, 1.0, 1.0];
        let corr = vec![1.0, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0];
        let serial = run_kernel(&betas, &ses, &loadings, &corr, &corr, 2, 3, 1, 1).unwrap();
        let parallel = run_kernel(&betas, &ses, &loadings, &corr, &corr, 2, 3, 1, 2).unwrap();
        assert_eq!(serial, parallel);
    }
}

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
    q_inverse: &'a [f64],
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

#[inline(always)]
fn symmetric_quadratic_form(inverse: &[f64], residual: &[f64]) -> Option<f64> {
    let size = residual.len();
    let mut value = 0.0;
    for row in 0..size {
        let residual_row = residual[row];
        value += inverse[row * size + row] * residual_row * residual_row;
        for col in 0..row {
            value += 2.0 * inverse[row * size + col] * residual_row * residual[col];
        }
    }
    value.is_finite().then_some(value)
}

fn compute_row_one_factor(
    row: usize,
    inputs: Inputs<'_>,
    beta_out: &mut [f64],
    se_out: &mut [f64],
    q_out: &mut f64,
    scratch: &mut Scratch,
) -> i32 {
    let k = inputs.traits;
    let mut normal = 0.0;
    let mut rhs = 0.0;

    for trait_index in 0..k {
        let beta = inputs.betas[row + inputs.n * trait_index];
        let se = inputs.ses[row + inputs.n * trait_index];
        let corr_diag = inputs.corr[trait_index * k + trait_index];
        let loading = inputs.loadings[trait_index];
        if !beta.is_finite()
            || !se.is_finite()
            || se <= 0.0
            || corr_diag <= 0.0
            || !loading.is_finite()
        {
            return STATUS_NON_FINITE;
        }

        let weight = 1.0 / (corr_diag * se * se);
        normal += loading * weight * loading;
        rhs += loading * weight * beta;
        scratch.weighted_loadings[trait_index] = loading / (corr_diag * se);
    }
    if !normal.is_finite() || normal <= 0.0 {
        return STATUS_SINGULAR;
    }

    let factor_beta = rhs / normal;
    if !factor_beta.is_finite() {
        return STATUS_NUMERICAL;
    }
    beta_out[0] = factor_beta;

    for trait_a in 0..k {
        let mut value = 0.0;
        for trait_b in 0..k {
            value += inputs.corr[trait_a * k + trait_b] * scratch.weighted_loadings[trait_b];
        }
        scratch.corr_weighted_loadings[trait_a] = value;
    }
    let mut meat = 0.0;
    for trait_index in 0..k {
        meat +=
            scratch.weighted_loadings[trait_index] * scratch.corr_weighted_loadings[trait_index];
    }
    let inverse = 1.0 / normal;
    let variance = inverse * meat * inverse;
    if !variance.is_finite() || variance < -1e-12 {
        return STATUS_NUMERICAL;
    }
    se_out[0] = variance.max(0.0).sqrt();

    for trait_index in 0..k {
        let beta = inputs.betas[row + inputs.n * trait_index];
        let se = inputs.ses[row + inputs.n * trait_index];
        scratch.residual[trait_index] = (beta - inputs.loadings[trait_index] * factor_beta) / se;
    }
    let Some(q) = symmetric_quadratic_form(inputs.q_inverse, &scratch.residual) else {
        return STATUS_NUMERICAL;
    };
    *q_out = q;

    STATUS_OK
}

fn compute_row_two_factors(
    row: usize,
    inputs: Inputs<'_>,
    beta_out: &mut [f64],
    se_out: &mut [f64],
    q_out: &mut f64,
    scratch: &mut Scratch,
) -> i32 {
    let k = inputs.traits;
    let mut normal_00 = 0.0;
    let mut normal_01 = 0.0;
    let mut normal_11 = 0.0;
    let mut rhs_0 = 0.0;
    let mut rhs_1 = 0.0;

    for trait_index in 0..k {
        let beta = inputs.betas[row + inputs.n * trait_index];
        let se = inputs.ses[row + inputs.n * trait_index];
        let corr_diag = inputs.corr[trait_index * k + trait_index];
        let loading_0 = inputs.loadings[trait_index];
        let loading_1 = inputs.loadings[trait_index + k];
        if !beta.is_finite()
            || !se.is_finite()
            || se <= 0.0
            || corr_diag <= 0.0
            || !loading_0.is_finite()
            || !loading_1.is_finite()
        {
            return STATUS_NON_FINITE;
        }

        let weight = 1.0 / (corr_diag * se * se);
        normal_00 += loading_0 * weight * loading_0;
        normal_01 += loading_0 * weight * loading_1;
        normal_11 += loading_1 * weight * loading_1;
        rhs_0 += loading_0 * weight * beta;
        rhs_1 += loading_1 * weight * beta;
        scratch.weighted_loadings[trait_index * 2] = loading_0 / (corr_diag * se);
        scratch.weighted_loadings[trait_index * 2 + 1] = loading_1 / (corr_diag * se);
    }

    if !normal_00.is_finite() || normal_00 <= 0.0 {
        return STATUS_SINGULAR;
    }
    let chol_00 = normal_00.sqrt();
    let chol_10 = normal_01 / chol_00;
    let chol_11_squared = normal_11 - chol_10 * chol_10;
    if !chol_10.is_finite() || !chol_11_squared.is_finite() || chol_11_squared <= 0.0 {
        return STATUS_SINGULAR;
    }
    let chol_11 = chol_11_squared.sqrt();

    let forward_0 = rhs_0 / chol_00;
    let forward_1 = (rhs_1 - chol_10 * forward_0) / chol_11;
    let factor_beta_1 = forward_1 / chol_11;
    let factor_beta_0 = (forward_0 - chol_10 * factor_beta_1) / chol_00;
    if !factor_beta_0.is_finite() || !factor_beta_1.is_finite() {
        return STATUS_NUMERICAL;
    }
    beta_out[0] = factor_beta_0;
    beta_out[1] = factor_beta_1;

    let inverse_chol_00 = 1.0 / chol_00;
    let inverse_11 = 1.0 / chol_11_squared;
    let inverse_01 = -chol_10 * inverse_chol_00 * inverse_11;
    let scaled_chol_10 = chol_10 * inverse_chol_00;
    let inverse_00 = 1.0 / normal_00 + scaled_chol_10 * scaled_chol_10 * inverse_11;
    if !inverse_00.is_finite() || !inverse_01.is_finite() || !inverse_11.is_finite() {
        return STATUS_NUMERICAL;
    }

    for trait_a in 0..k {
        let mut value_0 = 0.0;
        let mut value_1 = 0.0;
        for trait_b in 0..k {
            let corr_value = inputs.corr[trait_a * k + trait_b];
            value_0 += corr_value * scratch.weighted_loadings[trait_b * 2];
            value_1 += corr_value * scratch.weighted_loadings[trait_b * 2 + 1];
        }
        scratch.corr_weighted_loadings[trait_a * 2] = value_0;
        scratch.corr_weighted_loadings[trait_a * 2 + 1] = value_1;
    }

    let mut meat_00 = 0.0;
    let mut meat_01 = 0.0;
    let mut meat_10 = 0.0;
    let mut meat_11 = 0.0;
    for trait_index in 0..k {
        let weighted_0 = scratch.weighted_loadings[trait_index * 2];
        let weighted_1 = scratch.weighted_loadings[trait_index * 2 + 1];
        let corr_weighted_0 = scratch.corr_weighted_loadings[trait_index * 2];
        let corr_weighted_1 = scratch.corr_weighted_loadings[trait_index * 2 + 1];
        meat_00 += weighted_0 * corr_weighted_0;
        meat_01 += weighted_0 * corr_weighted_1;
        meat_10 += weighted_1 * corr_weighted_0;
        meat_11 += weighted_1 * corr_weighted_1;
    }

    let variance_0 = inverse_00 * inverse_00 * meat_00
        + inverse_00 * inverse_01 * (meat_01 + meat_10)
        + inverse_01 * inverse_01 * meat_11;
    let variance_1 = inverse_01 * inverse_01 * meat_00
        + inverse_01 * inverse_11 * (meat_01 + meat_10)
        + inverse_11 * inverse_11 * meat_11;
    if !variance_0.is_finite()
        || variance_0 < -1e-12
        || !variance_1.is_finite()
        || variance_1 < -1e-12
    {
        return STATUS_NUMERICAL;
    }
    se_out[0] = variance_0.max(0.0).sqrt();
    se_out[1] = variance_1.max(0.0).sqrt();

    for trait_index in 0..k {
        let fitted = inputs.loadings[trait_index] * factor_beta_0
            + inputs.loadings[trait_index + k] * factor_beta_1;
        let beta = inputs.betas[row + inputs.n * trait_index];
        let se = inputs.ses[row + inputs.n * trait_index];
        scratch.residual[trait_index] = (beta - fitted) / se;
    }
    let Some(q) = symmetric_quadratic_form(inputs.q_inverse, &scratch.residual) else {
        return STATUS_NUMERICAL;
    };
    *q_out = q;

    STATUS_OK
}

fn compute_row_generic(
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

    let Some(q) = symmetric_quadratic_form(inputs.q_inverse, &scratch.residual) else {
        return STATUS_NUMERICAL;
    };
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
        status_out[local_row] = match inputs.factors {
            1 => compute_row_one_factor(
                start_row + local_row,
                inputs,
                beta_row,
                se_row,
                &mut q_out[local_row],
                &mut scratch,
            ),
            2 => compute_row_two_factors(
                start_row + local_row,
                inputs,
                beta_row,
                se_row,
                &mut q_out[local_row],
                &mut scratch,
            ),
            _ => compute_row_generic(
                start_row + local_row,
                inputs,
                beta_row,
                se_row,
                &mut q_out[local_row],
                &mut scratch,
            ),
        };
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
    let mut q_inverse = vec![0.0; traits * traits];
    let mut q_work = vec![0.0; traits];
    for col in 0..traits {
        q_work.fill(0.0);
        q_work[col] = 1.0;
        lu_solve(&q_lu, &q_pivots, &mut q_work, traits);
        for row in 0..traits {
            q_inverse[row * traits + col] = q_work[row];
        }
    }
    for row in 0..traits {
        for col in 0..row {
            let symmetric = 0.5 * (q_inverse[row * traits + col] + q_inverse[col * traits + row]);
            q_inverse[row * traits + col] = symmetric;
            q_inverse[col * traits + row] = symmetric;
        }
    }

    let inputs = Inputs {
        betas,
        ses,
        loadings,
        corr: &corr,
        q_inverse: &q_inverse,
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

use rand_xoshiro::rand_core::{RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use std::f64::consts::PI;
use std::time::{Duration, Instant};

const HIT_MISS_MAX: usize = 10_000_000;
const HIT_MISS_MIN: usize = 10;
const HIT_MISS_MUL: usize = 10;
const HIT_MISS_ITER: usize = 100;

const INTEGRAL_MAX: usize = 1_000_000;
const INTEGRAL_MIN: usize = 10;
const INTEGRAL_MUL: usize = 10;
const INTEGRAL_ITER: usize = 100;

const I_VALUE : f64 = 0.602845;

const MONTE_CARLO_2D : usize = 1_000_000;
const QUADRATURE_2D : usize = 1_000;
const STRATIFICATION_2D : usize = 1_000;

const SEED: u64 = 7;

fn main() {
    let mut xoshiro = Xoshiro256PlusPlus::seed_from_u64(SEED);

    let mut iter = HIT_MISS_MIN;

    while iter <= HIT_MISS_MAX{

        let mut vec_pi: Vec<f64> = Vec::new();

        for _b in 0..HIT_MISS_ITER{

            let mut pi_buffer = 0u64;
            for _i in 0..iter {
                let x = normalize(xoshiro.next_u64());
                let y = normalize(xoshiro.next_u64());
                if x * x + y * y <= 1f64 {
                    pi_buffer += 1;
                }
            }
            vec_pi.push((pi_buffer as f64 / iter as f64) * 4f64);
        }

        let mut somme_err = 0f64;
        for i in &vec_pi{
            somme_err += f64::powi(i - PI,2);
        }

        let mean_squared_error = somme_err/iter as f64;

        println!("Monte Carlo MSE for {} iter is : {}",iter,mean_squared_error);
        //println!("{},{}",iter,mean_squared_error);
        //println!("pi : {}",vec_pi.get(0).unwrap());

        iter *= HIT_MISS_MUL;
    }

    //EXO 2
    let mut iter_integral = INTEGRAL_MIN;

    while iter_integral <= INTEGRAL_MAX{

        let mut vec_monte_carlo: Vec<f64> = Vec::new();
        let mut vec_quadrature: Vec<f64> = Vec::new();
        let mut vec_stratification: Vec<f64> = Vec::new();

        let mut monte_carlo_buffer = 0f64;
        let min = 0.1f64;
        let max = 0.9f64;

        for _a in 0..INTEGRAL_ITER {
        for _i in 0..iter_integral {
            monte_carlo_buffer += j_function(scale(xoshiro.next_u64(), min, max));
        }

        vec_monte_carlo.push(monte_carlo_buffer * (max - min) / iter_integral as f64);

        /*println!(
            "monte carlo : Integral of J is equal to : {}",
            monte_carlo_integral
        );*/

        let mut quadrature_buffer = 0f64;

        for i in 0..iter_integral {
            quadrature_buffer += j_function(
                min + (max - min) * (i as f64 / iter_integral as f64)
                    + (max - min) * (1f64 / iter_integral as f64),
            );
        }

        vec_quadrature.push(quadrature_buffer * (max - min) / iter_integral as f64);

        /*println!(
            "quadrature : Integral of J is equal to : {}",
            quadrature_integral
        );*/

        let mut stratification_buffer = 0f64;

        for i in 0..iter_integral {
            stratification_buffer += j_function(scale(
                xoshiro.next_u64(),
                min + (max - min) * (i as f64 / iter_integral as f64),
                min + (max - min) * ((i + 1) as f64 / iter_integral as f64),
            ));
        }

        vec_stratification.push(stratification_buffer * (max - min) / iter_integral as f64);
    }

    let mut somme_err_monte_carlo = 0f64;
    let mut somme_err_quadrature = 0f64;
    let mut somme_err_stratifiaction = 0f64;

        for i in &vec_monte_carlo{
            somme_err_monte_carlo += f64::powi(i - I_VALUE,2);
        }
        for i in &vec_quadrature{
            somme_err_quadrature += f64::powi(i - I_VALUE,2);
        }
        for i in &vec_stratification{
            somme_err_stratifiaction += f64::powi(i - I_VALUE,2);
        }

        let mse_monte_carlo = somme_err_monte_carlo/iter_integral as f64;
        let mse_quadrature = somme_err_quadrature/iter_integral as f64;
        let mse_stratifiaction = somme_err_stratifiaction/iter_integral as f64;

        println!("{},{},{},{}",iter_integral,mse_monte_carlo,mse_quadrature,mse_stratifiaction);



        /*println!(
            "stratification : Integral of J is equal to : {}",
            stratification_integral
        );*/

        iter_integral *= INTEGRAL_MUL;
    }

    // EXO 4

    let mut monte_carlo_buffer = 0f64;
    let min = 0f64;
    let max = 1f64;

    for _i in 0..MONTE_CARLO_2D {
        monte_carlo_buffer += i_function(
            scale(xoshiro.next_u64(), min, max),
            scale(xoshiro.next_u64(), min, max),
        );
    }

    let monte_carlo_integral = monte_carlo_buffer * (max - min) / MONTE_CARLO_2D as f64;

    println!(
        "monte carlo : Integral of I is equal to : {}",
        monte_carlo_integral
    );

    let mut quadrature_buffer = 0f64;

    for i in 0..QUADRATURE_2D {
        for j in 0..QUADRATURE_2D {
            quadrature_buffer += i_function(
                min + (max - min) * (i as f64 / QUADRATURE_2D as f64)
                    + (max - min) * (1f64 / QUADRATURE_2D as f64),
                min + (max - min) * (j as f64 / QUADRATURE_2D as f64)
                    + (max - min) * (1f64 / QUADRATURE_2D as f64),
            );
        }
    }

    let quadrature_integral =
        quadrature_buffer * (max - min) / (QUADRATURE_2D * QUADRATURE_2D) as f64;

    println!(
        "quadrature : Integral of I is equal to : {}",
        quadrature_integral
    );

    let mut stratification_buffer = 0f64;

    for i in 0..STRATIFICATION_2D {
        for j in 0..STRATIFICATION_2D {
            stratification_buffer += i_function(
                scale(
                    xoshiro.next_u64(),
                    min + (max - min) * (i as f64 / STRATIFICATION_2D as f64),
                    min + (max - min) * ((i + 1) as f64 / STRATIFICATION_2D as f64),
                ),
                scale(
                    xoshiro.next_u64(),
                    min + (max - min) * (j as f64 / STRATIFICATION_2D as f64),
                    min + (max - min) * ((j + 1) as f64 / STRATIFICATION_2D as f64),
                ),
            );
        }
    }

    let stratification_integral =
        stratification_buffer * (max - min) / (STRATIFICATION_2D * STRATIFICATION_2D) as f64;

    println!(
        "stratification : Integral of I is equal to : {}",
        stratification_integral
    );
}

fn scale(input: u64, min: f64, max: f64) -> f64 {
    return (input as f64 / u64::MAX as f64) * (max - min) + min;
}

fn normalize(input: u64) -> f64 {
    return input as f64 / u64::MAX as f64;
}

fn j_function(x: f64) -> f64 {
    return (1f64 / 2f64) * f64::exp(f64::asin(x)) * (x / (f64::sqrt(1f64 - (x * x))));
}

fn i_function(x: f64, y: f64) -> f64 {
    return x * y * f64::powi(f64::sin(1f64 / (x * y)), 2);
}

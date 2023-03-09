#include "sboxu_cpp_fp_lat.hpp"

using std::vector;
using std::complex;
using std::map;

int my_pow(int x, unsigned int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;
  
  int temp = my_pow(x, p/2);
  if (p%2 == 0) 
    return temp * temp;
  else
    return x * temp * temp;
}

int fp_n_scalar(int alpha, int x, int p, int m) {
    /***
    Returns alpha dot x, with x and alpha of the form x_{m-1}*p^{m-1} + ... + x_0
    ***/
    int sol = 0;
    int gen = 1;
    int temp = x;
    for (int i = 0; i < m; i++) {
        sol = sol + gen * ((alpha * (temp%p))%p);
        temp /= p;
        gen *= p;
    }
    return sol;
}

void fft(const vector<int>& pol, vector<vector<double>>& walsh_tab, fftw_plan pl, double *in, const int p, const int m, int a, int b) {
    /***
    Updates walsh_tab with the p coefficients outputted by the fft.
    ***/
    for (int i = 0; i < p; i++) {
        in[i] = pol[i];
    }

    fftw_execute(pl);
    double walsh_coef;
    // the output is of the form r0, r1, ... , r(p/2), i((n+1)/2-1), ... , i2, i1.
    // see fftw documentation.
    for (int i = 1; i <= (p-1)/2; i++) {
        walsh_coef = sqrt(in[i]*in[i] + in[p-i]*in[p-i]);
        walsh_tab[fp_n_scalar(i, a, p, m)][fp_n_scalar(i, b, p, m)] = walsh_coef > 1e-6 ? walsh_coef : 0.0;
        walsh_tab[fp_n_scalar((p-i), a, p, m)][fp_n_scalar((p-i), b, p, m)] = walsh_coef > 1e-6 ? walsh_coef : 0.0;
    }

    // p/2 is treated specially if p = 2.
    if (p == 2) {
        walsh_coef = sqrt(in[1]*in[1]);
        walsh_tab[a][b] =  walsh_coef > 1e-6 ? walsh_coef : 0.0;
    }
}

int my_dot_product(int a, int b, const int p, const int m) {
    /***
    If  a = a_{m-1} * p^(m-1) + ... + a_1 * p + a_0
    and b = b_{m-1} * p^(m-1) + ... + b_1 * p + b_0,
    Returns sum(a_i * b_i) mod p.
    ***/

    int sol = 0;
    int a_mod_p;
    int b_mod_p;
    for (int i = 0; i < m-1; i++) {
        a_mod_p = a%p;
        a = a/p;
        b_mod_p = b%p;
        b = b/p;
        sol += a_mod_p * b_mod_p;

    }
    a_mod_p = a%p;
    b_mod_p = b%p;
    sol += a_mod_p * b_mod_p;
    return sol % p;
}

void sum_of_polys(const int p, const vector<int>& array_of_polys, vector<int>& temp_poly, const int a_i, const int a_im1_to_1){
    /***
    Modifies temp_poly as the sum over x_i of the polynomials X^{a_i*x_i} * array_of_polys(poly of index a_im1_to_1 + x_i).
    This assumes that a_im1_to_1 is divisible by p.
    Time: O(p^2)
    ***/
    for (int x_i = 0; x_i < p; x_i++)
    {
        int prod = (a_i*x_i)%p;
        for (int j = 0; j < p - prod; j++) {
            temp_poly[prod + j] += array_of_polys[p*(a_im1_to_1 + x_i)+ j];
        }
        for (int j = p - prod; j < p; j++) {
            temp_poly[prod + j - p] += array_of_polys[p*(a_im1_to_1 + x_i)+ j];
        }
    }
}

void mid_sum_of_polys(const int p, const vector<int>& prev_result, vector<int>& array_of_polys, const int x_ip1, const int a_i, const int  a_im1_to_1, const int prev_size){
    /***
    Variant of sum_of_polys used in mid_compute.
    This time, takes two arrays of polys and updates one of them accordingly.
    Time: O(p^2)
    ***/
    for (int x_i = 0; x_i < p; x_i++)
    {
        for (int j = 0; j < p; j++) {
            array_of_polys[prev_size*a_i + p*(a_im1_to_1 + x_ip1) + (j + a_i*x_i)%p] += prev_result[p*(a_im1_to_1 + x_i)+ j];
        }
    }
}

void mid_compute(
    const int p, const int m, vector<int>& array_of_polys, const int array_size, const vector<int>& lut,
    const int b, const int a_0, const int x_mm1_to_i) 
{


    /***
    Recursive function that performs a bunch of middle computations
    and returns something that can be digested by the lat functions.
    Dynamic programming stuff. See doc for details.
    
    The lut is a function of F_{p^m}.
    array_to_modify is an array of size p^{m - depth} containing arrays of size p.
    At the start we always have that array_size = p^{m - depth + 1}.
    This function is used in the case m >= 2.
    x_to_depth is updated to be equal to, when called at the base case, p^{m-1}*x_{m-1} + ... + p*x_1.  

    First call will be mid_compute(p, m, p^m, lut, b, a_0)
    ***/
    int prev_size = array_size/p;
    vector<int> prev_result;
    int expo;

    if (prev_size == p) {
        // Base case
        for (int x_1 = 0; x_1 < p; x_1++) {
            for (int x_0 = 0; x_0 < p; x_0++) {
                expo = (a_0*x_0 - my_dot_product(b, lut[p*(p*x_mm1_to_i + x_1) + x_0], p, m) + p)%p;
                array_of_polys[p*x_1 + expo]++;
            }
        }

    }


    else {
        for (int x_ip1 = 0; x_ip1 < p; x_ip1++) {
            prev_result = vector<int>(prev_size);

            mid_compute(p, m, prev_result, prev_size, lut, b, a_0, p*x_mm1_to_i + x_ip1);
            for (int a_i = 0; a_i < p; a_i++) {
                for (int a_im1_to_1 = 0; a_im1_to_1 < prev_size/p; a_im1_to_1 += p) {
                    mid_sum_of_polys(p, prev_result, array_of_polys, x_ip1, a_i, a_im1_to_1, prev_size);
                }
            }

        }
    }
}

void lat_diagonal(const vector<int>& lut, vector<vector<double>>& walsh_tab, const int p, const int m, const int q, const int q_div_p, const int b, const fftw_plan pl, double* in) {
    // Now to compute the "polynomial" tied to all a's for a given b.

    if (m == 1) {
        // No need for intermediary computations. Compute for a, then fft.
        for (int a = 0; a < p; a++) {
            vector<int> temp_poly = vector<int>(q, 0);
            for (int x = 0; x < p; x++) {
                int expo = (a*x - (b*lut[x])%p + p)%p;
                temp_poly[expo]++;
            }
            fft(temp_poly, walsh_tab, pl, in, p, m, a, b);
        }
    }

    else {
        for (int a_0 = 0; a_0 < p; a_0++) {
            vector<int> array_of_polys = vector<int>(q);
            vector<int> temp_poly;
            mid_compute(p, m, array_of_polys, q, lut, b, a_0, 0);
            // iterate over all a's of fixed a_0
            for (int a_mm1 = 0; a_mm1 < p; a_mm1++) {
                for (int a_mm2_to_1 = 0; a_mm2_to_1 < q_div_p; a_mm2_to_1 += p) {
                    temp_poly = vector<int>(p);
                    sum_of_polys(p, array_of_polys, temp_poly, a_mm1, a_mm2_to_1);
                    fft(temp_poly, walsh_tab, pl, in, p, m, q_div_p * a_mm1 + a_mm2_to_1 + a_0, b);
                }
            }
        }
    }
}

vector<vector<double>> fpt_lat(const vector<int>& lut, const int p, const int m, const unsigned int n_threads) {

    fftw_init_threads();
    fftw_plan_with_nthreads(n_threads);

    const int q = my_pow(p, m);
    const int q_div_p = q/p;
    
    vector<vector<double>> walsh_tab(q, vector<double>(q, 0.0));
    walsh_tab[0][0] = q;

    vector<int> temp_poly;
    vector<int> array_of_polys; // This vector will contain q/p polynomials of degree p.

    // We will iterate over all b < p^m whose first nonzero coefficient is 1.
    // A way of doing that is by iterating over the part to the right of that 1.

    omp_set_num_threads(n_threads);
    #pragma omp parallel private(temp_poly, array_of_polys) firstprivate(lut, p, m, q, q_div_p)
    {
    // Next is fft stuff:
    // in-place dft of a real input to a hermitian-complex input.

    double* in;
    fftw_plan pl;
    #pragma omp critical
    {
    in = fftw_alloc_real(p);
    pl = fftw_plan_r2r_1d(p, in, in, FFTW_R2HC, FFTW_MEASURE);
    }

    int first_nonzero_index;
    int b_temp;

    #pragma omp for schedule(dynamic)
    for(int b_right = 0; b_right < q_div_p; b_right++) {
        b_temp = b_right;
        first_nonzero_index = 0;
        while (b_temp != 0) {
            b_temp /= p;
            first_nonzero_index++;
        }
        b_temp = my_pow(p, first_nonzero_index) + b_right;
        while (b_temp < q) {
            lat_diagonal(lut, walsh_tab, p, m, q, q_div_p, b_temp, pl, in);
            first_nonzero_index++;
            b_temp = b_right + my_pow(p, first_nonzero_index);
        }
    }

    #pragma omp critical
    {
    fftw_free(in);
    fftw_free(pl);
    }
    }

    return walsh_tab;
}

bool test_first_nonzero_is_one(const int a, const int p) {
    int temp = a;
    while (temp != 0 && temp%p == 0) {
        temp /= p;
    }
    return ((temp == 0) || (temp%p == 1));
}


void fft_column(const vector<complex<float>>& pol, vector<double>& walsh_column, const int a, const fftw_plan pl, fftw_complex* in, const int p, const int m) {
    
    for (int i = 0; i < p; i++) {
        in[i][0] = std::real(pol[i]);
        in[i][1] = std::imag(pol[i]);
    }
    
    fftw_execute(pl);
    double walsh_coef;
    
    for (int i = 0; i < p; i++) {
        walsh_coef = sqrt(in[i][0]*in[i][0] + in[i][1]*in[i][1]);
        walsh_column[fp_n_scalar(i, a, p, m)] = walsh_coef > 1e-6 ? walsh_coef : 0.0;
    }
}



vector<double> fpt_lat_column(const vector<int>& lut, const int p, const int m, const int b) {
    const int q = my_pow(p, m);
    const int q_div_p = q/p;
    vector<double> walsh_column = vector<double>(q, 0);

    vector<complex<float>> omega_power(p, complex<float>{0,0});
    for (int i = 0; i < p; i++) {
        float phase = 2 * M_PI * float(i)/p;
        omega_power[i].real(cos(phase));
        omega_power[i].imag(sin(phase));
    }

    if (m == 1) {   // No need for intermediary computations. Compute for b, then fft.

        fftw_complex* in;
        fftw_plan pl;
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * p);
        pl = fftw_plan_dft_1d(p, in, in, FFTW_BACKWARD, FFTW_MEASURE);

        vector<complex<float>> temp_poly = vector<complex<float>>(p);
        for (int x = 0; x < p; x++) {
            temp_poly[x] += omega_power[((-b*lut[x])%p + p)%p];
        }
        fft_column(temp_poly, walsh_column, 1, pl, in, p, m);
        fftw_free(in);
        fftw_free(pl);
    }

    else {
        // No ffts, as this would require many costly operations on complex numbers.
        for (int a_0 = 0; a_0 < p; a_0++) {
            vector<int> array_of_polys = vector<int>(q);
            vector<int> temp_poly;
            mid_compute(p, m, array_of_polys, q, lut, b, a_0, 0);
            for (int a_mm1 = 0; a_mm1 < p; a_mm1++) {
                for (int a_mm2_to_1 = 0; a_mm2_to_1 < q_div_p; a_mm2_to_1 += p) {
                    int a = q_div_p * a_mm1 + a_mm2_to_1 + a_0;
                    temp_poly = vector<int>(p);
                    sum_of_polys(p, array_of_polys, temp_poly, a_mm1, a_mm2_to_1);
                    complex<float> pol_eval = 0;
                    for (int index = 0; index <p; index++) {
                        pol_eval += float(temp_poly[index]) * omega_power[index];
                    }
                    double walsh_coeff = std::abs(pol_eval);
                    walsh_column[a] = walsh_coeff > 1e-6 ? walsh_coeff : 0;
                }
            }
        }
    }


    return walsh_column;    
}


vector<double> fpt_lat_row(const vector<int>& lut, const int p, const int m, const int a) {
    /***
    If m = 1, proceed as usual. If not, try to compute the column of the inverse function.
    Raise an exception otherwise as we are not aware of equivalently optimized algorithms.
    ***/
    const int q = my_pow(p, m);
    vector<int>inv_lut = vector<int>(q, -1);
    for (int x = 0; x < q; x++) {
        if (inv_lut[lut[x]] != -1) {
            throw std::invalid_argument("Not implemented for non-bijective functions. Try lat_column.");
        }
        inv_lut[lut[x]] = x;
    }
    return fpt_lat_column(inv_lut, p, m, a);
}

double max_fft(const vector<int>& pol, fftw_plan pl, double *in, const int p) {
    /***
    Returns the maximum squared norm of the FFT of pol's coefficients.
    ***/

    double sol = 0;

    for (int i = 0; i < p; i++) {
        in[i] = pol[i];
    }

    fftw_execute(pl);


    // the output is of the form r0, r1, ... , r(p/2), i((n+1)/2-1), ... , i2, i1.
    // see fftw documentation.
    for (int i = 1; i < (p + 1)/2; i++) {
        sol = std::max(sol,(in[i]*in[i] + in[p - i]*in[p - i]));
    }

    // p/2 is treated specially if p = 2.
    if (p == 2)
        sol = std::max(sol, in[p/2] * in[p/2]);

    return sol;
}

double max_lat_diagonal(const vector<int>& lut, const int p, const int m, const int q, const int q_div_p, const int b, const fftw_plan pl, double* in) {
    // Now to compute the "polynomial" tied to all a's for a given b and retrieve the max walsh coeff.
    double sol = 0;

    if (m == 1) {
        for (int a = 0; a < p; a++) {
            vector<int> temp_poly = vector<int>(q);
            for (int x = 0; x < p; x++) {
                int expo = (a*x - (b*lut[x])%p + p)%p;
                temp_poly[expo]++;
            }
            sol = std::max(sol, max_fft(temp_poly, pl, in, p));
        }
    }

    else {
        for (int a_0 = 0; a_0 < p; a_0++) {
            vector<int> array_of_polys = vector<int>(q);
            vector<int> temp_poly;
            mid_compute(p, m, array_of_polys, q, lut, b, a_0, 0);
            // iterate over all a's of fixed a_0
            for (int a_mm1 = 0; a_mm1 < p; a_mm1++) {
                for (int a_mm2_to_1 = 0; a_mm2_to_1 < q_div_p; a_mm2_to_1 += p) {
                    temp_poly = vector<int>(p);
                    sum_of_polys(p, array_of_polys, temp_poly, a_mm1, a_mm2_to_1);
                    sol = std::max(sol, max_fft(temp_poly, pl, in, p));
                }
            }
        }
    }

    return sol;    
}

double fpt_max_lat(const vector<int>& lut, const int p, const int m, const unsigned int n_threads) {
    /***
    Returns the linearity of the function whose lookup table is in lut.
    ***/

    fftw_init_threads();
    fftw_plan_with_nthreads(n_threads);

    const int q = my_pow(p, m);
    const int q_div_p = q/p;
    
    vector<int> temp_poly;
    vector<int> array_of_polys; // This vector will contain q/p polynomials of degree p.

    double sol = 0;

    omp_set_num_threads(n_threads);
    #pragma omp parallel private(temp_poly, array_of_polys) firstprivate(lut, p, m, q, q_div_p) reduction(max:sol)
    {

    double* in;
    fftw_plan pl;
    #pragma omp critical
    {
    in = fftw_alloc_real(p);
    pl = fftw_plan_r2r_1d(p, in, in, FFTW_R2HC, FFTW_MEASURE);
    }

    double thread_sol = 0;
    int first_nonzero_index;
    int b_temp;

    #pragma omp for schedule(dynamic)
    for(int b_right = 0; b_right < q_div_p; b_right++) {
        b_temp = b_right;
        first_nonzero_index = 0;
        while (b_temp != 0) {
            b_temp /= p;
            first_nonzero_index++;
        }
        b_temp = my_pow(p, first_nonzero_index) + b_right;
        while (b_temp < q) {
            thread_sol = std::max(thread_sol, max_lat_diagonal(lut, p, m, q, q_div_p, b_temp, pl, in));
            first_nonzero_index++;
            b_temp = b_right + my_pow(p, first_nonzero_index);
        }
    }

    #pragma omp critical
    {
    fftw_free(in);
    fftw_free(pl);
    }

    sol = std::max(sol, thread_sol);

    }

    fftw_cleanup_threads();

    return sqrt(sol);
}

void add_to_spectrum(map<double,int>& spectrum, const double key, const int value, const double epsilon) {
    // See if there is an element between key - epsilon and key + epsilon 
    auto low = spectrum.lower_bound(key - epsilon);
    auto high = spectrum.upper_bound(key + epsilon);
    if (low == high) {
        // Case where this is a new key
        spectrum[key] = value;
    }
    else {
        // Case where there already is a (key, int) in the spectrum, pointed to by low
        (*low).second += value;
    }
}

void spectrum_fft(const vector<int>& pol, map<double,int>& spectrum, fftw_plan pl, double *in, const int p, const double epsilon) {
    for (int i = 0; i < p; i++) {
        in[i] = pol[i];
    }

    fftw_execute(pl);

    // the output is of the form r0, r1, ... , r(p/2), i((n+1)/2-1), ... , i2, i1.
    // see fftw documentation.
    for (int i = 1; i < (p + 1)/2; i++) {
        add_to_spectrum(spectrum, sqrt(in[i]*in[i] + in[p - i]*in[p - i]), 2, epsilon);
    }

    // p/2 is treated specially if p = 2.
    if (p == 2)
        add_to_spectrum(spectrum, std::abs(in[p/2]), 1, epsilon);
}

void spectrum_lat_diagonal(const vector<int>& lut, map<double,int>& spectrum, const int p, const int m, const int q, const int q_div_p, const int b, const fftw_plan pl, double* in, const double epsilon) {
    if (m == 1) {
        for (int a = 0; a < p; a++) {
            vector<int> temp_poly = vector<int>(q);
            for (int x = 0; x < p; x++) {
                int expo = (a*x - (b*lut[x])%p + p)%p;
                temp_poly[expo]++;
            }
            spectrum_fft(temp_poly, spectrum, pl, in, p, epsilon);
        }
    }

    else {
        for (int a_0 = 0; a_0 < p; a_0++) {
            vector<int> array_of_polys = vector<int>(q);
            vector<int> temp_poly;
            mid_compute(p, m, array_of_polys, q, lut, b, a_0, 0);
            // iterate over all a's of fixed a_0
            for (int a_mm1 = 0; a_mm1 < p; a_mm1++) {
                for (int a_mm2_to_1 = 0; a_mm2_to_1 < q_div_p; a_mm2_to_1 += p) {
                    temp_poly = vector<int>(p);
                    sum_of_polys(p, array_of_polys, temp_poly, a_mm1, a_mm2_to_1);
                    spectrum_fft(temp_poly, spectrum, pl, in, p, epsilon);
                }
            }
        }
    }
}

map<double, int> fpt_walsh_spectrum(const vector<int>& lut, const int p, const int m, const double epsilon, const unsigned int n_threads) {
    /***
    Returns the linearity of the function whose lookup table is in lut.
    ***/

    fftw_init_threads();
    fftw_plan_with_nthreads(n_threads);

    const int q = my_pow(p, m);
    const int q_div_p = q/p;
    
    vector<int> temp_poly;
    vector<int> array_of_polys; // This vector will contain q/p polynomials of degree p.

    map<double,int> sol;
    sol[0.0] = (q - 1);    // We will skip the first column, so we count it now.
    sol[q] = 1;

    omp_set_num_threads(n_threads);
    #pragma omp parallel private(temp_poly, array_of_polys) firstprivate(lut, p, m, q, q_div_p)
    {

    map<double,int> thread_spectrum;
    double* in;
    fftw_plan pl;
    #pragma omp critical
    {
    in = fftw_alloc_real(p);
    pl = fftw_plan_r2r_1d(p, in, in, FFTW_R2HC, FFTW_MEASURE);
    }

    int first_nonzero_index;
    int b_temp;

    #pragma omp for schedule(dynamic)
    for(int b_right = 0; b_right < q_div_p; b_right++) {
        b_temp = b_right;
        first_nonzero_index = 0;
        while (b_temp != 0) {
            b_temp /= p;
            first_nonzero_index++;
        }
        b_temp = my_pow(p, first_nonzero_index) + b_right;
        while (b_temp < q) {
            spectrum_lat_diagonal(lut, thread_spectrum, p, m, q, q_div_p, b_temp, pl, in, epsilon);
            first_nonzero_index++;
            b_temp = b_right + my_pow(p, first_nonzero_index);
        }
    }

    #pragma omp critical
    {
    fftw_free(in);
    fftw_free(pl);
    }

    // Join all the thread solutions
    for (auto it = thread_spectrum.begin(); it != thread_spectrum.end(); it++) {
        double key = (*it).first;
        int value = (*it).second;
        #pragma omp critical
        {
        add_to_spectrum(sol, key, value, epsilon);
        }
    }

    }

    fftw_cleanup_threads();

    return sol;
}

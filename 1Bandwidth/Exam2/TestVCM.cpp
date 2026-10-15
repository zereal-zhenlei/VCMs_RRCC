#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]


///// Escanciano (2006)
double angle_xy(rowvec x, rowvec y) {
    double xnorm = norm(x);
    double ynorm = norm(y);
    if(xnorm == 0.0 || ynorm == 0.0) {
        return 0.0;
    } else {
        double cxy = norm_dot(x, y);
        if(cxy >= 1.0) cxy = 1.0;
        if(cxy <= -1.0) cxy = -1.0;
        return acos(cxy);
  }
}

//[[Rcpp::export]]
mat A_Esca(mat X) {
    const double pi = 3.14159265359;
    int n = X.n_rows, p = X.n_cols; 
    mat A = zeros(n, n), Xi, Xj, Xl;
    vec sups;
    IntegerVector index = seq(0, n - 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            Xi = X.row(i);
            Xj = X.row(j);
            sups = sapply(index, [&](double val) { Xl = X.row(val); 
                return pi - angle_xy(Xi - Xl, Xj - Xl); });
            A(i, j) = sum(sups);
            A(j, i) = A(i, j);
      }
    }
    return A * pow(pi, p / 2.0 - 1.0) / tgamma(p / 2.0 + 1.0);
}

// [[Rcpp::export]]
vec indexUvec(vec x,uvec inds){
    int n=x.size();
    vec newVec = zeros<vec>(n);
    for (int i = 0; i < n; i++){
        newVec(i) = x(inds(i));
    }
    return newVec;
}

// [[Rcpp::export]]
vec Weight_func_cpp(vec y, vec di){
    int n=y.size();
    vec eta = 1 - di, t = y;
    uvec index = sort_index(y);
    vec Y = indexUvec(y, index),Eta = indexUvec(eta, index), id = linspace(1, n, n), P = linspace(0, n  - 1, n);
    for(int i = 0; i < n; i++){
        P(i) = pow(1.0 - 1.0 / (id(n - 1 - i)), Eta(i));
    }
    vec Gn=linspace(0, n-1, n), wei =linspace(0, n-1, n);
    double Tn = 1;
    for(int i = 0; i < n; i++){
        Tn=1;
        for(int j = 0; j < n; j++){
        if(Y(j) <= t(i)) Tn = Tn * P(j);
        }
        if(Tn == 0) wei(i) = 0; //Gn(i) = 1e-6;//1.0 / n / n / n;
        else wei(i) = di(i) / Tn; //Gn(i)=Tn;//
    }
    return(wei);
}


//[[Rcpp::export]]
double Ker_Epa(double x, double h){
    double tmp = 0.0;
    tmp = (std::abs(x) <= h) * 3.0 * (1 - pow(x / h, 2.0)) / 4.0 / h; // std:: is crucial
    return tmp;
}

// // [[Rcpp::export]]
// mat U_dist_mat(vec U, vec u, double h) {
//     int n = U.size(), m = u.size();
//     mat a = zeros(n, m);
//     for(int i = 0; i < n; i++){
//         for(int j = 0; j < m; j++){
//         a(i, j) = (U(i) - u(j)) / h;
// //         a(j, i) = -a(i, j); // 非对称
//         }
//     }
//     return a;
// }

// [[Rcpp::export]]
mat U_dist_mat(vec U, vec u, double h) {
    int n = U.size(), m = u.size();
    mat A = zeros(n, m);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
        A(i, j) = (U(i) - u(j)) / h;
//         a(j, i) = -a(i, j); // 非对称
        }
    }
    return A;
}

// [[Rcpp::export]]
mat Ker_dist_mat(vec U, vec u, double h) {
    int n = U.size(), m = u.size();
    mat A = ones(n, m) * Ker_Epa(0.0, h); 
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            A(i, j) = Ker_Epa(U(i) - u(j), h);
            // A(j, i) = A(i, j); // 对称
        }
    }
    return A;
}

// [[Rcpp::export]]
double get_h(vec U, double h_con = 2.34, double n_con = 5.0) {
    int n = U.size();
    double h = h_con * pow(n, -1.0 / n_con) * stddev(U);
    return h;
}
// // [[Rcpp::export]]
// double get_h(vec U, double h_con = 2.34, double n_con = 5.0) {
//     int n = U.size();
//     const double qnorm4 = 0.67448975 * 2;
//     vec stds(3, fill::ones), quan;
//     stds(0) = stddev(U);
//     quan = NumericVector::create(0.25, 0.75);
//     vec Qr = quantile(U, quan);
//     stds(1) = (Qr(1) - Qr(0)) / qnorm4;
//     stds(2) = Qr(1);
//     double h = h_con * pow(n, -1.0 / n_con) * min(abs(stds));
//     return h;
// }

// vec Mammencpp(int n){
//     srand((unsigned)time(NULL)); 
//     double p = (sqrt(5.0) + 1.0) / (2.0 * sqrt(5.0));
//     vec v = zeros<vec>(n) - (sqrt(5.0) - 1.0) / 2.0;
//     vec b = randu<vec>(n);
//     for(int i = 0; i < n; i++){
//         if(b(i) > p) v(i) = (sqrt(5) + 1) / 2;
//     }
//     return(v); 
// }
// [[Rcpp::export]]
mat MammencppMat(int n, int B = 100) {
    mat V = zeros(n, B) - (sqrt(5.0) - 1.0) / 2.0;
    mat Boots = randu(n, B); //, distr_param(0, 1)
    double p = (sqrt(5.0) + 1.0) / (2.0 * sqrt(5.0));
    V = V % conv_to<mat>::from(Boots <= p);
    V.replace(0.0, (sqrt(5.0) + 1.0) / 2.0);
    return(V); 
}

// [[Rcpp::export]]
mat Boot_VcM(vec Y, mat X, vec T1, vec delta, vec U, vec alpha, double h_con, double n_con = 3.0, int B = 10) {
    int n = Y.size(), p = X.n_cols;
    vec wei, yhat1, resid1, resid1_w;
    double Tns = 0, TnBs = 0, h = get_h(U, h_con, n_con);
    mat TKer, Uu, XT, S1, W1, theta1, A_E, Ystars, result(2, 1);
    wei = Weight_func_cpp(T1, delta);
    TKer = Ker_dist_mat(U, U, h);
    Uu = U_dist_mat(U, U, h);
    XT = X; XT.insert_cols(p, T1);
    theta1 = zeros(n, 2 * (p + 1));
    for (int i = 0; i < n; i++) {
        S1 = XT;
        S1.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XT);
        W1 = diagmat(wei % TKer.col(i));
        theta1.row(i) = (pinv(S1.t() * W1 * S1) * S1.t() * W1 * Y).t();
    }   
    yhat1 = sum(XT % theta1.cols(0, p), 1);
    resid1 = Y - yhat1; 
    A_E = XT; A_E.insert_cols(p + 1, U); A_E = A_Esca(A_E);
    resid1_w  = wei % resid1; // + (1.0 - wei) * sum(resid1 % delta) / sum(delta);
    Tns = as_scalar(resid1_w.t() * A_E * resid1_w / n / n);

    vec yhat1_ast, resid1b, resid1b_w, Tn_ast(B, fill::zeros);
    // mat Ystars;
    Ystars = repmat(yhat1, 1, B) + repmat(resid1, 1, B) % MammencppMat(n, B);
    for(int b = 0; b < B; b++) {
        yhat1_ast = Ystars.col(b);
        for (int i = 0; i < n; i++) {
            S1 = XT;
            S1.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XT);
            W1 = diagmat(wei % TKer.col(i));
            theta1.row(i) = (pinv(S1.t() * W1 * S1) * S1.t() * W1 * yhat1_ast).t();
        }   
        resid1b = yhat1_ast - sum(XT % theta1.cols(0, p), 1);        
        resid1b_w  = wei % resid1b; // + (1 - wei) * sum(resid1b % delta) / sum(delta);
        Tn_ast(b) = as_scalar(resid1b_w.t() * A_E * resid1b_w / n / n);
   }
    // critical value judgment
    TnBs = as_scalar(1 * (Tns > as_scalar(quantile(Tn_ast, 1 - alpha))));
    result(0, 0) = Tns;
    result(1, 0) = TnBs;
    return result;
}

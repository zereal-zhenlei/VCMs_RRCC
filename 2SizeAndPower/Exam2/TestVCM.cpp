#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;

//[[Rcpp::plugins(openmp)]]
//[[Rcpp::depends(RcppArmadillo)]]

mat A_rand(mat X){
    int n = X.n_rows, m = X.n_cols;
    mat A = zeros(n, n), Xi, Xj, Xl;
    IntegerVector index = seq(0, n - 1);
    vec sups;
    for (int i = 0; i < n; i++){
        for (int j = 0; j <= i; j++){
            Xi = X.row(i);
            Xj = X.row(j);
            sups = sapply(index, [&](double val) { Xl = X.row(val); 
                return as_scalar(conv_to<vec>::from(Xi <= Xl).t() * conv_to<vec>::from(Xj <= Xl)); });
            A(i, j) = sum(sups) / m; // average for m and sum for l
        }
    }
    return A + A.t() - diagmat(A);
}

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


mat Indicator_weight(mat X, vec U){
    int n = X.n_rows;
    mat A = zeros(n, n);
    mat Xi, Xj;
    vec Xi_pr, Xj_pr;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            Xi = repmat(X.row(i), n, 1);
            Xj = repmat(X.row(j), n, 1);
            Xi_pr = prod(conv_to<mat>::from(Xi <= X), 1);
            Xj_pr = prod(conv_to<mat>::from(Xj <= X), 1);
            A(i, j) = mean(Xi_pr % (U[i] <= U) % Xj_pr % (U[j] <= U));
        }
    }
    return A;
}


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



double Ker_Epa(double x, double h){
    double tmp = 0.0;
    tmp = (std::abs(x) <= h) * 3.0 * (1 - pow(x / h, 2.0)) / 4.0 / h; // std:: is crucial
    return tmp;
}

mat X_dist_mat(vec X, double h) {
    int n = X.size();
    mat a = zeros(n, n);
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
        a(i, j) = (X(i) - X(j)) / h;
        a(j, i) = -a(i, j); // 非对称
        }
    }
    return a;
}

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

vec Mammencpp(int n){
    srand((unsigned)time(NULL)); 
    double p = (sqrt(5.0) + 1.0) / (2.0 * sqrt(5.0));
    vec v = zeros<vec>(n) - (sqrt(5.0) - 1.0) / 2.0;
    vec b = randu<vec>(n);
    for(int i = 0; i < n; i++){
        if(b(i) > p) v(i) = (sqrt(5) + 1) / 2;
    }
    return(v); 
}

mat MammencppMat(int n, int B = 100) {
    mat V = zeros(n, B) - (sqrt(5.0) - 1.0) / 2.0;
    mat Boots = randu(n, B); //, distr_param(0, 1)
    double p = (sqrt(5.0) + 1.0) / (2.0 * sqrt(5.0));
    V = V % conv_to<mat>::from(Boots <= p);
    V.replace(0.0, (sqrt(5.0) + 1.0) / 2.0);
    return(V); 
}

// [[Rcpp::export]]
double get_h(vec U, double h_con = 2.34, double n_con = 5.0) {
    int n = U.size();
    double h = h_con * pow(n, -1.0 / n_con) * stddev(U);
    return h;
}
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

// [[Rcpp::export]]
vec Test_VCM(vec Y, mat X, vec T1, vec delta, vec U, vec V, vec VE, vec VW, int B, vec alpha, double h_con = 2.34, double n_con = 3.0, int m = 50) {
    int stats = 9, n = Y.size(), p = X.n_cols;
    double h, hc;
    vec wei, YC, UC, TC, Tns(stats, fill::zeros), TnBs(stats, fill::zeros);
    mat TKer, TKerC, Uu, UuC, XV, XT, XC, XN, XE, XW, XI;
    mat S0, S1, S2, S3, S4, S5, S6, W0, W1, WC;
    mat theta0, theta1, theta2, theta3, theta4, theta5, theta6, varphiN, X_varphiN, resid1_w_t, resid1b_w_t, reBs(B, stats, fill::zeros);
    h = get_h(U, h_con, n_con);
    wei = Weight_func_cpp(T1, delta);
    TKer = Ker_dist_mat(U, U, h);
    Uu = U_dist_mat(U, U, h);
    // 0. Ideal /Oracle
    XV = X; XV.insert_cols(p, V);
    // 1. Proposed
    XT = X; XT.insert_cols(p, T1);
    // 2. CC
    uvec ind = find(delta == 1);
    YC = Y.rows(ind); 
    XC = X.rows(ind); 
    TC = T1.rows(ind);
    UC = U.rows(ind);
    hc = get_h(UC, h_con, n_con);
    TKerC = Ker_dist_mat(UC, UC, hc);
    UuC = U_dist_mat(UC, UC, hc);
    XC.insert_cols(p, TC);
    int n2 = YC.size();
    // 3. Naive
    XN = X; XN.insert_cols(p, T1); // XT
    // 4. Impute-CMI-exp
    XE = X; XE.insert_cols(p, VE);
    // 5. Impute-CMI-wei
    XW = X; XW.insert_cols(p, VW);
    // 6. Impute-Mean
    XI = X; XI.insert_cols(p, T1 % delta + mean(TC) * (1.0 - delta));

    theta0 = zeros(n, 2 * (p + 1));
    theta1 = zeros(n, 2 * (p + 1));
    theta2 = zeros(n2, 2 * (p + 1));
    theta3 = zeros(n, 2 * (p + 1));
    theta4 = zeros(n, 2 * (p + 1));
    theta5 = zeros(n, 2 * (p + 1));
    theta6 = zeros(n, 2 * (p + 1));
    for (int i = 0; i < n; i++) {
        S0 = XV;
        S1 = XT;
        S3 = XN;
        S4 = XE;
        S5 = XW;
        S6 = XI;
        S0.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XV);
        S1.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XT);
        S3.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XN);
        S4.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XE);
        S5.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XW);
        S6.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XI);
        W0 = diagmat(TKer.col(i));
        W1 = diagmat(wei % TKer.col(i));
        theta0.row(i) = (pinv(S0.t() * W0 * S0) * S0.t() * W0 * Y).t();
        theta1.row(i) = (pinv(S1.t() * W1 * S1) * S1.t() * W1 * Y).t();
        theta3.row(i) = (pinv(S3.t() * W0 * S3) * S3.t() * W0 * Y).t();
        theta4.row(i) = (pinv(S4.t() * W0 * S4) * S4.t() * W0 * Y).t();
        theta5.row(i) = (pinv(S5.t() * W0 * S5) * S5.t() * W0 * Y).t();
        theta6.row(i) = (pinv(S6.t() * W0 * S6) * S6.t() * W0 * Y).t();
    }   
    for (int i = 0; i < n2; i++) {
        S2 = XC;
        WC = diagmat(TKerC.col(i));
        S2.insert_cols(p + 1, repmat(UuC.col(i), 1, p + 1) % XC);
        theta2.row(i) = (pinv(S2.t() * WC * S2) * S2.t() * WC * YC).t();
    }  
    vec yhat0, yhat1, yhat2, yhat3, yhat4, yhat5, yhat6, resid0, resid1, resid2, resid3, resid4, resid5, resid6;
    yhat0 = sum(XV % theta0.cols(0, p), 1);
    yhat1 = sum(XT % theta1.cols(0, p), 1);
    yhat2 = sum(XC % theta2.cols(0, p), 1);
    yhat3 = sum(XN % theta3.cols(0, p), 1);
    yhat4 = sum(XE % theta4.cols(0, p), 1);
    yhat5 = sum(XW % theta5.cols(0, p), 1);
    yhat6 = sum(XI % theta6.cols(0, p), 1);
    resid0 = Y - yhat0;
    resid1 = Y - yhat1; 
    resid2 = YC - yhat2; 
    resid3 = Y - yhat3; 
    resid4 = Y - yhat4; 
    resid5 = Y - yhat5; 
    resid6 = Y - yhat6; 
    // vec wei2 = Gn_hats(T, delta, T);
    mat A0, A1, A2, A3, A4, A5, A6, A8, supsK(n, m);
    A0 = XV; A0.insert_cols(p + 1, U); A0 = A_Esca(A0);
    A1 = XT; A1.insert_cols(p + 1, U); A1 = A_Esca(A1);
    A2 = XC; A2.insert_cols(p + 1, UC); A2 = A_Esca(A2);
    // A3 = XN; A3.insert_cols(p + 1, U); A3 = A_Esca(A3);
    A3 = A1;
    A4 = XE; A4.insert_cols(p + 1, U); A4 = A_Esca(A4);
    A5 = XW; A5.insert_cols(p + 1, U); A5 = A_Esca(A5);
    A6 = XI; A6.insert_cols(p + 1, U); A6 = A_Esca(A6);
    // vec resid1_w  = wei % resid1 + (1.0 - wei) * sum(resid1 % delta) / sum(delta);
    vec resid1_w  = wei % resid1;// + (1.0 - wei) * sum(resid1 % delta) / n2;
    resid1_w_t = resid1_w.t();
    Tns(0) = as_scalar(resid0.t() * A0 * resid0 / n / n);
    Tns(1) = as_scalar(resid1_w_t * A1 * resid1_w / n / n);
    Tns(2) = as_scalar(resid2.t() * A2 * resid2 / n2 / n2);
    Tns(3) = as_scalar(resid3.t() * A3 * resid3 / n / n);
    Tns(4) = as_scalar(resid4.t() * A4 * resid4 / n / n);
    Tns(5) = as_scalar(resid5.t() * A5 * resid5 / n / n);
    Tns(6) = as_scalar(resid6.t() * A6 * resid6 / n / n);

    // Random asy
    /////////////## KS statistics ##/////////////
    varphiN = randn(p + 2, m); //, distr_param(0, 1)); 
    X_varphiN = XT; X_varphiN.insert_cols(p + 1, U); X_varphiN = X_varphiN * varphiN;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            supsK(i, j) = accu(abs(resid1_w_t * conv_to<mat>::from(X_varphiN <= as_scalar(X_varphiN(i, j)))));
        }
    }
    Tns(7) = supsK.max() / (m * sqrt(n));
    A8 = A_rand(X_varphiN);
    Tns(8) = as_scalar(resid1_w_t * A8 * resid1_w / n / n);

    vec yhat0_ast, yhat1_ast, yhat2_ast, yhat3_ast, yhat4_ast, yhat5_ast, yhat6_ast;
    // yhat0_ast = repmat(yhat0, 1, B) + repmat(resid0, 1, B) % MammencppMat(n, B);
    // yhat1_ast = repmat(yhat1, 1, B) + repmat(resid1, 1, B) % MammencppMat(n, B);
    // yhat2_ast = repmat(yhat2, 1, B) + repmat(resid2, 1, B) % MammencppMat(nc, B);
    // yhat3_ast = repmat(yhat3, 1, B) + repmat(resid3, 1, B) % MammencppMat(n, B);
    // yhat4_ast = repmat(yhat4, 1, B) + repmat(resid4, 1, B) % MammencppMat(n, B);
    // yhat5_ast = repmat(yhat5, 1, B) + repmat(resid5, 1, B) % MammencppMat(n, B);
    vec resid0b, resid1b, resid1b_w, resid2b, resid3b, resid4b, resid5b, resid6b;
    for(int b = 0; b < B; b++) {
        yhat0_ast = yhat0 + Mammencpp(n) % resid0;
        yhat1_ast = yhat1 + Mammencpp(n) % resid1;
        yhat2_ast = yhat2 + Mammencpp(n2) % resid2; //Mammencpp(n).rows(ind) % resid2;
        yhat3_ast = yhat3 + Mammencpp(n) % resid3;
        yhat4_ast = yhat4 + Mammencpp(n) % resid4;
        yhat5_ast = yhat5 + Mammencpp(n) % resid5;
        yhat6_ast = yhat6 + Mammencpp(n) % resid6;

        for (int i = 0; i < n; i++) {
            S0 = XV;
            S1 = XT;
            S3 = XN;
            S4 = XE;
            S5 = XW;
            S6 = XI;
            S0.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XV);
            S1.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XT);
            S3.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XN);
            S4.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XE);
            S5.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XW);
            S6.insert_cols(p + 1, repmat(Uu.col(i), 1, p + 1) % XI);
            W0 = diagmat(TKer.col(i));
            W1 = diagmat(wei % TKer.col(i));
            theta0.row(i) = (pinv(S0.t() * W0 * S0) * S0.t() * W0 * yhat0_ast).t();
            theta1.row(i) = (pinv(S1.t() * W1 * S1) * S1.t() * W1 * yhat1_ast).t();
            theta3.row(i) = (pinv(S3.t() * W0 * S3) * S3.t() * W0 * yhat3_ast).t();
            theta4.row(i) = (pinv(S4.t() * W0 * S4) * S4.t() * W0 * yhat4_ast).t();
            theta5.row(i) = (pinv(S5.t() * W0 * S5) * S5.t() * W0 * yhat5_ast).t();
            theta6.row(i) = (pinv(S6.t() * W0 * S6) * S6.t() * W0 * yhat6_ast).t();
        }   
        for (int i = 0; i < n2; i++) {
            S2 = XC;
            WC = diagmat(TKerC.col(i));
            S2.insert_cols(p + 1, repmat(UuC.col(i), 1, p + 1) % XC);
            theta2.row(i) = (pinv(S2.t() * WC * S2) * S2.t() * WC * yhat2_ast).t();
        }  
        resid0b = yhat0_ast - sum(XV % theta0.cols(0, p), 1);
        resid1b = yhat1_ast - sum(XT % theta1.cols(0, p), 1);
        resid2b = yhat2_ast - sum(XC % theta2.cols(0, p), 1);
        resid3b = yhat3_ast - sum(XN % theta3.cols(0, p), 1);
        resid4b = yhat4_ast - sum(XE % theta4.cols(0, p), 1);  
        resid5b = yhat5_ast - sum(XW % theta5.cols(0, p), 1);   
        resid6b = yhat6_ast - sum(XI % theta6.cols(0, p), 1);  
        resid1b_w  = wei % resid1b; // + (1 - wei) * sum(resid1b % delta) / n2;
        resid1b_w_t = resid1b_w.t();
        reBs(b, 0) = as_scalar(resid0b.t() * A0 * resid0b / n / n);
        reBs(b, 1) = as_scalar(resid1b_w_t * A1 * resid1b_w / n / n);
        reBs(b, 2) = as_scalar(resid2b.t() * A2 * resid2b / n2 / n2);
        reBs(b, 3) = as_scalar(resid3b.t() * A3 * resid3b / n / n);
        reBs(b, 4) = as_scalar(resid4b.t() * A4 * resid4b / n / n);
        reBs(b, 5) = as_scalar(resid5b.t() * A5 * resid5b / n / n);
        reBs(b, 6) = as_scalar(resid6b.t() * A6 * resid6b / n / n);
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < m; j++) {
                supsK(i, j) = accu(abs(resid1b_w_t * conv_to<mat>::from(X_varphiN <= as_scalar(X_varphiN(i, j)))));
            }
        }
        reBs(b, 7) = supsK.max() / (m * sqrt(n));
        reBs(b, 8)= as_scalar(resid1b_w_t * A8 * resid1b_w / n / n);
    }
    // critical value judgment
    for(int i = 0; i < stats; i++) {
        TnBs(i) = as_scalar(1.0 * (Tns(i) > as_scalar(quantile(reBs.col(i), 1.0 - alpha))));
    }
    return TnBs;
}

//     List result = List::create(
//         Named("Tns") = Tns,
//         Named("TnBs") = TnBs,
//         Named("XI") = XI

//     );
//     return result;
// }




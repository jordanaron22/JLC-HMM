//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double CalcCondMeanC(double mu1, double sig1, double mu2, double sig2, double bivar_corr, double obs2){
	double cond_mu1 = mu1 + bivar_corr*(sig1/sig2)*(obs2-mu2);
	return (cond_mu1);
}

// [[Rcpp::export]]
double CalcCondSigC(double sig1, double bivar_corr){
	double cond_sig1 = sig1*sqrt(1-pow(bivar_corr,2));
	return cond_sig1;
}

// [[Rcpp::export]]
NumericVector vectorEqBool(NumericVector vec, double lod) {
  NumericVector vecBool;
  for(NumericVector::iterator i = vec.begin(); i != vec.end(); ++i) {
    if (*i == lod) {
      vecBool.push_back(1);
    } else {
      vecBool.push_back(0);
    }
  }
  return vecBool;
}

// [[Rcpp::export]]
vec logClassificationCnonTobit(NumericVector act_obs, NumericVector light_obs, double mu_act, double sig_act, double mu_light, double sig_light,
                               double lod_act, double lod_light, double lambda_act, double lambda_light) {
  vec temp_act;
  vec temp_light;

  NumericVector vec_eq_act = vectorEqBool(act_obs, lod_act);
  NumericVector vec_eq_light = vectorEqBool(light_obs, lod_light);

  temp_act = ((log(1-lambda_act) + Rcpp::dnorm( act_obs, mu_act, sig_act, true )) * (1 - vec_eq_act)) + (log(lambda_act)*vec_eq_act);
  temp_light = ((log(1-lambda_light) + Rcpp::dnorm( light_obs, mu_light, sig_light, true )) * (1 - vec_eq_light)) + (log(lambda_light)*vec_eq_light);

  LogicalVector nan_vec_act = is_na(act_obs);
  LogicalVector nan_vec_light = is_na(light_obs);
  int n = act_obs.length();

  for (int i = 0; i < n; i++) {
    if(nan_vec_act[i] == true){
      temp_act[i] = 0;
    }

    if(nan_vec_light[i] == true){
      temp_light[i] = 0;
    }
  }

  return temp_act + temp_light;
}


// [[Rcpp::export]]
vec logClassificationCTobit(NumericVector act, NumericVector light, double mu_act, double sig_act, double mu_light, double sig_light,
                            double lod_act, double lod_light, double bivar_corr, double lintegral) {

	int n = act.length();
	vec log_dens_vec(n);

	LogicalVector nan_vec_act = is_na(act);
	LogicalVector nan_vec_light = is_na(light);

	for (int i = 0; i < n; i++) {
		double act_obs = act(i);
		double light_obs = light(i);

		// Observe both act and light
		if (nan_vec_act(i) != true && nan_vec_light(i) != true){
			// CASE 1
			if(act_obs > lod_act && light_obs > lod_light){

				double mu_light_cond = CalcCondMeanC(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs);
				double sig_light_cond = CalcCondSigC(sig_light,bivar_corr);

				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true ) + R::dnorm( light_obs, mu_light_cond, sig_light_cond, true ) ;
			}

			//CASE 2
			if(act_obs > lod_act && light_obs <= lod_light){

				double mu_light_cond = CalcCondMeanC(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs);
				double sig_light_cond = CalcCondSigC(sig_light,bivar_corr);

				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true ) + R::pnorm( light_obs, mu_light_cond, sig_light_cond, true, true ) ;
			}

			//CASE 3
			if(act_obs <= lod_act && light_obs > lod_light){

				double mu_act_cond = CalcCondMeanC(mu_act,sig_act,mu_light,sig_light,bivar_corr,light_obs);
				double sig_act_cond = CalcCondSigC(sig_act,bivar_corr);

				log_dens_vec(i) = R::pnorm( act_obs, mu_act_cond, sig_act_cond, true, true ) + R::dnorm( light_obs, mu_light, sig_light, true ) ;
			}

			//CASE 4
			if(act_obs <= lod_act && light_obs <= lod_light){
				log_dens_vec(i) = lintegral;
			}

			// Light Missing
		} else if (nan_vec_act(i) != true && nan_vec_light(i) == true){
			if(act_obs > lod_act){
				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true );
			} else {
				log_dens_vec(i) = R::pnorm( act_obs, mu_act, sig_act, true, true );
			}

			// Act missing
		} else if((nan_vec_act(i) == true && nan_vec_light(i) != true)){
			if(light_obs > lod_light){
				log_dens_vec(i) = R::dnorm( light_obs, mu_light, sig_light, true );
			} else {
				log_dens_vec(i) = R::pnorm( light_obs, mu_light, sig_light, true, true );
			}

			// Both Missing
		} else {
			log_dens_vec(i) = 0;
		}


	}

return log_dens_vec;

}

// [[Rcpp::export]]
vec logClassificationC(NumericVector act, NumericVector light, double mu_act, double sig_act, double mu_light, double sig_light,
                       double lod_act, double lod_light, double bivar_corr, double lintegral, double lambda_act, double lambda_light,
                       bool tobit) {

  vec log_class;

  if (tobit){
    log_class = logClassificationCTobit(act, light, mu_act, sig_act, mu_light, sig_light,
                                        lod_act, lod_light, bivar_corr, lintegral);
  } else {
    log_class = logClassificationCnonTobit(act, light, mu_act, sig_act, mu_light, sig_light,
                                           lod_act, lod_light, lambda_act, lambda_light);
  }
  return log_class;
}




/* This is from the seqHMM github*/
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExpC(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  LDOUBLE maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  LDOUBLE cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += EXPL(x(i) - maxv);
    }
  }

  return maxv + log1p(cumsum);
}

inline double LogSumExp2C(
    const double a,
    const double b
) {
  if (a == R_NegInf) {
    return b;
  }

  if (b == R_NegInf) {
    return a;
  }

  const double maximum =
    std::max(a, b);

  const double minimum =
    std::min(a, b);

  return maximum +
    std::log1p(
      std::exp(minimum - maximum)
    );
}

// [[Rcpp::export]]
mat ForwardIndC(const NumericVector& act_ind, const NumericVector& light_ind, NumericVector init, Rcpp::List tran_list,
                cube emit_act_week, cube emit_light_week,  cube emit_act_weekend, cube emit_light_weekend,
                int clust_i, double lod_act, double lod_light, cube corr_vec, vec beta_vec, double covar_risk, int event, double bline, double cbline, cube lintegral_mat, vec vcovar_vec,
                cube lambda_act_mat, cube lambda_light_mat, bool tobit, bool incl_surv){

	mat alpha( act_ind.length(), 2 );

  vec log_class_0_week = logClassificationC( act_ind, light_ind,
                                             emit_act_week(0,0,clust_i),
                                             emit_act_week(0,1,clust_i),
                                             emit_light_week(0,0,clust_i),
                                             emit_light_week(0,1,clust_i),
                                             lod_act, lod_light, corr_vec(clust_i,0,0), lintegral_mat(clust_i,0,0),
                                             lambda_act_mat(clust_i,0,0),lambda_light_mat(clust_i,0,0),tobit);

  vec log_class_1_week = logClassificationC( act_ind, light_ind,
                                             emit_act_week(1,0,clust_i),
                                             emit_act_week(1,1,clust_i),
                                             emit_light_week(1,0,clust_i),
                                             emit_light_week(1,1,clust_i),
                                             lod_act, lod_light, corr_vec(clust_i,1,0), lintegral_mat(clust_i,1,0),
                                             lambda_act_mat(clust_i,1,0),lambda_light_mat(clust_i,1,0),tobit);

  vec log_class_0_weekend = logClassificationC( act_ind, light_ind,
                                                emit_act_weekend(0,0,clust_i),
                                                emit_act_weekend(0,1,clust_i),
                                                emit_light_weekend(0,0,clust_i),
                                                emit_light_weekend(0,1,clust_i),
                                                lod_act, lod_light, corr_vec(clust_i,0,1), lintegral_mat(clust_i,0,1),
                                                lambda_act_mat(clust_i,0,1),lambda_light_mat(clust_i,0,1),tobit);

  vec log_class_1_weekend = logClassificationC( act_ind, light_ind,
                                                emit_act_weekend(1,0,clust_i),
                                                emit_act_weekend(1,1,clust_i),
                                                emit_light_weekend(1,0,clust_i),
                                                emit_light_weekend(1,1,clust_i),
                                                lod_act, lod_light, corr_vec(clust_i,1,1), lintegral_mat(clust_i,1,1),
                                                lambda_act_mat(clust_i,1,1),lambda_light_mat(clust_i,1,1),tobit);


	vec log_class_0 = (log_class_0_week % (1-vcovar_vec)) + (log_class_0_weekend % vcovar_vec);
	vec log_class_1 = (log_class_1_week % (1-vcovar_vec)) + (log_class_1_weekend % vcovar_vec);

	double surv_comp = 0.0;

	if (incl_surv){
	  if (event == 1){
	    surv_comp = log(bline) + beta_vec[clust_i]+covar_risk - (cbline * exp(beta_vec[clust_i]+(covar_risk)));
	  } else {
	    surv_comp =  -cbline * exp(beta_vec[clust_i]+covar_risk);
	  }
	}




	alpha(0,0) = log(init(0)) + log_class_0[0] + surv_comp;
	alpha(0,1) = log(init(1)) + log_class_1[0] + surv_comp;

	List tran_list_clust = tran_list[clust_i];

	for (int i = 1; i < act_ind.length(); i++) {
	  List tran_list_vcovar = tran_list_clust[vcovar_vec(i)];

	  NumericMatrix tran = tran_list_vcovar[i%96];

	  double fp_00 = alpha(i-1,0) + log(tran(0,0)) + log_class_0[i];

	  double fp_10 = alpha(i-1,1) + log(tran(1,0)) + log_class_0[i];

	  double fp_01 = alpha(i-1,0) + log(tran(0,1)) + log_class_1[i];

	  double fp_11 = alpha(i-1,1) + log(tran(1,1)) + log_class_1[i];

	  NumericVector fp_0 = NumericVector::create(fp_00,fp_10);
	  NumericVector fp_1 = NumericVector::create(fp_01,fp_11);

	  alpha(i,0) = logSumExpC(fp_0);
	  alpha(i,1) = logSumExpC(fp_1);
	}

	return(alpha);
}

// [[Rcpp::export]]
mat BackwardIndC(const NumericVector& act_ind, const NumericVector& light_ind, Rcpp::List tran_list,
                 cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend,
                 int clust_i, double lod_act, double lod_light, cube corr_vec, cube lintegral_mat, vec vcovar_vec,
                 cube lambda_act_mat, cube lambda_light_mat, bool tobit){

  int n = act_ind.length();
  mat beta( n, 2 );

  vec log_class_0_week = logClassificationC( act_ind, light_ind,
                                             emit_act_week(0,0,clust_i),
                                             emit_act_week(0,1,clust_i),
                                             emit_light_week(0,0,clust_i),
                                             emit_light_week(0,1,clust_i),
                                             lod_act, lod_light, corr_vec(clust_i,0,0), lintegral_mat(clust_i,0,0),
                                             lambda_act_mat(clust_i,0,0),lambda_light_mat(clust_i,0,0),tobit);

  vec log_class_1_week = logClassificationC( act_ind, light_ind,
                                             emit_act_week(1,0,clust_i),
                                             emit_act_week(1,1,clust_i),
                                             emit_light_week(1,0,clust_i),
                                             emit_light_week(1,1,clust_i),
                                             lod_act, lod_light, corr_vec(clust_i,1,0), lintegral_mat(clust_i,1,0),
                                             lambda_act_mat(clust_i,1,0),lambda_light_mat(clust_i,1,0),tobit);

  vec log_class_0_weekend = logClassificationC( act_ind, light_ind,
                                                emit_act_weekend(0,0,clust_i),
                                                emit_act_weekend(0,1,clust_i),
                                                emit_light_weekend(0,0,clust_i),
                                                emit_light_weekend(0,1,clust_i),
                                                lod_act, lod_light, corr_vec(clust_i,0,1), lintegral_mat(clust_i,0,1),
                                                lambda_act_mat(clust_i,0,1),lambda_light_mat(clust_i,0,1),tobit);

  vec log_class_1_weekend = logClassificationC( act_ind, light_ind,
                                                emit_act_weekend(1,0,clust_i),
                                                emit_act_weekend(1,1,clust_i),
                                                emit_light_weekend(1,0,clust_i),
                                                emit_light_weekend(1,1,clust_i),
                                                lod_act, lod_light, corr_vec(clust_i,1,1), lintegral_mat(clust_i,1,1),
                                                lambda_act_mat(clust_i,1,1),lambda_light_mat(clust_i,1,1),tobit);


	vec log_class_0 = (log_class_0_week % (1-vcovar_vec)) + (log_class_0_weekend % vcovar_vec);
	vec log_class_1 = (log_class_1_week % (1-vcovar_vec)) + (log_class_1_weekend % vcovar_vec);

  beta(n-1,0) = log(1);
  beta(n-1,1) = log(1);

  List tran_time_list = tran_list[clust_i];

  for (int i = n-2; i >= 0; i--) {

    List tran_list_vcovar = tran_time_list[vcovar_vec(i+1)];
    NumericMatrix tran = tran_list_vcovar[(i+1)%96];

    double bp_00 = log(tran(0,0)) + log_class_0[i+1] + beta(i+1,0);

    double bp_01 = log(tran(0,1)) + log_class_1[i+1] + beta(i+1,1);

    double bp_10 = log(tran(1,0)) + log_class_0[i+1] + beta(i+1,0);

    double bp_11 = log(tran(1,1)) + log_class_1[i+1] + beta(i+1,1);

    NumericVector bp_0 = NumericVector::create(bp_00,bp_01);
    NumericVector bp_1 = NumericVector::create(bp_10,bp_11);

    beta(i,0) = logSumExpC(bp_0);
    beta(i,1) = logSumExpC(bp_1);

  }

  return beta;

}

/*
// [[Rcpp::export]]
List ForwardC(const NumericMatrix& act, const NumericMatrix& light, NumericMatrix init, List tran_list,
              cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend,
              double lod_act, double lod_light, cube corr_mat, vec beta_vec, double beta_age,
              vec event_vec, vec bline_vec, vec cbline_vec, cube lintegral_mat, NumericVector sweights_vec,
              vec age_vec, mat vcovar_mat,cube lambda_act_mat, cube lambda_light_mat, bool tobit, bool incl_surv){

	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act_week.n_slices;
	List alpha_list(num_people);

	for (int ind = 0; ind < num_people; ind++) {
	  int age_ind = age_vec(ind);

		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		NumericVector light_ind = light.column(ind);
		double sweight = sweights_vec(ind);

		vec vcovar_vec = vcovar_mat.col(ind);

		for (int clust_i = 0; clust_i < num_re; clust_i++){
			NumericVector init_vec = init.row(clust_i);

			Cube1.slice(clust_i) = ForwardIndC(act_ind, light_ind, init_vec, tran_list, emit_act_week, emit_light_week,
               emit_act_weekend, emit_light_weekend,clust_i, lod_act, lod_light, corr_mat,
               beta_vec, beta_age,age_ind,event_vec[ind], bline_vec[ind], cbline_vec[ind], lintegral_mat,
               sweight, vcovar_vec,lambda_act_mat, lambda_light_mat, tobit, incl_surv);

		}

		alpha_list(ind) = Cube1;
	}
	return(alpha_list);
}
*/

// [[Rcpp::export]]
List BackwardC(const NumericMatrix& act, const NumericMatrix& light, List tran_list,
               cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend,
               double lod_act, double lod_light, cube corr_mat, cube lintegral_mat, mat vcovar_mat,
               cube lambda_act_mat, cube lambda_light_mat, bool tobit){

	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act_week.n_slices;
	List beta_list(num_people);

	for (int ind = 0; ind < num_people; ind++) {

		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		NumericVector light_ind = light.column(ind);

		vec vcovar_vec = vcovar_mat.col(ind);

		for (int clust_i = 0; clust_i < num_re; clust_i++){

			Cube1.slice(clust_i) = BackwardIndC(act_ind, light_ind, tran_list, emit_act_week, emit_light_week, emit_act_weekend,
               emit_light_weekend, clust_i, lod_act, lod_light, corr_mat, lintegral_mat, vcovar_vec,
               lambda_act_mat, lambda_light_mat, tobit);

		}

		beta_list(ind) = Cube1;
	}
	return(beta_list);
}

// [[Rcpp::export]]
List ForwardBackwardC(
    const NumericMatrix& act,
    const NumericMatrix& light,
    const NumericMatrix& init,
    const List& tran_list,
    const arma::cube& emit_act_week,
    const arma::cube& emit_light_week,
    const arma::cube& emit_act_weekend,
    const arma::cube& emit_light_weekend,
    const double lod_act,
    const double lod_light,
    const arma::cube& corr_mat,
    const NumericVector& beta_vec,
    const NumericVector& surv_covar_risk_vec,
    const NumericVector& event_vec,
    const NumericVector& bline_vec,
    const NumericVector& cbline_vec,
    const arma::cube& lintegral_mat,
    const arma::mat& vcovar_mat,
    const arma::cube& lambda_act_mat,
    const arma::cube& lambda_light_mat,
    const bool tobit,
    const bool incl_surv,
    const int vcovar_num,
    const int period_len
) {
  const int len =
    act.nrow();

  const int num_people =
    act.ncol();

  const int mix_num =
    emit_act_week.n_slices;

  if (len < 1) {
    stop("act must contain at least one time point");
  }

  if (
    light.nrow() != len ||
    light.ncol() != num_people
  ) {
    stop(
      "act and light dimensions do not match"
    );
  }

  if (
    vcovar_mat.n_rows !=
      static_cast<arma::uword>(len) ||
    vcovar_mat.n_cols !=
      static_cast<arma::uword>(num_people)
  ) {
    stop(
      "vcovar_mat dimensions do not match act"
    );
  }

  if (
    init.nrow() != mix_num ||
    init.ncol() != 2
  ) {
    stop(
      "init must have one row per class and two columns"
    );
  }

  if (tran_list.size() != mix_num) {
    stop(
      "tran_list must have one entry per latent class"
    );
  }

  if (
    beta_vec.size() != mix_num
  ) {
    stop(
      "beta_vec must contain one value per latent class"
    );
  }

  if (
    surv_covar_risk_vec.size() != num_people ||
    event_vec.size() != num_people ||
    bline_vec.size() != num_people ||
    cbline_vec.size() != num_people
  ) {
    stop(
      "Survival vectors must contain one value per participant"
    );
  }

  if (
    vcovar_num < 1 ||
    vcovar_num > 2
  ) {
    stop(
      "ForwardBackwardC currently supports one or two day types"
    );
  }

  if (period_len < 1) {
    stop(
      "period_len must be positive"
    );
  }

  if (
    emit_act_week.n_rows != 2 ||
    emit_act_week.n_cols != 2 ||
    emit_light_week.n_rows != 2 ||
    emit_light_week.n_cols != 2 ||
    emit_act_weekend.n_rows != 2 ||
    emit_act_weekend.n_cols != 2 ||
    emit_light_weekend.n_rows != 2 ||
    emit_light_weekend.n_cols != 2
  ) {
    stop(
      "Emission arrays must have dimensions 2 x 2 x class"
    );
  }

  if (
    corr_mat.n_rows !=
      static_cast<arma::uword>(mix_num) ||
    corr_mat.n_cols != 2 ||
    corr_mat.n_slices <
      static_cast<arma::uword>(vcovar_num)
  ) {
    stop(
      "corr_mat has unexpected dimensions"
    );
  }

  if (
    lintegral_mat.n_rows !=
      static_cast<arma::uword>(mix_num) ||
    lintegral_mat.n_cols != 2 ||
    lintegral_mat.n_slices <
      static_cast<arma::uword>(vcovar_num)
  ) {
    stop(
      "lintegral_mat has unexpected dimensions"
    );
  }

  if (
    lambda_act_mat.n_rows !=
      static_cast<arma::uword>(mix_num) ||
    lambda_act_mat.n_cols != 2 ||
    lambda_act_mat.n_slices <
      static_cast<arma::uword>(vcovar_num) ||
    lambda_light_mat.n_rows !=
      static_cast<arma::uword>(mix_num) ||
    lambda_light_mat.n_cols != 2 ||
    lambda_light_mat.n_slices <
      static_cast<arma::uword>(vcovar_num)
  ) {
    stop(
      "Lambda arrays have unexpected dimensions"
    );
  }

  List alpha_list(num_people);
  List beta_list(num_people);

  /*
   * If there is only one day type, the weekend
   * parameter slice is the same as the weekday slice.
   */
  const int weekend_slice =
    vcovar_num == 2 ? 1 : 0;

  for (int ind = 0;
       ind < num_people;
       ++ind) {

    NumericVector act_ind =
      act.column(ind);

    NumericVector light_ind =
      light.column(ind);

    arma::vec vcovar_vec =
      vcovar_mat.col(ind);

    /*
     * Validate and standardize the day-type values once
     * for this participant.
     */
    arma::vec day_type_vec(
      len,
      arma::fill::zeros
    );

    for (int time_ind = 0;
         time_ind < len;
         ++time_ind) {

      const double day_value =
        vcovar_vec[time_ind];

      if (
        NumericVector::is_na(day_value) ||
        !std::isfinite(day_value)
      ) {
        stop(
          "vcovar_mat contains a missing or nonfinite value"
        );
      }

      int day_type = 0;

      if (vcovar_num == 2) {
        const double rounded_value =
          std::round(day_value);

        if (
          std::abs(
            day_value - rounded_value
          ) > 1e-8
        ) {
          stop(
            "Day-type values must be coded as 0 or 1"
          );
        }

        day_type =
          static_cast<int>(
            rounded_value
          );
      }

      if (
        day_type < 0 ||
        day_type >= vcovar_num
      ) {
        stop(
          "Day-type values must be coded from 0 to vcovar_num - 1"
        );
      }

      day_type_vec[time_ind] =
        static_cast<double>(day_type);
    }

    arma::cube alpha_cube(
      len,
      2,
      mix_num,
      arma::fill::zeros
    );

    arma::cube beta_cube(
      len,
      2,
      mix_num,
      arma::fill::zeros
    );

    const double covar_risk =
      surv_covar_risk_vec[ind];

    for (int clust_i = 0;
         clust_i < mix_num;
         ++clust_i) {

      /*
       * Calculate the four emission vectors once.
       *
       * The same vectors are used by both Forward
       * and Backward.
       */

      arma::vec log_class_0_week =
        logClassificationC(
          act_ind,
          light_ind,
          emit_act_week(0, 0, clust_i),
          emit_act_week(0, 1, clust_i),
          emit_light_week(0, 0, clust_i),
          emit_light_week(0, 1, clust_i),
          lod_act,
          lod_light,
          corr_mat(clust_i, 0, 0),
          lintegral_mat(clust_i, 0, 0),
          lambda_act_mat(clust_i, 0, 0),
          lambda_light_mat(clust_i, 0, 0),
          tobit
        );

      arma::vec log_class_1_week =
        logClassificationC(
          act_ind,
          light_ind,
          emit_act_week(1, 0, clust_i),
          emit_act_week(1, 1, clust_i),
          emit_light_week(1, 0, clust_i),
          emit_light_week(1, 1, clust_i),
          lod_act,
          lod_light,
          corr_mat(clust_i, 1, 0),
          lintegral_mat(clust_i, 1, 0),
          lambda_act_mat(clust_i, 1, 0),
          lambda_light_mat(clust_i, 1, 0),
          tobit
        );

      arma::vec log_class_0_weekend =
        logClassificationC(
          act_ind,
          light_ind,
          emit_act_weekend(0, 0, clust_i),
          emit_act_weekend(0, 1, clust_i),
          emit_light_weekend(0, 0, clust_i),
          emit_light_weekend(0, 1, clust_i),
          lod_act,
          lod_light,
          corr_mat(
            clust_i,
            0,
            weekend_slice
          ),
          lintegral_mat(
            clust_i,
            0,
            weekend_slice
          ),
          lambda_act_mat(
            clust_i,
            0,
            weekend_slice
          ),
          lambda_light_mat(
            clust_i,
            0,
            weekend_slice
          ),
          tobit
        );

      arma::vec log_class_1_weekend =
        logClassificationC(
          act_ind,
          light_ind,
          emit_act_weekend(1, 0, clust_i),
          emit_act_weekend(1, 1, clust_i),
          emit_light_weekend(1, 0, clust_i),
          emit_light_weekend(1, 1, clust_i),
          lod_act,
          lod_light,
          corr_mat(
            clust_i,
            1,
            weekend_slice
          ),
          lintegral_mat(
            clust_i,
            1,
            weekend_slice
          ),
          lambda_act_mat(
            clust_i,
            1,
            weekend_slice
          ),
          lambda_light_mat(
            clust_i,
            1,
            weekend_slice
          ),
          tobit
        );

      arma::vec log_class_0 =
        log_class_0_week;

      arma::vec log_class_1 =
        log_class_1_week;

      if (vcovar_num == 2) {
        log_class_0 =
          (
            log_class_0_week %
            (1.0 - day_type_vec)
          ) +
          (
            log_class_0_weekend %
            day_type_vec
          );

        log_class_1 =
          (
            log_class_1_week %
            (1.0 - day_type_vec)
          ) +
          (
            log_class_1_weekend %
            day_type_vec
          );
      }

      arma::mat alpha_mat(
        len,
        2,
        arma::fill::zeros
      );

      arma::mat beta_mat(
        len,
        2,
        arma::fill::zeros
      );

      /*
       * Survival contribution is class-specific but
       * constant across longitudinal time.
       *
       * Preserve the current implementation by adding
       * it only to the Forward initialization.
       */

      double surv_comp = 0.0;

      if (incl_surv) {
        const double linear_predictor =
          beta_vec[clust_i] +
          covar_risk;

        if (event_vec[ind] == 1.0) {
          surv_comp =
            std::log(bline_vec[ind]) +
            linear_predictor -
            cbline_vec[ind] *
            std::exp(linear_predictor);
        } else {
          surv_comp =
            -cbline_vec[ind] *
            std::exp(linear_predictor);
        }
      }

      alpha_mat(0, 0) =
        std::log(init(clust_i, 0)) +
        log_class_0[0] +
        surv_comp;

      alpha_mat(0, 1) =
        std::log(init(clust_i, 1)) +
        log_class_1[0] +
        surv_comp;

      List tran_list_clust =
        tran_list[clust_i];

      if (
        tran_list_clust.size() <
        vcovar_num
      ) {
        stop(
          "tran_list contains too few day types"
        );
      }

      /*
       * Forward recursion.
       */

      for (int time_ind = 1;
           time_ind < len;
           ++time_ind) {

        const int day_type =
          static_cast<int>(
            day_type_vec[time_ind]
          );

        List tran_list_vcovar =
          tran_list_clust[day_type];

        const int transition_index =
          time_ind % period_len;

        if (
          transition_index >=
          tran_list_vcovar.size()
        ) {
          stop(
            "Transition index is outside tran_list"
          );
        }

        NumericMatrix tran =
          tran_list_vcovar[
            transition_index
          ];

        const double fp_00 =
          alpha_mat(time_ind - 1, 0) +
          std::log(tran(0, 0)) +
          log_class_0[time_ind];

        const double fp_10 =
          alpha_mat(time_ind - 1, 1) +
          std::log(tran(1, 0)) +
          log_class_0[time_ind];

        const double fp_01 =
          alpha_mat(time_ind - 1, 0) +
          std::log(tran(0, 1)) +
          log_class_1[time_ind];

        const double fp_11 =
          alpha_mat(time_ind - 1, 1) +
          std::log(tran(1, 1)) +
          log_class_1[time_ind];

        alpha_mat(time_ind, 0) =
          LogSumExp2C(
            fp_00,
            fp_10
          );

        alpha_mat(time_ind, 1) =
          LogSumExp2C(
            fp_01,
            fp_11
          );
      }

      /*
       * Backward recursion using the same emission
       * vectors calculated above.
       */

      beta_mat(len - 1, 0) = 0.0;
      beta_mat(len - 1, 1) = 0.0;

      for (int time_ind = len - 2;
           time_ind >= 0;
           --time_ind) {

        const int next_time =
          time_ind + 1;

        const int day_type =
          static_cast<int>(
            day_type_vec[next_time]
          );

        List tran_list_vcovar =
          tran_list_clust[day_type];

        const int transition_index =
          next_time % period_len;

        if (
          transition_index >=
          tran_list_vcovar.size()
        ) {
          stop(
            "Transition index is outside tran_list"
          );
        }

        NumericMatrix tran =
          tran_list_vcovar[
            transition_index
          ];

        const double bp_00 =
          std::log(tran(0, 0)) +
          log_class_0[next_time] +
          beta_mat(next_time, 0);

        const double bp_01 =
          std::log(tran(0, 1)) +
          log_class_1[next_time] +
          beta_mat(next_time, 1);

        const double bp_10 =
          std::log(tran(1, 0)) +
          log_class_0[next_time] +
          beta_mat(next_time, 0);

        const double bp_11 =
          std::log(tran(1, 1)) +
          log_class_1[next_time] +
          beta_mat(next_time, 1);

        beta_mat(time_ind, 0) =
          LogSumExp2C(
            bp_00,
            bp_01
          );

        beta_mat(time_ind, 1) =
          LogSumExp2C(
            bp_10,
            bp_11
          );
      }

      alpha_cube.slice(clust_i) =
        alpha_mat;

      beta_cube.slice(clust_i) =
        beta_mat;
    }

    alpha_list[ind] =
      wrap(alpha_cube);

    beta_list[ind] =
      wrap(beta_cube);
  }

  return List::create(
    _["alpha"] = alpha_list,
    _["beta"] = beta_list
  );
}

// [[Rcpp::export]]
mat CalcTranHelperC(int init_state, int new_state, NumericMatrix act, NumericMatrix light, List tran_list_mat,
                    cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend,
                    NumericVector ind_like_vec, List alpha, List beta, double lod_act, double lod_light,
                    cube corr_mat, cube lintegral_mat, mat pi_l, int clust_i, mat vcovar_mat,
                    cube lambda_act_mat, cube lambda_light_mat, bool tobit){
  int num_people = act.ncol();
  int len = act.nrow();

  mat tran_vals_re_mat( len-1, num_people );

  for (int ind = 0; ind < num_people; ind++) {

    vec vcovar_vec = vcovar_mat.col(ind);

    List tran_list_clust = tran_list_mat[clust_i];
    mat tran_mat = tran_list_clust(0);

    if (tran_list_clust.size() > 1){
      mat tran_mat_week = tran_list_clust(0);
      mat tran_mat_weekend = tran_list_clust(1);


      tran_mat_week.each_col() %= (1-vcovar_vec);
      tran_mat_weekend.each_col() %= vcovar_vec;
      tran_mat = tran_mat_week + tran_mat_weekend;
    }



    //0,0->0 & 1,0->1 & 0,1->2 & 1,1->3
    int tran_vec_ind = init_state + (new_state * 2);

    arma::cube alpha_ind = alpha(ind);
    arma::cube beta_ind = beta(ind);
    double likelihood = ind_like_vec(ind);

    NumericMatrix act_ind = act( Range(1,len-1) , Range(ind,ind) );
    NumericVector act_ind_m1 = act_ind.column(0);

    NumericMatrix light_ind = light( Range(1,len-1) , Range(ind,ind) );
    NumericVector light_ind_m1 = light_ind.column(0);

    vec class_vec_week = logClassificationC( act_ind, light_ind,
                                             emit_act_week(new_state,0,clust_i),
                                             emit_act_week(new_state,1,clust_i),
                                             emit_light_week(new_state,0,clust_i),
                                             emit_light_week(new_state,1,clust_i),
                                             lod_act, lod_light, corr_mat(clust_i,new_state,0), lintegral_mat(clust_i,new_state,0),
                                             lambda_act_mat(clust_i,new_state,0),lambda_light_mat(clust_i,new_state,0),tobit);

    vec class_vec_weekend = logClassificationC( act_ind, light_ind,
                                                emit_act_weekend(new_state,0,clust_i),
                                                emit_act_weekend(new_state,1,clust_i),
                                                emit_light_weekend(new_state,0,clust_i),
                                                emit_light_weekend(new_state,1,clust_i),
                                                lod_act, lod_light, corr_mat(clust_i,new_state,1), lintegral_mat(clust_i,new_state,1),
                                                lambda_act_mat(clust_i,new_state,1),lambda_light_mat(clust_i,new_state,1),tobit);

    vec class_vec = (class_vec_week % (1-vcovar_vec)) + (class_vec_weekend % vcovar_vec);

    vec alpha_ind_slice = alpha_ind(span(0,len-2),span(init_state,init_state),span(clust_i,clust_i));
    vec beta_ind_slice = beta_ind(span(1,len-1),span(new_state,new_state),span(clust_i,clust_i));
    vec tran_vec_slice = tran_mat.col(tran_vec_ind);

    vec temp = alpha_ind_slice + beta_ind_slice + log(tran_vec_slice) + log(pi_l(ind,clust_i)) + class_vec - likelihood;
    vec tran_vals_re_ind = arma::exp(temp);

    tran_vals_re_mat.col(ind) = tran_vals_re_ind;

  }


  return tran_vals_re_mat;
}

inline R_xlen_t ArrayIndex3C(
    const int i,
    const int j,
    const int k,
    const int dim1,
    const int dim2
) {
  // Index for an R array with dimensions:
  // dim1 x dim2 x dim3
  return
    static_cast<R_xlen_t>(i) +
    static_cast<R_xlen_t>(dim1) *
      (
        static_cast<R_xlen_t>(j) +
        static_cast<R_xlen_t>(dim2) *
          static_cast<R_xlen_t>(k)
      );
}


inline R_xlen_t ArrayIndex4C(
    const int i,
    const int j,
    const int k,
    const int l,
    const int dim1,
    const int dim2,
    const int dim3
) {
  // Index for an R array with dimensions:
  // dim1 x dim2 x dim3 x dim4
  return
    static_cast<R_xlen_t>(i) +
    static_cast<R_xlen_t>(dim1) *
      (
        static_cast<R_xlen_t>(j) +
        static_cast<R_xlen_t>(dim2) *
          (
            static_cast<R_xlen_t>(k) +
            static_cast<R_xlen_t>(dim3) *
              static_cast<R_xlen_t>(l)
          )
      );
}

inline double ExpTransitionPosteriorC(
    const double log_value
) {
  if (log_value == R_NegInf) {
    return 0.0;
  }

  const double value = std::exp(log_value);

  if (!std::isfinite(value)) {
    stop(
      "Nonfinite posterior transition probability"
    );
  }

  return value;
}

// [[Rcpp::export]]
List CalcTranGradHessFastC(
    const NumericMatrix& act,
    const NumericMatrix& light,
    const List& tran_list_mat,
    const arma::cube& emit_act_week,
    const arma::cube& emit_light_week,
    const arma::cube& emit_act_weekend,
    const arma::cube& emit_light_weekend,
    const NumericVector& ind_like_vec,
    const List& alpha,
    const List& beta,
    const double lod_act,
    const double lod_light,
    const arma::cube& corr_mat,
    const arma::cube& lintegral_mat,
    const arma::mat& pi_l,
    const NumericMatrix& vcovar_mat,
    const NumericVector& sweights_vec,
    const arma::cube& lambda_act_mat,
    const arma::cube& lambda_light_mat,
    const bool tobit,
    const int vcovar_num,
    const double period_len
) {
  const int len = act.nrow();
  const int num_people = act.ncol();
  const int mix_num = pi_l.n_cols;

  if (light.nrow() != len ||
      light.ncol() != num_people) {
    stop("act and light dimensions do not match");
  }

  if (vcovar_mat.nrow() != len ||
      vcovar_mat.ncol() != num_people) {
    stop(
      "vcovar_mat must have the same dimensions as act"
    );
  }

  if (alpha.size() != num_people ||
      beta.size() != num_people) {
    stop(
      "alpha and beta must contain one entry per participant"
    );
  }

  if (ind_like_vec.size() != num_people ||
    sweights_vec.size() != num_people) {
    stop(
      "ind_like_vec and sweights_vec must contain "
      "one entry per participant"
    );
  }

  if (vcovar_num < 1 || vcovar_num > 2) {
    stop(
      "CalcTranGradHessFastC currently supports "
      "one or two transition day types"
    );
  }

  /*
   * Output dimensions:
   *
   * gradient:
   *   6 parameters x class x day type
   *
   * Hessian:
   *   6 parameters x 6 parameters x class x day type
   *
   * Parameter order:
   *   0: state-0 intercept
   *   1: state-0 cosine
   *   2: state-0 sine
   *   3: state-1 intercept
   *   4: state-1 cosine
   *   5: state-1 sine
   */

  NumericVector grad_array(
    static_cast<R_xlen_t>(6) *
      mix_num *
      vcovar_num,
    0.0
  );

  NumericVector hess_array(
    static_cast<R_xlen_t>(6) *
      6 *
      mix_num *
      vcovar_num,
    0.0
  );

  grad_array.attr("dim") =
    IntegerVector::create(
      6,
      mix_num,
      vcovar_num
    );

  hess_array.attr("dim") =
    IntegerVector::create(
      6,
      6,
      mix_num,
      vcovar_num
    );

  /*
   * The first transition corresponds to original R time 2.
   *
   * Therefore:
   *   C++ transition index 0 -> R time 2
   *   C++ transition index 1 -> R time 3
   *   ...
   */

  std::vector<double> cos_vec(len - 1);
  std::vector<double> sin_vec(len - 1);

  for (int t = 0; t < len - 1; ++t) {
    const double original_time =
      static_cast<double>(t + 2);

    const double angle =
      2.0 *
      arma::datum::pi *
      original_time /
      period_len;

    cos_vec[t] = std::cos(angle);
    sin_vec[t] = std::sin(angle);
  }

  /*
   * Participant loop.
   */

  for (int ind = 0; ind < num_people; ++ind) {
    NumericVector alpha_ind = alpha[ind];
    NumericVector beta_ind = beta[ind];

    IntegerVector alpha_dim =
      alpha_ind.attr("dim");

    IntegerVector beta_dim =
      beta_ind.attr("dim");

    if (
      alpha_dim.size() != 3 ||
      beta_dim.size() != 3 ||
      alpha_dim[0] != len ||
      beta_dim[0] != len ||
      alpha_dim[1] != 2 ||
      beta_dim[1] != 2 ||
      alpha_dim[2] != mix_num ||
      beta_dim[2] != mix_num
    ) {
      stop(
        "Unexpected alpha or beta dimensions"
      );
    }

    const double likelihood =
      ind_like_vec[ind];

    const double sample_weight =
      sweights_vec[ind];

    if (
      !std::isfinite(likelihood) ||
      !std::isfinite(sample_weight)
    ) {
      stop(
        "Nonfinite individual likelihood or sample weight"
      );
    }

    /*
     * Emission data for times 2 through len.
     * These vectors are reused across all latent classes.
     */

    NumericVector act_transition(len - 1);
    NumericVector light_transition(len - 1);

    for (int t = 0; t < len - 1; ++t) {
      act_transition[t] =
        act(t + 1, ind);

      light_transition[t] =
        light(t + 1, ind);
    }

    /*
     * Latent-class loop.
     */

    for (int clust_i = 0;
         clust_i < mix_num;
         ++clust_i) {

      const double pi_value =
        pi_l(ind, clust_i);

      const double log_pi =
        pi_value > 0.0
          ? std::log(pi_value)
          : R_NegInf;

      /*
       * Retrieve the transition-probability matrices once
       * for this class.
       *
       * Each matrix has:
       *   rows = times 2 through len
       *   columns = p00, p10, p01, p11
       */

      List tran_list_clust =
        tran_list_mat[clust_i];

      if (tran_list_clust.size() < vcovar_num) {
        stop(
          "tran_list_mat contains too few day types"
        );
      }

      std::vector<NumericMatrix> tran_mats;
      tran_mats.reserve(vcovar_num);

      for (int vcovar_ind = 0;
           vcovar_ind < vcovar_num;
           ++vcovar_ind) {
        NumericMatrix tran_mat =
          tran_list_clust[vcovar_ind];

        if (
          tran_mat.nrow() != len - 1 ||
          tran_mat.ncol() != 4
        ) {
          stop(
            "Each transition matrix must have "
            "len - 1 rows and four columns"
          );
        }

        tran_mats.push_back(tran_mat);
      }

      /*
       * Calculate the four emission vectors only once:
       *
       * state 0, weekday
       * state 1, weekday
       * state 0, weekend
       * state 1, weekend
       *
       * The old CalcTranHelperC recalculates destination-state
       * emissions separately for each initial state.
       */

      arma::vec log_class_0_week =
        logClassificationC(
          act_transition,
          light_transition,
          emit_act_week(0, 0, clust_i),
          emit_act_week(0, 1, clust_i),
          emit_light_week(0, 0, clust_i),
          emit_light_week(0, 1, clust_i),
          lod_act,
          lod_light,
          corr_mat(clust_i, 0, 0),
          lintegral_mat(clust_i, 0, 0),
          lambda_act_mat(clust_i, 0, 0),
          lambda_light_mat(clust_i, 0, 0),
          tobit
        );

      arma::vec log_class_1_week =
        logClassificationC(
          act_transition,
          light_transition,
          emit_act_week(1, 0, clust_i),
          emit_act_week(1, 1, clust_i),
          emit_light_week(1, 0, clust_i),
          emit_light_week(1, 1, clust_i),
          lod_act,
          lod_light,
          corr_mat(clust_i, 1, 0),
          lintegral_mat(clust_i, 1, 0),
          lambda_act_mat(clust_i, 1, 0),
          lambda_light_mat(clust_i, 1, 0),
          tobit
        );

      arma::vec log_class_0_weekend =
        log_class_0_week;

      arma::vec log_class_1_weekend =
        log_class_1_week;

      if (vcovar_num == 2) {
        log_class_0_weekend =
          logClassificationC(
            act_transition,
            light_transition,
            emit_act_weekend(0, 0, clust_i),
            emit_act_weekend(0, 1, clust_i),
            emit_light_weekend(0, 0, clust_i),
            emit_light_weekend(0, 1, clust_i),
            lod_act,
            lod_light,
            corr_mat(clust_i, 0, 1),
            lintegral_mat(clust_i, 0, 1),
            lambda_act_mat(clust_i, 0, 1),
            lambda_light_mat(clust_i, 0, 1),
            tobit
          );

        log_class_1_weekend =
          logClassificationC(
            act_transition,
            light_transition,
            emit_act_weekend(1, 0, clust_i),
            emit_act_weekend(1, 1, clust_i),
            emit_light_weekend(1, 0, clust_i),
            emit_light_weekend(1, 1, clust_i),
            lod_act,
            lod_light,
            corr_mat(clust_i, 1, 1),
            lintegral_mat(clust_i, 1, 1),
            lambda_act_mat(clust_i, 1, 1),
            lambda_light_mat(clust_i, 1, 1),
            tobit
          );
      }

      /*
       * Time loop.
       */

      for (int t = 0; t < len - 1; ++t) {
        double day_value =
          vcovar_mat(t + 1, ind);

        if (
          NumericVector::is_na(day_value) ||
          !std::isfinite(day_value)
        ) {
          stop(
            "vcovar_mat contains a missing transition day type"
          );
        }

        int day_type = 0;

        if (vcovar_num == 2) {
          day_type =
            static_cast<int>(
              std::round(day_value)
            );
        }

        if (
          day_type < 0 ||
          day_type >= vcovar_num
        ) {
          stop(
            "Transition day-type values must be "
            "coded from 0 to vcovar_num - 1"
          );
        }

        const NumericMatrix& tran_mat =
          tran_mats[day_type];

        /*
         * Column order created by Params2TranVectorT:
         *
         * column 0 = p00
         * column 1 = p10
         * column 2 = p01
         * column 3 = p11
         */

        const double p00 =
          tran_mat(t, 0);

        const double p10 =
          tran_mat(t, 1);

        const double p01 =
          tran_mat(t, 2);

        const double p11 =
          tran_mat(t, 3);

        const double emit_state0 =
          day_type == 0
            ? log_class_0_week[t]
            : log_class_0_weekend[t];

        const double emit_state1 =
          day_type == 0
            ? log_class_1_week[t]
            : log_class_1_weekend[t];

        const R_xlen_t alpha_0_index =
          ArrayIndex3C(
            t,
            0,
            clust_i,
            len,
            2
          );

        const R_xlen_t alpha_1_index =
          ArrayIndex3C(
            t,
            1,
            clust_i,
            len,
            2
          );

        const R_xlen_t beta_0_index =
          ArrayIndex3C(
            t + 1,
            0,
            clust_i,
            len,
            2
          );

        const R_xlen_t beta_1_index =
          ArrayIndex3C(
            t + 1,
            1,
            clust_i,
            len,
            2
          );

        /*
         * Posterior transition probabilities:
         *
         * xi00: state 0 -> state 0
         * xi01: state 0 -> state 1
         * xi10: state 1 -> state 0
         * xi11: state 1 -> state 1
         */

        const double xi00 =
          sample_weight *
          ExpTransitionPosteriorC(
            alpha_ind[alpha_0_index] +
            std::log(p00) +
            emit_state0 +
            beta_ind[beta_0_index] +
            log_pi -
            likelihood
          );

        const double xi01 =
          sample_weight *
          ExpTransitionPosteriorC(
            alpha_ind[alpha_0_index] +
            std::log(p01) +
            emit_state1 +
            beta_ind[beta_1_index] +
            log_pi -
            likelihood
          );

        const double xi10 =
          sample_weight *
          ExpTransitionPosteriorC(
            alpha_ind[alpha_1_index] +
            std::log(p10) +
            emit_state0 +
            beta_ind[beta_0_index] +
            log_pi -
            likelihood
          );

        const double xi11 =
          sample_weight *
          ExpTransitionPosteriorC(
            alpha_ind[alpha_1_index] +
            std::log(p11) +
            emit_state1 +
            beta_ind[beta_1_index] +
            log_pi -
            likelihood
          );

        /*
         * Derivatives for transitions leaving state 0.
         *
         * p01 is modeled with a logistic regression:
         *
         * d log(p00) / d eta = -p01
         * d log(p01) / d eta =  p00
         *
         * Both second derivatives are:
         *
         * -p00 * p01
         */

        const double score_state0 =
          -xi00 * p01 +
           xi01 * p00;

        const double curvature_state0 =
          -(xi00 + xi01) *
          p00 *
          p01;

        /*
         * Derivatives for transitions leaving state 1.
         *
         * p10 is modeled with a logistic regression:
         *
         * d log(p10) / d eta =  p11
         * d log(p11) / d eta = -p10
         */

        const double score_state1 =
           xi10 * p11 -
           xi11 * p10;

        const double curvature_state1 =
          -(xi10 + xi11) *
          p10 *
          p11;

        const double design[3] = {
          1.0,
          cos_vec[t],
          sin_vec[t]
        };

        /*
         * Add state-0 gradient and Hessian block.
         */

        for (int row = 0; row < 3; ++row) {
          const int parameter_row = row;

          const R_xlen_t grad_index =
            ArrayIndex3C(
              parameter_row,
              clust_i,
              day_type,
              6,
              mix_num
            );

          grad_array[grad_index] +=
            score_state0 *
            design[row];

          for (int col = 0; col < 3; ++col) {
            const int parameter_col = col;

            const R_xlen_t hess_index =
              ArrayIndex4C(
                parameter_row,
                parameter_col,
                clust_i,
                day_type,
                6,
                6,
                mix_num
              );

            hess_array[hess_index] +=
              curvature_state0 *
              design[row] *
              design[col];
          }
        }

        /*
         * Add state-1 gradient and Hessian block.
         */

        for (int row = 0; row < 3; ++row) {
          const int parameter_row =
            row + 3;

          const R_xlen_t grad_index =
            ArrayIndex3C(
              parameter_row,
              clust_i,
              day_type,
              6,
              mix_num
            );

          grad_array[grad_index] +=
            score_state1 *
            design[row];

          for (int col = 0; col < 3; ++col) {
            const int parameter_col =
              col + 3;

            const R_xlen_t hess_index =
              ArrayIndex4C(
                parameter_row,
                parameter_col,
                clust_i,
                day_type,
                6,
                6,
                mix_num
              );

            hess_array[hess_index] +=
              curvature_state1 *
              design[row] *
              design[col];
          }
        }
      }
    }
  }

  return List::create(
    _["grad_array"] = grad_array,
    _["hess_array"] = hess_array
  );
}



// [[Rcpp::export]]
mat ForwardIndAltC(vec decoded_ind, rowvec init, Rcpp::List tran_list, int clust_i, vec vcovar_vec){

  mat alpha( decoded_ind.n_elem, 2 );

  alpha(0,0) = log(init(0));
  alpha(0,1) = log(init(1));

  List tran_list_clust = tran_list[clust_i];

  double log_class_0;
  double log_class_1;

  for (int i = 1; i < decoded_ind.n_elem; i++) {
    List tran_list_vcovar = tran_list_clust[vcovar_vec(i)];

    NumericMatrix tran = tran_list_vcovar[i%96];


    log_class_0 = 0;
    log_class_1 = 0;

    if (decoded_ind(i) == 0){
      log_class_0 = 0;
      log_class_1 = -9999999;
    }

    if (decoded_ind(i) == 1){
      log_class_0 = -9999999;
      log_class_1 = 0;
    }

    double fp_00 = alpha(i-1,0) + log(tran(0,0)) + log_class_0;

    double fp_10 = alpha(i-1,1) + log(tran(1,0)) + log_class_0;

    double fp_01 = alpha(i-1,0) + log(tran(0,1)) + log_class_1;

    double fp_11 = alpha(i-1,1) + log(tran(1,1)) + log_class_1;

    NumericVector fp_0 = NumericVector::create(fp_00,fp_10);
    NumericVector fp_1 = NumericVector::create(fp_01,fp_11);

    alpha(i,0) = logSumExpC(fp_0);
    alpha(i,1) = logSumExpC(fp_1);
  }

  return(alpha);
}




// [[Rcpp::export]]
double weightedLogClassificationCTobitC(
    const NumericVector& act,
    const NumericVector& light,
    const NumericVector& weights,
    const double mu_act,
    const double sig_act,
    const double mu_light,
    const double sig_light,
    const double lod_act,
    const double lod_light,
    const double bivar_corr,
    const double lintegral
) {
  const R_xlen_t n = act.size();

  if (light.size() != n || weights.size() != n) {
    stop("act, light, and weights must have the same length");
  }

  if (
    sig_act <= 0.0 ||
    sig_light <= 0.0 ||
    std::abs(bivar_corr) >= 1.0
  ) {
    return R_PosInf;
  }

  // These quantities do not vary across observations.
  const double one_minus_rho_squared =
    1.0 - bivar_corr * bivar_corr;

  const double sqrt_one_minus_rho_squared =
    std::sqrt(one_minus_rho_squared);

  const double sig_light_cond =
    sig_light * sqrt_one_minus_rho_squared;

  const double sig_act_cond =
    sig_act * sqrt_one_minus_rho_squared;

  const double light_cond_slope =
    bivar_corr * sig_light / sig_act;

  const double act_cond_slope =
    bivar_corr * sig_act / sig_light;

  double negative_weighted_loglike = 0.0;

  for (R_xlen_t i = 0; i < n; ++i) {
    const double weight = weights[i];

    // Posterior probabilities can underflow to exactly zero.
    // Such observations contribute exactly zero to the objective.
    if (weight == 0.0) {
      continue;
    }

    if (!R_FINITE(weight)) {
      return NA_REAL;
    }

    const double act_obs = act[i];
    const double light_obs = light[i];

    const bool act_missing =
      NumericVector::is_na(act_obs);

    const bool light_missing =
      NumericVector::is_na(light_obs);

    double log_density = 0.0;

    // Both measurements observed.
    if (!act_missing && !light_missing) {
      if (act_obs > lod_act) {
        if (light_obs > lod_light) {
          // Case 1: both above their limits of detection.
          const double mu_light_cond =
            mu_light +
            light_cond_slope * (act_obs - mu_act);

          log_density =
            R::dnorm(
              act_obs,
              mu_act,
              sig_act,
              true
            ) +
            R::dnorm(
              light_obs,
              mu_light_cond,
              sig_light_cond,
              true
            );

        } else {
          // Case 2: activity observed, light censored.
          const double mu_light_cond =
            mu_light +
            light_cond_slope * (act_obs - mu_act);

          log_density =
            R::dnorm(
              act_obs,
              mu_act,
              sig_act,
              true
            ) +
            R::pnorm(
              light_obs,
              mu_light_cond,
              sig_light_cond,
              true,
              true
            );
        }

      } else {
        if (light_obs > lod_light) {
          // Case 3: activity censored, light observed.
          const double mu_act_cond =
            mu_act +
            act_cond_slope * (light_obs - mu_light);

          log_density =
            R::pnorm(
              act_obs,
              mu_act_cond,
              sig_act_cond,
              true,
              true
            ) +
            R::dnorm(
              light_obs,
              mu_light,
              sig_light,
              true
            );

        } else {
          // Case 4: both censored.
          log_density = lintegral;
        }
      }

    } else if (!act_missing && light_missing) {
      // Light missing; use marginal activity distribution.
      if (act_obs > lod_act) {
        log_density =
          R::dnorm(
            act_obs,
            mu_act,
            sig_act,
            true
          );
      } else {
        log_density =
          R::pnorm(
            act_obs,
            mu_act,
            sig_act,
            true,
            true
          );
      }

    } else if (act_missing && !light_missing) {
      // Activity missing; use marginal light distribution.
      if (light_obs > lod_light) {
        log_density =
          R::dnorm(
            light_obs,
            mu_light,
            sig_light,
            true
          );
      } else {
        log_density =
          R::pnorm(
            light_obs,
            mu_light,
            sig_light,
            true,
            true
          );
      }

    } else {
      // Both missing: contribution is zero.
      log_density = 0.0;
    }

    // Preserve the current R behavior.
    if (log_density == R_NegInf) {
      log_density = -9999.0;
    }

    if (NumericVector::is_na(log_density)) {
      return NA_REAL;
    }

    negative_weighted_loglike -=
      weight * log_density;
  }

  return negative_weighted_loglike;
}


inline R_xlen_t ConvertCaseIndex(
    const int one_based_index,
    const R_xlen_t vector_length
) {
  if (IntegerVector::is_na(one_based_index)) {
    stop("Case index contains NA");
  }

  const R_xlen_t zero_based_index =
    static_cast<R_xlen_t>(one_based_index - 1);

  if (
    zero_based_index < 0 ||
    zero_based_index >= vector_length
  ) {
    stop("Case index is outside the observation vector");
  }

  return zero_based_index;
}


inline double WeightedNormalLogLikeFromStats(
    const NumericVector& stats,
    const double mu,
    const double sigma
) {
  if (stats.size() != 3) {
    stop("Univariate sufficient-statistic vector must have length 3");
  }

  if (!(sigma > 0.0) || !R_FINITE(sigma)) {
    return R_NegInf;
  }

  const long double sum_w = stats[0];

  if (sum_w == 0.0L) {
    return 0.0;
  }

  const long double sum_w_x = stats[1];
  const long double sum_w_x2 = stats[2];
  const long double mu_ld = mu;
  const long double sigma_ld = sigma;

  const long double centered_sum =
    sum_w_x2 -
    2.0L * mu_ld * sum_w_x +
    mu_ld * mu_ld * sum_w;

  const long double log_2pi =
    std::log(2.0L * arma::datum::pi);

  const long double result =
    -0.5L * sum_w *
      (log_2pi + 2.0L * std::log(sigma_ld)) -
    centered_sum /
      (2.0L * sigma_ld * sigma_ld);

  return static_cast<double>(result);
}


inline double WeightedBivariateLogLikeFromStats(
    const NumericVector& stats,
    const double mu_act,
    const double sig_act,
    const double mu_light,
    const double sig_light,
    const double rho
) {
  if (stats.size() != 6) {
    stop("Bivariate sufficient-statistic vector must have length 6");
  }

  if (
    !(sig_act > 0.0) ||
    !(sig_light > 0.0) ||
    std::abs(rho) >= 1.0
  ) {
    return R_NegInf;
  }

  const long double sum_w = stats[0];

  if (sum_w == 0.0L) {
    return 0.0;
  }

  const long double sum_w_act = stats[1];
  const long double sum_w_light = stats[2];
  const long double sum_w_act2 = stats[3];
  const long double sum_w_light2 = stats[4];
  const long double sum_w_act_light = stats[5];

  const long double mu_a = mu_act;
  const long double mu_l = mu_light;
  const long double sig_a = sig_act;
  const long double sig_l = sig_light;
  const long double rho_ld = rho;

  const long double centered_act2 =
    sum_w_act2 -
    2.0L * mu_a * sum_w_act +
    mu_a * mu_a * sum_w;

  const long double centered_light2 =
    sum_w_light2 -
    2.0L * mu_l * sum_w_light +
    mu_l * mu_l * sum_w;

  const long double centered_cross =
    sum_w_act_light -
    mu_l * sum_w_act -
    mu_a * sum_w_light +
    mu_a * mu_l * sum_w;

  const long double one_minus_rho2 =
    1.0L - rho_ld * rho_ld;

  const long double quadratic =
    centered_act2 / (sig_a * sig_a) -
    2.0L * rho_ld * centered_cross /
      (sig_a * sig_l) +
    centered_light2 / (sig_l * sig_l);

  const long double log_2pi =
    std::log(2.0L * arma::datum::pi);

  const long double result =
    -sum_w *
      (
        log_2pi +
        std::log(sig_a) +
        std::log(sig_l) +
        0.5L * std::log(one_minus_rho2)
      ) -
    quadratic / (2.0L * one_minus_rho2);

  return static_cast<double>(result);
}


// [[Rcpp::export]]
List buildTobitCaseCacheC(
    const NumericVector& act,
    const NumericVector& light,
    const NumericVector& weights,
    const List& case_indices
) {
  const R_xlen_t n = act.size();

  if (light.size() != n || weights.size() != n) {
    stop("act, light, and weights must have equal lengths");
  }

  IntegerVector case1_index = case_indices["case1"];
  IntegerVector case2_index = case_indices["case2"];
  IntegerVector case3_index = case_indices["case3"];
  IntegerVector case4_index = case_indices["case4"];
  IntegerVector case5_index = case_indices["case5"];
  IntegerVector case6_index = case_indices["case6"];
  IntegerVector case7_index = case_indices["case7"];
  IntegerVector case8_index = case_indices["case8"];

  IntegerVector case2_group =
  case_indices["case2_group"];

  NumericVector case2_unique_act =
    case_indices["case2_unique_act"];

  if (case2_group.size() != case2_index.size()) {
    stop(
      "Case 2 group mapping must have one value "
      "per Case 2 observation"
    );
  }

  if (
    case2_index.size() > 0 &&
    case2_unique_act.size() == 0
  ) {
    stop(
      "Case 2 observations exist but no unique "
      "activity values were supplied"
    );
  }

  long double total_weight = 0.0L;

  for (R_xlen_t i = 0; i < n; ++i) {
    if (!R_FINITE(weights[i])) {
      stop("Nonfinite emission weight");
    }

    total_weight += weights[i];
  }

  // Case 1: bivariate normal sufficient statistics.
  long double c1_w = 0.0L;
  long double c1_wa = 0.0L;
  long double c1_wl = 0.0L;
  long double c1_waa = 0.0L;
  long double c1_wll = 0.0L;
  long double c1_wal = 0.0L;

  for (R_xlen_t j = 0; j < case1_index.size(); ++j) {
    const R_xlen_t i =
      ConvertCaseIndex(case1_index[j], n);

    const long double w = weights[i];
    const long double a = act[i];
    const long double l = light[i];

    c1_w += w;
    c1_wa += w * a;
    c1_wl += w * l;
    c1_waa += w * a * a;
    c1_wll += w * l * l;
    c1_wal += w * a * l;
  }

  
  
  // Case 2: summarize the marginal activity density
  // and aggregate posterior weights by unique activity value.
  long double c2_w = 0.0L;
  long double c2_wa = 0.0L;
  long double c2_waa = 0.0L;

  const R_xlen_t case2_group_num =
    case2_unique_act.size();

  std::vector<long double>
    case2_group_weight_accumulator(
      static_cast<std::size_t>(case2_group_num),
      0.0L
    );

  for (R_xlen_t j = 0;
      j < case2_index.size();
      ++j) {

    const R_xlen_t i =
      ConvertCaseIndex(
        case2_index[j],
        n
      );

    const int group_one_based =
      case2_group[j];

    if (
      IntegerVector::is_na(group_one_based) ||
      group_one_based < 1 ||
      static_cast<R_xlen_t>(group_one_based) >
        case2_group_num
    ) {
      stop(
        "Case 2 group number is outside "
        "the unique-activity vector"
      );
    }

    const R_xlen_t group_zero_based =
      static_cast<R_xlen_t>(
        group_one_based - 1
      );

    const long double w =
      weights[i];

    const long double a =
      act[i];

    case2_group_weight_accumulator[
      static_cast<std::size_t>(
        group_zero_based
      )
    ] += w;

    c2_w += w;
    c2_wa += w * a;
    c2_waa += w * a * a;
  }

  NumericVector case2_group_weights(
    case2_group_num
  );

  for (R_xlen_t group_ind = 0;
      group_ind < case2_group_num;
      ++group_ind) {

    case2_group_weights[group_ind] =
      static_cast<double>(
        case2_group_weight_accumulator[
          static_cast<std::size_t>(
            group_ind
          )
        ]
      );
  }


  // Case 3: summarize marginal light density,
  // retain weights for conditional activity CDF.
  long double c3_w = 0.0L;
  long double c3_wl = 0.0L;
  long double c3_wll = 0.0L;

  NumericVector case3_weights(case3_index.size());

  for (R_xlen_t j = 0; j < case3_index.size(); ++j) {
    const R_xlen_t i =
      ConvertCaseIndex(case3_index[j], n);

    const long double w = weights[i];
    const long double l = light[i];

    case3_weights[j] = static_cast<double>(w);

    c3_w += w;
    c3_wl += w * l;
    c3_wll += w * l * l;
  }

  // Case 4: every observation has the same log probability.
  long double case4_sum_w = 0.0L;

  for (R_xlen_t j = 0; j < case4_index.size(); ++j) {
    const R_xlen_t i =
      ConvertCaseIndex(case4_index[j], n);

    case4_sum_w += weights[i];
  }

  // Case 5: marginal activity density.
  long double c5_w = 0.0L;
  long double c5_wa = 0.0L;
  long double c5_waa = 0.0L;

  for (R_xlen_t j = 0; j < case5_index.size(); ++j) {
    const R_xlen_t i =
      ConvertCaseIndex(case5_index[j], n);

    const long double w = weights[i];
    const long double a = act[i];

    c5_w += w;
    c5_wa += w * a;
    c5_waa += w * a * a;
  }

  // Case 6: activity CDF must still be evaluated observation by observation.
  NumericVector case6_weights(case6_index.size());

  for (R_xlen_t j = 0; j < case6_index.size(); ++j) {
    const R_xlen_t i =
      ConvertCaseIndex(case6_index[j], n);

    case6_weights[j] = weights[i];
  }

  // Case 7: marginal light density.
  long double c7_w = 0.0L;
  long double c7_wl = 0.0L;
  long double c7_wll = 0.0L;

  for (R_xlen_t j = 0; j < case7_index.size(); ++j) {
    const R_xlen_t i =
      ConvertCaseIndex(case7_index[j], n);

    const long double w = weights[i];
    const long double l = light[i];

    c7_w += w;
    c7_wl += w * l;
    c7_wll += w * l * l;
  }

  // Case 8: light CDF must still be evaluated observation by observation.
  NumericVector case8_weights(case8_index.size());

  for (R_xlen_t j = 0; j < case8_index.size(); ++j) {
    const R_xlen_t i =
      ConvertCaseIndex(case8_index[j], n);

    case8_weights[j] = weights[i];
  }

  return List::create(
    _["case1_stats"] = NumericVector::create(
      static_cast<double>(c1_w),
      static_cast<double>(c1_wa),
      static_cast<double>(c1_wl),
      static_cast<double>(c1_waa),
      static_cast<double>(c1_wll),
      static_cast<double>(c1_wal)
    ),

    _["case2_stats"] = NumericVector::create(
      static_cast<double>(c2_w),
      static_cast<double>(c2_wa),
      static_cast<double>(c2_waa)
    ),

    _["case2_group_weights"] =case2_group_weights,

    _["case3_stats"] = NumericVector::create(
      static_cast<double>(c3_w),
      static_cast<double>(c3_wl),
      static_cast<double>(c3_wll)
    ),

    _["case3_weights"] = case3_weights,

    _["case4_sum_w"] =
      static_cast<double>(case4_sum_w),

    _["case5_stats"] = NumericVector::create(
      static_cast<double>(c5_w),
      static_cast<double>(c5_wa),
      static_cast<double>(c5_waa)
    ),

    _["case6_weights"] = case6_weights,

    _["case7_stats"] = NumericVector::create(
      static_cast<double>(c7_w),
      static_cast<double>(c7_wl),
      static_cast<double>(c7_wll)
    ),

    _["case8_weights"] = case8_weights,

    _["total_weight"] =
      static_cast<double>(total_weight)
  );
}

// [[Rcpp::export]]
double weightedTobitObjectiveByCaseC(
    const NumericVector& act,
    const NumericVector& light,
    const List& case_indices,
    const List& case_cache,
    const double mu_act,
    const double sig_act,
    const double mu_light,
    const double sig_light,
    const double lod_act,
    const double lod_light,
    const double bivar_corr,
    const double lintegral
) {
  const R_xlen_t n = act.size();

  if (light.size() != n) {
    stop("act and light must have equal lengths");
  }

  if (
    !(sig_act > 0.0) ||
    !(sig_light > 0.0) ||
    std::abs(bivar_corr) >= 1.0
  ) {
    return R_PosInf;
  }

  NumericVector case1_stats =
    case_cache["case1_stats"];

  NumericVector case2_stats =
  case_cache["case2_stats"];

  NumericVector case2_group_weights =
    case_cache["case2_group_weights"];

  NumericVector case2_unique_act =
    case_indices["case2_unique_act"];

  const double case2_light_value =
    as<double>(
      case_indices["case2_light_value"]
    );

  NumericVector case3_stats =
    case_cache["case3_stats"];

  NumericVector case3_weights =
    case_cache["case3_weights"];

  const double case4_sum_w =
    as<double>(case_cache["case4_sum_w"]);

  NumericVector case5_stats =
    case_cache["case5_stats"];

  NumericVector case6_weights =
    case_cache["case6_weights"];

  NumericVector case7_stats =
    case_cache["case7_stats"];

  NumericVector case8_weights =
    case_cache["case8_weights"];

  IntegerVector case3_index =
    case_indices["case3"];

  IntegerVector case6_index =
    case_indices["case6"];

  IntegerVector case8_index =
    case_indices["case8"];

  if (
    case2_unique_act.size() !=
      case2_group_weights.size() ||
    case3_index.size() !=
      case3_weights.size() ||
    case6_index.size() !=
      case6_weights.size() ||
    case8_index.size() !=
      case8_weights.size()
  ) {
    stop(
      "Case data and cached weights do not match"
    );
  }

  if (
    case2_group_weights.size() > 0 &&
    !R_FINITE(case2_light_value)
  ) {
    stop(
      "Case 2 censored-light value is nonfinite"
    );
  } 
    


  double log_like = 0.0;

  // Case 1: complete bivariate normal density.
  log_like += WeightedBivariateLogLikeFromStats(
    case1_stats,
    mu_act,
    sig_act,
    mu_light,
    sig_light,
    bivar_corr
  );

  // Case 2: marginal activity density from sufficient statistics.
  log_like += WeightedNormalLogLikeFromStats(
    case2_stats,
    mu_act,
    sig_act
  );

  const double one_minus_rho2 =
    1.0 - bivar_corr * bivar_corr;

  const double sqrt_one_minus_rho2 =
    std::sqrt(one_minus_rho2);

  const double sig_light_cond =
    sig_light * sqrt_one_minus_rho2;

  const double sig_act_cond =
    sig_act * sqrt_one_minus_rho2;

  const double light_cond_slope =
    bivar_corr * sig_light / sig_act;

  const double act_cond_slope =
    bivar_corr * sig_act / sig_light;

  // Case 2: conditional light CDF.
  //
  // All observations with the same activity value
  // have the same conditional probability. Their
  // posterior weights have already been summed.
  for (R_xlen_t group_ind = 0;
      group_ind < case2_group_weights.size();
      ++group_ind) {

    const double w =
      case2_group_weights[group_ind];

    if (w == 0.0) {
      continue;
    }

    const double act_obs =
      case2_unique_act[group_ind];

    const double mu_light_cond =
      mu_light +
      light_cond_slope *
        (act_obs - mu_act);

    double log_probability =
      R::pnorm(
        case2_light_value,
        mu_light_cond,
        sig_light_cond,
        true,
        true
      );

    if (log_probability == R_NegInf) {
      log_probability = -9999.0;
    }

    if (!R_FINITE(log_probability)) {
      return R_PosInf;
    }

    log_like +=
      w * log_probability;
  }

  // Case 3: marginal light density from sufficient statistics.
  log_like += WeightedNormalLogLikeFromStats(
    case3_stats,
    mu_light,
    sig_light
  );

  // Case 3: conditional activity CDF.
  for (R_xlen_t j = 0; j < case3_index.size(); ++j) {
    const double w = case3_weights[j];

    if (w == 0.0) {
      continue;
    }

    const R_xlen_t i =
      ConvertCaseIndex(case3_index[j], n);

    const double act_obs = act[i];
    const double light_obs = light[i];

    const double mu_act_cond =
      mu_act +
      act_cond_slope * (light_obs - mu_light);

    double log_probability =
      R::pnorm(
        act_obs,
        mu_act_cond,
        sig_act_cond,
        true,
        true
      );

    if (log_probability == R_NegInf) {
      log_probability = -9999.0;
    }

    if (!R_FINITE(log_probability)) {
      return R_PosInf;
    }

    log_like += w * log_probability;
  }

  // Case 4: both censored.
  log_like += case4_sum_w * lintegral;

  // Case 5: marginal activity density.
  log_like += WeightedNormalLogLikeFromStats(
    case5_stats,
    mu_act,
    sig_act
  );

  // Case 6: marginal activity CDF.
  for (R_xlen_t j = 0; j < case6_index.size(); ++j) {
    const double w = case6_weights[j];

    if (w == 0.0) {
      continue;
    }

    const R_xlen_t i =
      ConvertCaseIndex(case6_index[j], n);

    double log_probability =
      R::pnorm(
        act[i],
        mu_act,
        sig_act,
        true,
        true
      );

    if (log_probability == R_NegInf) {
      log_probability = -9999.0;
    }

    if (!R_FINITE(log_probability)) {
      return R_PosInf;
    }

    log_like += w * log_probability;
  }

  // Case 7: marginal light density.
  log_like += WeightedNormalLogLikeFromStats(
    case7_stats,
    mu_light,
    sig_light
  );

  // Case 8: marginal light CDF.
  for (R_xlen_t j = 0; j < case8_index.size(); ++j) {
    const double w = case8_weights[j];

    if (w == 0.0) {
      continue;
    }

    const R_xlen_t i =
      ConvertCaseIndex(case8_index[j], n);

    double log_probability =
      R::pnorm(
        light[i],
        mu_light,
        sig_light,
        true,
        true
      );

    if (log_probability == R_NegInf) {
      log_probability = -9999.0;
    }

    if (!R_FINITE(log_probability)) {
      return R_PosInf;
    }

    log_like += w * log_probability;
  }

  return -log_like;
}
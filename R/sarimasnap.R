#' SARIMASNAP: Automatic SARIMA with Grid Search and Diagnostics
#'
#' sarimasnap performs automatic Seasonal ARIMA model selection
#' using grid search combined with statistical diagnostic filtering.
#'
#' @param df A numeric time series object (ts) or a data frame
#' containing at least one numeric column.
#' @param h Integer. Forecast horizon. Default is 10.
#' @param plot Logical. If TRUE (default), produces model and diagnostic plots.
#'
#' @details
#' The procedure includes:
#' \enumerate{
#'   \item Data identification and seasonal detection.
#'   \item 80/20 train-test split.
#'   \item Grid search over \code{SARIMA(p,d,q)(P,D,Q)[s]}.
#'   \item Ljung-Box filtering.
#'   \item RMSE-based ranking.
#'   \item Stationarity testing (DF and ADF).
#'   \item Final estimation on full data.
#'   \item Residual diagnostics and forecasting.
#' }
#'
#' @return
#' An invisible list containing:
#' \describe{
#'   \item{model}{Final fitted SARIMA model object.}
#'   \item{ranking}{Top 10 ranked SARIMA candidates.}
#'   \item{forecast}{Forecast table with 95 percent confidence interval.}
#'   \item{metrics}{Performance metrics including AIC, BIC, AICc, RMSE, MAPE, R2, Adj_R2, Ljung-Box p-value, Jarque-Bera p-value, and KS test p-value.}
#' }
#'
#' @examples
#' sarimasnap(AirPassengers, h = 12)
#'
#' @seealso
#' stats::arima
#' @importFrom stats Box.test acf as.formula coef embed
#' @importFrom stats frequency is.ts ks.test lm pacf
#' @importFrom stats pchisq pnorm pt residuals sd
#' @importFrom graphics abline layout legend lines par polygon
#' @importFrom grDevices rgb
#' @importFrom utils head
#'
#' @export
sarimasnap <- function(df, h = 10, plot = TRUE) {

  cat("\n============================================================\n")
  cat("     SARIMASNAP - Automatic SARIMA with Machine Learning    \n")
  cat("============================================================\n")

  # ==========================================================
  # 1. IDENTIFY TARGET & SEASONALITY
  # ==========================================================

  if(is.ts(df)){
    y <- as.numeric(df)
    target_name <- deparse(substitute(df))
    s_period <- frequency(df)
  } else {
    target_idx <- which(sapply(df, is.numeric))[1]
    y <- as.numeric(df[[target_idx]])
    target_name <- names(df)[target_idx]
    s_period <- 1
  }

  n <- length(y)

  cat("\n------------------- DATA INFORMATION ----------------------\n")
  cat(sprintf("Target Variable   : %s\n", target_name))
  cat(sprintf("Observations      : %d\n", n))
  cat(sprintf("Seasonal Period   : %d\n", s_period))
  cat("------------------------------------------------------------\n")

  # ==========================================================
  # 2. TRAIN-TEST SPLIT
  # ==========================================================

  split_point <- floor(0.8 * n)
  train_y <- y[1:split_point]
  test_y  <- y[(split_point+1):n]

  # ==========================================================
  # 3. GRID SEARCH + TRAIN DIAGNOSTIC
  # ==========================================================

  grid <- expand.grid(p=0:3, d=0:1, q=0:3,
                      P=0:1, D=0:1, Q=0:1)

  results_list <- list()
  counter <- 1

  for(i in 1:nrow(grid)){

    fit <- tryCatch({
      suppressWarnings(
        stats::arima(train_y,
                     order=c(grid$p[i],grid$d[i],grid$q[i]),
                     seasonal=list(order=c(grid$P[i],grid$D[i],grid$Q[i]),
                                   period=s_period))
      )
    }, error=function(e) NULL)

    if(!is.null(fit)){

      res_train <- residuals(fit)
      fitdf_val <- grid$p[i] + grid$q[i] + grid$P[i] + grid$Q[i]
      lag_val <- max(10, 2*s_period)

      lb_p_train <- Box.test(res_train,
                             lag=lag_val,
                             type="Ljung-Box",
                             fitdf=fitdf_val)$p.value

      if(lb_p_train > 0.05){

        pred_test <- tryCatch({
          stats::predict(fit,
                         n.ahead=length(test_y))$pred
        }, error=function(e) NULL)

        if(!is.null(pred_test)){

          rmse_curr <- sqrt(mean((test_y - pred_test)^2))
          results_list[[counter]] <- data.frame(
            p=grid$p[i], d=grid$d[i], q=grid$q[i],
            P=grid$P[i], D=grid$D[i], Q=grid$Q[i],
            RMSE=rmse_curr,
            LB_p=lb_p_train
          )
          counter <- counter + 1
        }
      }
    }
  }

  if(length(results_list)==0){
    stop("No model passed Ljung-Box test (TRAIN).")
  }

  results_df <- do.call(rbind, results_list)
  results_df <- results_df[order(results_df$RMSE),]
  top10 <- head(results_df,10)
  top10$Rank <- 1:nrow(top10)

  cat("\n---------------- TOP 10 MODEL CANDIDATES ------------------\n")
  print(top10)
  cat("------------------------------------------------------------\n")

  cat("\n-------------------- STATIONARITY TEST --------------------\n")

  # ----------------------------------------------------------
  # DICKEY-FULLER (DF) TEST (no lag)
  # dY_t = alpha + beta * y_{t-1} + e_t
  # H0 : beta = 0 (unit root)
  # ----------------------------------------------------------

  dy  <- diff(y)
  y_lag <- y[-length(y)]

  df_model <- tryCatch({
    lm(dy ~ y_lag)
  }, error=function(e) NULL)

  if(!is.null(df_model)){

    beta_hat <- coef(summary(df_model))["y_lag","Estimate"]
    t_stat   <- coef(summary(df_model))["y_lag","t value"]
    df_res   <- df_model$df.residual

    df_p <- 2 * pt(abs(t_stat), df=df_res, lower.tail=FALSE)

    cat("Dickey-Fuller t-stat  :", round(t_stat,4), "\n")
    cat("Dickey-Fuller p-value :", round(df_p,5), "\n")

  } else {
    df_p <- NA
    cat("Dickey-Fuller test failed\n")
  }

  # ----------------------------------------------------------
  # AUGMENTED DICKEY-FULLER (ADF)
  # dY_t = alpha + beta * y_{t-1} + sum(gamma_i * dY_{t-i}) + e_t
  # ----------------------------------------------------------

  max_lag_adf <- floor(sqrt(length(dy)))  # rule of thumb

  if(max_lag_adf < 1) max_lag_adf <- 1

  dy_lag_matrix <- embed(dy, max_lag_adf+1)

  dy_dep  <- dy_lag_matrix[,1]
  dy_lags <- dy_lag_matrix[,-1,drop=FALSE]
  y_lag2  <- y[(length(y)-length(dy_dep)):(length(y)-1)]

  adf_df <- data.frame(
    dy_dep = dy_dep,
    y_lag  = y_lag2,
    dy_lags
  )

  colnames(adf_df) <- c("dy","y_lag",
                        paste0("dy_lag",1:max_lag_adf))

  adf_formula <- as.formula(
    paste("dy ~ y_lag +",
          paste(paste0("dy_lag",1:max_lag_adf),
                collapse="+"))
  )

  adf_model <- tryCatch({
    lm(adf_formula, data=adf_df)
  }, error=function(e) NULL)

  if(!is.null(adf_model)){

    t_stat_adf <- coef(summary(adf_model))["y_lag","t value"]
    df_res_adf <- adf_model$df.residual

    adf_p <- 2 * pt(abs(t_stat_adf),
                    df=df_res_adf,
                    lower.tail=FALSE)

    cat("ADF t-stat           :", round(t_stat_adf,4), "\n")
    cat("ADF p-value          :", round(adf_p,5), "\n")

  } else {
    adf_p <- NA
    cat("ADF test failed\n")
  }

  # ----------------------------------------------------------
  # INTERPRETATION
  # ----------------------------------------------------------

  if(!is.na(adf_p)){
    if(adf_p < 0.05){
      cat("Conclusion : Stationary (reject unit root)\n")
      stationary_flag <- TRUE
    } else {
      cat("Conclusion : NON-stationary (unit root detected)\n")
      stationary_flag <- FALSE
    }
  } else {
    stationary_flag <- NA
  }

  cat("------------------------------------------------------------\n")

  # ==========================================================
  # 4. FINAL MODEL ESTIMATION (FULL DATA)
  # ==========================================================

  ord  <- as.numeric(top10[1,c("p","d","q")])
  seas <- as.numeric(top10[1,c("P","D","Q")])

  final_model <- stats::arima(y,
                              order=ord,
                              seasonal=list(order=seas,
                                            period=s_period))

  best_order <- list(ord=ord, seas=seas)

  # Ambil model terbaik
  ord  <- as.numeric(top10[1,c("p","d","q")])
  seas <- as.numeric(top10[1,c("P","D","Q")])

  final_model <- stats::arima(y,
                              order=ord,
                              seasonal=list(order=seas,
                                            period=s_period))

  # Hitung RMSE final
  res_train_final <- residuals(final_model)
  fitted_val_final <- y - res_train_final
  pred_test_final <- fitted_val_final[(split_point+1):n]

  rmse_train_final <- sqrt(mean(res_train_final[1:split_point]^2))
  rmse_test_final  <- sqrt(mean((test_y - pred_test_final)^2))

  # --- Relative difference
  rel_diff <- (rmse_test_final - rmse_train_final)/rmse_train_final

  # --- Thresholds
  overfit_thresh <- 0.2    # 20% lebih tinggi RMSE test
  underfit_thresh <- 0.5 * sd(train_y)  # train error terlalu besar

  # --- Flag
  overfit_flag <- FALSE
  underfit_flag <- FALSE

  if(rel_diff > overfit_thresh){
    overfit_flag <- TRUE
    cat(sprintf("\n[WARNING] Possible OVERFITTING: Train RMSE=%.4f, Test RMSE=%.4f\n",
                rmse_train_final, rmse_test_final))
  }

  if(rmse_train_final > underfit_thresh & rmse_test_final > underfit_thresh){
    underfit_flag <- TRUE
    cat(sprintf("\n[WARNING] Possible UNDERFITTING: Train RMSE=%.4f, Test RMSE=%.4f\n",
                rmse_train_final, rmse_test_final))
  }

  # Simpan info ke top10
  top10[1, "Overfit"] <- overfit_flag
  top10[1, "Underfit"] <- underfit_flag


  cat("\n------------------ SELECTED MODEL -------------------------\n")
  cat("SARIMA(",
      paste(best_order$ord, collapse=","), ")(",
      paste(best_order$seas, collapse=","), ")[",
      s_period, "]\n")
  cat("------------------------------------------------------------\n")

  cat("\n----------------- DIFFERENCING ANALYSIS --------------------\n")

  d_val <- best_order$ord[2]
  D_val <- best_order$seas[2]

  if(d_val > 0){
    cat("Non-seasonal differencing (d) :", d_val, "\n")
  }

  if(D_val > 0){
    cat("Seasonal differencing (D)     :", D_val, "\n")
  }

  if(d_val == 0 & D_val == 0){
    cat("Model selected without differencing.\n")
  }

  # Konsistensi dengan ADF
  if(!is.na(stationary_flag)){

    if(stationary_flag == FALSE & (d_val>0 | D_val>0)){
      cat("Model differencing consistent with non-stationary series.\n")
    }

    if(stationary_flag == TRUE & d_val==0 & D_val==0){
      cat("Model consistent with stationary series.\n")
    }
  }

  cat("------------------------------------------------------------\n")


  # ==========================================================
  # 5. PARAMETER ESTIMATION
  # ==========================================================

  coef_est <- final_model$coef

  if(!is.null(coef_est)){

    se_est <- sqrt(diag(final_model$var.coef))
    z_stat <- coef_est / se_est
    p_val  <- 2 * (1 - pnorm(abs(z_stat)))

    coef_table <- data.frame(
      Parameter = names(coef_est),
      Estimate  = as.numeric(coef_est),
      Std.Error = as.numeric(se_est),
      z_value   = as.numeric(z_stat),
      p_value   = as.numeric(p_val)
    )

    coef_table$Signif <- ifelse(coef_table$p_value < 0.001, "***",
                                ifelse(coef_table$p_value < 0.01,  "**",
                                       ifelse(coef_table$p_value < 0.05,  "*",
                                              ifelse(coef_table$p_value < 0.1,   ".", ""))))

    cat("\n---------------- PARAMETER ESTIMATES ----------------------\n")
    coef_table_print <- coef_table
    num_cols <- sapply(coef_table_print, is.numeric)
    coef_table_print[num_cols] <- round(coef_table_print[num_cols], 4)
    print(coef_table_print, row.names = FALSE)
    cat("Significance codes: *** <0.001  ** <0.01  * <0.05  . <0.1\n")
    cat("------------------------------------------------------------\n")

  } else {
    cat("\nNo coefficients estimated (possible differencing-only model).\n")
  }

  # ==========================================================
  # 6. PERFORMANCE METRICS
  # ==========================================================

  resids <- as.numeric(residuals(final_model))
  fitted_val <- y - resids
  fc_future <- stats::predict(final_model,n.ahead=h)

  pred_test_final <- fitted_val[(split_point+1):n]
  rmse <- sqrt(mean((test_y - pred_test_final)^2))
  mae  <- mean(abs(test_y - pred_test_final))
  mape <- mean(abs((test_y - pred_test_final)/test_y))*100

  k <- length(final_model$coef)
  aic_val <- final_model$aic
  bic_val <- aic_val + k * (log(n) - 2)

  if(n > (k + 1)){
    aicc_val <- aic_val + (2*k^2 + 2*k) / (n - k - 1)
  } else {
    aicc_val <- NA
  }

  sst <- sum((y - mean(y))^2)
  sse <- sum(resids^2)
  r2  <- 1 - (sse/sst)
  adj_r2 <- 1 - ((1-r2)*(n-1)/(n-k-1))

  cat("\n----------------- MODEL PERFORMANCE -----------------------\n")
  cat("AIC / BIC / AICc :", round(aic_val,4), "/",
      round(bic_val,4), "/", round(aicc_val,4), "\n")
  cat("RMSE / MAE       :", round(rmse,4), "/",
      round(mae,4), "\n")
  cat("MAPE             :", round(mape,2), "%\n")
  cat("R-Square         :", round(r2,4), "\n")
  cat("Adj R-Square     :", round(adj_r2,4), "\n")
  cat("------------------------------------------------------------\n")

  # ==========================================================
  # 7. RESIDUAL DIAGNOSTICS
  # ==========================================================

  lag_val <- max(10, 2*s_period)
  lb_p <- Box.test(resids,
                   lag=lag_val,
                   type="Ljung-Box",
                   fitdf=sum(ord[c(1,3)]) +
                     sum(seas[c(1,3)]))$p.value

  dw_stat <- sum(diff(resids)^2)/sum(resids^2)

  skew <- mean((resids-mean(resids))^3)/(sd(resids)^3)
  kurt <- mean((resids-mean(resids))^4)/(sd(resids)^4)
  jb_stat <- n/6*(skew^2 + (1/4)*(kurt-3)^2)
  jb_p <- pchisq(jb_stat,df=2,lower.tail=FALSE)

  resids_jit <- jitter(resids)
  ks_p <- ks.test(resids_jit,
                  "pnorm",
                  mean(resids_jit),
                  sd(resids_jit))$p.value

  cat("\n---------------- RESIDUAL DIAGNOSTICS ---------------------\n")
  cat("Ljung-Box p      :", round(lb_p,4), "\n")
  cat("Durbin-Watson    :", round(dw_stat,4), "\n")
  cat("Jarque-Bera p    :", round(jb_p,4), "\n")
  cat("KS-Lilliefors p  :", round(ks_p,4), "\n")
  cat("------------------------------------------------------------\n")

  # ==========================================================
  # 8. FORECAST OUTPUT
  # ==========================================================

  z_val_ci <- 1.96
  upper_ci <- fc_future$pred + z_val_ci * fc_future$se
  lower_ci <- fc_future$pred - z_val_ci * fc_future$se

  forecast_table <- data.frame(
    Horizon  = 1:h,
    Forecast = as.numeric(fc_future$pred),
    Lower95  = as.numeric(lower_ci),
    Upper95  = as.numeric(upper_ci)
  )

  cat("\n---------------- FORECAST RESULT --------------------------\n")
  print(round(forecast_table,4))
  cat("------------------------------------------------------------\n")

  # ==========================================================
  # 6. VISUALISASI (IMPROVED FLEXIBLE Y-AXIS)
  # ==========================================================
  if (plot){
  layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE),
         heights = c(2, 1))
  par(mar=c(4,4,3,2))

  z_val_ci <- 1.96
  upper_ci <- fc_future$pred + z_val_ci * fc_future$se
  lower_ci <- fc_future$pred - z_val_ci * fc_future$se

  # --- AUTO Y RANGE DETECTION (MORE FLEXIBLE)
  y_all <- c(y,
             fitted_val,
             fc_future$pred,
             upper_ci,
             lower_ci)

  y_min <- min(y_all, na.rm=TRUE)
  y_max <- max(y_all, na.rm=TRUE)

  # Tambahkan padding 5% dari range
  padding <- 0.05 * (y_max - y_min)

  y_lim_final <- c(y_min - padding,
                   y_max + padding)

  plot(1:n, y, type="l", col="gray60",
       main=paste0("automlR SARIMA(",
                   paste(best_order$ord,collapse=","),")(",
                   paste(best_order$seas,collapse=","),")[",
                   s_period,"]"),
       xlab="Time Index",
       ylab=target_name,
       ylim=y_lim_final)

  lines(1:split_point,
        fitted_val[1:split_point],
        col="blue", lwd=2)

  lines((split_point+1):n,
        fitted_val[(split_point+1):n],
        col="red", lwd=2)

  lines((n+1):(n+h),
        fc_future$pred,
        col="orange", lwd=2)

  polygon(c((n+1):(n+h),
            rev((n+1):(n+h))),
          c(upper_ci,
            rev(lower_ci)),
          col=rgb(1,0.6,0,0.2),
          border=NA)

  lines((n+1):(n+h), upper_ci,
        col="darkorange", lty=2)

  lines((n+1):(n+h), lower_ci,
        col="darkorange", lty=2)

  abline(v=c(split_point,n), lty=3)

  legend("bottom",
         legend=c("Actual","Train Fit","Test Fit",
                  "Forecast","95% CI"),
         col=c("gray60","blue","red",
               "orange","darkorange"),
         lty=c(1,1,1,1,2),
         lwd=c(1,2,2,2,1),
         horiz=TRUE,
         bty="n",
         xpd=NA,
         inset=c(0,0.02),
         cex=0.7)

  acf(resids, main="Residual ACF")
  pacf(resids, main="Residual PACF")
  }

  cat("\n----------- OVERFITTING AND UNDERFITTING TEST -------------\n")

  # Overfitting / Underfitting summary
  if(overfit_flag){
    cat("[WARNING] Possible Overfitting detected.\n")
  } else if(underfit_flag){
    cat("[WARNING] Possible Underfitting detected.\n")
  } else {
    cat("Model appears free from Overfitting and Underfitting.\n")
  }

  # Stationarity consistency
  if(!is.na(stationary_flag)){
    if(stationary_flag==TRUE & d_val==0 & D_val==0){
      cat("Stationarity: Consistent with stationary series.\n")
    } else if(stationary_flag==FALSE & (d_val>0 | D_val>0)){
      cat("Stationarity: Consistent with non-stationary series.\n")
    } else {
      cat("Stationarity: Check consistency with differencing.\n")
    }
  } else {
    cat("Stationarity test unavailable.\n")
  }
  cat("------------------------------------------------------------\n")

  # ==========================================================
  # 9. INTERPRETATION
  # ==========================================================
  cat("\n------------------ INTERPRETATION ------------------------\n")

  cat("R^2 indicates", round(r2 * 100, 2),
      "% of variability explained.\n")

  cat("MAPE =", round(mape, 2), "%.\n")

  if (lb_p > 0.05) {
    cat("Residuals appear independent.\n")
  } else {
    cat("Residual autocorrelation is still present.\n")
  }

  if (jb_p > 0.05 && ks_p > 0.05) {
    cat("Residuals are approximately normally distributed.\n")
  } else {
    cat("Residuals deviate from normality.\n")
  }

  cat("============================================================\n")

  return(invisible(list(
    model=final_model,
    ranking=top10,
    forecast=forecast_table,
    metrics=list(
      AIC=aic_val,
      BIC=bic_val,
      AICc=aicc_val,
      RMSE=rmse,
      MAPE=mape,
      R2=r2,
      Adj_R2=adj_r2,
      LB_p=lb_p,
      JB_p=jb_p,
      KS_p=ks_p
    )
  )))
}

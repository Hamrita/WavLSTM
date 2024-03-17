##########################################
##########################################

#  load packages

library(crypto2)
library(ggplot2)
library(waveslim)
library(forecast)

# load data (BTC)

BTC=crypto_history(limit = 1)
BTC=BTC[,c("timestamp", "open","high", "low", "close", "volume")]
dd=as.Date(1:NROW(BTC), origin="2013-04-28")
BTC=cbind(Date=dd,BTC[,-1])
BTC_C=BTC[,c("Date", "close")]

# Graph of close price

g1=ggplot(BTC, aes(x=Date, y=close))+
  geom_line(col="blue",linewidth=0.7)
g1


################################################
################################################

# Split data into train and test data

split_ratio=0.85
n_tr=floor(NROW(BTC)*split_ratio)
n_tst=NROW(BTC)-n_tr

#############################################
#     ARIMA modelling
############################################

ar=auto.arima(BTC_C[1:n_tr,2],ic="aic")
ar1=Arima(BTC_C[1:n_tr,2], c(0,1,1) )
ar_pred=ar1$fitted
#ar_For=Arima(BTC_C[-(1:n_tr),2], c(0,1,1) )$fitted
ar_forecast=forecast(ar1,n_tst)$mean

# accuracy for arima model
Accuracy=function(actual, predicted){
  MAE=Metrics::mae(actual, predicted)
  MAPE=Metrics::mape(actual, predicted)
  RMSE=Metrics::rmse(actual, predicted)
  RAE=Metrics::rae(actual, predicted)
  data.frame(MAE=MAE,MAPE=MAPE,RMSE=RMSE,RAE=RAE)
}
acc_arima_tr=Accuracy(BTC_C[(1:n_tr ),2], ar_pred)
acc_arima_tst=Accuracy(BTC_C[-(1:n_tr ),2], ar_forecast)

acc_arima=rbind(Train=acc_arima_tr,Test=acc_arima_tst)

############################################
# Wavelet arima modelling
############################################

# Wavelet decomposition

# wavelet additive decomposition (MRA)
WD=function(x, wf = "la8", J = 4, method = "modwt", boundary = "periodic"){
  wd=mra(as.vector(x), wf=wf, J=J, method=method, boundary = boundary)
  coeff=matrix(unlist(wd), nc=(J+1))
  colnames(coeff)=names(wd)
  coeff
}

# select best wavelet filter and best level of decomposition respect to RMSE
wavSelect=function(x){
  wf=c("haar", "d4", "d6", "d8", "d16", "la8", "la16", "la8")
  J=2:floor(log(length(x),2))
  boundary=c("periodic","reflection")
  param=expand.grid(wf=wf,J=J,boundary=boundary)
  rmse=NULL
  for(i in 1:nrow(param)){
    dec=WD(x, wf=as.character(param$wf[i]),J=param$J[i],
           boundary=as.character(param$boundary[i]))
    rmse[i]=sqrt(mean((x-rowSums(dec))^2))
  }
  i_min=which.min(rmse)
  wf_best=param[i_min,]
  return(list(RMSE=rmse, wf_best=wf_best))
}

# Select the best wavelet filter and level of 
# decomposition respect to the RMSE

best=wavSelect(BTC$close)
J=best$wf_best[,2]
wf=as.character(best$wf_best[,1])
boundary=as.character(best$wf_best[,3])

WavDecomopsition=WD(BTC$close,wf=wf,J=J,boundary=boundary)
WWD=as.vector(WavDecomopsition)
Scale=rep(colnames(WavDecomopsition),each=NROW(WavDecomopsition))
dd=rep(BTC$Date, length(colnames(WavDecomopsition)))
DecData=data.frame(Date=dd,Level=WWD,Scale)

g2=ggplot(DecData, aes(Date, Level))+
  geom_line()+facet_grid(rows=vars(Scale))
g2


wav_tr=WaveletArima::WaveletFittingarma(BTC_C[1:n_tr,2],wf,J,
  MaxARParam = 5,MaxMAParam = 5, NForecast = n_tst)
acc_wav_tr=Accuracy(BTC_C[1:n_tr,2],wav_tr$FinalPrediction )
acc_wav_tst=Accuracy(BTC_C[-(1:n_tr),2],wav_tr$Finalforecast )

acc_wav=rbind(Train=acc_wav_tr, Test=acc_wav_tst )

##############################################
#  Wavelet ANN modelling (best: hidden=2)
##############################################

#  Tuning: hidden and nonseasLag
hidden=1:3
llag=1:10
parms=expand.grid(hidden,llag)
nn=NROW(parms)
rrmse=NULL
for(i in 1:3){
  for(j in 1:15){
  ww=WaveletANN::WaveletFittingann(BTC_C[1:n_tr,2],
      J,wf,nonseaslag = j,hidden = i,
      NForecast = n_tst)
  rr=Metrics::rmse(BTC_C[1:n_tr,2],ww$FinalPrediction )
  rrmse=c(rrmse,rr)
}}
idx=which.min(rrmse)
wavANN3_tr=WaveletANN::WaveletFittingann(BTC_C[1:n_tr,2],J,wf,
 hidden=3, NForecast = n_tst, nonseaslag = 10)
acc_wavANN3_tr=Accuracy(BTC_C[1:n_tr,2],wavANN3_tr$FinalPrediction )
acc_wavANN3_tst=Accuracy(BTC_C[-(1:n_tr),2],wavANN3_tr$Finalforecast )

acc_wavANN3=rbind(Train=acc_wavANN3_tr, Test=acc_wavANN3_tst )

##############################################
##############################################
tr=WavDecomopsition[1:n_tr,]
tst=WavDecomopsition[-(1:n_tr),]

# AR models for predicting details
D_predict=NULL
D_forecast=NULL
for(i in 1:(NCOL(tr))){
  fit=auto.arima(tr[,i])
  pred=fit$fitted
  forcst=forecast(fit, h=n_tst)
  D_predict=cbind(D_predict,pred)
  D_forecast=cbind(D_forecast,as.matrix(forcst$mean))
}
DD_forecast=rowSums(D_forecast)
DD_predict=rowSums(D_predict)
result_tr=list(fost=DD_forecast,pred=DD_predict)



ar_d1=auto.arima(tr[,"D1"], max.order = 10)
fit1
order1 <- arimaorder(ar_d1)
model_arima1<-arima(WavDecomopsition[,"D1"],order=c(order1[1], order1[2], order1[3]))
pred_arima1 <-ar_d1$fitted
forecast_arima1 <- data.frame(predict(ar_d1,n.ahead=n_tst))
forecast_arima1 <-forecast_arima1$pred


grid_search_results <- ts_lstm_x_tuning(
  X_train, y_train, X_val, y_val,
  embedded_colnames, custom_loss, early_stopping,
  n_lag = 2, # desired lag value
  lstm_units_list = c(32, 64,128),
  learning_rate_list = c(0.001, 0.01),
  batch_size_list = c(32, 64,128),
  dropout_list = c(0.1,0.2,0.3),
  l1_reg_list = c(0.001),
  l2_reg_list = c(0.001),
  n_iter = 10,
  n_verbose = 0 # or 1
)




#' @title Wavelet Based LSTM Model
#'
#' @param ts Time Series Data
#' @param MLag Maximum Lags
#' @param split_ratio Training and Testing Split
#' @param wlevels Wavelet Levels
#' @param epochs Number of epochs
#' @param LSTM_unit Number of unit in LSTM layer
#' @import caret dplyr caretForecast tseries stats wavelets TSLSTM
#' @return
#' \itemize{
#'   \item Train_actual: Actual train series
#'   \item Test_actual: Actual test series
#'   \item Train_fitted: Fitted train series
#'   \item Test_predicted: Predicted test series
#'   }
#' @export
#'
#' @examples
#' \donttest{
#'y<-rnorm(100,mean=100,sd=50)
#'WTSLSTM<-WaveletLSTM(ts=y)
#'}
#' @references
#' Paul, R.K. and Garai, S. (2021). Performance comparison of wavelets-based machine learning technique for forecasting agricultural commodity prices, Soft Computing, 25(20), 12857-12873

WaveletLSTM<-function(ts,MLag=12,split_ratio=0.8,wlevels=3,epochs=25,LSTM_unit=20){
  SigLags<-NULL
  SigLags<-function(Data,MLag){
    ts<-as.ts(na.omit(Data))
    adf1<-adf.test(na.omit(ts))
    if (adf1$p.value>0.05){
      ts<-ts
    } else {
      ts<-diff(ts)
    }
    adf2<-adf.test(ts)
    if (adf2$p.value>0.05){
      ts<-ts
    } else {
      ts<-diff(ts)
    }
    
    CorrRes<-NULL
    for (i in 1:MLag) {
      # i=1
      ts_y<-dplyr::lag(as.vector(ts), i)
      t<-cor.test(ts,ts_y)
      corr_res<-cbind(Corr=t$statistic,p_value=t$p.value)
      CorrRes<-rbind(CorrRes,corr_res)
    }
    rownames(CorrRes)<-seq(1:MLag)
    Sig_lags<-rownames(subset(CorrRes,CorrRes[,2]<=0.05))
    maxlag<-max(as.numeric(Sig_lags))
    return(list(Result=as.data.frame(CorrRes),SigLags=as.numeric(Sig_lags),MaxSigLag=maxlag))
  }
  ntest<-round(length(ts)*(1-split_ratio), digits = 0)
  Split1 <- caretForecast::split_ts(as.ts(ts), test_size = ntest)
  train_data1 <- Split1$train
  test_data1 <- Split1$test
  Wvlevels<-wlevels
  mraout <- wavelets::modwt(as.vector(ts), filter="haar", n.levels=Wvlevels)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  ts_fitted<-NULL
  ts_foreast<-NULL
  model_par=NULL
  for (j in 1:ncol(WaveletSeries)) {
    w<-as.ts(WaveletSeries[,j])
    maxl<-SigLags(Data=w,MLag = MLag)$MaxSigLag
    model<-TSLSTM::ts.lstm(ts=w,xreg = NULL,tsLag=maxl,xregLag = 0,LSTMUnits=LSTM_unit, Epochs=epochs,SplitRatio =split_ratio)
    model_par<-rbind(model_par,model$Param)
    ts_fitted<-model$TrainFittedValue
    ts_foreast<-model$TestPredictedValue
  }
  
  trainf <- apply(ts_fitted,1,sum)
  testf <- apply(ts_foreast,1,sum)
  return(list(Train_actual=train_data1,Test_actual=test_data1,Train_fitted=trainf,
              test_pred=testf))
}

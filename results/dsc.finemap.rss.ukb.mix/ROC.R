library(dplyr)
library(scam)
## Functions
roc_data = function(d1, cutoff = c(pip_cutoff, 1), connect_org = FALSE) {
  grid = 200
  ttv = seq(1:grid)/grid
  ttv = ttv[which(ttv>=cutoff[1] & ttv<=cutoff[2])]
  # see SuSiE-Manuscript issue 2
  d1 = d1[order(d1[,1]), ]
  end = tail(d1[which(d1[,2] == 0),][,1],1)
  # ttv = c(ttv[-length(ttv)], min(ttv[length(ttv)], end))
  ttv = ttv[-length(ttv)]
  # end of issue 2
  rst1 = t(sapply(ttv, function(x) c(sum(d1[,2][d1[,1]>=x]), length(d1[,2][d1[,1]>=x]))))
  rst1 = cbind(rst1, sum(d1[,2]), sum(1-d1[,2]))
  # colnames(rst1) = c('true_positive', 'total_positive', 'total_signal','total_null')
  if (connect_org) {
    # connect to origin
    last_row = tail(rst1, 1)
    rst1 = rbind(rst1, c(last_row[1], last_row[2]-1, last_row[3]), c(0.001,0.001,last_row[3]))
  }
  rst1 = as.data.frame(rst1)
  colnames(rst1) = c('true_positive', 'total_positive', 'total_signal', 'total_null')
  if (connect_org) {
    rst2 = as.data.frame(cbind(rst1$true_positive / rst1$total_positive, rst1$true_positive / rst1$total_signal,  c(ttv, ttv[length(ttv)], 1)))
  } else {
    rst2 = as.data.frame(cbind(rst1$true_positive / rst1$total_positive, rst1$true_positive / rst1$total_signal,  ttv, 
                               (rst1$total_positive - rst1$true_positive)/rst1$total_null))
  }
  colnames(rst2) = c('Precision', 'Recall', 'Threshold', 'FPR')
  return(list(counts = rst1, rates = rst2))
}


create_chunks = function(item, n) {
  splitted = suppressWarnings(split(item, 1:n))
  return(c(splitted[[1]], splitted[[length(splitted)]][length(splitted[[length(splitted)]])]))
}
make_smooth = function(x,y,subset=chunks, sm = smooth) {
  if (sm) {
    if (subset < length(x) && subset > 0) {
      x = create_chunks(x, subset)
      y = create_chunks(y, subset)
    }
    dat = data.frame(cbind(x,y))
    colnames(dat) = c('x','y')
    y=predict(scam(y ~ s(x, bs = "mpi"), data = dat))
  }
  return(list(x=x,y=y))
}
add_text = function(thresholds, x, y, threshold, color, delta = 0.015, y.delta=0, print.text = T) {
  idx = which(thresholds == threshold)
  if(print.text){
    text(x[idx] - delta, y[idx] +y.delta, labels = threshold, col = color)
  }
  points(x[idx],y[idx], col = color, cex = 1.1, lwd=2.5)
}

## parameters
estimate_resid = c(TRUE, FALSE)
all.comb = expand.grid(estimate_resid)
colnames(all.comb) = c('estimate_resid')

pip_cutoff = 0.3

chunks = 0
smooth = FALSE
colors = c('cyan', '#7A68A6', '#348ABD', '#A60628', '#FF0000', '#188487', '#E2A233',
           '#A9A9A9', '#000000', '#FF00FF', '#FFD700', '#ADFF2F', '#00FFFF')

for (case in 1:nrow(all.comb)){
  estimate_resid = all.comb[case, 'estimate_resid']
  input = paste0('susie_rss_ukb_mix_pip_extraction/susie_rss_ukb_mix_pip_ER',
                 estimate_resid,'.rds')
  output2_name = paste0('susie_rss_ukb_mix_roc/susie_rss_ukb_mix_pip_ER',
                        estimate_resid, '_roc_2')
  output2 = paste0(output2_name,'.rds')
  
  print("Computing ROC data ...")
  result = readRDS(input)
  susierss_lm_insample = roc_data(do.call(rbind, lapply(1:length(result), function(i) cbind(result[[i]]$susierss_lm_insample, result[[i]]$is_signal))))
  susierss_mix_insample = roc_data(do.call(rbind, lapply(1:length(result), function(i) cbind(result[[i]]$susierss_mix_insample, result[[i]]$is_signal))))
  susierss_lm_insamplez = roc_data(do.call(rbind, lapply(1:length(result), function(i) cbind(result[[i]]$susierss_lm_insamplez, result[[i]]$is_signal))))
  susierss_mix_insamplez = roc_data(do.call(rbind, lapply(1:length(result), function(i) cbind(result[[i]]$susierss_mix_insamplez, result[[i]]$is_signal))))
  
  saveRDS(list(susierss_lm_insample = susierss_lm_insample, susierss_mix_insample = susierss_mix_insample,
               susierss_lm_insamplez = susierss_lm_insamplez, susierss_mix_insamplez = susierss_mix_insamplez), output2)
  
  ## plots
  d1 = readRDS(output2)
  pdf(paste0(output2_name,"_power.pdf"), width=5, height=5, pointsize=15)
  yy = make_smooth(1 - d1$susierss_lm_insample$rates$Precision, d1$susierss_lm_insample$rates$Recall)
  plot(yy$x, yy$y, t="l", col='black', ylab = "Power", xlab ="False Discovery Rate", main = "2 effect variables",
       bty='l', lwd = 2, xlim = c(0,0.5), ylim = c(0,0.5))
  add_text(d1$susierss_lm_insample$rates$Threshold, yy$x, yy$y, 0.95, 'black', 0.01, 0.01)
  
  ## power vs FDR
  
  yy = make_smooth(1-d1$susierss_mix_insample$rates$Precision, d1$susierss_mix_insample$rates$Recall)
  lines(yy$x, yy$y, col = colors[5], lwd = 2, xlim = c(0,0.5), ylim = c(0,0.5))
  add_text(d1$susierss_mix_insample$rates$Threshold, yy$x, yy$y, 0.95, colors[5], 0.00, -0.02)

  yy = make_smooth(1-d1$susierss_lm_insamplez$rates$Precision, d1$susierss_lm_insamplez$rates$Recall)
  lines(yy$x, yy$y, col = 'green2', lwd = 2, xlim = c(0,0.5), ylim = c(0,0.5))
  add_text(d1$susierss_lm_insamplez$rates$Threshold, yy$x, yy$y, 0.95, 'green2', 0.00, -0.02)

  yy = make_smooth(1-d1$susierss_mix_insamplez$rates$Precision, d1$susierss_mix_insamplez$rates$Recall)
  lines(yy$x, yy$y, col = colors[7], lwd = 2, xlim = c(0,0.5), ylim = c(0,0.5))
  add_text(d1$susierss_mix_insamplez$rates$Threshold, yy$x, yy$y, 0.95, colors[7], 0.00, -0.02)

  legend("topleft", legend=c("lm insample", "mix insample", "lm insample addz", "mix insample addz"),
         col=c('black', colors[5], 'green2', colors[7]), lty=1, cex=0.8)
  dev.off()
  system(paste0("convert -flatten -density 120 ", output2_name, '_power.pdf', " ", output2_name, '_power.png'))
  
  ## TPR vs FPR
  yy = make_smooth(d1$susierss_lm_insample$rates$FPR, d1$susierss_lm_insample$rates$Recall)
  pdf(paste0(output2_name,"_roc.pdf"), width=5, height=5, pointsize=15)
  plot(yy$x, yy$y, t="l", col='black', ylab = "True Discovery Rate", xlab ="False Positive Rate", main = "2 effect variables",
       bty='l', lwd = 2, xlim = c(0,0.005), ylim = c(0,0.3))
  add_text(d1$susierss_lm_insample$rates$Threshold, yy$x, yy$y, 0.95, 'black', -0.015, print.text = F)
  
  yy = make_smooth(d1$susierss_mix_insample$rates$FPR, d1$susierss_mix_insample$rates$Recall)
  lines(yy$x, yy$y, col = colors[5], lwd = 2, xlim = c(0,0.005), ylim = c(0,0.3))
  add_text(d1$susierss_mix_insample$rates$Threshold, yy$x, yy$y, 0.95, colors[5], -0.015, print.text = F)
  
  yy = make_smooth(d1$susierss_lm_insamplez$rates$FPR, d1$susierss_lm_insamplez$rates$Recall)
  lines(yy$x, yy$y, col = 'green2', lwd = 2, xlim = c(0,0.005), ylim = c(0,0.3))
  add_text(d1$susierss_lm_insamplez$rates$Threshold, yy$x, yy$y, 0.95, 'green2', -0.015, print.text = F)
  
  yy = make_smooth(d1$susierss_mix_insamplez$rates$FPR, d1$susierss_mix_insamplez$rates$Recall)
  lines(yy$x, yy$y, col = colors[7], lwd = 2, xlim = c(0,0.005), ylim = c(0,0.3))
  add_text(d1$susierss_mix_insamplez$rates$Threshold, yy$x, yy$y, 0.95, colors[7], -0.015, print.text = F)

  legend("topleft", legend=c("lm insample", "mix insample", "lm insample addz", "mix insample addz"),
         col=c('black', colors[5], 'green2', colors[7]), lty=1, cex=0.8)
  dev.off()
  system(paste0("convert -flatten -density 120 ", output2_name, '_roc.pdf', " ", output2_name, '_roc.png'))
}

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
ldmethod = c('in_sample_Z', 'ref_sample_Z')
all.comb = expand.grid(ldmethod)
colnames(all.comb) = c('ldmethod')

pip_cutoff = 0.3

chunks = 0
smooth = FALSE
colors = c('cyan', '#7A68A6', '#348ABD', '#A60628', '#FF0000', '#188487', '#E2A233',
           '#A9A9A9', '#000000', '#FF00FF', '#FFD700', '#ADFF2F', '#00FFFF')

for (case in 1:nrow(all.comb)){
  ldmethod = all.comb[case, 'ldmethod']
  input = paste0('susie_rss_ukb_weight_pip_extraction/susie_rss_ukb_weight_pip_ld',ldmethod,'.rds')
  output_name = paste0('susie_rss_ukb_weight_roc/susie_rss_ukb_weight_ld',ldmethod)
  output = paste0(output_name,'.rds')
  
  print("Computing ROC data ...")
  result = readRDS(input)
  rss_z0 = roc_data(cbind(result$RSS_z0, result$is_signal))
  rss_z001 = roc_data(cbind(result$RSS_z001, result$is_signal))
  rss_z002 = roc_data(cbind(result$RSS_z002, result$is_signal))
  rss_z005 = roc_data(cbind(result$RSS_z005, result$is_signal))
  rss_z01 = roc_data(cbind(result$RSS_z01, result$is_signal))
  susie = roc_data(cbind(result$susie, result$is_signal))
  
  saveRDS(list(susie = susie, rss_z0 = rss_z0, rss_z001 = rss_z001, rss_z002 = rss_z002,
               rss_z005 = rss_z005, rss_z01 = rss_z01), output)
  
  ## plots
  d = readRDS(output)
  pdf(paste0(output_name,"_power.pdf"), width=5, height=5, pointsize=15)
  par(xpd=T, mar=par()$mar+c(0,0,0,4))
  yy = make_smooth(1 - d$susie$rates$Precision, d$susie$rates$Recall)
  plot(yy$x, yy$y, t="l", col='black', ylab = "Power", xlab ="False Discovery Rate", main = "2 effect variables",
       bty='L', lwd = 2, xlim = c(0,0.4), ylim = c(0,0.4))
  add_text(d$susie$rates$Threshold, yy$x, yy$y, 0.95, 'black', 0.01, 0.01, print.text = F)
  
  ## power vs FDR
  yy = make_smooth(1-d$rss_z0$rates$Precision, d$rss_z0$rates$Recall)
  lines(yy$x, yy$y, col = colors[5], lwd = 2, xlim = c(0,0.4), ylim = c(0,0.4))
  add_text(d$rss_z0$rates$Threshold, yy$x, yy$y, 0.95, colors[5], 0.00, -0.02, print.text = F)
  
  yy = make_smooth(1-d$rss_z001$rates$Precision, d$rss_z001$rates$Recall)
  lines(yy$x, yy$y, col = colors[2], lwd = 2, xlim = c(0,0.4), ylim = c(0,0.4))
  add_text(d$rss_z001$rates$Threshold, yy$x, yy$y, 0.95, colors[2], 0.00, -0.02, print.text = F)
  
  yy = make_smooth(1-d$rss_z002$rates$Precision, d$rss_z002$rates$Recall)
  lines(yy$x, yy$y, col = colors[3], lwd = 2, xlim = c(0,0.4), ylim = c(0,0.4))
  add_text(d$rss_z002$rates$Threshold, yy$x, yy$y, 0.95, colors[3], 0.00, -0.02, print.text = F)
  
  yy = make_smooth(1-d$rss_z005$rates$Precision, d$rss_z005$rates$Recall)
  lines(yy$x, yy$y, col = colors[4], lwd = 2, xlim = c(0,0.4), ylim = c(0,0.4))
  add_text(d$rss_z005$rates$Threshold, yy$x, yy$y, 0.95, colors[4], 0.00, -0.02, print.text = F)
  
  yy = make_smooth(1-d$rss_z01$rates$Precision, d$rss_z01$rates$Recall)
  lines(yy$x, yy$y, col = colors[6], lwd = 2, xlim = c(0,0.4), ylim = c(0,0.4))
  add_text(d$rss_z01$rates$Threshold, yy$x, yy$y, 0.95, colors[6], 0.00, -0.02, print.text = F)
  
  legend("topright", inset=c(-0.6,0), 
         legend=c("SuSiE", "w = 0", "w = 0.001", "w = 0.002", "w = 0.005", "w = 0.01"),
         col=c('black', colors[5], colors[2], colors[3], colors[4], colors[6]), 
         lty=1, cex=0.8)
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  system(paste0("convert -flatten -density 120 ", output_name, '_power.pdf', " ", output_name, '_power.png'))
  
  ## TPR vs FPR
  yy = make_smooth(d$susie$rates$FPR, d$susie$rates$Recall)
  pdf(paste0(output_name,"_roc.pdf"), width=5, height=5, pointsize=15)
  par(xpd=T, mar=par()$mar+c(0,0,0,4))
  plot(yy$x, yy$y, t="l", col='black', ylab = "True Discovery Rate", xlab ="False Positive Rate", main = "2 effect variables",
       bty='l', lwd = 2, xlim = c(0,0.0005), ylim = c(0,0.4))
  add_text(d$susie$rates$Threshold, yy$x, yy$y, 0.95, 'black', -0.015, print.text = F)
  
  yy = make_smooth(d$rss_z0$rates$FPR, d$rss_z0$rates$Recall)
  lines(yy$x, yy$y, col = colors[5], lwd = 2, xlim = c(0,0.0005), ylim = c(0,0.4))
  add_text(d$rss_z0$rates$Threshold, yy$x, yy$y, 0.95, colors[5], -0.015, print.text = F)
  
  yy = make_smooth(d$rss_z001$rates$FPR, d$rss_z001$rates$Recall)
  lines(yy$x, yy$y, col = colors[2], lwd = 2, xlim = c(0,0.0005), ylim = c(0,0.4))
  add_text(d$rss_z001$rates$Threshold, yy$x, yy$y, 0.95, colors[2], -0.015, print.text = F)
  
  yy = make_smooth(d$rss_z002$rates$FPR, d$rss_z002$rates$Recall)
  lines(yy$x, yy$y, col = colors[3], lwd = 2, xlim = c(0,0.0005), ylim = c(0,0.4))
  add_text(d$rss_z002$rates$Threshold, yy$x, yy$y, 0.95, colors[3], -0.015, print.text = F)
  
  yy = make_smooth(d$rss_z005$rates$FPR, d$rss_z005$rates$Recall)
  lines(yy$x, yy$y, col = colors[4], lwd = 2, xlim = c(0,0.0005), ylim = c(0,0.4))
  add_text(d$rss_z005$rates$Threshold, yy$x, yy$y, 0.95, colors[4], -0.015, print.text = F)
  
  yy = make_smooth(d$rss_z01$rates$FPR, d$rss_z01$rates$Recall)
  lines(yy$x, yy$y, col = colors[6], lwd = 2, xlim = c(0,0.0005), ylim = c(0,0.4))
  add_text(d$rss_z01$rates$Threshold, yy$x, yy$y, 0.95, colors[6], -0.015, print.text = F)
  
  legend("topright", inset=c(-0.6,0), 
         legend=c("SuSiE", "w = 0", "w = 0.001", "w = 0.002", "w = 0.005", "w = 0.01"),
         col=c('black', colors[5], colors[2], colors[3], colors[4], colors[6]), 
         lty=1, cex=0.8)
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  system(paste0("convert -flatten -density 120 ", output_name, '_roc.pdf', " ", output_name, '_roc.png'))
}

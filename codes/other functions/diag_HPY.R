# read diagnostics 
library(R.matlab)
library(tidyverse)
library(coda)

data = readMat('/Users/felpo/Dropbox/J100_diagnostics_10000.mat')

# there are in total ncol(M_par) parametri
# nrow(M_par) = Nruns * samples each Run
# Nruns is the number of additional samples of the algorithm (addsample in matlab)
# Samples each Run (N_iter in matlab)
N_runs = 50; N_iter = 10000
M_par = data$M.parametri.storage %>% 
  array(.,dim = c(N_iter, ncol(data$M.parametri.storage),N_runs))
dim(M_par)
M_par %>% head()

M_par_last = data$M.parametri
dim(M_par_last)
plot(M_par_last[,2],ty = "l")
hist(M_par_last[,10])

# dim(M_par)
(M_par[,,50] == M_par_last) %>% mean()

# ogni 10 
plot(M_par_last[seq(1,10000,by = 5),10],ty = "l")
M_par_last[,10] %>% unique() %>% length()
M_par_last[,10] %>% table()

# running mean 
cum_sum = M_par_last[,10] %>% cumsum()
plot(cum_sum/1:length(cum_sum),type = "l")

# Effective sampe size
apply(M_par_last, 2, effectiveSize) 
apply(M_par_last, 2, effectiveSize) %>% mean()
apply(M_par_last,2, function(x) x %>% unique() %>% length())

# Effective sampe size for all the Parameters
eff_size = apply(M_par, c(2,3), effectiveSize) 
eff_size %>% dim()
eff_size_mean = apply(eff_size, 2, mean)
hist(eff_size_mean)
eff_size %>% mean()
plot(1:50,eff_size_mean, ty = "p", c)

# Effective sampe size for all the Parameters
uniq_length = apply(M_par, c(2,3), function(x) unique(x) %>% length()) 
hist(uniq_length)
df = tibble(x = uniq_length %>% as.vector())
ggplot(df, aes(x = x)) + geom_histogram() + xlab("") + ggtitle("Histogram of Unique Particles")+
  theme(plot.title = element_text(hjust = 0.5, size = 13))

# PLOT ACF
library(latex2exp)
par(mfrow = c(2,2))
acf(M_par[,101,50], lag.max = 10, main = TeX('$\\sigma$'))
acf(M_par[,1,50], lag.max = 10, main = TeX('$\\theta$'))
acf(M_par[,102,50], lag.max = 10, main = TeX('$\\sigma_1$'))
acf(M_par[,2,50], lag.max = 10, main = TeX('$\\theta_1$'))
dev.off()

# PLOT ACF GGPLOT
library(forecast)
library(tidyverse)
sig_plot = forecast::Acf(M_par[,101,50])  %>% autoplot +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\sigma$'))
theta_plot = forecast::Acf(M_par[,1,50])  %>% autoplot +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\theta$'))
sig1_plot = forecast::Acf(M_par[,102,50])  %>% autoplot+
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\sigma_1$'))
theta1_plot = forecast::Acf(M_par[,2,50])  %>% autoplot +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\theta_1$'))
library(gridExtra);library(grid);library(lattice)
grid.arrange(sig_plot, theta_plot, sig1_plot, theta1_plot, ncol=2,
             top="Autocorrelation Plots")


# TRACEPLOT 
library(latex2exp)
par(mfrow = c(2,2))
plot(M_par[seq(1,10000,by = 5),101,50], main = TeX('$\\sigma$'),type = "l")
plot(M_par[seq(1,10000,by = 5),1,50], main = TeX('$\\theta$'),type = "l")
plot(M_par[seq(1,10000,by = 5),102,50], main = TeX('$\\sigma_1$'),type = "l")
plot(M_par[seq(1,10000,by = 5),2,50], main = TeX('$\\theta_1$'),type = "l")

# TRACEPLOTS GGPLOT
df = tibble(sig = M_par[seq(1,10000,by = 5),101,50],
            theta = M_par[seq(1,10000,by = 5),1,50],
            sig1 = M_par[seq(1,10000,by = 5),102,50],
            theta1 = M_par[seq(1,10000,by = 5),2,50],
            row_n = 1:length(sig))
sig_plot = ggplot(df) + geom_line(aes(x = row_n, y = sig))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\sigma$')) + xlab("")  + ylab("")+
  geom_hline(yintercept = mean(df$sig), color = "red", linetype ="dashed")
theta_plot = ggplot(df) + geom_line(aes(x = row_n, y = theta))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\theta$')) + xlab("")  + ylab("")+
  geom_hline(yintercept = mean(df$theta), color = "red", linetype ="dashed")
sig1_plot = ggplot(df) + geom_line(aes(x = row_n, y = sig1))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\sigma_1$')) + xlab("")  + ylab("")+
  geom_hline(yintercept = mean(df$sig1), color = "red", linetype ="dashed")
theta1_plot = ggplot(df) + geom_line(aes(x = row_n, y = theta1))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\theta_1$')) + xlab("")  + ylab("")+
  geom_hline(yintercept = mean(df$theta1), color = "red", linetype ="dashed")
grid.arrange(sig_plot, theta_plot, sig1_plot, theta1_plot, ncol=2,
             top="Trace Plots")

# Running mean ggplo
df = tibble(sig = cumsum(M_par[seq(1,10000,by = 5),101,50])/1:2000,
            theta = cumsum(M_par[seq(1,10000,by = 5),1,50])/1:2000,
            sig1 = cumsum(M_par[seq(1,10000,by = 5),102,50])/1:2000,
            theta1 = cumsum(M_par[seq(1,10000,by = 5),2,50])/1:2000,
            row_n = 1:length(sig))
sig_plot = ggplot(df) + geom_line(aes(x = row_n, y = sig))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\sigma$')) + xlab("")  + ylab("") +
  geom_hline(yintercept = mean(df$sig), color = "red", linetype ="dashed") 
theta_plot = ggplot(df) + geom_line(aes(x = row_n, y = theta))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\theta$')) + xlab("")  + ylab("")+
  geom_hline(yintercept = mean(df$theta), color = "red", linetype ="dashed") 
sig1_plot = ggplot(df) + geom_line(aes(x = row_n, y = sig1))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\sigma_1$')) + xlab("")  + ylab("")+
  geom_hline(yintercept = mean(df$sig1), color = "red", linetype ="dashed") 
theta1_plot = ggplot(df) + geom_line(aes(x = row_n, y = theta1))  +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) +
  ggtitle(TeX('$\\theta_1$')) + xlab("")  + ylab("")+
  geom_hline(yintercept = mean(df$theta1), color = "red", linetype ="dashed") 
grid.arrange(sig_plot, theta_plot, sig1_plot, theta1_plot, ncol=2,
             top="Running Mean Plots")



apply(M_par[29,,],2, function(x) x %>% unique() %>% length())
apply(M_par[29,,], 2, effectiveSize)
plot(M_par[,29,29],ty = "l")
acf(M_par[,29,5], lag.max = 50)
acf(M_par_last[,12], lag.max = 50)



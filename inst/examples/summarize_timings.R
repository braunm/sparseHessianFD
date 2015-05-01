library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

theme_set(theme_bw())

data.file <- "inst/examples/timings3.Rdata"
load(data.file)

tab <- mutate(res, ms=bench.time/1000000) %>%
  dcast(N+k+T+ord+bench.rep+ncolors~bench.expr, value.var="ms")  %>%
  mutate(nvars=N*k+k, hess_f=hess/f, hess_df=hess/df,
         hess_k=hess/k, hess_ncol = hess/ncolors, total=setup+hess,
         total_f=total/f, total_df = total/df, df/nvars)

tab2 <- gather(tab, stat, ms, c(setup:hess,hess_f:total_df)) %>%
  group_by(N, k, T, ord, stat, nvars) %>%
  summarize(mean=mean(ms), sd=sd(ms), lower=mean-1.96*sd, upper=mean+1.96*sd)


D2 <- subset(data.frame(tab2), stat %in% c("f", "df", "hess",
                                           "colors", "setup","hess_df"))
D2$stat <- revalue(D2$stat, c("f"="Function", "df"="Gradient", "hess"="Hessian",
                              "colors"="Partitioning",
                              "setup"="Initialization",
                              "hess_df"="Hessian/Gradient"))
D2$stat <- relevel(D2$stat, "Hessian/Gradient") %>%
  relevel("Initialization") %>%
  relevel("Partitioning") %>%
  relevel("Hessian") %>%
  relevel("Gradient") %>%
  relevel("Function")



P2 <- ggplot(D2, aes(x=N,y=mean, color=as.factor(k), linetype=as.factor(k))) %>%
  + geom_line(size=.4) %>%
  + scale_x_continuous("Number of heterogeneous units") %>%
  + scale_y_continuous("Computation time (milliseconds)") %>%
  + guides(color=guide_legend("k"), linetype=guide_legend("k")) %>%
  + facet_wrap(~stat, scales="free") %>%
  + theme(text=element_text(size=8))


pdf(file="vignettes/timings.pdf", width=6.5, height=4)
print(P2)
dev.off()


## P3 <- ggplot(filter(D2,k==4), aes(x=N)) %>%
##   + geom_line(aes(y=mean), linetype=1) %>%
##   + geom_line(aes(y=upper), linetype=3) %>%
##   + geom_line(aes(y=lower), linetype=3) %>%
##   + facet_wrap(~stat, scales="free", nrow=2, ncol=3)





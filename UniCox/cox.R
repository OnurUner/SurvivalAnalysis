library(R.matlab)
library("survival")

mat <- readMat('C:\\r workspace\\MultiSconES\\data\\breast_dataset.mat')
X_mat <- mat$data
cencoring <- mat$cencoring + 1
column.names <- lapply(c(mat$x.columns), function(x){paste('A', trimws(x), sep = '')})

rownames(X_mat) <- c(mat$x.rows)
colnames(X_mat) <- column.names
breast <- as.data.frame(X_mat)
Y <- mat$labels
breast$cencoring <- c(cencoring)
breast$time <- c(Y[,6])

genes = c(column.names)
res.cox <- coxph(Surv(time, cencoring) ~ genes[1], data = breast)

covariates <- c(column.names)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, cencoring)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = breast)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                       "p.value")
                         res <- c(p.value)
                         names(res) <- c("p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- c(univ_results)
lapply(res, function(x) write.table( data.frame(x), 'cox_pvalues.csv'  , append= T, sep=',' ))
write.csv(res, file = "C:\\r workspace\\MultiSconES\\data\\unicox_scores.csv",row.names=FALSE, col.names = FALSE)
#res <- t(as.data.frame(univ_results, check.names = FALSE))
#as.data.frame(res)

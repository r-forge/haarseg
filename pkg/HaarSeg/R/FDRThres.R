FDRThres <- function(x, q, sdev) {
	M = length(x);
	if (M < 2) { 
		T = 0;
	} else {
		m = (1:M) / M;
		sortedX = sort(abs(x),decreasing = TRUE);	
		p = 2*(1 - pnorm(sortedX, sd = sdev));
		k = which(p <= m*q);
		k = k[length(k)];
		if (length(k) == 0) {
			T = sortedX[1] + 1e-16;  #2^-52 is like MATLAB "eps"
		} else {
			T = sortedX[k];
		}
	}
}#FDRThres


##############################################################################
# HISTORY:
# 2007-12-17 [HB]
# o Created.
############################################################################## 

SegmentByPeaks <- function(data, peaks, weights = 0) {
	st = c(1,peaks);
	ed = c(peaks-1,length(data));

	segs = data;
	for (k in 1:length(st)) {
		if (length(weights) > 1) {
			segs[st[k]:ed[k]] = sum(weights[st[k]:ed[k]] * data[st[k]:ed[k]]) / sum(weights[st[k]:ed[k]]);
		} else {
			segs[st[k]:ed[k]] = mean(data[st[k]:ed[k]]);
		}
	}#for k
	return (segs);
}#SegmentByPeaks 


##############################################################################
# HISTORY:
# 2007-12-17 [HB]
# o Created.
############################################################################## 

real.data = c(rep.int(0,2000),rep.int(1,100),rep.int(0,2000));
noisy.data = real.data + rnorm(length(real.data),sd = 0.2);
plot(noisy.data)

# using default parameters
seg.data = haarSeg(noisy.data);

#segments result table: segment start index | segment size | segment value
print(seg.data$SegmentsTable)

# the complete segmented signal
lines(seg.data$Segmented, col="red", lwd=3)

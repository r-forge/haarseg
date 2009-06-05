%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  haarSeg.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{haarSeg}
\alias{haarSeg}


\title{Performs segmentation according to the HaarSeg algorithm}

\description{
  Performs segmentation according to the HaarSeg algorithm.
   HaarSeg segmentation is based on detecting local maxima in the wavelet
   domain, using Haar wavelet.  The main algorithm parameter is breaksFdrQ,
   which controls the sensitivity of the segmentation result.
   This function includes several optional extentions, supporting the use
   of weights (also known as quality of measurments) and raw measurments.
   We recommend using both extentions where possible, as it greatly
   improves the segmentation result.
   Raw red / green measurments are used to detect low value probes,
   which are more sensitive to noise.
}

\usage{haarSeg(I, W=vector(), rawI=vector(), chromPos=matrix(c(1, length(I)), nrow = 1, ncol = 2), breaksFdrQ=0.001, haarStartLevel=1, haarEndLevel=5)}

\arguments{
 \item{I}{a single array of log(R/G) measurements, sorted according to
   their genomic location.}
 \item{W}{Weight matrix, corresponding to quality of measurment.
   Insert \eqn{1/(\sigma^2)} as weights if your platform output
   \eqn{\sigma} as the quality of measurment. W must have the same
   size as I.}
 \item{rawI}{The mininum between the raw red and raw green measurment
   (before applying log ratio, but after any background reduction
   and/or normalization).
   rawI is used for the non-stationary variance compensation.
   rawI must have the same size as I.}
 \item{chromPos}{A matrix of two columns. The first column is the start
   index of each chromosome. The second column is the end index of each
   chromosome.}
 \item{breaksFdrQ}{The FDR q parameter. This value should lie between
   0 and 0.5. The smaller this value is, the less sensitive the
   segmentation result will be.
   For example, we will detect less breaks in the segmentation result when
   using Q = 1e-4, compared to the amounts of breaks when using Q = 1e-3.
   Common used values are 1e-2, 1e-3, 1e-4. Default value is 1e-3.}
 \item{haarStartLevel}{The detail subband from which we start to detect
   peaks. The higher this value is, the less sensitive we are to short
   segments. The default is value is 1, corresponding to segments of 2
   probes.}
 \item{haarEndLevel}{The detail subband until which we use to detect
   peaks. The higher this value is, the more sensitive we are to large
   trends in the data. This value DOES NOT indicate the largest possible
   segment that can be detected.  The default is value is 5,
   corresponding to step of 32 probes in each direction.}
}

\value{
 A \code{\link[base]{list}} containing two elements:
  \item{SegmentsTable}{Segments result table:
        (segment start index, segment size, segment value)}
  \item{Segmented}{The complete segmented signal (same size as I).}
}

\examples{
real.data = c(rep.int(0,2000),rep.int(1,100),rep.int(0,2000));
noisy.data = real.data + rnorm(length(real.data),sd = 0.2);
plot(noisy.data)

# using default parameters
seg.data = haarSeg(noisy.data);

#segments result table: segment start index | segment size | segment value
print(seg.data$SegmentsTable)

# the complete segmented signal
lines(seg.data$Segmented, col="red", lwd=3)
}

\author{Erez Ben-Yaacov}



\keyword{iteration}
\keyword{logic}

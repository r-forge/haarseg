###########################################################################/**
# @RdocFunction haarSeg
#
# @title "Performs segmentation according to the HaarSeg algorithm"
#
# \description{
#   @get "title". 
#    HaarSeg segmentation is based on detecting local maxima in the wavelet 
#    domain, using Haar wavelet.  The main algorithm parameter is breaksFdrQ,
#    which controls the sensitivity of the segmentation result. 
#    This function includes several optional extentions, supporting the use
#    of weights (also known as quality of measurments) and raw measurments.
#    We recommend using both extentions where possible, as it greatly 
#    improves the segmentation result. 
#    Raw red / green measurments are used to detect low value probes, 
#    which are more sensitive to noise.
# }
#
# @synopsis
#
# \arguments{
#  \item{I}{a single array of log(R/G) measurements, sorted according to
#    their genomic location.}
#  \item{W}{Weight matrix, corresponding to quality of measurment. 
#    Insert \eqn{1/(\sigma^2)} as weights if your platform output 
#    \eqn{\sigma} as the quality of measurment. W must have the same 
#    size as I.}
#  \item{rawI}{The mininum between the raw red and raw green measurment 
#    (before applying log ratio, but after any background reduction 
#    and/or normalization).  
#    rawI is used for the non-stationary variance compensation. 
#    rawI must have the same size as I.}
#  \item{chromPos}{A matrix of two columns. The first column is the start
#    index of each chromosome. The second column is the end index of each
#    chromosome.}
#  \item{breaksFdrQ}{The FDR q parameter. This value should lie between
#    0 and 0.5. The smaller this value is, the less sensitive the 
#    segmentation result will be.
#    For example, we will detect less breaks in the segmentation result when 
#    using Q = 1e-4, compared to the amounts of breaks when using Q = 1e-3. 
#    Common used values are 1e-2, 1e-3, 1e-4. Default value is 1e-3.}
#  \item{haarStartLevel}{The detail subband from which we start to detect
#    peaks. The higher this value is, the less sensitive we are to short 
#    segments. The default is value is 1, corresponding to segments of 2 
#    probes.}
#  \item{haarEndLevel}{The detail subband until which we use to detect
#    peaks. The higher this value is, the more sensitive we are to large
#    trends in the data. This value DOES NOT indicate the largest possible
#    segment that can be detected.  The default is value is 5, 
#    corresponding to step of 32 probes in each direction.}
# }
#
# \value{
#  A @list containing two elements:
#   \item{SegmentsTable}{Segments result table: 
#         (segment start index, segment size, segment value)}
#   \item{Segmented}{The complete segmented signal (same size as I).}
# }
#
# @examples "../incl/haarSeg.Rex"
#
# \author{Erez Ben-Yaacov}
#
# @keyword iteration
# @keyword logic
#*/###########################################################################
haarSeg <- function(I, 
        W = vector(),
        rawI = vector(), 
        chromPos = matrix(c(1,length(I)), nrow=1, ncol=2),
        breaksFdrQ = 0.001,        
        haarStartLevel = 1,
        haarEndLevel = 5) 
{
  ProbeNum = length(I);
  weightsFlag = length(W);
  nsvFlag = length(rawI);
  
  if (nsvFlag) {
    # non stationary variance empirical threshold set to 50
    NSV_TH = 50;
    varMask = (rawI < NSV_TH);
  }
  
  S = I;
  allSt = vector();
  allSize = vector();
  allVal = vector();
  CFun = .C("rConvAndPeak", 
            as.double(I), 
            as.integer(ProbeNum), 
            as.integer(1), 
            convResult = double(ProbeNum), 
            peakLoc = integer(ProbeNum),
            PACKAGE="HaarSeg");
  diffI = CFun$convResult;

  if (nsvFlag) {
    pulseSize = 2;
    CFun = .C("rPulseConv",
              as.double(varMask),
              as.integer(ProbeNum),
              as.integer(pulseSize),
              as.double(1/pulseSize),
              res = double(ProbeNum),
              PACKAGE="HaarSeg");
    diffMask = (CFun$res >= 0.5);
  
    peakSigmaEst = median(abs(diffI[!diffMask])) / 0.6745;
    noisySigmaEst = median(abs(diffI[diffMask])) / 0.6745;
    
    if (is.na(peakSigmaEst)) {  peakSigmaEst = 0; }
    if (is.na(noisySigmaEst)) {  noisySigmaEst = 0; }  
  } else {
    peakSigmaEst = median(abs(diffI)) / 0.6745;
  }
  
  # segmentation is done on each chromosome seperatly
  for (chr in 1:nrow(chromPos)) {
    y = I[chromPos[chr,1]:chromPos[chr,2]];
    if (nsvFlag) {
      yVarMask = varMask[chromPos[chr,1]:chromPos[chr,2]];
    }

    if (weightsFlag) {
      wei = W[chromPos[chr,1]:chromPos[chr,2]];
    }

    uniPeakLoc = as.integer(-1);
    for (level in haarStartLevel:haarEndLevel) {
      stepHalfSize = 2^(level);
      if (weightsFlag) {
        CFun = .C("rWConvAndPeak",
                  as.double(y),
                  as.double(wei),
                  as.integer(length(y)),
                  as.integer(stepHalfSize),
                  convResult = double(length(y)),
                  peakLoc = integer(length(y)),
                  PACKAGE="HaarSeg");
      } else {
        CFun = .C("rConvAndPeak",
                  as.double(y),
                  as.integer(length(y)),
                  as.integer(stepHalfSize),
                  convResult = double(length(y)),
                  peakLoc = integer(length(y)),
                  PACKAGE="HaarSeg");
      }

      convRes = CFun$convResult;
      peakLocForC = CFun$peakLoc;
      peakLoc = peakLocForC[1:match(-1,peakLocForC)-1]+1;
  
      if (nsvFlag) {
        pulseSize = 2*stepHalfSize;
        CFun = .C("rPulseConv",
                   as.double(yVarMask),
                   as.integer(length(yVarMask)),
                   as.integer(pulseSize),
                   as.double(1/pulseSize),
                   res = double(length(yVarMask)),
                   PACKAGE="HaarSeg");
        convMask = as.double(CFun$res >= 0.5);
        sigmaEst = (1-convMask)*peakSigmaEst + convMask*noisySigmaEst;
        T = FDRThres(convRes[peakLoc] / sigmaEst[peakLoc], breaksFdrQ, 1);
      } else {
        T = FDRThres(convRes[peakLoc], breaksFdrQ, peakSigmaEst);
      }
  
      unifyWin = as.integer(2^(level - 1));
      tmpPeakLoc = uniPeakLoc;
  
      if (nsvFlag) {
        convRes = convRes / sigmaEst;
      }

      CThres <- .C("rThresAndUnify", 
                    as.double(convRes), 
                    as.integer(length(y)), 
                    peakLocForC,
                    tmpPeakLoc,
                    as.double(T),
                    as.integer(unifyWin),
                    uniPeakLoc = integer(length(y)),
                    PACKAGE="HaarSeg");
      uniPeakLoc = CThres$uniPeakLoc;
    } # for (level ...)

    breakpoints = uniPeakLoc[1:match(-1,uniPeakLoc)-1] + 1;
    
    if (weightsFlag) {
      segs = SegmentByPeaks(y, breakpoints, wei);
    } else {
      segs = SegmentByPeaks(y, breakpoints);
    }
      
    dsegs = which(diff(segs) != 0);
    segSt = c(1,dsegs + 1);
    segEd = c(dsegs,length(segs));
    segSize = segEd - segSt + 1;
    allSt = c(allSt,(segSt + chromPos[chr,1] - 1));
    allSize = c(allSize,segSize);
    allVal = c(allVal,segs[segSt]);
    S[chromPos[chr,1]:chromPos[chr,2]] = segs;
  } # for (chr ...)
  
  segTable = matrix(c(allSt,allSize,allVal), nrow=length(allSt), ncol=3);
  return(list(SegmentsTable = segTable, Segmented = S));
  } # haarSeg()


##############################################################################
# HISTORY:
# 2009-01-06 [EBY]
# o Updated Rdoc comments.
# 2008-12-17 [HB]
# o Rename from HaarSeg() to haarSeg() to follow RCC naming conventions.
# o Created.
############################################################################## 

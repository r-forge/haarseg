#include <R.h>
#include "HaarSeg.h"

void rConvAndPeak(const double * signal,
                  const int * signalSize,
                  const int * stepHalfSize,
                  double * convResult,
                  int * peakLoc) {
    HaarConv(signal, NULL, *signalSize, *stepHalfSize, convResult);
    FindLocalPeaks(convResult, *signalSize, peakLoc);                     
}//rConvAndPeak

void rWConvAndPeak(const double * signal,
                  const double * weight,
                  const int * signalSize,
                  const int * stepHalfSize,
                  double * convResult,
                  int * peakLoc) {
    HaarConv(signal, weight, *signalSize, *stepHalfSize, convResult);
    FindLocalPeaks(convResult, *signalSize, peakLoc);                     
}//rWConvAndPeak


void rThresAndUnify(const double * addon,
                    const int * signalSize,
                    int * addonPeaks,
                    const int * basePeaks,
                    const double * threshold,
                    const int * windowSize,
                    int * uniPeaks) {
    HardThreshold(addon, *threshold, addonPeaks);
    UnifyLevels(basePeaks, addonPeaks, *windowSize, *signalSize, uniPeaks);
}//rThresAndUnify

void rAdjustBreaks(const double * signal,
		   const int * signalSize,
		   const int * peakLoc,
		   int * newPeakLoc) {
    AdjustBreaks(signal, *signalSize, peakLoc, newPeakLoc);
}//rAdjustBreaks

void rPulseConv(const double * signal,
		 	  const int * signalSize, 
		 	  const int * pulseSize,
        const double * pulseHeight, 
		 	  double * result) {
   PulseConv(signal, *signalSize, *pulseSize, *pulseHeight, result);
}//rPulseConv

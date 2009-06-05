/* 
 *    Copyright (C) 2008  Erez Ben-Yaacov
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    <http://www.gnu.org/licenses/>
 *
 */

#include "HaarSeg.h"

#define OK 0
#define ERROR -1
#define NOT_VALID -1
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define abs(X)  ((X) > 0 ? (X) : -(X))
#define MAXINT 2147483647

/*
 * HaarConv : convolve haar wavelet function with a signal, 
 * applying circular padding to the signal.
 * supports weights when weight pointer is not NULL.
 */
 int HaarConv(const double * signal,
              const double * weight,
		 	  int signalSize, 
		 	  int stepHalfSize, 
		 	  double * result)
 {
	int k;
	int highEnd, lowEnd;
	double stepNorm;
    double lowWeightSum;
    double highWeightSum;
    double lowSquareSum;
    double highSquareSum;
    double lowNonNormed;
    double highNonNormed;
    double totalNorm;

	if (stepHalfSize > signalSize) {
		return ERROR; /* TODO: handle this endcase */
	}
    result[0] = 0;
    if (weight != NULL) {
        /* init weight sums */
        highWeightSum = 0;
        highSquareSum = 0;
        highNonNormed = 0;
        for (k = 0; k < stepHalfSize; k++) {
            highWeightSum += weight[k];
            highSquareSum += weight[k]*weight[k]; 
            highNonNormed += weight[k]*signal[k];
        }
	/* circular padding */
        lowWeightSum = highWeightSum; 
        lowSquareSum = highSquareSum;
        lowNonNormed = -highNonNormed;
    }/*if (weight != NULL) */
    for (k = 1; k < signalSize; k++) {
        highEnd = k + stepHalfSize - 1;
        if (highEnd >= signalSize) {
            highEnd = signalSize - 1 - (highEnd - signalSize);
        }
        lowEnd = k - stepHalfSize - 1;
        if (lowEnd < 0) {
            lowEnd = - lowEnd - 1; 
        }
        if (weight != NULL) {
            lowNonNormed += signal[lowEnd]*weight[lowEnd] - signal[k-1]*weight[k-1];
            highNonNormed += signal[highEnd]*weight[highEnd] - signal[k-1]*weight[k-1]; 
            lowWeightSum += weight[k-1] - weight[lowEnd];
            highWeightSum += weight[highEnd] - weight[k-1];
            lowSquareSum += weight[k-1]*weight[k-1] - weight[lowEnd]*weight[lowEnd];
            highSquareSum += weight[highEnd]*weight[highEnd] - weight[k-1]*weight[k-1];
            result[k] = (lowNonNormed / lowWeightSum + highNonNormed / highWeightSum) * sqrt(stepHalfSize/2);            
/*            totalNorm = lowSquareSum / (lowWeightSum*lowWeightSum) + highSquareSum / (highWeightSum*highWeightSum); */
/*            result[k] = (lowNonNormed / lowWeightSum + highNonNormed / highWeightSum) / sqrt(totalNorm); */
        }/*if (weight != NULL) */
        else {
            result[k] = result[k-1] + signal[highEnd] + signal[lowEnd] - 2*signal[k-1];
        }
    }/* for k */
    
    if (weight == NULL) {
        stepNorm = sqrt((double)(2*stepHalfSize));
        for (k = 1; k < signalSize; k++) {
            result[k] /= stepNorm;
        }
    }

	return OK;
 }/* int HaarConv */
 
/*
 * FindLocalPeaks: find local maxima on positive values,
 * and local minima on negative values.
 * First and last index are never considered extramum.
 */
 int FindLocalPeaks(const double * signal, int signalSize, int * peakLoc)
 {
	 int k,j;
	 int maxSuspect, minSuspect;
	 int peakLocInd;
	 
	 maxSuspect = NOT_VALID;
	 minSuspect = NOT_VALID;
	 peakLocInd = 0;
	 for (k = 1; k < signalSize-1; k++) {
		 if (signal[k] > 0) {
			 if ((signal[k] > signal[k-1]) && (signal[k] > signal[k+1])) {
				 peakLoc[peakLocInd] = k;
				 peakLocInd++;
			 }
			 else if ((signal[k] > signal[k-1]) && (signal[k] == signal[k+1])) {
				 maxSuspect = k;
			 }
			 else if (signal[k] == signal[k-1] && (signal[k] > signal[k+1])) {
				 if (maxSuspect != NOT_VALID) {
                     peakLoc[peakLocInd] = maxSuspect;
                     peakLocInd++;
                     /*
					 for (j = maxSuspect; j <= k; j++) {
						 peakLoc[peakLocInd] = j;
						 peakLocInd++;
					 }/* for j */
					 maxSuspect = NOT_VALID;					 
				 }
			 }
			 else if ((signal[k] == signal[k-1]) && (signal[k] < signal[k+1])) {
				 maxSuspect = NOT_VALID;
			 }
		 }/* if (signal[k] > 0) */
		 else if (signal[k] < 0) {
			 if ((signal[k] < signal[k-1]) && (signal[k] < signal[k+1])) {
				 peakLoc[peakLocInd] = k;
				 peakLocInd++;
			 }
			 else if ((signal[k] < signal[k-1]) && (signal[k] == signal[k+1])) {
				 minSuspect = k;
			 }
			 else if ((signal[k] == signal[k-1]) && (signal[k] < signal[k+1])) {
				 if(minSuspect != NOT_VALID) {
                     peakLoc[peakLocInd] = minSuspect;
                     peakLocInd++;
                     /*
					 for (j = minSuspect; j <= k; j++) {
						 peakLoc[peakLocInd] = j;
						 peakLocInd++;
					 }/* for j */
					 minSuspect = NOT_VALID;					 
				 }
			 }
			 else if ((signal[k] == signal[k-1]) && (signal[k] > signal[k+1])) {
				 minSuspect = NOT_VALID;
			 }
		 }/* else if (signal[k] < 0) */
	 }/* for k */
	 
	 peakLoc[peakLocInd] = NOT_VALID;
	 
	 return OK;
 }/* int FindLocalPeaks */
 
 /*
  * HardThreshold: Apply hard thresholding
  */
 int HardThreshold(const double * signal, 
		 		   double threshold,
		 		   int * peakLoc) 
 {
	 int k,l;
	 
     k = 0;
     l = 0;
     while (peakLoc[k] != NOT_VALID) {
		 if ((signal[peakLoc[k]] >= threshold) || (signal[peakLoc[k]] <= -threshold)) {
			 /* peak is over the threshold */
             peakLoc[l] = peakLoc[k];
             l++;
		 }
         k++;
     }
     peakLoc[l] = NOT_VALID;
     
	 return OK;
 }/* int HardThreshold */
 
 /*
  * UnifyLevels: Unify several decomposition levels
  */
 int UnifyLevels(const int * baseLevel,
		 		 const int * addonLevel,
		 		 int windowSize,
		 		 int signalSize,
		 		 int * joinedLevel) 
 {
	 int baseInd,addonInd,joinedInd;
	 	 
	 baseInd = 0;
	 addonInd = 0;
	 joinedInd = 0;
     
	 /* going over all base */
	 while (baseLevel[baseInd] != NOT_VALID) {
		 while ((addonLevel[addonInd] != NOT_VALID) && 
				(addonLevel[addonInd] <= (baseLevel[baseInd] + windowSize))) {
			 if (addonLevel[addonInd] < (baseLevel[baseInd] - windowSize)) {
				 joinedLevel[joinedInd] = addonLevel[addonInd];
				 joinedInd++;
			 }
			 addonInd++;
		 }/* while ((addonLevel[addonInd] ... */
		 joinedLevel[joinedInd] = baseLevel[baseInd];
		 joinedInd++;
		 baseInd++;
	 }/* while (baseLevel[baseInd] */
	 
	 /* insert remaining indexes in addon to joined */
	 while (addonLevel[addonInd] != NOT_VALID) {
		 joinedLevel[joinedInd] = addonLevel[addonInd];
		 joinedInd++;
		 addonInd++;
	 }
	 joinedLevel[joinedInd] = NOT_VALID;
	 
	 return OK;
 }/* int UnifyLevels */
 
 /*
  * CopyLocVec: copy source index vector to target index vector
  */
 int CopyLocVec(const int * source, int * target) {
     int k;
     k = 0;
     while (source[k] != NOT_VALID) {
         target[k] = source[k];
         k++;
     }
/*      mexPrintf("CopyLocVec: copied %d elements\n",k); */
     target[k] = NOT_VALID;
     
     return OK;
 }/* int CopyLocVec */
 
/*
 *  AdjustBreaks: improving localization of breaks by using a suboptimal,
 *  linear complexity procedure. We try to move each break 1 sample 
 *  left/right, choosing the offset which leads to minimum data error. 	
 */
 int AdjustBreaks(const double * signal,
                  int signalSize,
                  const int * peakLoc,
                  int * newPeakLoc) {
    int k,m,p;
    int n1,n2;
    int bestOffset;
    double s1, s2, ss1, ss2;
    double score, bestScore;
    
    k = 0;
    while (peakLoc[k] != NOT_VALID) {
        newPeakLoc[k] = peakLoc[k];
        k++;
    }
    newPeakLoc[k] = NOT_VALID;
    
    k = 0;
    n1 = 0;
    n2 = 0;
    while (newPeakLoc[k] != NOT_VALID) {
        /* calculating width of segments around the breakpoint */
        if (k == 0) {
            n1 = newPeakLoc[k];
        }
        else {
            n1 = newPeakLoc[k] - newPeakLoc[k-1];
        }
        if (newPeakLoc[k+1] == NOT_VALID) {
            n2 = signalSize - newPeakLoc[k];
        }
        else {
            n2 = newPeakLoc[k+1] - newPeakLoc[k];
        }
        
        /* finding the best offset for current breakpoint, trying only 1 sample offset */
        bestScore = MAXINT;
        bestOffset = 0;
        for (p = -1; p <= 1; p++) {
            /* pointless to try and remove single sample segments */
            if ((n1 == 1) && (p == -1)) {
                continue;
            }
            if ((n2 == 1) && (p == 1)) {
                continue;
            }
                
            s1 = 0;
            for (m = (newPeakLoc[k] - n1); m <= (newPeakLoc[k] + p - 1); m++) {
                s1 += signal[m];
            }
            s1 = s1 / (n1 + p);
            s2 = 0;
            for (m = (newPeakLoc[k] + p); m <= (newPeakLoc[k] + n2 - 1); m++) {
                s2 += signal[m];
            }
            s2 = s2 / (n2 - p);
            
            ss1 = 0;
            for (m = (newPeakLoc[k] - n1); m <= (newPeakLoc[k] + p - 1); m++) {
                ss1 += (signal[m] - s1)*(signal[m] - s1);
            }
            ss2 = 0;
            for (m = (newPeakLoc[k] + p); m <= (newPeakLoc[k] + n2 - 1); m++) {
                ss2 += (signal[m] - s2)*(signal[m] - s2);
            }
            score = ss1 + ss2;
            if (score < bestScore) {
                bestScore = score;
                bestOffset = p;
            }
        }/* for p */
        newPeakLoc[k] += bestOffset; 
        k++;
    }/* while newPeakLoc */
    
    return OK;
 }/* int AdjustBreaks */
 
/*
 * StepConv : convolve a pulse function with a signal, 
 * applying circular padding to the signal.
 */
 int PulseConv(const double * signal,
		 	  int signalSize, 
		 	  int pulseSize,
        double pulseHeight, 
		 	  double * result)
 {
	int k, n, tail, head;

	if (pulseSize > signalSize) {
		return ERROR; /* TODO: handle this endcase */
	}
  /* circular padding init */
  result[0] = 0;
  for (k = 0; k < ((pulseSize + 1)/2); k++) {
    result[0] += signal[k];
  }
  for (k = 0; k < (pulseSize/2); k++) {
    result[0] += signal[k];
  }  
  result[0] *= pulseHeight;  
  n = 1;
  for (k = (pulseSize/2); k < signalSize + (pulseSize/2) - 1; k++) {
        tail = k - pulseSize;
        if (tail < 0) {
            tail = -tail - 1; 
        }
        head = k;
        if (head >= signalSize) {
          head = signalSize - 1 - (head - signalSize);
        }
        result[n] = result[n-1] + ((signal[head] - signal[tail]) * pulseHeight);
        n++;
    }/* for k */
    
	return OK;
 }/* int PulseConv */

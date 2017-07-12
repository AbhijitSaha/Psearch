;$Id: //depot/Release/ENVI50_IDL82/idl/idldir/lib/moment.pro#1 $
;
; Copyright (c) 1994-2012, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       MOMENT
;
; PURPOSE:
;       This function computes the mean, variance, skewness and kurtosis
;       of an N-element vector. IF the vector contains N identical elements, 
;       the mean and variance are computed; the skewness and kurtosis are 
;       not defined and are returned as IEEE NaN (Not a Number). Optionally, 
;       the mean absolute deviation and standard deviation may also be 
;       computed. The returned result is a 4-element vector containing the
;       mean, variance, skewness and kurtosis of the input vector X.
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = Moment(X)
;
; INPUTS:
;       X:      An N-element vector of type integer, float or double.
;
; KEYWORD PARAMETERS:
;       DIMENSION: Set this keyword to a scalar indicating the dimension
;         across which to calculate the moment. If this keyword is not
;         present or is zero, then the moment is computed across all
;         dimensions of the input array.
;         
;         If this keyword is present, then the moment is only calculated
;         only across a single dimension. In this case, the results for
;         the KURTOSIS, MDEV, MEAN, SDEV, SKEWNESS, and VARIANCE keywords
;         will be arrays with one less dimension than the input.
;         The result will be an array where the supplied dimension has been
;         removed, but will have a trailing dimension of 4, with the 4
;         slices corresponding to the mean, variance, skewness, and kurtosis.
;
;         For example, if X has dimensions [L,M,N], and DIMENSION is set to 2,
;         then the Result will have dimensions [L,N,4]. Result[*,*,0] will
;         contain the mean, Result[*,*,1] will contain the variance, etc.
;          
;       DOUBLE: IF set to a non-zero value, computations are done in
;               double precision arithmetic.
;
;       KURTOSIS: Set this keyword to a named variable that will contain
;                 the kurtosis of X.
;
;       MAXMOMENT:
;               Use this keyword to limit the number of moments:
;               Maxmoment = 1  Calculate only the mean.
;               Maxmoment = 2  Calculate the mean and variance.
;               Maxmoment = 3  Calculate the mean, variance, and skewness.
;               Maxmoment = 4  Calculate the mean, variance, skewness,
;                              and kurtosis (the default).
;       Note: If X only contains 1 element, then MAXMOMENT is ignored,
;             and MOMENT returns a Mean equal to X, and NaN for all
;             other output values.
;
;       MDEV: Set this keyword to a named variable that will contain
;             the mean absolute deviation of X.
;
;       MEAN: Set this keyword to a named variable that will contain
;                 the mean of X.
;
;       NAN:    Treat NaN elements as missing data.
;               (Due to problems with IEEE support on some platforms,
;                infinite data will be treated as missing as well. )
;
;       SDEV:   Set this keyword to a named variable that will contain
;               the standard deviation of X.
;
;       SKEWNESS: Set this keyword to a named variable that will contain
;                 the skewness of X.
;
;       VARIANCE: Set this keyword to a named variable that will contain
;                 the variance of X.
;
; EXAMPLE:
;       Define the N-element vector of sample data.
;         x = [65, 63, 67, 64, 68, 62, 70, 66, 68, 67, 69, 71, 66, 65, 70]
;       Compute the mean, variance, skewness and kurtosis.
;         result = moment(x)
;       The result should be the 4-element vector: 
;       [66.7333, 7.06667, -0.0942851, -1.18258]
;  
;
; PROCEDURE:
;       MOMENT computes the first four "moments" about the mean of an N-element
;       vector of sample data. The computational formulas are given in the IDL 
;       Reference Guide. 
;
; REFERENCE:
;       APPLIED STATISTICS (third edition)
;       J. Neter, W. Wasserman, G.A. Whitmore
;       ISBN 0-205-10328-6
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, August 1994
;       Modified:    GGS, RSI, September 1995
;                    Added DOUBLE keyword. 
;                    Added checking for N identical elements. 
;                    Added support for IEEE NaN (Not a Number).
;                    Modified variance formula.
;       Modified:    GGS, RSI, April 1996
;                    Modified keyword checking and use of double precision. 
;                    GSL, RSI, August 1997
;                    Added Maxmoment keyword.
;       Modified:    Wed Jan 28 13:28:07 1998, Scott Lett, RSI Added
;                    NAN keyword.
;       Modified: CT, VIS, Jan 2009: Don't throw error for 1-element input,
;                 instead just return [Value, NaN, NaN, NaN].
;       CT, Dec 2009: (from J. Bailin) Added DIMENSION keyword. Also added
;                     MEAN, SKEWNESS, KURTOSIS, VARIANCE keywords.
;       CT, Dec 2010: Fix NAN keyword to filter out both NaN and Infinity
;             when calculating the count of good values.
;-
FUNCTION Moment, X, DIMENSION=dimensionIn, $
    Double = Double, MEAN=mean, Mdev = Mdev, Sdev = Sdev, $
    SKEWNESS=skew, KURTOSIS=kurt, VARIANCE=var, $
    Maxmoment = MaxmomentIn, NaN = nan
    
  compile_opt idl2, hidden, logical_predicate
  ON_ERROR, 2
  
  ndims = SIZE(X, /N_DIMENSION)
  if (KEYWORD_SET(dimensionIn) && (N_ELEMENTS(dimensionIn) gt 1 || $
    dimensionIn[0] lt 1 || dimensionIn[0] gt ndims)) then begin
    MESSAGE, 'Illegal keyword value for DIMENSION.'
  endif
  
  if (~KEYWORD_SET(dimensionIn) || ndims le 1) then begin
    if (KEYWORD_SET(nan)) then begin
      nX = TOTAL(FINITE(x), /INTEGER)
    endif else begin
      nX = N_ELEMENTS(X)
    endelse
    dimension = 0
    fulldimens = [nX]
    Xdimens = [nX]
  endif else begin
    dimension = dimensionIn[0]
    fulldimens = SIZE(x, /DIMENSIONS)
    Xdimens = fulldimens
    Xdimens[dimension-1] = 1
    nResult = PRODUCT(Xdimens, /INTEGER)
    nX = fulldimens[dimension-1]
    if (KEYWORD_SET(nan)) then begin
      nX = TOTAL(FINITE(x), dimension, /INTEGER)
    endif else begin
      nX = fulldimens[dimension-1]
    endelse
  endelse
  
  Maxmoment = Keyword_Set(MaxmomentIn) ? MaxmomentIn[0] : 4
  
  if (N_ELEMENTS(x) lt 2) then Maxmoment = 1
  
  ;If the DOUBLE keyword is not set then the internal precision and
  ;result are identical to the type of input.
  type = SIZE(x, /TYPE)
  IF N_ELEMENTS(Double) EQ 0 THEN $
    Double = type EQ 5 || type EQ 9
    
  Mean = TOTAL(X, dimension, DOUBLE = Double, NAN=nan) / nX
  
  Var  = !VALUES.F_NAN
  Skew = !VALUES.F_NAN
  Kurt = !VALUES.F_NAN
  Mdev = !VALUES.F_NAN
  Sdev = !VALUES.F_NAN
  
  IF Maxmoment GT 1 THEN BEGIN    ; Calculate higher moments.
    Resid = dimension ? (X - REBIN(REFORM(Mean,Xdimens),fulldimens)) : (X - Mean)
    
    ;   Var = TOTAL(Resid^2, Double = Double) / (nX-1.0);Simple formula
    
    ; Numerically-stable "two-pass" formula, which offers less
    ; round-off error. Page 613, Numerical Recipes in C.
    Var = (TOTAL(Resid^2, dimension, DOUBLE = Double, NAN=nan) - $
      (TOTAL(Resid, dimension, DOUBLE = Double, NAN=nan)^2)/nX)/(nX-1.0)
      
    ;Mean absolute deviation (returned through the Mdev keyword).
    if arg_present(Mdev) then $
      Mdev = TOTAL(ABS(Resid), dimension, DOUBLE = Double, NAN=nan) / nX
      
    ; Standard deviation (returned through the Sdev keyword).
    Sdev = SQRT(Var)
    
    ; When dimension is specified, don't worry about the Sdev test, just accept the NaNs
    if (dimension gt 0 || sdev ne 0) then begin
      if maxmoment gt 2 then $
        Skew = TOTAL(Resid^3, dimension, DOUBLE = Double, NAN=nan) / (nX * Sdev ^ 3)
        
      ; The "-3" term makes the kurtosis value zero for normal distributions.
      ; Positive values of the kurtosis (lepto-kurtic) indicate pointed or
      ; peaked distributions; Negative values (platy-kurtic) indicate flat-
      ; tened or non-peaked distributions.
      if maxmoment gt 3 then $
        Kurt = TOTAL(Resid^4, dimension, DOUBLE = Double, NAN=nan) / (nX * Sdev ^ 4) - 3.0
    endif
  endif
  
  if ~dimension then RETURN, [Mean, Var, Skew, Kurt]
  
  rdims = SIZE(Mean, /DIMENSIONS)
  result = KEYWORD_SET(double) ? DBLARR([rdims, 4], /NOZERO) : $
    FLTARR([rdims, 4], /NOZERO)
    
  ; Jump through some REFORM hoops to avoid making copies of the arrays.
  Mean = REFORM(Mean, nResult, /OVERWRITE)
  if (N_ELEMENTS(Var) gt 1) then Var = REFORM(Var, nResult, /OVERWRITE)
  if (N_ELEMENTS(Skew) gt 1) then Skew = REFORM(Skew, nResult, /OVERWRITE)
  if (N_ELEMENTS(Kurt) gt 1) then Kurt = REFORM(Kurt, nResult, /OVERWRITE)
  
  result[0] = Mean
  result[nResult] = Var
  result[2*nResult] = Skew
  result[3*nResult] = Kurt
  
  Mean = REFORM(Mean, rdims, /OVERWRITE)
  if (N_ELEMENTS(Var) gt 1) then Var = REFORM(Var, rdims, /OVERWRITE)
  if (N_ELEMENTS(Skew) gt 1) then Skew = REFORM(Skew, rdims, /OVERWRITE)
  if (N_ELEMENTS(Kurt) gt 1) then Kurt = REFORM(Kurt, rdims, /OVERWRITE)
  
  return, result
  
end

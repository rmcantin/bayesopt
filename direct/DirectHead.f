C+-----------------------------------------------------------------------+
C| Program       : DIRectHead.f                                          |
C| Last modified : 07-25-2007                                            |
C| Written by    : Ruben Martinez-Cantin                                 |
C| Simplified Call function for DIRECT.                                  |
C+-----------------------------------------------------------------------+

      SUBROUTINE Directhead(fcn, x, n, fmin, l, u, Ierror, maxf, maxT,
     +                      iidata)

      IMPLICIT NONE

C+----------------------------------------------------------------------+
C| DIRECT parameters.                                                   |
C+----------------------------------------------------------------------+
      INTEGER maxf, maxT, alg
      PARAMETER (alg = 1)
      
      Double Precision eps, fglobal, fglper, volper, sigmaper
      PARAMETER (eps = 1.D-4)
      PARAMETER (fglobal = -1.D100)
      PARAMETER (fglper = 0.D0)
      PARAMETER (volper = 0.D0)
      PARAMETER (sigmaper = 0.D0)

C+-----------------------------------------------------------------------+
C| EXTERNAL Variables.                                                   |
C+-----------------------------------------------------------------------+
      EXTERNAL fcn
      INTEGER n, Ierror
      DOUBLE PRECISION x(n),fmin,l(n),u(n)

C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize  
      Parameter (iisize = 3)
      Parameter (idsize = 3)
      Parameter (icsize = 3)    
      INTEGER iidata(iisize)
      Double Precision ddata(idsize)
      Character*40 cdata(icsize)    

      INTEGER logfile

C+----------------------------------------------------------------------+
C|  Define and open the logfile and the resultfile, where we store the  |
C|  results of the run. We save the problem number, the number of       |
C|  function evaluations needed, the time used by DIRECT, the percent   |
C|  error and a flag signaling if we used DIRECT or DIRECT-l in the     |
C|  resultfile.                                                         |
C+----------------------------------------------------------------------+

C+----------------------------------------------------------------------+
C|  It seems that DIRECT modify the value of maxf. Do not assign the    |
C|  value as a constant parameter!                                      |
C+----------------------------------------------------------------------+
C      maxf = 300
C      maxT = 1000
      logfile    = 2
      open(logfile, file='direct.out')

C+----------------------------------------------------------------------+
C| Call the optimization method.                                        |
C+----------------------------------------------------------------------+
        CALL DIRect(fcn, x, n, eps, maxf, maxT,
     +              fmin, l, u, alg, Ierror, logfile, 
     +              fglobal, fglper, volper, sigmaper,
     +              iidata, iisize, ddata, idsize, cdata, icsize)

      close(logfile)

      end

c  irirtam.for
c ---------------------------------------------------------------
c  IRI subroutines READIRTAMCOF, FOUT1, and GAMMA2 that read the
c  IRTAM coefficients and calculate NmF2, hmF2, and B0 for the
c  Real-Time IRI. The process of assimilating ionsonde and
c  Digisonde data into a Real-Time IRI is described in: Galkin et
c  al.,Radio Sci., 47, RS0L07, doi:10.1029/2011RS004952, 2012.
c ---------------------------------------------------------------
C  2016.01 06/17/16 First release
C  2016.02 11/03/17 Added B1 input
C  2016.03 12/20/19 Correct directory and name for IRTAM coeffs
C
c ---------------------------------------------------------------
c

		SUBROUTINE READIRTAMCOF_2020(ISEL,IDATE,IHHMM,MFF,FF)
C -----------------------------------------------------------
C Finds and reads IRTAM coefficients for foF2, hmF2, B0, and
C B1 for date IDATE (yyyymmdd) and time HOURUT (decimal hours,
C Universal Time):
C ISEL parameter   filename
C   0    foF2    IRTAM_foF2_COEFFS_yyyymmdd_hhmm.asc
C   1    hmF2    IRTAM_hmF2_COEFFS_yyyymmdd_hhmm.asc
C   2    B0      IRTAM_B0in_COEFFS_yyyymmdd_hhmm.asc
C   3    B1      IRTAM_B1in_COEFFS_yyyymmdd_hhmm.asc
c The coefficient array is stored in FF(MFF).
C D. Bilitza, Jun 17.2016
C Added B1 input      D. Bilitza, Nov 3, 2017
C -----------------------------------------------------------

c		CHARACTER	FILNAM*28
		integer		iuccir
		CHARACTER*100	FILNAM
		CHARACTER*120	LINE
		CHARACTER*12    INAME
		DIMENSION 	    FF(MFF)
        konsol=6

        IUCCIR=10
		iname='foF2_COEFFS_'
		if(isel.gt.0) iname='hmF2_COEFFS_'
		if(isel.gt.1) iname='B0in_COEFFS_'
		if(isel.gt.2) iname='B1in_COEFFS_'

c read foF2 coefficients

        WRITE(FILNAM,104) iname,idate,ihhmm
ccc104     FORMAT('/home/bilitza/tango_home/IRI_data/iri/',
ccc     &     'iri_dev/IRTAM_dir/IRTAM_',A12,I8,'_',I4.4,'.ASC')
104     FORMAT('IRTAM_',A12,I8,'_',I4.4,'.ASC')
        OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=201,
     &          FORM='FORMATTED')
     	print*,mff,filnam

C skip header with comments
c	    do 1,13 READ(iuccir,1289) LINE
4686	READ(iuccir,1289) LINE
1289	Format(A120)
		if(LINE(1:12).NE."# END_HEADER") goto 4686

c read coefficients
        READ(iuccir,4689) FF
4699    FORMAT(E16.8)
4689    FORMAT(4E16.8)
4690    CLOSE(10)

		goto 300

201     WRITE(6,203) FILNAM
203     FORMAT(1X////,
     &    ' The file ',A100,' is not in your directory.')

300		CONTINUE
		RETURN
		END
C
C
      real function FOUT1_2020(XMODIP,XLATI,XLONGI,UT,TOV,FF0)
c--------------------------------------------------------------
C CALCULATES CRITICAL FREQUENCY FOF2/MHZ USING SUBROUTINE GAMMA2.
C XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
C LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME
C (DEC. HOURS), FF0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
C D.BILITZA,JULY 85.
C Modified to accept IRTAM coefficients FFO and TOV
C D. Bilitza, Jun 17.2016
c--------------------------------------------------------------
      DIMENSION FF0(1064)
      INTEGER QF(9)
      DATA QF/11,11,8,4,1,0,0,0,0/
      fout1_2020=GAMMA2_2020(XMODIP,XLATI,XLONGI,UT,TOV,6,QF,9,76,13,106
     &4,FF0)
      RETURN
      END
C
C
        REAL FUNCTION GAMMA2_2020(SMODIP,SLAT,SLONG,HOUR,TOV,
     &                          IHARM,NQ,K1,M,MM,M3,SFE)
C---------------------------------------------------------------
C CALCULATES GAMMA2=FOF2 OR M3000 USING CCIR NUMERICAL MAP
C COEFFICIENTS SFE(M3) FOR MODIFIED DIP LATITUDE (SMODIP/DEG)
C GEOGRAPHIC LATITUDE (SLAT/DEG) AND LONGITUDE (SLONG/DEG)
C AND UNIVERSIAL TIME (HOUR/DECIMAL HOURS). IHARM IS THE MAXIMUM
C NUMBER OF HARMONICS USED FOR DESCRIBING DIURNAL VARIATION.
C NQ(K1) IS AN INTEGER ARRAY GIVING THE HIGHEST DEGREES IN
C LATITUDE FOR EACH LONGITUDE HARMONIC WHERE K1 GIVES THE NUMBER
C OF LONGITUDE HARMONICS. M IS THE NUMBER OF COEFFICIENTS FOR
C DESCRIBING VARIATIONS WITH SMODIP, SLAT, AND SLONG. MM IS THE
C NUMBER OF COEFFICIENTS FOR THE FOURIER TIME SERIES DESCRIBING
C VARIATIONS WITH UT.
C M=1+NQ(1)+2*[NQ(2)+1]+2*[NQ(3)+1]+... , MM=2*IHARM+1, M3=M*MM
C SHEIKH,4.3.77.
C Modified to accept IRTAM coefficients SFE(M3) and TOV
C D. Bilitza, Jun 17.2016
C---------------------------------------------------------------
      REAL*8 C(12),S(12),COEF(100),SUM
      DIMENSION NQ(K1),XSINX(13),SFE(M3)
      COMMON/CONST/UMR,PI

      XMIN=(HOUR-TOV)*60.+720
      HOU=(15.0*HOUR-180.0)*UMR
      S(1)=SIN(HOU)
      C(1)=COS(HOU)

      DO 250 I=2,IHARM
        C(I)=C(1)*C(I-1)-S(1)*S(I-1)
        S(I)=C(1)*S(I-1)+S(1)*C(I-1)
250     CONTINUE

      MMM=M*MM
      DO 300 I=1,M
        MI=(I-1)*MM
        COEF(I)=SFE(MI+1)
        DO 301 J=1,IHARM
          COEF(I)=COEF(I)+SFE(MI+2*J)*S(J)+SFE(MI+2*J+1)*C(J)
301       CONTINUE
        COEF(I)=COEF(I)+SFE(MMM+I)*XMIN
300     CONTINUE

      SUM=COEF(1)
      SS=SIN(SMODIP*UMR)
      S3=SS
      XSINX(1)=1.0
      INDEX=NQ(1)

      DO 350 J=1,INDEX
        SUM=SUM+COEF(1+J)*SS
        XSINX(J+1)=SS
        SS=SS*S3
350     CONTINUE

      XSINX(NQ(1)+2)=SS
      NP=NQ(1)+1
      SS=COS(SLAT*UMR)
      S3=SS

      DO 400 J=2,K1
        S0=SLONG*(J-1.)*UMR
        S1=COS(S0)
        S2=SIN(S0)
        INDEX=NQ(J)+1
        DO 450 L=1,INDEX
          NP=NP+1
          SUM=SUM+COEF(NP)*XSINX(L)*SS*S1
          NP=NP+1
          SUM=SUM+COEF(NP)*XSINX(L)*SS*S2
450       CONTINUE
        SS=SS*S3
400     CONTINUE

      gamma2_2020=SUM

      RETURN
      END
C
C

c  irirtam-test.for
c ---------------------------------------------------------------
c  Test program for the IRI subroutines READIRTAMCOF, FOUT1,
c  and GAMMA2 that read the IRTAM coefficients and calculate
c  NmF2, hmF2, and B0 for the Real-Time IRI.
c ---------------------------------------------------------------
C  2016.01 06/17/16 First release
C  2017.01 11/03/17 Added B1 input
c ---------------------------------------------------------------
c
		Real ff2(1064),fh2(1064),fb0(1064),fb1(1064)
        COMMON/CONST/UMR,PI

        pi=ATAN(1.0)*4.
        UMR=pi/180.

		idate=20160523
		alati=-11.95
		along=283.13
		ryear=2016.6
c        call igrf_dip(alati,along,ryear,300.0,dec,dip,xmagbr,xmodip)
        xmodip=-0.37

		hourut=0.0
		tov=12
		itovhm=1200

		call READIRTAMCOF_2020(0,idate,itovhm,1064,ff2)
		call READIRTAMCOF_2020(1,idate,itovhm,1064,fh2)
		call READIRTAMCOF_2020(2,idate,itovhm,1064,fb0)
c		call READIRTAMCOF(3,idate,itovhm,1064,fb1)

1122	fof2rt=FOUT1_2020(XMODIP,ALATI,ALONG,HOURUT,TOV,ff2)
		hmf2rt=FOUT1_2020(XMODIP,ALATI,ALONG,HOURUT,TOV,fh2)
		B0rt=FOUT1_2020(XMODIP,ALATI,ALONG,HOURUT,TOV,fb0)
c		B1rt=FOUT1(XMODIP,ALATI,ALONG,HOURUT,TOV,fb1)

		print*,idate,hourut,itovhm,fof2rt,hmf2rt,B0rt
c		print*,idate,hourut,itovhm,fof2rt,hmf2rt,B0rt,B1rt
		hourut=hourut+1.0
		if(hourut.lt.18.0) goto 1122
		stop
		end
c

!
!     Parameters NINFI must be updated after changes (for parallelization)
!
!     NOTE: Reals and logicals should appear at the end.
!
      INTEGER NINFI
      PARAMETER (NINFI = 333)
      INTEGER MULD2H, NRHF,NROHF,NVIR, NFRO,                               &
            NISH,NASH,NSSH,NOCC,NORB,NBAS,                                 &
            NNORB,NNBAS, N2ORB,N2BAS,                                      &
            IISH,IASH,ISSH,IOCC,IORB,IBAS,                                 &
            IIISH,IIASH,IIORB,IIBAS,I2ORB,I2BAS,                           &
            ICMO, NSYM,                                                    &
            NISHT,NASHT,NSSHT,NOCCT,NORBT,NBAST,NCMOT,NRHFT,NVIRT,         &
            N2ISHX,NNASHX,N2ASHX,NNASHY,NNOCCX,N2OCCX,                     &
            NNORBT,NNORBX,N2ORBT,N2ORBX,NNBAST,N2BAST,NNBASX,N2BASX,       &
            NNRHFT,NNRHFX,N2RHFT,N2RHFX,NNVIRT,NNVIRX,N2VIRT,N2VIRX,       &
            NAS1,NAS2,NAS3,NNOCCT,N2OCCT,                                  &
            NAS1T,NAS2T,NAS3T
      COMMON /PCM_INFORB/ MULD2H(8,8), NRHF(8),NROHF(8),NVIR(8),NFRO(8),   &
            NISH(8),NASH(8),NSSH(8),NOCC(8),NORB(8),NBAS(8),               &
            NNORB(8),NNBAS(8), N2ORB(8),N2BAS(8),                          &
            IISH(8),IASH(8),ISSH(8),IOCC(8),IORB(8),IBAS(8),               &
            IIISH(8),IIASH(8),IIORB(8),IIBAS(8),I2ORB(8),I2BAS(8),         &
            ICMO(8), NSYM,                                                 &
            NISHT,NASHT,NSSHT,NOCCT,NORBT,NBAST,NCMOT,NRHFT,NVIRT,         &
            N2ISHX,NNASHX,N2ASHX,NNASHY,NNOCCX,N2OCCX,                     &
            NNORBT,NNORBX,N2ORBT,N2ORBX,NNBAST,N2BAST,NNBASX,N2BASX,       &
            NNRHFT,NNRHFX,N2RHFT,N2RHFX,NNVIRT,NNVIRX,N2VIRT,N2VIRX,       &
            NAS1(8),NAS2(8),NAS3(8),NNOCCT,N2OCCT,                         &
            NAS1T,NAS2T,NAS3T
!     MXSSYM = maximum number of "super symmetries"
      INTEGER MXSSYM
      PARAMETER ( MXSSYM = 100 )
      INTEGER NSSYM, NORBSS, IORBSS, NINFSS, MXDGSS
      COMMON /PCM_INFOSS/ NSSYM, NORBSS(MXSSYM), IORBSS(MXSSYM),           &
            NINFSS(MXSSYM,3), MXDGSS

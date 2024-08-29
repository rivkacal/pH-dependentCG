!***********************************************************************
! Common parameters and variables for DualGo.f                         *
!***********************************************************************

        IMPLICIT NONE

! These are settings
!      	integer startt, N, stepstop, thermtime, ThermInt, neartime, WO,
!     Q WOT
	integer(KIND=8)stepstop
      	integer startt, N, thermtime, ThermInt, neartime, WO,
     Q WOT, DEPROTCHARGEFREQ, NHIS, Naromatic, Nhydrophob,PROTCHARGEFREQ
	real percent1, tau, T, pchange, RsTau,gamma, pH
	character (LEN=60) initval, output, finalpx1,
     Q Conf, symtype


	character(LEN=20) PDB, atom
	parameter(PDB="(A6,I5,A4,TR1,A3,I6,3F8.3)")
	parameter(atom='ATOM  ')

	character(LEN=45)FMTX
        parameter(FMTX="(I5,I4,A4,A3,4F8.3)")

! LD parameters
	real c_e, c_i

! Pi is Pi
       real Pi
       parameter (pi = 3.141592653589793)
       integer Nmax, i, j, k, l,JN,MC,
     Q KEtimes

       integer writecount
	character(LEN = 15) list,list1, list2
	parameter(list="(3F8.3,F6.2)")
	parameter(list1="(I5,4F8.3)")
	parameter(list2="(I5,3F8.3)")
       real  X, Y, Z, ms, Vx, Vy, Vz,
     Q Fx, Fy, Fz

! This limits the size of the simulation.
       parameter (Nmax=20000)
       dimension X(Nmax), Y(Nmax), Z(Nmax), ms(Nmax),
     Q Vx(Nmax), Vy(Nmax), Vz(Nmax),
     Q Fx(Nmax), Fy(Nmax), Fz(Nmax)


       real temp, rand, temprand, xrandom, yrandom, zrandom,
     Q KEaveall, KEdev2


! These store information about the chain
	integer ResSeq,MDT,ChainIndex,GroupIndex,BeadIndex
	character(LEN=4) AtType
	character(LEN=3) ResID
	dimension ResSeq(Nmax)
	dimension GroupIndex(Nmax)
	dimension BeadIndex(Nmax)
	dimension ChainIndex(Nmax)
	dimension AtType(Nmax)
	dimension ResID(Nmax)
	integer CLmax
	parameter (clmax = 14)
	integer ChainLength
	dimension ChainLength(ClMax)


! XT, YT, and ZT are temp arrays to be used in one routine at a time
       real XT, YT, ZT
       dimension XT(Nmax), YT(Nmax), ZT(Nmax)
       integer AN, ANo

! these are dummy variables used in several routines
      real XIJ,YIJ,ZIJ,XKJ,YKJ, ZKJ, XKL,YKL,ZKL,DX,DY,
     + DZ, GX,GY,GZ,CT,CPHI,SPHI,Z1, Z2,FXI,FYI,FZI,
     + FXJ,FYJ,FZJ, FXK,FYK,FZK,FXL,FYL,FZL,DF, CT0, CT1,CT2
! end of dummy variables.

! These varaibles are used for indexing the atom types
	integer nAT, ATn, AT, ATbyType
	character (LEN = 2) ATindex
	parameter (nAT = 1)
	dimension ATindex(nAT), ATn(nAT), AT(Nmax),
     Q ATbyType(nAT, Nmax)
	parameter (ATindex =(/'CA'/))

! End of index varaibles

! Energy Variables
	real E, ET, P1, KE

! End of Energy Variables


! Variables used for contacts LJ
	integer NC, trip, tripi
	real PeakH
	parameter (peakH = 0.5)
	integer maxCon, CO, PairNum,
     Q noPairNum, N1, N2, IC, JC, Conts
	real sigma, epsC, epsC1, epsC2, ConCut, Shift
	parameter (ConCut = 1.2)
	parameter (maxCon =300*Nmax)
	dimension sigma(maxCon), epsC(Maxcon), epsC1(Maxcon), epsC2(Maxcon)
	dimension IC(maxcon), JC(maxcon), trip(3,maxcon)

! these are the non native contact pairs' variables
	integer NNoPairNum, NPairNum, INC,JNC, NNC, NNCmax,NNCt,NCset
	parameter (NNCmax = (Nmax-4)**2/2)
	real NCSigma, NNCsigma
	dimension  INC(NNCmax),JNC(NNCmax), 
     Q NCSigma(NNCmax), NNCsigma(NNCmax),NCset(NNCmax)
! End of non native

! these are the ellipsoid repulsions pairs' variables
	integer IEllipsoid, JEllipsoid, KEllipsoid,
     Q	        ellipsoidRepulsionsNum,ellipsoidRepulsionsMax,
     Q          currentEllipsoidRepulsions, 
     Q          currentEllipsoidRepulsionsNum
	parameter (ellipsoidRepulsionsMax = (Nmax-4)**2/2)
	real ellipsoidSigma, ellipsoidCoeff
	dimension  IEllipsoid(ellipsoidRepulsionsMax),
     Q             JEllipsoid(ellipsoidRepulsionsMax),
     Q             KEllipsoid(ellipsoidRepulsionsMax),
     Q		   ellipsoidSigma(ellipsoidRepulsionsMax),
     Q             ellipsoidCoeff(ellipsoidRepulsionsMax),
     Q             currentEllipsoidRepulsions(ellipsoidRepulsionsMax)
! End of ellipsoid repulsions variables

! these are the electrostatics variables
	integer esAtomsNum, maxEsAtoms, esFirstAtomIndex, esSecondAtomIndex,
     Q          esPairsNum, esMinBeadDistance, esMinNeighbor
	real esCharge, esEnergyCutoff, ionicStrength,ionicRadius,
     Q       solventDensity, esDistanceCutoff, DebyeHuckelPotentials,
     Q       DebyeHuckelForces, esCutoffMax, HisCharge, tempChargeByIndex
        parameter (esCutoffMax = 40000.0) 
	parameter (maxEsAtoms = Nmax)
	dimension esFirstAtomIndex(maxEsAtoms**2)
	dimension esSecondAtomIndex(maxEsAtoms**2)
        dimension esCharge(maxEsAtoms**2)
        dimension HisCharge(maxESAtoms)
        dimension DebyeHuckelPotentials(8000000)
        dimension DebyeHuckelForces(8000000)
        dimension tempChargeByIndex(maxESAtoms)
! End of electrostatics  	

! Variables for bonded pairs
	integer NBmax
	parameter (NBmax=Nmax*2)
	integer Ib1, Ib2, nBA, intB, nPROT
	dimension Ib1(NBmax), Ib2(NBmax)
	real Rb, RbT, bK, RBC
	dimension  Rb(NBmax), RbT(NBmax),
     Q  bK(NBmax), RBC(NBmax)
! End of varaible for bonded pairs

! Variables for bond angles
	integer NTmax
	parameter (NTmax =Nmax*2)
	real ANTT, TK, intT, ANTC
	dimension ANTT(Ntmax), TK(Ntmax), ANTC(Ntmax)
	integer IT, JT, KT, nTA
	dimension IT(Ntmax), JT(Ntmax), KT(Ntmax)
! End of varaibles for Bond Angles

! Variables for Phi angles
	integer npmax
	parameter (npmax =NMAX*3)
	integer nPA, IP, JP, KP, LP, intP
	real APT, PK, AP0, AP1, Z10, Z20, Z11, Z22, Z12, Dums, 
     Q DFLIM, DF1, DF0, DR1, DR2, DR3, DR4, DR5, DR6, DRX, DRY, DRZ
	dimension APT(Npmax), IP(Npmax), JP(Npmax), KP(Npmax),
     Q LP(Npmax), PK(Npmax)

	real GAMC1, GAMC3, GAMS1, GAMS3, DihAng
	dimension GAMC1(Npmax), GAMC3(Npmax), GAMS1(Npmax),
     Q DihAng(Npmax), GAMS3(Npmax)

      real GMUL, TM24,TM06,tenm3, zero,one,two,four,six,twelve,ellipsoidRepulsionsNum ftem,
     Q S, COSNP, SINNP, COSNP3, SINNP3, DC1, DC2, DC3,
     Q DC4, DC5, DC6, EPW
	real cosarray, sinarray
	integer refinephi
	parameter (refinephi=10000)
	dimension cosarray(3*refinephi), sinarray(3*refinephi)

! end of varaibles for Phi Angles.
	
! Variables for Chiral angles
	integer nChiralMax
	parameter (nChiralMax =NMAX*3)
	integer nChirals, Ichiral, Jchiral, Kchiral, Lchiral
        real chiralCoeff, chiralValue
	dimension chiralValue(Npmax), Ichiral(Npmax), Jchiral(Npmax),
     Q             Kchiral(Npmax), Lchiral(Npmax), chiralCoeff(Npmax)
! end of variables for chiral angles

! file name variables
	character(LEN=200) Trajectory, EnergyTot,EnergyTerm,ContactFile,
     Q TemperatureVtime, ThreeBodyFile, TrajDist, allContactsFile, 
     Q ContactRangesFile, ContactRangesTwoBodyFile, 
     Q ContactRangesEnergyFile

! variables for contact ranges option
	integer rangesNumber,rangeContactsNumber,rangeContacts, maxRanges,
     Q          TwoBodyRanges
	real    EtotalRanges
        parameter (maxRanges = 50)
	dimension rangeContactsNumber(maxCon)
	dimension rangeContacts(maxRanges,maxCon)
	dimension TwoBodyRanges(maxRanges)
	dimension EtotalRanges(maxRanges)

	

! conditional flags for additional execution settings 
        character(LEN=25) hasStaticAtoms, confineInBox,
     Q                    useElectrostatics, useDebyeHuckel, useChirals,
     Q                    useEllipsoidRepulsions, esCutoffType,
     Q                    compensateElectrostaticContacts,
     Q			  writeAllContacts, useDHEnergyTable,
     Q			  writeContactsRanges, changeHisCharge,
     Q                changeHydroEps


! addtional parameters, used only if conditional flag is set 

        integer DynamicAtomRange(1000)
	integer DynLength
	integer numDyn
	integer useESCutoff, useDHTable
        real boxMin,boxMax, boxCoeff
	dimension boxMin(3),boxMax(3)
        real deConstant, screeningFactor, saltCoefficient

    ! changed during simulation
       integer maxHisNeigh
       integer HisResID
       dimension HisResID(450)

       COMMON /char/ initval, output, finalpx1, 
     Q conf, AtType, ResID,Trajectory, EnergyTot,EnergyTerm,ContactFile,
     Q TemperatureVtime, ThreeBodyFile, symtype, hasStaticAtoms,
     Q confineInBox, useElectrostatics, useDebyeHuckel, useChirals,
     Q useEllipsoidRepulsions, TrajDist, esCutoffType, 
     Q writeAllContacts, allContactsFile,useDHEnergyTable,
     Q writeContactsRanges,ContactRangesFile, ContactRangesTwoBodyFile,
     Q ContactRangesEnergyFile, changeHisCharge, changeHydroEps

       COMMON /real/ X, Y, Z, ms, Vx,
     Q Vy, Vz, Fx, Fy, Fz,rand, temprand, xrandom, yrandom, zrandom,
     Q KEaveall, KEdev2, tau, rstau,gamma,c_i, c_e, T,
     Q pchange, temp, Rb, RbT, RBC, bK, ANTC, ANTT, APT, 
     Q GAMS1, GAMS3, GAMC1, GAMC3, PK, sigma, shift,
     Q epsC, epsC1, epsC2, NCSigma,NNCsigma, esCharge, DihAng,
     Q xt, yt, zt, cosarray, sinarray,boxMin,boxMax,boxCoeff,
     Q deConstant, screeningFactor, saltCoefficient,
     Q chiralCoeff, chiralValue, ellipsoidSigma, ellipsoidCoeff,
     Q esEnergyCutoff, ionicStrength, ionicRadius, solventDensity,
     Q esDistanceCutoff, DebyeHuckelPotentials, DebyeHuckelForces,
     Q EtotalRanges, tempChargeByIndex, HisCharge, pH
     
     
       COMMON /int/ KEtimes, writecount,WO,WOT,NNC,NNCt,NCset,
     Q startt, N, stepstop, thermtime, ThermInt, neartime,
     Q AN, ANo, ATn, AT, ATbyType, Ib1, Ib2, trip, tripi, DEPROTCHARGEFREQ,
     Q nBA, IT, JT, KT, nTA, TK, intT, intB, intP, IP, JP, KP,PROTCHARGEFREQ,
     Q  LP, nPA, IC, JC, N1, N2, PairNum, NoPairNum, INC,JNC,
     Q NPairNum, NNoPairNum, ResSeq, chainlength, MDT, NC, nPROT,
     Q DynamicAtomRange, esAtomsNum,esFirstAtomIndex,esSecondAtomIndex,
     Q esPairsNum,nChirals, Ichiral, Jchiral, Kchiral, Lchiral,
     Q IEllipsoid, JEllipsoid, KEllipsoid,ellipsoidRepulsionsNum,
     Q currentEllipsoidRepulsions, currentEllipsoidRepulsionsNum,
     Q esMinBeadDistance, BeadIndex,GroupIndex,ChainIndex,
     Q numDyn,DynLength, useESCutoff,useDHTable,rangesNumber,
     Q rangeContactsNumber,rangeContacts,TwoBodyRanges, NHIS,maxHisNeigh,
     Q HisResID

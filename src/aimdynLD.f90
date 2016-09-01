subroutine aimdynLD( ntime,nsmll,nmsex,nmage,nmhiv,nmart, &  !extents
     popt,poph,pop1549,poph1549,popa,yrz, &  !records
     U, L, tstep, &                                    !initial state
     FINC, &                                             !final incidence
     hivtarget,hivtarget2,arttarget,arttarget2, & !HIV/ART target data
     nbmv, nbfv, srv, &                           !births and sex ratios
     migrates, mortrates, &                       !migration and mortality
     ar, cdinit, pp, mh, ma, &      !birth data, hiv SR, AR, pp, h mort, a mort
     tbp, record, DD, &                           !TB parameters, record
     stoch, coarse, &                                        !flag for stochastic
     indat, LL &                                             !data and likelihood for return
     )
  use randomxtra
  implicit none
  ! input parameter
  integer, intent(in) :: ntime,nsmll,nmsex,nmage,nmhiv,nmart
  double precision, intent(inout), dimension(ntime) :: popt,poph,pop1549,poph1549,popa,yrz !records
  double precision, intent(inout), dimension(nmsex,nmage,nmhiv,nmart) :: U, L
  double precision, intent(in) :: tstep
  double precision, intent(out), dimension(nmsex,nmage,nmhiv,nmart) :: FINC !final incidence
  double precision, intent(in), dimension(nsmll) :: hivtarget, hivtarget2, arttarget, arttarget2 !array inputs for ART/HIV
  double precision, intent(in), dimension(nsmll) :: nbmv, nbfv, srv !array inputs for births, hivSR
  double precision, intent(in), dimension(nsmll,nmsex,nmage) :: migrates, mortrates !migr/mort inputs
  double precision, intent(in) :: ar(nmage,nmsex), cdinit(nmsex,nmage,nmhiv)
  double precision, intent(in), dimension(nmsex,nmage,nmhiv,nmart) :: pp, mh, ma
  double precision, intent(in) :: tbp(23)  !TB parameters
  double precision, intent(inout), dimension(ntime,20) :: record !records
  double precision, intent(inout), dimension(nmsex,nmage,nmhiv,nmart) :: DD !for outputting both at end
  integer, intent(in) :: stoch
  ! working variables
  double precision, dimension(nmage) :: noma, nofa, fracf, fracm, fracuv, fraclv
  double precision, dimension(nmsex,nmage,nmhiv,nmart) :: X, dX, dI, mu, hivinu, hivinl, mr, ap, P1, P2, &
       artin,artinu,artinl, fracu, IRR
  integer :: i=1, j=1, k=1, ii=1, jj=1, npy
  integer ::  i1549(35), io15(66)
  double precision :: htarget, htarget2, atarget, atarget2, htgt, atgt, newI, newA !working scalars
  double precision :: tp, p1549, yr, nbm, nbf, SR, foi0  !working scalars
  double precision :: ltbinit(nmage)
  double precision, dimension(nmsex,nmage,nmhiv,nmart) :: DU, vDD, vDU, vnotes, vunnotes,vdeaths, vdreg !TBD
  double precision :: cdrg,vrg  !TB globals
  double precision, dimension(8) :: rho !a vector of HIV IRR by CD4 cats
  double precision, dimension(4) :: hra !a vector of ART HRs by time on ART
  real, dimension(nmsex,nmage,nmhiv,nmart) :: CFRd , CFRu, SSd, SSu !arrays of CFRs, durations
  integer, intent(in) :: coarse
  integer :: flf
  integer, intent(in), dimension(nsmll,4) :: indat !in data
  double precision, intent(inout) :: LL            !log-likelihood value

  ! logging!!
  logical :: quiet
  ! end declarations
  integer*4 today(3), now(3)    !for logging
  quiet = .TRUE.               !logging on/off
  if (.NOT. quiet) then
     call idate(today)   ! today(1)=day, (2)=month, (3)=year
     call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
     open(1,file='aimdynlog.txt',status='replace')
     write(1,*) 'opening log file...'
     write(1, 1000 )  today(2), today(1), today(3), now
1000 format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',i2.2, ':', i2.2, ':', i2.2 )
     write(1,*) 'parameters received'
     ! write(1,*) parms
  end if
  ! end logging start

  ! initialize RNG
  if (stoch .EQ. 1) call rndstart()

  ! safety zeroing
  ! r1
  noma = 0.0
  nofa = 0.0
  fracf = 0.0
  fracm = 0.0
  fracuv = 0.0
  fraclv = 0.0
  ! r2
  X = U + L
  dX = 0.0
  dI = 0.0                        !disaggregated hiv incidence
  mu = 0.0                        !mortality
  hivinu = 0.0                     !for disaggregatign hiv incidence
  hivinl = 0.0                     !for disaggregatign hiv incidence
  mr = 0.0                        !migration rates
  ap = 0.0                        !art aging rates
  ap(:,:,2:nmhiv,2) = 2.0
  ap(:,:,2:nmhiv,3) = 4.0                          ! ART aging rates (out)
  P1 = 0.0
  P2 = 0.0
  artin = 0.0
  artinu = 0.0
  artinl = 0.0
  ! dA = 0.0
  fracu = 0.0
  ! r3
  i=1
  j=0                           !this may actually have been necessary
  k=1
  ii=1
  jj=1
  npy = ceiling(1.0/tstep) !number per year
  ! r4
  i1549 = (/ (I, I=16,50) /)                     !indices relevant to 15-49
  io15 = (/ (I, I=16,nmage) /)                     !indices relevant to 15-
  ! r5
  htarget = 0.0
  htarget2 = 0.0
  atarget = 0.0
  atarget2 = 0.0
  htgt = 0.0
  atgt = 0.0
  newI = 0.0
  newA = 0.0
  ! r6
  tp = 0.0
  p1549 = 0.0
  yr = 0.0
  nbm = 0.0
  nbf = 0.0
  SR = 0.0
  foi0 = tbp(1)! 1.0e-2                    !force of infection
  ltbinit = (/ (exp(-foi0*I), I=0,(nmage-1)) /)
  flf = 1

  if( coarse .eq. 1) flf = npy  !flag for some calculations below
  
  do ii = 1, nmhiv
     do jj = 1, nmart
        fracu(1,:,ii,jj) = ltbinit
        fracu(2,:,ii,jj) = ltbinit
     enddo
  enddo



  U = X * fracu
  L = X * (1.0-fracu)
  if( .not.quiet ) write(1,*) 'X=',sum(X)
  if( .not.quiet ) write(1,*) 'U=',sum(U)
  if( .not.quiet ) write(1,*) 'L=',sum(L)
  
  ! new initialisation
  DU = ceiling(tbp(1) * sum(X) / tbp(2)) * (1-fracu) * X/ sum((1-fracu)*X) !down as DT in R code
  DD = ceiling( tbp(6)*exp(-1.0/tbp(12))* DU / (tbp(6)*exp(-1/tbp(12)) + (1-tbp(6))*exp(-1/tbp(11))) )
  DU = floor(DU-DD)
  vDD = 0.0
  vDU = 0.0
  vnotes = 0.0
  vunnotes = 0.0
  vdeaths = 0.0
  vdreg = 0.0
  ! other globals for change or efficiency
  cdrg = tbp(6)
  vrg = tbp(13)

  ! set up the HIV/TB parameters
  rho = (/ 0.0, 250.0, 575.0, 700.0, 775.0, 850.0, 925.0, 975.0 /) !CD4 midpoints' decrement
  rho = rho * tbp(14) * 1e-2                                       !use parameter
  rho = exp( rho )              !actual IRR
  hra = (/ 0.0, 0.25, 0.5, 1.0 /) !geometric progression over: none, <6mo, 12mo, >12mo 
  hra = tbp(15) ** hra            !the HRs by time on ART
  IRR = 1.0
  do ii = 1,nmhiv
     IRR(:,:,ii,:) = IRR(:,:,ii,:) * rho(ii)
  enddo
  do ii = 1,nmart
     IRR(:,:,:,ii) = IRR(:,:,:,ii) * hra(ii)
  enddo

  ! CFRs
  CFRu(:,:,:,1) = tbp(16)                !HIV+ve no ART, off  
  CFRd(:,:,:,1) = tbp(17)                !HIV+ve no ART, on
  CFRu(:,:,:,2:3) = tbp(18)              !HIV+ve ART < 1y, off  
  CFRd(:,:,:,2:3) = tbp(19)              !HIV+ve ART < 1y, on
  CFRu(:,:,:,4) = tbp(20)                !HIV+ve ART > 1y, off  
  CFRd(:,:,:,4) = tbp(21)                !HIV+ve ART > 1y, on
  CFRu(:,:,1,:) = tbp(9)                 !HIV-ve, off  (overwriting)
  CFRd(:,:,1,:) = tbp(10)                !HIV-ve, on

  ! duration/survivals
  SSu = tbp(11)                 !HIV-ve as default
  SSd = tbp(12)
  SSu(:,:,2:4,:) = tbp(22)                 !HIV+ve
  SSd(:,:,2:4,:) = tbp(23)
  ! calculate
  SSu = exp(-tstep*flf/SSu)
  SSd = exp(-tstep*flf/SSd)
  SSu = 1 - SSu                 !now survival function
  SSd = 1 - SSd
    
  if( .not.quiet ) write(1,*) '---> IRR = ',rho
  if( .not.quiet ) write(1,*) '---> HR = ',hra
  
  if( .not.quiet ) write(1,*) 'DD=',sum(DD)
  if( .not.quiet ) write(1,*) 'DU=',sum(DU)
  
  ! time loop
  do i = 1, ntime

     if ( MOD(i-1,npy) .eq. 0 ) then  !update data once a year
        k = 1                         !reset within-year count
        j = j + 1                     !how many updates done
        yr = floor(yrz(i))

        if( .not.quiet ) write(1,*) 'i,j,k,yr=', i,j,k,yr
        ! data updates
        ! external inputs from getcountryparms
        ! HIV & ART targets
        htarget = hivtarget(j)
        htarget2 = hivtarget2(j)
        atarget = arttarget(j)
        atarget2 = arttarget2(j)

        ! births
        nbm = nbmv(j)
        nbf = nbfv(j)
        ! HIV sex disaggregation
        SR = srv(j)
        ! migration rates
        mr(1,:,1,1) = migrates(j,1,:) !M
        mr(2,:,1,1) = migrates(j,2,:) !F

        ! mortality rates (check fortran syntax when porting)
        do ii = 1,nmhiv
           do jj = 1,nmart
              mu(1,:,ii,jj) = mortrates(j,1,:) !M
              mu(2,:,ii,jj) = mortrates(j,2,:) !F
           enddo
        enddo

        ! next 4 are state dependent
        noma = sum(sum(X(1,:,:,:),dim=3),dim=2)
        nofa = sum(sum(X(2,:,:,:),dim=3),dim=2)
        tp = sum(X)
        p1549 = sum(X(:,i1549,:,:))

        ! u fracs
        fracuv = sum(sum(sum(U(:,:,:,:),dim=1),dim=3),dim=2) !U by age
        fraclv = sum(sum(sum(L(:,:,:,:),dim=1),dim=3),dim=2) !L by age
        fracuv = fracuv / (fracuv+fraclv+1e-10)
        fraclv = 1.0 - fracuv

        ! introduce hivinu/l
        !  recalculate
        fracm = ar(:,1) * noma / (sum(ar(:,1) * noma) + 1e-10)
        fracm = (1-SR) * fracm / sum(fracm)
        fracf = ar(:,2) * nofa / (sum(ar(:,2) * nofa) + 1e-10)
        fracf = SR * fracf / sum(fracf)
        do ii=1, nmhiv
           hivinu(1,:,ii,1) = fracuv * fracm * cdinit(1,:,ii)
           hivinu(2,:,ii,1) = fracuv * fracf * cdinit(2,:,ii)
           hivinl(1,:,ii,1) = fraclv * fracm * cdinit(1,:,ii)
           hivinl(2,:,ii,1) = fraclv * fracf * cdinit(2,:,ii)
        enddo
        hivinu(1,:,1,1) = -fracuv * fracm
        hivinu(2,:,1,1) = -fracuv * fracf
        hivinl(1,:,1,1) = -fraclv * fracm
        hivinl(2,:,1,1) = -fraclv * fracf
     endif                      !end data update

     !  ------- TB infection ---------
     dX = U * (1-exp(-foi0*tstep))
     U = U - dX
     L = L + dX
     ! state updates
     ! ------- background mortality ---------
     U = U * mu                     !
     L = L * mu                     !
     ! ! ------- migration ---------
     U = U + mr * tstep * fracu
     L = L + mr * tstep * (1.0-fracu)
     ! -------- aging ---------
     dX = U * tstep
     U = U * (1-tstep)
     U(:,2:nmage,:,:) = U(:,2:nmage,:,:) + dX(:,1:(nmage-1),:,:)
     dX = L * tstep
     L = L * (1-tstep)
     L(:,2:nmage,:,:) = L(:,2:nmage,:,:) + dX(:,1:(nmage-1),:,:)
     ! ------ new births --------
     U(1,1,1,1) = U(1,1,1,1) + tstep * nbm ! new boys
     U(2,1,1,1) = U(2,1,1,1) + tstep * nbf ! new girls
     !     ## child sector...

     ! ------- new infections ------
     ! interpolate target & relax towards
     htgt = htarget2*k/npy + htarget*(npy-k)/npy ! interpolate
     newI = max(real( htgt*sum(X(:,i1549,:,:)) - sum(X(:,i1549,2:nmhiv,:)) ),0.0) * (1-exp(-20*tstep))
     U = U + newI * hivinu
     L = L + newI * hivinl

     !--------- starting ART.--------
     ! data to target
     atgt = atarget2*k/npy + atarget*(npy-k)/npy  ! interpolate
     if (atgt .gt. 0 ) then
        newA = max(real(atgt*sum(X(:,io15,2:nmhiv,:)) - sum(X(:,io15,:,2:4))),0.0) * (1-exp(-20*tstep))
        ! start pattern
        P1 = X
        P1(:,:,1:2,:) = 0              ! not eligible if HIV- or >500
        P1(:,:,:,2:4) = 0            ! not eligible if on ART
        P1 = P1 / sum(P1)
        P2 = X * mh              ! proportional to HIV mortality
        P2 = P2 / sum(P2)
        artin = -0.5 * (P1 + P2)    ! coming out of 1 = HIV
        artin(:,:,:,2) = -artin(:,:,:,1)  ! going into 2=ART
        do ii=1,nmhiv
           do jj=1,nmart
              artinu(1,:,ii,jj) = fracuv * artin(1,:,ii,jj)
              artinu(2,:,ii,jj) = fracuv * artin(2,:,ii,jj)
              artinl(1,:,ii,jj) = fraclv * artin(1,:,ii,jj)
              artinl(2,:,ii,jj) = fraclv * artin(2,:,ii,jj)
           enddo
        enddo
        ! NB updating this here avoids having negative folk in low CD4
        U = U + newA * artinu
        L = L + newA * artinl
     endif

     ! interpolator within year
     k = k+1                     ! iterate interpolator
     ! ------- progression ----------
     !  HIV
     dX = U * tstep * pp                      ! 
     U = U - dX                       ! out from the old
     U(:,:,3:nmhiv,:) = U(:,:,3:nmhiv,:) + dX(:,:,2:(nmhiv-1),:) !in with the new
     dX = L * tstep * pp                      ! 
     L = L - dX                       ! out from the old
     L(:,:,3:nmhiv,:) = L(:,:,3:nmhiv,:) + dX(:,:,2:(nmhiv-1),:) !in with the new
     ! ART
     dX = U * tstep * ap
     U = U - dX
     U(:,:,:,3:4) = U(:,:,:,3:4) +  dX(:,:,:,2:3)
     dX = L * tstep * ap
     L = L - dX
     L(:,:,:,3:4) = L(:,:,:,3:4) +  dX(:,:,:,2:3)
     ! -------- HIV mortality ----------
     U = U * (1 - tstep * mh)
     L = L * (1 - tstep * mh)
     ! -------- ART mortality ----------
     U = U * (1 - tstep * ma)
     L = L * (1 - tstep * ma)

     ! total
     X = U + L

          
     ! TB update
     if (stoch .EQ. 1) then
        call stochstep()
     else
        call appstep()
     endif

     ! recording
     popt(i) = sum(X)                 ! population total
     poph(i) = sum(X(:,:,2:nmhiv,:))         ! hiv total
     poph1549(i) = sum(X(:,i1549,2:nmhiv,:)) ! hiv in 15-49
     pop1549(i) = sum(X(:,i1549,:,:))   ! 15-49 population
     popa(i) = sum(X(:,:,:,2:4))       ! ART
  enddo                         !end time loop

  DD = DD + DU                  !total for output

  if (stoch .EQ. 1) call rndend()

  if (.NOT. quiet) then
     write(1,*) 'closing log file...'
     close(1)
  end if

  ! --------------- end of main programme --------------
contains
  ! ------------- subroutines and functions...

  ! test subroutine
  subroutine testsr()
    implicit none
    i = i + 1
  end subroutine testsr

  ! deterministic approximate update
  subroutine appstep()
    implicit none
    ! tbp = foi,bet,v,pr,eps,cdr0,dcdr,theta,CFRu,CFRd,TSu,TSd,VR
    double precision :: bet,v,pr,eps,cdr0,dcdr,theta,TSu,TSd,VR
    double precision :: cdr,cdr2,vcdr,VRL,vVR,VR2
    double precision, dimension(nmsex,nmage,nmhiv,nmart) :: TBI,lr,dDU,dDD,notes,unnotes,deaths,dreg
    double precision, dimension(2) :: cdrp,vrp
    double precision :: ffc
    
    ! set
    bet=tbp(2);v=tbp(3);pr=tbp(4);eps=tbp(5);cdr0=tbp(6);dcdr=tbp(7);theta=tbp(8)
    VR=tbp(13)

    !! --- means/variances for random detections ---
    cdrp = momentslogitN( logit(cdr0) + i*tstep*dcdr, sqrt(i*tstep)/theta )
    cdr = cdrp(1)
    vcdr = cdrp(2)
    cdr2 = cdr**2 + vcdr

    VRp = momentslogitN( logit(VR) + i*tstep*0, sqrt(i*tstep)/theta )
    VRL = VRp(1)
    vVR = VRp(2)
    VR2 = VR**2 + vVR

    ! make prate and lr etc
    TBI = ((U + v * L) * pr * foi0 + eps * L) * tstep ! poisson rate
    TBI = TBI * IRR                                   !IRR for incidence
    lr = (U + v * L) * pr * (bet/sum(X)) * tstep     !bit proportional to prev

    !! --- mean increments ---
    dDD = TBI*cdr               !new detects
    dDU = TBI - dDD                  !new non-detects
    DD = DD + dDD
    DU = DU + dDU
    notes = (DD-dDD/2) * SSd
    unnotes = (DU-dDU/2) * SSu
    DD = DD - notes
    DU = DU - unnotes
    deaths = unnotes*CFRu
    deaths = deaths + notes*CFRd
    dreg = deaths*VRL

    ! variance approximation
    ffc = (1 - exp(-i * tstep / 20.0)) !scale-up
    ffc = ffc * sqrt(10.0 / theta)    !scaling with theta
    vnotes = (notes / 4) * ffc
    vDU = (DU /2) * ffc         !these are really sds
    VDD = (DD /2) * ffc
    vdreg = (dreg * 3/4) * ffc

    foi0 = bet * sum(DD+DU) / sum(X) !update the FOI outside!

    ! record
    record(i,1) = i * tstep
    record(i,2) = cdr
    record(i,3) = sum(DD)
    record(i,4) = sum(DU)
    record(i,5) = VR
    record(i,6) = sum(notes)
    record(i,7) = sum(unnotes)
    record(i,8) = sum(deaths)
    record(i,9) = sum(dreg)
    record(i,10) = sum(TBI)
    record(i,11) = foi0
    record(i,12) = sum(X)
    record(i,13) = sum(vDD) ! 
    record(i,14) = sum(vDU) ! 
    record(i,15) = sum(vnotes)  !
    record(i,16) = 1e5*(sum(vDD(:,16:nmage,:,:) + vDU(:,16:nmage,:,:)))/sum(X(:,16:nmage,:,:))  !variance for per capita prev
    record(i,17) = sum(vdeaths) !
    record(i,18) = sum(vdreg) ! 
    record(i,19) = vcdr
    record(i,20) = 1e5*(sum(DD(:,16:nmage,:,:) + DU(:,16:nmage,:,:)))/sum(X(:,16:nmage,:,:))  !prev per 100k
    if (i .EQ. ntime) then
       FINC = TBI*flf !record final incidence
       if( .not.quiet ) write(1,*) 'INC HIV% = ', 1e2*sum(FINC(:,:,2:4,:))/sum(FINC)
       if( .not.quiet ) write(1,*) 'PRV HIV% = ', 1e2*sum(DU(:,:,2:4,:)+DD(:,:,2:4,:))/sum(DU+DD)
    endif
    
    ! likeliehood
    if ( MOD(i-1,npy) .eq. 0 ) then  !update data once a year
       
       ! notifications
       if (indat(j,1) .ge. 0 ) LL = LL - DBLE((indat(j,1) - npy*record(i,6))**2) &
            / (2*record(i,15)**2+1e-2) 
       ! vital registrations
       if (indat(j,2) .ge. 0 ) LL = LL - DBLE((indat(j,2) - npy*record(i,9))**2) &
            / (2*record(i,18)**2+1e-2)
       ! prevalence
       if (indat(j,3) .ge. 0 ) LL = LL - 1e3*DBLE((indat(j,3) - record(i,20))**2) &
            / (1*2*record(i,16)**2 + 2*indat(j,4)**2 + 1e-2) !


       if( .not.quiet .and. (indat(j,1) .ge. 0) ) write(1,*) 'LLN = ', &
            DBLE((indat(j,1) - npy*record(i,6))**2) ! / (2*record(i,15)**2+1e-2)
       if( .not.quiet .and. (indat(j,2) .ge. 0) ) write(1,*) 'LLM = ', &
            DBLE((indat(j,2) - npy*record(i,9))**2) ! / (2*record(i,18)**2+1e-2)
       if( .not.quiet .and. (indat(j,3) .ge. 0) ) write(1,*) 'LLP = ', &
            DBLE((indat(j,3) - record(i,20))**2) ! / (2*record(i,16)**2 + 2*indat(j,4)**2 + 1e-2) !

    endif

  end subroutine appstep

  ! stochastic update
  subroutine stochstep()
    implicit none
    ! tbp = foi,bet,v,pr,eps,cdr0,dcdr,theta,CFRu,CFRd,TSu,TSd,VR
    double precision :: bet,v,pr,eps,cdr0,dcdr,theta,TSu,TSd,VR
    integer, dimension(nmsex,nmage,nmhiv,nmart) :: TBI,dDU,dDD,notes,unnotes,deaths,dreg
    integer, dimension(4) :: dms
    integer :: ili              
    ili = i
    
    if ( (coarse .ne. 1) .or. (MOD(i-1,npy) .eq. 0) ) then  !update data once a year
       ! set
       bet=tbp(2);v=tbp(3);pr=tbp(4);eps=tbp(5);cdr0=tbp(6);dcdr=tbp(7);theta=tbp(8)
       VR=tbp(13)
       dms = (/nmsex,nmage,nmhiv,nmart/) !dimensions

       cdrg = logit(cdrg) + tstep*flf * dcdr + sqrt(tstep*flf) * normrnd() /theta !random walk
       cdrg = exp(cdrg)
       cdrg = cdrg / (1.0 + cdrg)
       vrg =  logit(vrg) + sqrt(tstep*flf) * normrnd() /theta !random walk
       vrg = exp(vrg)
       vrg = vrg / (1.0 + vrg)

       ! attention int/dble etc
       TBI = rporay4( REAL(((U + v * L) * pr * foi0 + eps * L) * tstep * flf * IRR) , dms )
       dDD = rbiray4s(TBI,REAL(cdrg),dms)   ! #new detects
       dDU = TBI - dDD                  ! #new non-detects
       DD = DD + dDD
       DU = DU + dDU
       notes = rbiray4( NINT(DD-dDD/2.0), SSd, dms)
       unnotes = rbiray4( NINT(DU-dDU/2.0), SSu, dms)
       DD = DD - notes
       DU = DU - unnotes
       ! next shouldn't be needed but to avoid any <0
       DD = MAX(DD,0.0)
       DU = MAX(DU,0.0)
       deaths = rbiray4(unnotes,CFRu,dms)
       deaths = deaths + rbiray4(notes,CFRd,dms)
       dreg = rbiray4s(deaths,REAL(vrg),dms)

       ! need to make sure cdr and vr are globals and updated, and not clashing with appstep
       foi0 = bet * sum(DD+DU) / sum(X) !update the FOI outside!

       ! record
       if( coarse .eq. 1) ili = j
       record(ili,1) = (ili-1) * tstep * flf
       record(ili,2) = cdrg
       record(ili,3) = sum(DD)
       record(ili,4) = sum(DU)
       record(ili,5) = VRg
       record(ili,6) = sum(notes)
       record(ili,7) = sum(unnotes)
       record(ili,8) = sum(deaths)
       record(ili,9) = sum(dreg)
       record(ili,10) = sum(TBI)
       record(ili,11) = foi0
       record(ili,12) = sum(X)
       record(ili,13) = 0.0
       record(ili,14) = 0.0
       record(ili,15) = 0.0
       record(ili,16) = 0.0
       record(ili,17) = 0.0
       record(ili,18) = 0.0
       record(ili,19) = 0.0
       record(ili,20) = 1e5*(sum(DD(:,16:nmage,:,:) + DU(:,16:nmage,:,:)))/sum(X(:,16:nmage,:,:))       !
    endif
    if (i .EQ. ntime) FINC = REAL(TBI) * flf !record final incidence
  end subroutine stochstep



  ! ------- functions below -------

  ! --- some functions for Gauss-Hermite evaluation of logistic normal moments
  ! logistic distrution shifted and scaled
  function plogis( ngh, xz, mu, sig )
    implicit none
    integer, intent(in) :: ngh
    double precision :: plogis(ngh)
    double precision, intent(in) :: xz(ngh), mu, sig
    plogis = xz * sig * sqrt(2.0) + mu
    plogis = 1.0 / (1 + exp(-plogis))
  end function plogis

  ! approximately calculate the moments of logistic normal
  function momentslogitN( mu, sig )
    implicit none
    double precision :: momentslogitN(2)
    double precision, intent(in) :: mu, sig
    double precision :: mn, vr
    double precision, parameter, dimension(20) :: ghxz = &
         (/ 5.3874809,  4.6036824,  3.9447640,  3.3478546,  2.7888061,  2.2549740, &
         1.7385377,  1.2340762,  0.7374737,  0.2453407, -0.2453407, -0.7374737, &
         -1.2340762, -1.7385377, -2.2549740, -2.7888061, -3.3478546, -3.9447640, &
         -4.6036824, -5.3874809 /) !points
    double precision, parameter, dimension(20) :: ghwts = &
         (/ 1.257801e-13, 2.482062e-10, 6.127490e-08, 4.402121e-06, 1.288263e-04, &
         1.830103e-03, 1.399784e-02, 6.150637e-02, 1.617393e-01, 2.607931e-01, &
         2.607931e-01, 1.617393e-01, 6.150637e-02, 1.399784e-02, 1.830103e-03, &
         1.288263e-04, 4.402121e-06, 6.127490e-08, 2.482062e-10, 1.257801e-13 /) !wts
    mn = sum( ghwts * plogis(20, ghxz, mu, sig) )
    vr = sum( ghwts * (plogis(20, ghxz, mu, sig) - mn)**2 )
    momentslogitN = (/mn,vr/)
  end function momentslogitN
  ! --- end functions for  evaluation of logistic normal moments

  ! logit function
  function logit(x)
    implicit none
    double precision, intent(in) :: x
    double precision :: logit
    logit = log(x/(1.0-x))
  end function logit

end subroutine aimdynLD

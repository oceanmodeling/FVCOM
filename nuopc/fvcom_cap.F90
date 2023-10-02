!>
!! @mainpage FVCOM NUOPC Cap
!! @author Jianhua Qi (jqi@umassd.edu) basded on Saeed Moghimi (moghimis@gmail.com)
!! @date 05/25/2020 Original documentation
!------------------------------------------------------

#define MULTIPROCESSOR
#define SPHERICAL
#define WET_DRY
#define TWO_D_MODEL          !___disabled
#define WAVE_ROLLER___disabled

module fvcom_cap

  !-----------------------------------------------------------------------------
  ! FVCOM Component.
  !-----------------------------------------------------------------------------

  use mpi
  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices,       &
    model_label_SetClock    => label_SetClock,    &
    model_label_CheckImport => label_CheckImport, &    
    model_label_Advance     => label_Advance,     &
    model_label_Finalize    => label_Finalize

!JQI  USE MESSENGER, ONLY  : MPI_COMM_FVCOM, UPDATER

  use NUOPC_FVCOM, only : NUOPC_FVCOM_init
  use NUOPC_FVCOM, only : NUOPC_FVCOM_run
  use NUOPC_FVCOM, only : NUOPC_FVCOM_Final
  use LIMS       , only : MGL,MT,M,NGL,NT,N,KB,KBM1,MYID,NPROCS    
  use ALL_VARS   , only : NTVE,NBVE,EL,U,V,VX,YC
  use ALL_VARS   , only : WAVESTRX_2D, WAVESTRY_2D,WAVESTRX_3D,WAVESTRY_3D,nv
  use ALL_VARS   , only : U_STOKES_2D, V_STOKES_2D,U_STOKES_3D,V_STOKES_3D
  use ALL_VARS   , only : PAR, INTSTEP_SECONDS
  use ALL_VARS   , only : UUWIND, VVWIND
  USE ALL_VARS   , ONLY : E2N2D
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
  USE ALL_VARS   , ONLY : PA_AIR
# endif
  use MOD_PAR    , only : NBN,BN_MLT,BN_LOC,BNC,NODE_MATCH
  use MOD_PAR    , only : NC,EC,AEXCHANGE
  use MOD_PAR    , only : NMAP,MPI_F,MPI_FVCOM_GROUP,NGID_X,ACOLLECT
  use MOD_PAR    , only : EMAP, EGID_X
  use LIMS       , only : MSRID
  use mod_prec   , only : SP

  use mod_wd     , only : ISWETC,ISWETCT
  use mod_spherical, only : TPI,DEG2RAD,DELTUY
  use mod_driver , only : syear, smonth, sday, shour, sminute, ssecond
  use mod_driver , only : eyear, emonth, eday, ehour, eminute, esecond

  use fvcom_mod, only: NUOPC4WAV,NUOPC4MET

  use fvcom_mod, only: meshdata
  use fvcom_mod, only: create_parallel_esmf_mesh_from_meshdata
  use fvcom_mod, only: extract_parallel_data_from_mesh
  use fvcom_mod, only: eliminate_ghosts
  
  implicit none

  private
  public SetServices
  public fvcom_name

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: unit
    logical           :: assoc    ! is the farrayPtr associated with internal data
    logical           :: connected
    real(ESMF_KIND_R8), dimension(:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToFvcom_num = 0
  type (fld_list_type) :: fldsToFvcom(fldsMax)
  integer :: fldsFrFvcom_num = 0
  type (fld_list_type) :: fldsFrFvcom(fldsMax)

  type(meshdata),save  :: mdataIn, mdataOut
  
  character(len=2048):: info
  integer :: dbrc     ! temporary debug rc value

  !to test field halo update.
!JQI  type (ESMF_RouteHandle), save :: ATM_HaloRouteHandel
  
  !real(ESMF_KIND_R8)      :: WaveCouplingIntervalSec, WindCouplingIntervalSec  !in seconds
  !type(ESMF_TimeInterval) :: WaveCouplingInterval, WindCouplingInterval

!JQI2021  logical, save            :: first_exchange = .true.
!JQI2021  integer, save            :: iunit_log = 10000
    logical :: wave_forcing, meteo_forcing, surge_forcing


!JQI2021  real,parameter :: wave_force_limmit = 0.05

! The FVCOM domain name
  character(len=20) :: fvcom_name

  ! config options
  type(ESMF_MeshLoc) :: meshloc
  logical :: dbug
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  !> NUOPC SetService method is the only public entry point.
  !! SetServices registers all of the user-provided subroutines
  !! in the module with the NUOPC layer.
  !!
  subroutine SetServices(model, rc)
  
    implicit none
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    character(len=*),parameter   :: subname='(fvcom_cap:SetServices)'
    ! Local variables
    integer                      :: num,i
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    !Assume no need to change clock settings
    ! attach specializing method(s)
!    call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
!      specRoutine=SetClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    !comment out for now to avoid over writing NOUPC check import
!    call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
!      specPhaseLabel="RunPhase1", specRoutine=CheckImport, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
      specRoutine=FVCOM_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!!!JQI    call FVCOM_FieldsSetup()
!
!    do num = 1,fldsToFvcom_num
!        print *,  "fldsToFvcom_num  ", fldsToFvcom_num
!        print *,  "fldsToFvcom(num)%shortname  ", fldsToFvcom(num)%shortname
!        print *,  "fldsToFvcom(num)%stdname  ", fldsToFvcom(num)%stdname
!     write(info,*) subname,'fldsToFvcom(num)%stdname  ', fldsToFvcom(num)%stdname
!!     write(info,*) subname,"fldsToFvcom(num)%shortname  ", fldsToFvcom(num)%shortname
!   end do
!
!   do num = 1,fldsFrFvcom_num
!     !print *,  "fldsFrFvcom_num  ", fldsFrFvcom_num
!     !print *,  "fldsFrFvcom(num)%shortname  ", fldsFrFvcom(num)%shortname
!     !print *,  "fldsFrFvcom(num)%stdname  ", fldsFrFvcom(num)%stdname
!     write(info,*) subname,'fldsFrFvcom(num)%stdname  ', fldsFrFvcom(num)%stdname
!   end do
!
    
    write(info,*) subname,' --- fvcom SetServices completed --- '
    !print *,      subname,' --- fvcom SetServices completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
  end subroutine SetServices
  
  !-----------------------------------------------------------------------------
  !> First initialize subroutine called by NUOPC.  The purpose
  !! is to set which version of the Initialize Phase Definition (IPD)
  !! to use.
  !!
  !! For this FVCOM cap, we are using IPDv01.
  !!
  !! @param model an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    
    implicit none
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                :: vm
    type(ESMF_Clock)             :: dclock
    type(ESMF_Time)              :: currTime, stopTime
    integer                      :: esmf_comm,fvcom_comm,ierr
    logical                      :: isPresent, isSet
    character(len=ESMF_MAXSTR)   :: message, cvalue
    character(len=*),parameter   :: subname='(fvcom_cap:InitializeP1)'

    rc = ESMF_SUCCESS

    write(info,*) subname,' --- initialization phase 1 begin --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    ! query attributes
    ! fvcom name
    call NUOPC_CompAttributeGet(model, name='case_name', value=cvalue, &
      isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (isPresent .and. isSet) then
      fvcom_name = trim(cvalue)
    else
      fvcom_name = 'sci'
    end if
    call ESMF_LogWrite(trim(subname)//': fvcom_name = '//trim(fvcom_name), ESMF_LOGMSG_INFO)

    ! mesh location for fields
    call NUOPC_CompAttributeGet(model, name='meshloc', value=cvalue, & 
      isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (isPresent .and. isSet) then
       if (trim(cvalue) == 'node') then
          meshloc = ESMF_MESHLOC_NODE
       else
          meshloc = ESMF_MESHLOC_ELEMENT
       end if
    else
       cvalue = 'node'
       meshloc = ESMF_MESHLOC_NODE
    end if
    write(message, '(A)') trim(subname)//' meshloc is set to '//trim(cvalue)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

    ! debug option
    call NUOPC_CompAttributeGet(model, name='dbug', value=cvalue, &
      isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dbug = .false.
    if (isPresent .and. isSet) then
       if (trim(cvalue) .eq. '.true.' .or. trim(cvalue) .eq. 'true') dbug = .true.
    end if
    write(message, '(A,L)') trim(subname)//' debug flag is set to ', dbug
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

    ! query the driver for its clock
    call NUOPC_ModelGet(model, driverClock=dclock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! query start and stop time 
    call ESMF_ClockGet(dclock, currTime=currTime, stopTime=stoptime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(currTime, yy=syear, mm=smonth, dd=sday, h=shour, &
      m=sminute, s=ssecond, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(stopTime, yy=eyear, mm=emonth, dd=eday, h=ehour, &
      m=eminute, s=esecond, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! details Get current ESMF VM.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! details Get MPI_communicator from ESMF VM.
    call ESMF_VMGet(vm, mpiCommunicator=esmf_comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call MPI_Comm_dup(esmf_comm, fvcom_comm, ierr)
    !fvcom_comm = esmf_comm
    !print*, 'MPI_3 ',esmf_comm,fvcom_comm
    ! Duplicate the MPI communicator not to interfere with ESMF communications.
    ! The duplicate MPI communicator can be used in any MPI call in the user
    ! code. Here the MPI_Barrier() routine is called.
    call MPI_Barrier(fvcom_comm, ierr)
    !Initialize fvcom before setting up fields
    
    !NUOPC4MET = .true.
    !NUOPC4WAV = .true.
    
    call NUOPC_FVCOM_init(fvcom_comm,fvcom_name)

    call FVCOM_FieldsSetup()

    !WTIMINC = 3*3600.0
    !RSTIMINC = adc_cpl_int +  adc_cpl_num / adc_cpl_den
    !print *,   'WTIMINC > ',WTIMINC,'  RSTIMINC > ',RSTIMINC

    call FVCOM_AdvertiseFields(importState, fldsToFvcom_num, fldsToFvcom, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call FVCOM_AdvertiseFields(exportState, fldsFrFvcom_num, fldsFrFvcom, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- initialization phase 1 completed --- '
    !print *,      subname,' --- initialization phase 1 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine InitializeP1


  !> Advertises a set of fields in an ESMF_State object by calling
  !! NUOPC_Advertise in a loop.
  !!
  !! @param state the ESMF_State object in which to advertise the fields
  !! @param nfields number of fields
  !! @param field_defs an array of fld_list_type listing the fields to advertise
  !! @param rc return code
  subroutine FVCOM_AdvertiseFields(state, nfields, field_defs, rc)
  
    implicit none
    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(fvcom_cap:FVCOM_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields
      !print *, 'Advertise: '//trim(field_defs(i)%stdname)//'---'//trim(field_defs(i)%shortname)
      call ESMF_LogWrite('Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    !print *,      subname,' --- IN   --- '
    write(info,*) subname,' --- completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine FVCOM_AdvertiseFields
  !
  !----------------------------------------------------------------------------------
  subroutine FVCOM_FieldsSetup
  
    implicit none
    integer                    :: rc,k
    character*(ESMF_MAXSTR)    :: tmpName, tmpShortName
    character(len=*),parameter :: subname='(fvcom_cap:FVCOM_FieldsSetup)'
    logical :: hasEntry

    !--------- import fields to Sea FVCOM -------------
    !TODO: Consider moving these lines to driver to avoid doing it in both CAPS

    write(tmpName,'(a)') 'significant_wave_height'
    write(tmpShortName,'(a)') 'wavhss'

    hasEntry = NUOPC_FieldDictionaryHasEntry(trim(tmpName),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(.not. hasEntry) call NUOPC_FieldDictionaryAddEntry(trim(tmpName), "mx", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=trim(tmpName), shortname=trim(tmpShortName))

    write(tmpName,'(a)') 'average_wave_length'
    write(tmpShortName,'(a)') 'wavlen'

    hasEntry = NUOPC_FieldDictionaryHasEntry(trim(tmpName),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(.not. hasEntry) call NUOPC_FieldDictionaryAddEntry(trim(tmpName), "mx", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=trim(tmpName), shortname=trim(tmpShortName))

    write(tmpName,'(a)') 'average_wave_direction'
    write(tmpShortName,'(a)') 'wavdir'

    hasEntry = NUOPC_FieldDictionaryHasEntry(trim(tmpName),rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(.not. hasEntry) call NUOPC_FieldDictionaryAddEntry(trim(tmpName), "mx", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=trim(tmpName), shortname=trim(tmpShortName))

    !--------- import fields from atm to Fvcom -------------
    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=   &
    "air_pressure_at_sea_level", shortname= "pmsl" )
    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom, stdname=   &
    "inst_zonal_wind_height10m", shortname= "izwh10m" )
    call fld_list_add(num=fldsToFvcom_num, fldlist=fldsToFvcom,            &
    stdname="inst_merid_wind_height10m" , shortname= "imwh10m" )

    !--------- export fields from Sea Fvcom -------------
    call fld_list_add(num=fldsFrFvcom_num, fldlist=fldsFrFvcom, &
         stdname="ocean_mask", shortname= "omask" )
    call fld_list_add(num=fldsFrFvcom_num, fldlist=fldsFrFvcom, &
         stdname="sea_surface_height_above_sea_level", shortname= "seahgt" )
    call fld_list_add(num=fldsFrFvcom_num, fldlist=fldsFrFvcom, &
         stdname="surface_eastward_sea_water_velocity", shortname= "uucurr" )
    call fld_list_add(num=fldsFrFvcom_num, fldlist=fldsFrFvcom, &
         stdname="surface_northward_sea_water_velocity",shortname= "vvcurr" )

!NEMS hycod standard names>>>
!https://esgf.esrl.noaa.gov/projects/couplednems/coupling_fields
!ocn_current_zonal
!ocncz
!m s-1	Ocean current X component.	 	 	 	
!ocn_current_merid
!ocncm
!m s-1	Ocean current Y component.	 	

   !
    write(info,*) subname,' --- Passed--- '
    !print *,      subname,' --- Passed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine FVCOM_FieldsSetup

  !---------------------------------------------------------------------------------

  subroutine fld_list_add(num, fldlist, stdname, data, shortname, unit)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    implicit none
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    real(ESMF_KIND_R8), dimension(:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname
    character(len=*),    intent(in),optional :: unit

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(fvcom_cap:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc)
      return
    endif

    fldlist(num)%stdname = trim(stdname)
    if (present(shortname)) then
      fldlist(num)%shortname = trim(shortname)
    else
      fldlist(num)%shortname = trim(stdname)
    endif

    if (present(data)) then
      fldlist(num)%assoc     = .true.
      fldlist(num)%farrayPtr => data
    else
      fldlist(num)%assoc     = .false.
    endif

    if (present(unit)) then
      fldlist(num)%unit      = unit
    endif


    write(info,*) subname,' --- Passed--- '
    !print *,      subname,' --- Passed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine fld_list_add

  !-----------------------------------------------------------------------------
  !> Called by NUOPC to realize import and export fields.

  !! The fields to import and export are stored in the fldsToFvcom and fldsFrFvcom
  !! arrays, respectively.  Each field entry includes the standard name,
  !! information about whether the field's grid will be provided by the cap,
  !! and optionally a pointer to the field's data array.  Currently, all fields
  !! are defined on the same mesh defined by the cap.
  !! The fields are created by calling fvcom_cap::fvcom_XXXXXXXXXXXXXXXXXXX.
  !!
  !! @param model an ESMF_GridComp object
  !! @param importState an ESMF_State object for import fields
  !! @param exportState an ESMF_State object for export fields
  !! @param clock an ESMF_Clock object
  !! @param rc return code

  subroutine InitializeP2(model, importState, exportState, clock, rc)
    
    implicit none
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field)        :: field
    !Saeed added
    type(meshdata)               :: mdata
    type(ESMF_Mesh)              :: ModelMesh,meshIn,meshOut
    type(ESMF_VM)                :: vm
    real(ESMF_KIND_R8), pointer  :: dataPtr_mask(:)
    integer                      :: i, localPet, petCount
    character(len=*),parameter   :: subname='(fvcom_cap:RealizeFieldsProvidingGrid)'

    rc = ESMF_SUCCESS

    !> \details Get current ESMF VM.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Get query local pet information for handeling global node information
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    ! call ESMF_VMPrint(vm, rc=rc)

    !! Assign VM to mesh data type.
    mdata%vm = vm
    !print *,localPet,"< LOCAL pet, ADC ..1.............................................. >> "
    ! create a Mesh object for Fields
    call extract_parallel_data_from_mesh(mdata, localPet)
    ! keep only non-ghost elements, required for CMEPS coupling
    call eliminate_ghosts(mdata, localPet, dbug)
!        print *,"ADC ..2.............................................. >> "
    call create_parallel_esmf_mesh_from_meshdata(mdata,ModelMesh )
!        print *,"ADC ..3.............................................. >> "
    !
    
    if (dbug) then  
       call ESMF_MeshWrite(ModelMesh, filename="fvcom_mesh", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    end if

!    print *,"FVCOM >> "
!    print *,"NumNd", mdata%NumNd
!    print *,"NumOwnedNd", mdata%NumOwnedNd
!    print *,"NumEl", mdata%NumEl
!    print *,"NumND_per_El", mdata%NumND_per_El


    meshIn  = ModelMesh ! for now out same as in
    meshOut = meshIn

    mdataIn  = mdata
    mdataOut = mdata

!    print *,"FVCOM ..4.............................................. >> "
    call FVCOM_RealizeFields(importState, meshIn , mdata, fldsToFvcom_num, fldsToFvcom, "FVCOM import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !
    call FVCOM_RealizeFields(exportState, meshOut, mdata, fldsFrFvcom_num, fldsFrFvcom, "FVCOM export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!
    ! initialize export field omask (required for coupling through CMEPS)
    ! assuming that mask is not changing in time
    ! we could move this to Advance call and try to use ISWETN (at node points),
    ! ISWETC (at element cells)
    do i = 1, fldsFrFvcom_num
       if (trim(fldsFrFvcom(i)%shortname) == 'omask' .and. fldsFrFvcom(i)%connected) then
          call State_getFldPtr_(ST=exportState, fldname='omask', fldptr=dataPtr_mask, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          dataPtr_mask = 0.0
       end if
    end do

    write(info,*) subname,' --- initialization phase 2 completed --- '
    print *,      subname,' --- initialization phase 2 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
  end subroutine InitializeP2
  

  !> Adds a set of fields to an ESMF_State object.  Each field is wrapped
  !! in an ESMF_Field object.  Memory is either allocated by ESMF or
  !! an existing FVCOM pointer is referenced.
  !!
  !! @param state the ESMF_State object to add fields to
  !! @param grid the ESMF_Grid object on which to define the fields
  !! @param nfields number of fields
  !! @param field_defs array of fld_list_type indicating the fields to add
  !! @param tag used to output to the log
  !! @param rc return code
  subroutine FVCOM_RealizeFields(state, mesh, mdata, nfields, field_defs, tag, rc)

    implicit none
    type(ESMF_State), intent(inout)    :: state
    type(ESMF_Mesh), intent(in)        :: mesh
    type(meshdata)                     :: mdata
    integer, intent(in)                :: nfields
    type(fld_list_type), intent(inout) :: field_defs(:)
    character(len=*), intent(in)       :: tag
    integer, intent(inout)             :: rc


    type(ESMF_Field)                   :: field
    type(ESMF_DistGrid)                :: nodeDistgrid
    type(ESMF_RouteHandle)             :: haloHandle
    integer                            :: i
    character(len=*),parameter         :: subname='(fvcom_cap:FVCOM_RealizeFields)'

    rc = ESMF_SUCCESS
    
    ! Get node DistGrid from the Mesh.
    call ESMF_MeshGet(mesh, nodalDistgrid=nodeDistgrid, rc=rc)
    
    ! Create an ESMF Array with a halo region from a node DistGrid.

    do i = 1, nfields
      field = ESMF_FieldCreate(name=field_defs(i)%shortname, mesh=mesh, &
         typekind=ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)

      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
  
      ! Create the RouteHandle for the halo communication.
      !call ESMF_FieldHaloStore(field,routehandle=haloHandle, rc=rc)
      !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !  line=__LINE__, &
      !  file=__FILE__)) &
      !  return  ! bail out

      ! Do the halo communication.
!!JQI      call ESMF_FieldHalo(field, routehandle=haloHandle, rc=rc)  

        if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

          call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
            ESMF_LOGMSG_INFO, &
            line=__LINE__, &
            file=__FILE__, &
            rc=dbrc)

          !print *,      subname,' --- Connected --- '
          field_defs(i)%connected = .true.
        else
          call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
            ESMF_LOGMSG_INFO, &
            line=__LINE__, &
            file=__FILE__, &
            rc=dbrc)
          ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
          !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
          ! remove a not connected Field from State
          call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          !print *,      subname," Field ", field_defs(i)%stdname ,' --- Not-Connected --- '
          field_defs(i)%connected = .false.
        end if
    end do

    ! After its last use the RouteHandle can be released.
!    call ESMF_FieldHaloRelease(haloHandle, rc=rc)
	
    ! The Field can now be destroyed.
!    call ESMF_FieldDestroy(field, rc=rc)
	
    ! The Array can now be destroyed.
!    call ESMF_ArrayDestroy(array, rc=rc)
	
    write(info,*) subname,' --- OUT--- '
    !print *,      subname,' --- OUT --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
  end subroutine FVCOM_RealizeFields
  !-----------------------------------------------------------------------------

!  subroutine SetClock_mine_not_active(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc

!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ADCTimeStep

!    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
!    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! initialize internal clock
    ! - on entry, the component clock is a copy of the parent clock
    ! - the parent clock is on the slow timescale atm timesteps
    ! - reset the component clock to have a timeStep that is for fvcom-wav of the parent
    !   -> timesteps
    
    !call ESMF_TimeIntervalSet(ADCTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
    !TODO: use nint !!?
!    call ESMF_TimeIntervalSet(ADCTimeStep, s= adc_cpl_int , rc=rc) ! 5 minute steps
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call NUOPC_CompSetClock(model, clock, ADCTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
      
!    print *, "ADC Timeinterval1 = "
!    call ESMF_TimeIntervalPrint(ADCTimeStep, options="string", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!        line=__LINE__, &
!        file=__FILE__)) &
!        return  ! bail out    

!  end subroutine
  !-----------------------------------------------------------------------------
  ! From CICE model uses same clock as parent gridComp
!  subroutine SetClock_not_active(model, rc)
!    type(ESMF_GridComp)  :: model
!    integer, intent(out) :: rc
!    
!    ! local variables
!    type(ESMF_Clock)              :: clock
!    type(ESMF_TimeInterval)       :: ADCTimeStep, timestep
!    character(len=*),parameter  :: subname='(fvcom_cap:SetClock)'

!    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
!    call ESMF_GridCompGet(model, clock=clock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    !call ESMF_TimeIntervalSet(ADCTimeStep, s=     adc_cpl_int, sN=adc_cpl_num, sD=adc_cpl_den, rc=rc) ! 5 minute steps
    ! tcraig: dt is the cice thermodynamic timestep in seconds
!   call ESMF_TimeIntervalSet(timestep, s=adc_cpl_int, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call ESMF_ClockSet(clock, timestep=timestep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
!    call ESMF_TimeIntervalSet(ADCTimeStep, s=adc_cpl_int, rc=rc) 
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call NUOPC_CompSetClock(model, clock, ADCTimeStep, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
!  end subroutine


  !-----------------------------------------------------------------------------

  !> Called by NUOPC to advance the FVCOM model a single timestep >>>>>
  !! <<<<<<<<<  TODO: check! this is not what we want!!!.
  !!
  !! This subroutine copies field data out of the cap import state and into the
  !! model internal arrays.  Then it calls FVCOM_Run to make a NN timesteps.
  !! Finally, it copies the updated arrays into the cap export state.
  !!
  !! @param model an ESMF_GridComp object
  !! @param rc return code
  !-----------------------------------------------------------------------------
  subroutine ModelAdvance(model, rc)

    implicit none

    ! input variables
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_State)           :: importState, exportState
    type(ESMF_Time)            :: currTime
    type(ESMF_TimeInterval)    :: timeStep

    real(ESMF_KIND_R8), ALLOCATABLE, TARGET :: tmp(:)

    ! exports
    real(ESMF_KIND_R8), pointer:: dataPtr_zeta(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_velx(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_vely(:)

    !imports
    real(ESMF_KIND_R8), pointer:: dataPtr_wavhs(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_wavlen(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_wavdir(:)

    real(ESMF_KIND_R8), pointer:: dataPtr_uwnd(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_vwnd(:)
    real(ESMF_KIND_R8), pointer:: dataPtr_pres(:)

    type(ESMF_StateItem_Flag)  :: itemType
    type(ESMF_Mesh)            :: mesh
    type(ESMF_Field)           :: lfield
    character(len=128)         :: fldname,timeStr
    integer                    :: i1,num
    integer                    :: ITIME_BGN_FVCOM, ITIME_END_FVCOM
    integer                    :: nCplFVCOM
    real(ESMF_KIND_R8)         :: timeStepAbs
    integer :: ss,ssN,ssD
    integer :: j, k,ierr
    character*(ESMF_MAXSTR)    :: tmpShortName

    real(ESMF_KIND_R8), allocatable :: unode(:),vnode(:)
    
    real(sp), allocatable :: TMP_WHS(:),TMP_WLEN(:),TMP_WDIR(:)
    real(sp), allocatable :: TMP1_WHS(:),TMP1_WLEN(:),TMP1_WDIR(:)
    real(sp), allocatable :: WHS(:),WLEN(:),WDIR(:)

    real(sp), allocatable :: TMP_WVNX(:),TMP_WVNY(:),TMP_PRN(:)
    real(sp), allocatable :: TMP1_WVNX(:),TMP1_WVNY(:),TMP1_PRN(:)
    real(sp), parameter   :: missing_sp = -999.9_sp
    character(len=*), parameter :: subname='(fvcom_cap:ModelAdvance)'
    ! ----------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS
    dbrc = ESMF_SUCCESS

    ! Query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
       exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

    ! Write current time 
    call ESMF_ClockPrint(clock, options="currTime", &
       preString="------>Advancing FVCOM from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

    ! Query clock
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

    ! Write next time
    call ESMF_TimePrint(currTime + timeStep, &
       preString="------------FVCOM---------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, &
       file=__FILE__)) &
       return  ! bail out

    ! Query time interval
    call ESMF_TimeIntervalGet(timeStep, s=ss,sN=ssN,sD=ssD,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Calculate time to run model 
    timeStepAbs = real(ss)+real(ssN/ssD)
    if (mod(timeStepAbs, IntStep_Seconds) == 0) then
       nCplFVCOM = nint(timeStepAbs/IntStep_Seconds)
    else
       stop
    endif

    ! Query time string that will be used to write states
    call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !-----------------------------------------
    ! Import from WAV
    !-----------------------------------------

    ! Check for import fields
    wave_forcing= .false. !.true.
    !do num = 1,fldsToFvcom_num
    !   if (fldsToFvcom(num)%shortname == 'wavhss') wave_forcing = wave_forcing .and. fldsToFvcom(num)%connected
    !   if (fldsToFvcom(num)%shortname == 'wavlen') wave_forcing = wave_forcing .and. fldsToFvcom(num)%connected
    !   if (fldsToFvcom(num)%shortname == 'wavdir') wave_forcing = wave_forcing .and. fldsToFvcom(num)%connected
    !end do

    ! Get and fill imported wave model related fields
    if (wave_forcing) then
       ! Significant wave height
       call State_getFldPtr(ST=importState,fldname='wavhss',fldptr=dataPtr_wavhs,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

       write(info,'(A,3g14.7,i8)') trim(subname)//' import wavhss ', &
          minval(dataPtr_wavhs), maxval(dataPtr_wavhs), sum(dataPtr_wavhs), size(dataPtr_wavhs)
       call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO)
        
       ! Wave length
       call State_getFldPtr(ST=importState,fldname='wavlen',fldptr=dataPtr_wavlen,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

       write(info,'(A,3g14.7,i8)') trim(subname)//' import wavlen ', &
          minval(dataPtr_wavlen), maxval(dataPtr_wavlen), sum(dataPtr_wavlen), size(dataPtr_wavlen)
       call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO)

       ! Wave direction
       call State_getFldPtr(ST=importState,fldname='wavdir',fldptr=dataPtr_wavdir,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

       write(info,'(A,3g14.7,i8)') trim(subname)//' import wavdir ', &
          minval(dataPtr_wavdir), maxval(dataPtr_wavdir), sum(dataPtr_wavdir), size(dataPtr_wavdir)
       call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO)

       ! Fill internal data arrays
       if (par) then
          ! Allocate temporary arrays
          if (.not. allocated(tmp_whs )) allocate(tmp_whs (0:ngl))
          if (.not. allocated(tmp_wlen)) allocate(tmp_wlen(0:ngl))
          if (.not. allocated(tmp_wdir)) allocate(tmp_wdir(0:ngl))
          tmp_whs (:) = missing_sp
          tmp_wlen(:) = missing_sp
          tmp_wdir(:) = missing_sp

          if (.not. allocated(tmp1_whs )) allocate(tmp1_whs (0:nt))
          if (.not. allocated(tmp1_wlen)) allocate(tmp1_wlen(0:nt))
          if (.not. allocated(tmp1_wdir)) allocate(tmp1_wdir(0:nt))
          tmp1_whs(:)    = missing_sp
          tmp1_whs(1:n)  = dataptr_wavhs (1:n)
          tmp1_wlen(:)   = missing_sp
          tmp1_wlen(1:n) = dataptr_wavlen(1:n)
          tmp1_wdir(:)   = missing_sp
          tmp1_wdir(1:n) = dataptr_wavdir(1:n)

          ! Collect data from processors
          call acollect(myid, msrid, nprocs, emap, tmp1_whs , tmp_whs )
          call acollect(myid, msrid, nprocs, emap, tmp1_wlen, tmp_wdir)
          call acollect(myid, msrid, nprocs, emap, tmp1_wdir, tmp_wlen)

          ! Send the entire array to everyone
          call mpi_bcast(tmp_whs , (ngl+1), mpi_f, 0, mpi_fvcom_group, ierr)
          call mpi_bcast(tmp_wlen, (ngl+1), mpi_f, 0, mpi_fvcom_group, ierr)
          call mpi_bcast(tmp_wdir, (ngl+1), mpi_f, 0, mpi_fvcom_group, ierr)

          ! Mask data and fix range
          ! TODO: Check this is actually needed or not
          where(abs(tmp_whs)  .gt. 1e6) tmp_whs  = 1e-10
          where(abs(tmp_wlen) .gt. 1e6) tmp_wlen = 1e-10
          where(abs(tmp_wdir) .gt. 1e6) tmp_wdir = 1e-10

          where(tmp_whs  >  1e3) tmp_whs  =  1e3
          where(tmp_wlen >  1e3) tmp_wlen =  1e3
          where(tmp_wdir >  1e3) tmp_wdir =  1e3
          where(tmp_whs  < -1e3) tmp_whs  = -1e3
          where(tmp_wlen < -1e3) tmp_wlen = -1e3
          where(tmp_wdir < -1e3) tmp_wdir = -1e3

          ! Fill internal data structures that is used by the model
          do i1 = 1, nt
             whs (i1) = tmp_whs (egid_x(i1))
             wlen(i1) = tmp_wlen(egid_x(i1))
             wdir(i1) = tmp_wdir(egid_x(i1))
          end do
       end if

       ! Compute wave related properties
       call radiation_stress(whs,wlen,wdir)
    else
       nuopc4wav = .false.
       write(info,*) subname,' --- no wave forcing exchange / waves are not all connected --- '
       call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    end if        

    !-----------------------------------------
    ! Import from ATM
    !-----------------------------------------

    ! Check for import fields
    meteo_forcing= .true.
    do num = 1,fldsToFvcom_num
      if (fldsToFvcom(num)%shortname == 'pmsl')    meteo_forcing = meteo_forcing .and. fldsToFvcom(num)%connected
      if (fldsToFvcom(num)%shortname == 'imwh10m') meteo_forcing = meteo_forcing .and. fldsToFvcom(num)%connected
      if (fldsToFvcom(num)%shortname == 'izwh10m') meteo_forcing = meteo_forcing .and. fldsToFvcom(num)%connected
    end do
    
    !Get and fill imported forcing fields
    if (meteo_forcing) then
       ! Mean sea level pressure
       call State_getFldPtr(ST=importState,fldname='pmsl',fldptr=dataPtr_pres,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

       write(info,'(A,3g14.7,i8)') trim(subname)//' import pmsl ', &
          minval(dataPtr_pres), maxval(dataPtr_pres), sum(dataPtr_pres), size(dataPtr_pres)
       call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO)

       ! Meridional wind
       call State_getFldPtr(ST=importState,fldname='imwh10m',fldptr=dataPtr_vwnd,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

       write(info,'(A,3g14.7,i8)') trim(subname)//' import imwh10m ', &
          minval(dataPtr_vwnd), maxval(dataPtr_vwnd), sum(dataPtr_vwnd), size(dataPtr_vwnd)
       call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO)

       ! Zonal wind
       call State_getFldPtr(ST=importState,fldname='izwh10m',fldptr=dataPtr_uwnd,rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
       write(info,'(A,3g14.7,i8)') trim(subname)//' import izwh10m ', &
          minval(dataPtr_uwnd), maxval(dataPtr_uwnd), sum(dataPtr_uwnd), size(dataPtr_uwnd)
       call ESMF_LogWrite(trim(info), ESMF_LOGMSG_INFO)
        
       ! Fill internal data arrays
       if (par) then
          ! Allocate temporary arrays
          if (.not. allocated(tmp_wvnx)) allocate(tmp_wvnx(0:ngl))
          if (.not. allocated(tmp_wvny)) allocate(tmp_wvny(0:ngl))
          if (.not. allocated(tmp_prn )) allocate(tmp_prn (0:ngl))
          tmp_wvnx(:) = missing_sp
          tmp_wvny(:) = missing_sp
          tmp_prn (:) = missing_sp

          if (.not. allocated(tmp1_wvnx)) allocate(tmp1_wvnx(0:nt))
          if (.not. allocated(tmp1_wvny)) allocate(tmp1_wvny(0:nt))
          if (.not. allocated(tmp1_prn )) allocate(tmp1_prn (0:nt))
          tmp1_wvnx(:)   = missing_sp
          tmp1_wvnx(1:n) = dataptr_uwnd(1:n)
          tmp1_wvny(:)   = missing_sp
          tmp1_wvny(1:n) = dataptr_vwnd(1:n)
          tmp1_prn(:)    = missing_sp
          tmp1_prn(1:n)  = dataptr_pres(1:n)

          ! Collect data from processors
          call acollect(myid, msrid, nprocs, emap, tmp1_wvnx, tmp_wvnx)
          call acollect(myid, msrid, nprocs, emap, tmp1_wvny, tmp_wvny)
          call acollect(myid, msrid, nprocs, emap, tmp1_prn , tmp_prn )

          ! Send the entire array to everyone
          call mpi_bcast(tmp_wvnx, (ngl+1), mpi_f, 0, mpi_fvcom_group, ierr)
          call mpi_bcast(tmp_wvny, (ngl+1), mpi_f, 0, mpi_fvcom_group, ierr)
          call mpi_bcast(tmp_prn , (ngl+1), mpi_f, 0, mpi_fvcom_group, ierr)

          ! Fill internal data structures that is used by the model
          do i1 = 1, nt
             uuwind(i1) = tmp_wvnx(egid_x(i1))
             vvwind(i1) = tmp_wvny(egid_x(i1))
          end do
# if defined (AIR_PRESSURE) || (HEATING_CALCULATED)
          ! Fill with the average of surrounding elements
          call e2n2d(tmp_prn, pa_air)
# endif
       end if
    else
       nuopc4met = .false.
       write(info,*) subname,' --- no meteo forcing exchange / atm feilds are not all connected --- '
       call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! dump import state into VTK file
    if (dbug) then
       call StateWriteVTK(importState, 'state_imp_'//trim(timeStr), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    end if

    !-----------------------------------------
    ! Run model
    !-----------------------------------------

    call NUOPC_FVCOM_Run(nCplFVCOM)

    !-----------------------------------------
    ! Export to other components
    !-----------------------------------------

    ! Check for export fields
    surge_forcing= .true.
    do num = 1,fldsFrFvcom_num
       if (fldsFrFvcom(num)%shortname == 'seahgt') surge_forcing = surge_forcing .and. fldsFrFvcom(num)%connected
       if (fldsFrFvcom(num)%shortname == 'uucurr') surge_forcing = surge_forcing .and. fldsFrFvcom(num)%connected
       if (fldsFrFvcom(num)%shortname == 'vvcurr') surge_forcing = surge_forcing .and. fldsFrFvcom(num)%connected
    end do

    if (surge_forcing) then
    
      allocate(unode(mdataOut%NumOwnedNd)); unode = 0.0
      allocate(vnode(mdataOut%NumOwnedNd)); vnode = 0.0
      
      do i1 = 1,mdataOut%NumOwnedNd
        unode(i1) = 0.0
        vnode(i1) = 0.0
	do j = 1, ntve(i1)
	  unode(i1) = unode(i1)+u(nbve(i1,j),1)
	  vnode(i1) = vnode(i1)+v(nbve(i1,j),1)
	end do
	unode(i1) = unode(i1)/ntve(i1)  
	vnode(i1) = vnode(i1)/ntve(i1)  
      end do
      
      !-----------------------------------------
      !   EXPORT
      !-----------------------------------------
      allocate (dataPtr_zeta(mdataOut%NumOwnedNd))

      !pack and send exported fields
      
      allocate (tmp(mdataOut%NumOwnedNd))

      ! >>>>> PACK and send ZETA
      call State_getFldPtr_(ST=exportState,fldname='seahgt',fldptr=dataPtr_zeta, &
!JQI        rc=rc,dump=.true.,timeStr=timeStr)
        rc=rc,dump=.false.,timeStr=timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!      print*,"INSIDE FVCOM CAP ",size(dataPtr_zeta)
      !fill only owned nodes for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd_NoHalo
        tmp(i1) = EL(mdataOut%owned_to_present_nodes(i1))
      end do
      do i1 = 1, mdataOut%NumOwnedNd_Halo
        tmp(i1+mdataOut%NumOwnedNd_NoHalo) = EL(mdataOut%owned_to_present_halo_nodes(i1))
      end do

      !assign to field
      dataPtr_zeta = tmp
      
      deallocate(tmp)
      
      allocate(tmp(mdataOut%NumOwnedNd))
      allocate(dataPtr_velx(mdataOut%NumOwnedNd))
      !----------------------------------------
      ! >>>>> PACK and send VELX
      call State_getFldPtr(ST=exportState,fldname='uucurr',fldptr=dataPtr_velx,rc=rc)
!JQI      call State_getFldPtr_(ST=exportState,fldname='velx',fldptr=dataPtr_velx, &
!JQI        rc=rc,dump=.true.,timeStr=timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill elements for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd_NoHalo
        tmp(i1) = UNODE(mdataOut%owned_to_present_nodes(i1))
      end do
      do i1 = 1, mdataOut%NumOwnedNd_Halo
        tmp(i1+mdataOut%NumOwnedNd_NoHalo) = UNODE(mdataOut%owned_to_present_halo_nodes(i1))
      end do
      !assign to field
      dataPtr_velx = tmp
      deallocate(tmp)
      
      allocate(tmp(mdataOut%NumOwnedNd))
      allocate(dataPtr_vely(mdataOut%NumOwnedNd))
      !----------------------------------------
      ! >>>>> PACK and send VELY
      call State_getFldPtr(ST=exportState,fldname='vvcurr',fldptr=dataPtr_vely,rc=rc)
!JQI      call State_getFldPtr_(ST=exportState,fldname='vely',fldptr=dataPtr_vely, &
!JQI        rc=rc,dump=.true.,timeStr=timeStr)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      !fill elements for tmp vector
      do i1 = 1, mdataOut%NumOwnedNd_NoHalo
        tmp(i1) = VNODE(mdataOut%owned_to_present_nodes(i1))
      end do
      do i1 = 1, mdataOut%NumOwnedNd_Halo
        tmp(i1+mdataOut%NumOwnedNd_NoHalo) = VNODE(mdataOut%owned_to_present_halo_nodes(i1))
      end do
      !assign to field
      dataPtr_vely = tmp
      deallocate(tmp)

      deallocate(unode,vnode)
    else
      write(info,*) subname,' --- no surge forcing for wave. 1way coupled WW3 -> FVCOM  ---'
      call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
      !print *, info
    end if

    ! dump export state
    if (dbug) then
       call StateWriteVTK(exportState, 'state_exp_'//trim(timeStr), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
    end if

!    deallocate(unode)
!    deallocate(vnode)
      
  end subroutine ModelAdvance

!-----------------------------------------------------------
  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  subroutine State_GetFldPtr_(ST, fldname, fldptr, rc, dump,timeStr)
    
    implicit none
    type(ESMF_State), intent(in) :: ST
    type(ESMF_RouteHandle)       :: haloHandle
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:)
    integer, intent(out), optional :: rc
    logical, intent(in), optional  :: dump
    character(len=128),intent(inout), optional :: timeStr
    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(fvcom_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (present(rc)) rc = lrc

    !call ESMF_FieldHaloStore(lfield, routehandle = haloHandle, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !Halo update
    !call ESMF_FieldHalo(lfield, routehandle = haloHandle, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !!
    !!TODO: this should not be here. It should finalize once
!!JQI    call ESMF_FieldHaloRelease (routehandle = haloHandle, rc=lrc)
!!JQI    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 
    !if (present(rc)) rc = lrc
    !write(info,*) ' --- ATM  halo routehandel in work >>>>>  ---'
    !call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    if (present(dump)) then
      if (dump) then
        if (.not. present(timeStr)) timeStr="_"
        call ESMF_FieldWrite(lfield, &
          fileName='field_ocn_'//trim(fldname)//trim(timeStr)//'.nc', &
          rc=rc,overwrite=.true.)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if
  end subroutine State_GetFldPtr_

  !> Retrieve a pointer to a field's data array from inside an ESMF_State object.
  !!
  !! @param ST the ESMF_State object
  !! @param fldname name of the fields
  !! @param fldptr pointer to 1D array
  !! @param rc return code
  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    
    implicit none
    type(ESMF_State), intent(in) :: ST
    type(ESMF_RouteHandle)       :: haloHandle
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(FVCOM:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !call ESMF_FieldHaloStore(lfield, routehandle = haloHandle, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    !Halo update
    !call ESMF_FieldHalo(lfield,routehandle = haloHandle, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    !!
    !!TODO: this should not be here. It should finalize once
    !call ESMF_FieldHaloRelease (routehandle = haloHandle, rc=lrc)
    !if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    if (present(rc)) rc = lrc
  end subroutine State_GetFldPtr

  subroutine CheckImport_not_comp(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc
    
    ! This is the routine that enforces correct time stamps on import Fields
    
    ! local variables
    type(ESMF_Clock)        :: driverClock
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_State)        :: importState
    type(ESMF_Field)        :: field
    logical                 :: atCorrectTime

    rc = ESMF_SUCCESS
    return
    
!    ! query Component for the driverClock
!    call NUOPC_ModelGet(model, driverClock=driverClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
!    ! get the start time and current time out of the clock
!    call ESMF_ClockGet(driverClock, startTime=startTime, &
!      currTime=currTime, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    
   
  end subroutine

  !-----------------------------------------------------------------------------
  


  !-----------------------------------------------------------------------------
  !> Called by NUOPC at the end of the run to clean up.  The cap does
  !! this simply by calling FVCOM_Final.
  !!
  !! @param model the ESMF_GridComp object
  !! @param rc return code
  subroutine FVCOM_model_finalize(model, rc)

    implicit none  
    !input arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    !local variables
    type(ESMF_Clock)     :: clock
    type(ESMF_Time)      :: currTime
    character(len=*),parameter  :: subname='(fvcom_cap:fvcom_model_finalize)'

    rc = ESMF_SUCCESS

    write(info,*) subname,' --- finalize called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    call NUOPC_FVCOM_Final()

    write(info,*) subname,' --- finalize completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine FVCOM_model_finalize

!==============================================================================|
   SUBROUTINE RADIATION_STRESS(WHS,WLEN,WDIR)
!------------------------------------------------------------------------------|
   
   USE MOD_PREC
   USE ALL_VARS, ONLY : ISBCE,DLTXC,DLTYC,IEC,IENODE,NE,RAMP,DEG2RAD,GRAV_N,D1,ZZ,D,DZ,DZ1,H1
   IMPLICIT NONE

   REAL(SP), DIMENSION(:) :: WHS,WLEN,WDIR
   REAL(SP), ALLOCATABLE     :: SXX(:,:),SXY(:,:),SYY(:,:)
   REAL(SP), DIMENSION(0:NT,KB) :: PSXXPX,PSXYPX,PSXYPY,PSYYPY
   REAL(SP), DIMENSION(0:NT,KB) :: PSPXPZ,PSPYPZ

   REAL(SP), DIMENSION(0:MT) :: WAVE_NUMBER,WAVE_NUMBER_X,WAVE_NUMBER_Y,SIN_DIR,COS_DIR
   REAL(SP), DIMENSION(0:MT) :: WAVE_ENERGY,KD,WAVE_C
   REAL(SP), DIMENSION(0:MT) :: O_WAVE_NUMBER
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH,O_2SINH
   REAL(SP) :: EXFLUX

   INTEGER  :: I,K,IA,IB,J1,J2
   REAL(SP) :: FSS,FCS,FSC,FCC
   REAL(SP) :: CFF1,CFF2,CFF3,CFF4,CFF5,CFF6,FAC2,sum3dsxx
   REAL(SP) :: SXXIJ,SXYIJ,SYYIJ,DIJ

   REAL(SP) :: XTMP,XTMP1
   REAL(SP), ALLOCATABLE :: TPZDIST(:,:)


   real(sp), allocatable :: sxx_gl(:),sxy_gl(:),syy_gl(:)
   real(sp), allocatable :: sxx_tmp(:),sxy_tmp(:),syy_tmp(:)
   real(sp), allocatable :: U_STOKES_3D_TMP(:,:), V_STOKES_3D_TMP(:,:)

!==============================================================================|
   REAL(SP), PARAMETER :: KDMAX = 3.0_SP ! Based on MELLOR(2015). For old code: KDMAX = 5.0_SP
   REAL(SP), PARAMETER :: eps1 = 1E-14_SP
   REAL(SP), PARAMETER :: WAVE_LENGTH_MIN = 0.01_SP
   REAL(SP), PARAMETER :: PI = 3.1415926
   
   IF(wave_forcing)THEN
!---------------Jianzhong----------------------------
   IF(.NOT.ALLOCATED(SXX)) ALLOCATE(SXX(0:MT,KB))
   IF(.NOT.ALLOCATED(SXY)) ALLOCATE(SXY(0:MT,KB))
   IF(.NOT.ALLOCATED(SYY)) ALLOCATE(SYY(0:MT,KB))
   IF(.NOT.ALLOCATED(U_STOKES_3D_TMP)) ALLOCATE(U_STOKES_3D_TMP(0:MT,KB))
   IF(.NOT.ALLOCATED(V_STOKES_3D_TMP)) ALLOCATE(V_STOKES_3D_TMP(0:MT,KB))
!----------------------------------------------------
   ALLOCATE(TPZDIST(0:NT,KB));       TPZDIST     = 0.0_SP

   WAVE_NUMBER   = 0.0_SP   ;WAVE_NUMBER_X = 0.0_SP   ;WAVE_NUMBER_Y = 0.0_SP
   O_COSH        = 0.0_SP   ;O_SINH        = 0.0_SP   ;O_2SINH       = 0.0_SP
   O_WAVE_NUMBER = 0.0_SP

!
!  Compute wave numbers and wave energy.
!
   DO I=1,MT
    WAVE_NUMBER(I) = 2.0_SP*PI/MAX(WLEN(I),WAVE_LENGTH_MIN)
   END DO 
   O_WAVE_NUMBER = 1.0_SP/WAVE_NUMBER
!JQI   SIN_DIR       = SIN(WDIR*DEG2RAD)
!JQI   COS_DIR       = COS(WDIR*DEG2RAD)
   SIN_DIR       = SIN(WDIR)
   COS_DIR       = COS(WDIR)
   WAVE_NUMBER_X = WAVE_NUMBER*COS_DIR
   WAVE_NUMBER_Y = WAVE_NUMBER*SIN_DIR
   WAVE_ENERGY   = 0.0625_SP*GRAV_N*WHS*WHS
!
!  Compute wave celerity and phase velocity.
!
   DO I=1,MT
!     KD(I) = MIN(WAVE_NUMBER(I)*D(I)+eps1,KDMAX)
!JQI error     KD(I) = WAVE_NUMBER(I)*D1(I)+eps1
     KD(I) = WAVE_NUMBER(I)*D(I)+eps1
   END DO 
   
   WHERE(KD <= KDMAX) 
    WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))

    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)
    O_2SINH = 1.0_SP/SINH(2.0_SP*KD)
   ELSEWHERE
    WAVE_C = SQRT(GRAV_N*O_WAVE_NUMBER*TANH(KD))
   END WHERE

#if defined(WAVE_ROLLER)
   OROLLER = 0.0_SP;GAMW = 0.0_SP; ROLLA = 0.0_SP
   DO I=1,MT
     GAMW(I) = MIN(D(I)/(HSC1(I)+eps1),5.0_SP) 
     DO K=1,KBM1
        OROLLER(I)=OROLLER(I)+D(I)*DZ(I,K)*(1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4))
     END DO
     OROLLER(I)=1.0_SP/(OROLLER(I)+eps1)
     ROLLA(I)=0.0424*HSC1(I)*QB1(I)*WLEN(I)
   END DO
#endif
   
!----------INITIALIZE STRESS ARRAY ----------------------------------------------!

   SXX    = 0.0_SP   ;SXY    = 0.0_SP   ;SYY    = 0.0_SP
!   SXXA   = 0.0_SP   ;SXYA   = 0.0_SP   ;SYYA   = 0.0_SP
   PSXXPX = 0.0_SP   ;PSXYPX = 0.0_SP   ;PSXYPY = 0.0_SP   ;PSYYPY = 0.0_SP
   PSPXPZ = 0.0_SP   ;PSPYPZ = 0.0_SP

   DO I=1,M
     sum3dsxx=0
    IF(KD(I) <= KDMAX)THEN 
     DO K=1,KBM1
       FAC2 = 1.0_SP+ZZ(I,K)
       FCC  = COSH(KD(I)*FAC2)*O_COSH(I)
       FCS  = COSH(KD(I)*FAC2)*O_SINH(I)
       FSC  = SINH(KD(I)*FAC2)*O_COSH(I)
       FSS  = SINH(KD(I)*FAC2)*O_SINH(I)

       CFF1 = WAVE_NUMBER(I)*WAVE_ENERGY(I)
!       CFF4 = CFF1*FCS*FSS
       CFF4 = CFF1*FSC*FSS
       CFF6 = CFF1*FCS*FCC
       CFF5 = CFF1*FCS*FCC*O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
#if defined(WAVE_ROLLER)
       CFF3 = 1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4)
       CFF3 = CFF3*OROLLER(I)*ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4 + &
                  + CFF6 + CFF3*COS_DIR(I)*COS_DIR(I)
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4 + &
                  + CFF6 + CFF3*SIN_DIR(I)*SIN_DIR(I)
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)      + &
                  CFF3*SIN_DIR(I)*COS_DIR(I)
#else
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)+CFF6-CFF4
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)+CFF6-CFF4
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)
#endif       
       TPZDIST(I,K) = FCC*FSS
     END DO  
    ELSE
     DO K=1,KBM1
       FAC2 = ZZ(I,K)
       FCC  = EXP(KD(I)*FAC2)
       FCS  = FCC
       FSC  = FCC
       FSS  = FCC

       CFF1 = WAVE_NUMBER(I)*WAVE_ENERGY(I)
!       CFF4 = CFF1*FCS*FSS
       CFF4 = CFF1*FSC*FSS
       CFF6 = CFF1*FCS*FCC
       CFF5 = CFF1*FCS*FCC*O_WAVE_NUMBER(I)*O_WAVE_NUMBER(I)
#if defined(WAVE_ROLLER)
       CFF3 = 1.0_SP-TANH((2.0_SP*ZZ(I,K)*GAMW(I))**4)
       CFF3 = CFF3*OROLLER(I)*ROLLA(I)/(WLEN(I)+eps1)*WAVE_C(I)**2
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)-CFF4 + &
                  + CFF6 + CFF3*COS_DIR(I)*COS_DIR(I)
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)-CFF4 + &
                  + CFF6 + CFF3*SIN_DIR(I)*SIN_DIR(I)
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)      + &
                  CFF3*SIN_DIR(I)*COS_DIR(I)
#else
       SXX(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_X(I)+CFF6-CFF4
       SYY(I,K) = CFF5*WAVE_NUMBER_Y(I)*WAVE_NUMBER_Y(I)+CFF6-CFF4
       SXY(I,K) = CFF5*WAVE_NUMBER_X(I)*WAVE_NUMBER_Y(I)
#endif       
       TPZDIST(I,K) = FCC*FSS
     END DO  
    END IF
   END DO  

#if defined (MULTIPROCESSOR)
   IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,KB,MYID,NPROCS,SXX,SXY,SYY)
   IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SXX,SXY,SYY)   !Jianzhong
#endif

!JQI       allocate(sxx_gl(0:MGL)); sxx_gl = 0.0_sp
!JQI       allocate(sxy_gl(0:MGL)); sxy_gl = 0.0_sp
!JQI       allocate(syy_gl(0:MGL)); syy_gl = 0.0_sp
!JQI       allocate(sxx_tmp(0:MT)); sxx_tmp = 0.0_sp
!JQI       allocate(sxy_tmp(0:MT)); sxy_tmp = 0.0_sp
!JQI       allocate(syy_tmp(0:MT)); syy_tmp = 0.0_sp

 !JQI      sxx_tmp(:) = sxx(:,1)
!JQI       sxy_tmp(:) = sxy(:,1)
!JQI       syy_tmp(:) = syy(:,1)

!JQI       CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,sxx_tmp,sxx_gl)
!JQI       CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,sxy_tmp,sxy_gl)
!JQI       CALL ACOLLECT(MYID,MSRID,NPROCS,NMAP,syy_tmp,syy_gl)
!JQI       do i=1,mgl
!JQI         write(700+myid,*) i,sxx_gl(i),sxy_gl(i),syy_gl(i)
!JQI       end do

!JQI       deallocate(sxx_gl)
!JQI       deallocate(sxy_gl)
!JQI       deallocate(syy_gl)
!JQI       deallocate(sxx_tmp)
!JQI       deallocate(sxy_tmp)
!JQI       deallocate(syy_tmp)
!JQI

   DO I=1,NE
     IA=IEC(I,1)
     IB=IEC(I,2)
     J1=IENODE(I,1)
     J2=IENODE(I,2)

#    if defined (WET_DRY)
     IF(ISWETC(IA) == 1 .OR. ISWETC(IB) == 1)THEN
#    endif

     DO K=1,KBM1
       SXXIJ=0.5_SP*(SXX(J1,K)+SXX(J2,K))
       SXYIJ=0.5_SP*(SXY(J1,K)+SXY(J2,K))
       SYYIJ=0.5_SP*(SYY(J1,K)+SYY(J2,K))
       DIJ = 0.5_SP*(D(J1)*DZ(J1,K)+D(J2)*DZ(J2,K))

#      if defined (SPHERICAL)
       !for spherical coordinator and domain across 360^o          
       XTMP  = VX(J2)*TPI-VX(J1)*TPI
       XTMP1 = VX(J2)-VX(J1)
       IF(XTMP1 >  180.0_SP)THEN
         XTMP = -360.0_SP*TPI+XTMP
       ELSE IF(XTMP1 < -180.0_SP)THEN
         XTMP =  360.0_SP*TPI+XTMP
       END IF

       EXFLUX       = DIJ*SXXIJ*DLTYC(I)
       PSXXPX(IA,K) = PSXXPX(IA,K) - EXFLUX
       PSXXPX(IB,K) = PSXXPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SXYIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSXYPY(IA,K) = PSXYPY(IA,K) + EXFLUX
       EXFLUX       = DIJ*SXYIJ*XTMP*COS(DEG2RAD*YC(IB))
       PSXYPY(IB,K) = PSXYPY(IB,K) - EXFLUX

       EXFLUX     = DIJ*SXYIJ*DLTYC(I)
       PSXYPX(IA,K) = PSXYPX(IA,K) - EXFLUX
       PSXYPX(IB,K) = PSXYPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SYYIJ*XTMP*COS(DEG2RAD*YC(IA))
       PSYYPY(IA,K) = PSYYPY(IA,K) + EXFLUX
       EXFLUX       = DIJ*SYYIJ*XTMP*COS(DEG2RAD*YC(IB))
       PSYYPY(IB,K) = PSYYPY(IB,K) - EXFLUX

#      else
       EXFLUX       = DIJ*SXXIJ*DLTYC(I)
       PSXXPX(IA,K) = PSXXPX(IA,K) - EXFLUX
       PSXXPX(IB,K) = PSXXPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SXYIJ*DLTXC(I)
       PSXYPY(IA,K) = PSXYPY(IA,K) + EXFLUX
       PSXYPY(IB,K) = PSXYPY(IB,K) - EXFLUX

       EXFLUX       = DIJ*SXYIJ*DLTYC(I)
       PSXYPX(IA,K) = PSXYPX(IA,K) - EXFLUX
       PSXYPX(IB,K) = PSXYPX(IB,K) + EXFLUX

       EXFLUX       = DIJ*SYYIJ*DLTXC(I)
       PSYYPY(IA,K) = PSYYPY(IA,K) + EXFLUX
       PSYYPY(IB,K) = PSYYPY(IB,K) - EXFLUX
#      endif     
     END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO

   WAVESTRX_3D = 0.0_SP
   WAVESTRY_3D = 0.0_SP
   
   CALL RADIATION_STRESS_Z(WAVE_ENERGY,KD,KDMAX,PSPXPZ,PSPYPZ)

   WAVESTRX_3D = PSXXPX + PSXYPY - PSPXPZ
   WAVESTRY_3D = PSXYPX + PSYYPY - PSPYPZ

!qxu set rediation stress limit to 200 Pa 01/19/2021
     WAVESTRX_3D = max(min(WAVESTRX_3D,200.0_SP),-200.0_SP)
     WAVESTRY_3D = max(min(WAVESTRY_3D,200.0_SP),-200.0_SP)
!qxu}

#  if defined (WET_DRY)
   DO I = 1,NT
     WAVESTRX_3D(I,:) = WAVESTRX_3D(I,:)*ISWETC(I)
     WAVESTRY_3D(I,:) = WAVESTRY_3D(I,:)*ISWETC(I)

!JQI
     IF(H1(I) <= 0.0_sp)THEN
       WAVESTRX_3D(I,:) = 0.0_SP
       WAVESTRY_3D(I,:) = 0.0_SP
     END IF
!JQI    
     IF(ISBCE(I) == 2)THEN
       WAVESTRX_3D(I,:) = 0.0_SP
       WAVESTRY_3D(I,:) = 0.0_SP
     END IF
   END DO  
#  endif
#  if !defined (TWO_D_MODEL)
   WAVESTRX_2D(:) = 0.0_SP; WAVESTRY_2D(:) = 0.0_SP
   DO I = 1,NT
     DO K=1,KBM1
        WAVESTRX_2D(I) = WAVESTRX_2D(I)+WAVESTRX_3D(I,K)
        WAVESTRY_2D(I) = WAVESTRY_2D(I)+WAVESTRY_3D(I,K)
     END DO
   END DO
#  else
!JQI202306   CALL RADIATION_STRESS_2D
#  endif
   WAVESTRX_2D = WAVESTRX_2D*RAMP
   WAVESTRY_2D = WAVESTRY_2D*RAMP
   WAVESTRX_3D = WAVESTRX_3D*RAMP
   WAVESTRY_3D = WAVESTRY_3D*RAMP

#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,WAVESTRX_3D,WAVESTRY_3D) 
#  endif
!Calculate stokes velocity
   U_STOKES_3D_TMP = 0.0_SP; U_STOKES_3D = 0.0_SP; U_STOKES_2D = 0.0_SP
   V_STOKES_3D_TMP = 0.0_SP; V_STOKES_3D = 0.0_SP; V_STOKES_2D = 0.0_SP
   DO I=1,M
    IF(KD(I) <= KDMAX)THEN
     DO K=1,KBM1
        FAC2 = 1.0_SP+ZZ(I,K)
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*COSH(2*KD(I)*FAC2)/SINH(2*KD(I))*WAVE_ENERGY(I)
# endif
        U_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_Y(I)
     ENDDO
    ELSE
     DO K=1,KBM1
        FAC2 = ZZ(I,K)
# if defined (WAVE_ROLLER)
        CFF2=2/WAVE_C(I)*EXP(KD(I)*FAC2)*(WAVE_ENERGY(I)+D(I)*GRAV_N(I)*ROLLA(I)/(WLEN(I)+eps1))
# else
        CFF2=2/WAVE_C(I)*EXP(KD(I)*FAC2)*WAVE_ENERGY(I)
# endif
        U_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_X(I)
        V_STOKES_3D_TMP(I,K)=CFF2*WAVE_NUMBER_Y(I)
     ENDDO
    END IF 
   ENDDO
   DO I=1,NT
     DO K=1,KBM1
       U_STOKES_3D(I,K)=(U_STOKES_3D_TMP(NV(I,1),K)+U_STOKES_3D_TMP(NV(I,2),K)+U_STOKES_3D_TMP(NV(I,3),K))/3.0_SP
       V_STOKES_3D(I,K)=(V_STOKES_3D_TMP(NV(I,1),K)+V_STOKES_3D_TMP(NV(I,2),K)+V_STOKES_3D_TMP(NV(I,3),K))/3.0_SP
       U_STOKES_2D(I)=U_STOKES_2D(I)+U_STOKES_3D(I,K)*DZ1(I,K)
       V_STOKES_2D(I)=V_STOKES_2D(I)+V_STOKES_3D(I,K)*DZ1(I,K)
     ENDDO
   ENDDO
#  if defined(MULTIPROCESSOR)
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U_STOKES_3D,V_STOKES_3D) 
   IF(PAR)CALL AEXCHANGE(EC,MYID,NPROCS,U_STOKES_2D,V_STOKES_2D) 
#  endif

   END IF

   RETURN
   END SUBROUTINE RADIATION_STRESS  
!==============================================================================|

   SUBROUTINE RADIATION_STRESS_Z(WAVE_ENERGY,KD,KDMAX,PSPXPZ,PSPYPZ) 

!==============================================================================|
   USE MOD_PREC
   USE ALL_VARS, ONLY : Z,Z1,VX,VY,N2E2D

   IMPLICIT NONE
   REAL(SP), INTENT(IN)  :: WAVE_ENERGY(0:MT),KD(0:MT),KDMAX
   REAL(SP), INTENT(OUT) :: PSPXPZ(0:NT,KB),PSPYPZ(0:NT,KB)
   REAL(SP)              :: SPX(KB),SPY(KB)

   REAL(SP) :: WAVE_ENERGY1(0:NT), KD1(0:NT)
   REAL(SP), DIMENSION(0:MT) :: O_COSH,O_SINH
   INTEGER  :: I,K,J,J1,J2,I1,I2,I3
   REAL(SP) :: FSS1,FCS1,FSC1,FCC1
   REAL(SP) :: CFF1,CFF2,FAC1,FAC2,FAC3
   REAL(SP) :: WEIJ,KDIJ,DIJ,SIJ
#  if defined (SPHERICAL)
   REAL(SP) :: XTMP,XTMP1
#  endif
!==============================================================================|

!----------INITIALIZE ARRAYS---------------------------------------------------!
    CALL N2E2D(WAVE_ENERGY,WAVE_ENERGY1)
    CALL N2E2D(KD,KD1)
    PSPXPZ  = 0.0_SP   ;PSPYPZ  = 0.0_SP  
    O_COSH  = 1.0_SP/COSH(KD)
    O_SINH  = 1.0_SP/SINH(KD)

   DO I = 1, N
     SPX = 0.0_SP; SPY = 0.0_SP
#    if defined (WET_DRY)
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#    endif
       I1=NV(I,1);I2=NV(I,2);I3=NV(I,3)
       IF(KD1(I) <= KDMAX)THEN
       
       DO K=1,KBM1
! Calculate some coefficients
         FAC1 = 1.0_SP+Z(I1,K) 
         FAC2 = 1.0_SP+Z(I2,K)
         FAC3 = 1.0_SP+Z(I3,K)
         FCC1 = (COSH(KD(I1)*FAC1)*O_COSH(I1)+COSH(KD(I2)*FAC2)*O_COSH(I2)+COSH(KD(I3)*FAC3)*O_COSH(I3))/3
         FCS1 = (COSH(KD(I1)*FAC1)*O_SINH(I1)+COSH(KD(I2)*FAC2)*O_SINH(I2)+COSH(KD(I3)*FAC3)*O_SINH(I3))/3
         FSC1 = (SINH(KD(I1)*FAC1)*O_COSH(I1)+SINH(KD(I2)*FAC2)*O_COSH(I2)+SINH(KD(I3)*FAC3)*O_COSH(I3))/3
         FSS1 = (SINH(KD(I1)*FAC1)*O_SINH(I1)+SINH(KD(I2)*FAC2)*O_SINH(I2)+SINH(KD(I3)*FAC3)*O_SINH(I3))/3
         CFF1=(FCC1-FSS1)*(FSS1*0.5_SP)
         CFF2=(FCC1-FSS1)*(FCS1*(1+Z1(I,K))*WAVE_ENERGY1(I)-WAVE_ENERGY1(I)*FSS1/TANH(KD1(I)))
         DO J = 1, 3
           J1=J+1-INT((J+1)/4)*3
           J2=J+2-INT((J+2)/4)*3
           WEIJ=0.5_SP*(WAVE_ENERGY(NV(I,J1))+WAVE_ENERGY(NV(I,J2)))*CFF1
           KDIJ=0.5_SP*(KD(NV(I,J1))+KD(NV(I,J2)))*CFF2
           SIJ=WEIJ+KDIJ
#          if defined (SPHERICAL)
           SPX(K)=SPX(K)-DELTUY(I,J)*SIJ
#          else
           SPX(K)=SPX(K)-(VY(NV(I,J2))-VY(NV(I,J1)))*SIJ
#          endif

#          if defined (SPHERICAL)
           XTMP  = VX(NV(I,J2))*TPI-VX(NV(I,J1))*TPI
           XTMP1 = VX(NV(I,J2))-VX(NV(I,J1))
           IF(XTMP1 >  180.0_SP)THEN
             XTMP = -360.0_SP*TPI+XTMP
           ELSE IF(XTMP1 < -180.0_SP)THEN
             XTMP =  360.0_SP*TPI+XTMP
           END IF  

           SPY(K)=SPY(K)+XTMP*COS(DEG2RAD*YC(I))*SIJ
#          else
           SPY(K)=SPY(K)+(VX(NV(I,J2))-VX(NV(I,J1)))*SIJ
#          endif
         END DO
       END DO
       
       ELSE
       
       DO K=1,KBM1
! Calculate some coefficients
         FAC1 = Z(I1,K) 
         FAC2 = Z(I2,K)
         FAC3 = Z(I3,K)
         FCC1 = (EXP(KD(I1)*FAC1)+EXP(KD(I2)*FAC2)+EXP(KD(I3)*FAC3))/3.0_SP
	 FCS1 = FCC1
         FSC1 = FCC1
         FSS1 = FCC1
	 
         CFF1=(FCC1-FSS1)*(FSS1*0.5_SP)
         CFF2=(FCC1-FSS1)*(FCS1*(1+Z1(I,K))*WAVE_ENERGY1(I)-WAVE_ENERGY1(I)*FSS1)
         DO J = 1, 3
           J1=J+1-INT((J+1)/4)*3
           J2=J+2-INT((J+2)/4)*3
           WEIJ=0.5_SP*(WAVE_ENERGY(NV(I,J1))+WAVE_ENERGY(NV(I,J2)))*CFF1
           KDIJ=0.5_SP*(KD(NV(I,J1))+KD(NV(I,J2)))*CFF2
           SIJ=WEIJ+KDIJ
#          if defined (SPHERICAL)
           SPX(K)=SPX(K)-DELTUY(I,J)*SIJ
#          else
           SPX(K)=SPX(K)-(VY(NV(I,J2))-VY(NV(I,J1)))*SIJ
#          endif

#          if defined (SPHERICAL)
           XTMP  = VX(NV(I,J2))*TPI-VX(NV(I,J1))*TPI
           XTMP1 = VX(NV(I,J2))-VX(NV(I,J1))
           IF(XTMP1 >  180.0_SP)THEN
             XTMP = -360.0_SP*TPI+XTMP
           ELSE IF(XTMP1 < -180.0_SP)THEN
             XTMP =  360.0_SP*TPI+XTMP
           END IF  

           SPY(K)=SPY(K)+XTMP*COS(DEG2RAD*YC(I))*SIJ
#          else
           SPY(K)=SPY(K)+(VX(NV(I,J2))-VX(NV(I,J1)))*SIJ
#          endif
         END DO
       END DO
       
       END IF
       
       DO K = 1,KBM1
         PSPXPZ(I,K) = SPX(K)-SPX(K+1) 
	 PSPYPZ(I,K) = SPY(K)-SPY(K+1)
       END DO
#    if defined (WET_DRY)
     END IF
#    endif
   END DO

   RETURN
   END SUBROUTINE RADIATION_STRESS_Z
!==============================================================================|

  subroutine StateWriteVTK(state, prefix, rc)

    ! input arguments
    type(ESMF_State), intent(in) :: state
    character(len=*), intent(in) :: prefix
    integer, intent(out), optional :: rc

    ! local variables
    integer :: i, itemCount
    type(ESMF_Field) :: field
    character(ESMF_MAXSTR), allocatable :: itemNameList(:)
    character(len=*),parameter  :: subname='(fvcom_cap:StateWriteVTK)'

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    ! get number of fields in the state
    call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    ! get item names
    if (.not. allocated(itemNameList)) allocate(itemNameList(itemCount))

    call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    ! loop over fields and write them
    do i = 1, itemCount
       ! get field
       call ESMF_StateGet(state, itemName=trim(itemNameList(i)), field=field, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, &
           file=__FILE__)) &
           return  ! bail out

       ! write it
       call ESMF_FieldWriteVTK(field, trim(prefix)//'_'//trim(itemNameList(i)), rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, &
           file=__FILE__)) &
           return  ! bail out
    end do

    ! clean temporary variables
    if (allocated(itemNameList)) deallocate(itemNameList)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine StateWriteVTK

end module

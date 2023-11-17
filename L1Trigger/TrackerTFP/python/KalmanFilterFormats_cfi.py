# configuration of internal KF variable bases which can be shifted by powers of 2 w.r.t. KF output track parameter
# TrackerTFPProducer_params.PrintKFDebug printouts unused MSB for each variable, so that one could consider decreasing the basseshift by that amount
# numerical instabillity (negative C00, C11, C22, C33) requires smaller baseshifts of related variables (rx, Sxx, Kxx, Rxx, invRxx)
# if a variable overflows an Exception will be thrown and the corresponding baseshift needs to be increased.
import FWCore.ParameterSet.Config as cms

TrackTriggerKalmanFilterFormats_params = cms.PSet (

  BaseShiftx0           = cms.int32(  -4 ),
  BaseShiftx1           = cms.int32(  -9 ),
  BaseShiftx2           = cms.int32(  -4 ),
  BaseShiftx3           = cms.int32(  -5 ),
  BaseShiftv0           = cms.int32(  -3 ),
  BaseShiftv1           = cms.int32(   9 ),
  BaseShiftr0           = cms.int32(  -9 ),
  BaseShiftr1           = cms.int32(  -2 ),
  BaseShiftS00          = cms.int32(   1 ),
  BaseShiftS01          = cms.int32(  -5 ),
  BaseShiftS12          = cms.int32(   2 ),
  BaseShiftS13          = cms.int32(   1+1 ),
  BaseShiftK00          = cms.int32( -17 ),
  BaseShiftK10          = cms.int32( -22 ),
  BaseShiftK21          = cms.int32( -23 ),
  BaseShiftK31          = cms.int32( -23 ),
  BaseShiftR00          = cms.int32(  -2 ),
  BaseShiftR11          = cms.int32(   9 ),
  BaseShiftInvR00Approx = cms.int32( -29 ),
  BaseShiftInvR11Approx = cms.int32( -40 ),
  BaseShiftInvR00Cor    = cms.int32( -15 ),
  BaseShiftInvR11Cor    = cms.int32( -15 ),
  BaseShiftInvR00       = cms.int32( -26 ),
  BaseShiftInvR11       = cms.int32( -33 ),
  BaseShiftC00          = cms.int32(   5 ),
  BaseShiftC01          = cms.int32(  -1 ),
  BaseShiftC11          = cms.int32(  -6 ),
  BaseShiftC22          = cms.int32(   1 ),
  BaseShiftC23          = cms.int32(   0 ),
  BaseShiftC33          = cms.int32(   1 ),
  BaseShiftr02          = cms.int32(  -2 ),
  BaseShiftr12          = cms.int32(  11 ),
  BaseShiftchi20        = cms.int32( -14+9 ),
  BaseShiftchi21        = cms.int32( -14+9 )

)
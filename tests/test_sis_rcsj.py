import numpy as np
from shs.physics.josephson_rcsj import *
from shs.utils.constants import SUPERCONDUCTING_FLUX_QUANTUM

def test_ambegaokar_baratoff_is_positive_below_tc_and_zero_at_tc():
    assert ambegaokar_baratoff_critical_current(4.2,9.2,20)>0
    assert ambegaokar_baratoff_critical_current(9.2,9.2,20)==0

def test_voltage_phase_relation_and_current_components_close():
    j=RCSJJunction(2e-6,20,2e-13); phase=.4; rate=3e9; accel=7e17
    p=junction_components(j,phase,rate,accel)
    assert np.isclose(p["voltage_V"], SUPERCONDUCTING_FLUX_QUANTUM*rate/(2*np.pi))
    total=p["supercurrent_A"]+p["resistive_current_A"]+p["capacitive_current_A"]
    _,calculated=current_biased_rhs(j,phase,rate,total)
    assert np.isclose(calculated,accel)

def test_half_flux_quantum_suppresses_uniform_junction_critical_current():
    j=RCSJJunction(2e-6,20,0,magnetic_flux_Wb=SUPERCONDUCTING_FLUX_QUANTUM)
    assert abs(j.effective_critical_current_A)<1e-20

def test_nearby_vortex_changes_extended_junction_interference():
    plain,_,_=extended_junction_interference(2e-6,0,0,2e-7,[])
    vortex,_,_=extended_junction_interference(2e-6,0,0,2e-7,[{"x_m":-2e-7,"y_m":2e-7,"charge":1,"electrode":"left"}])
    assert abs(vortex)<abs(plain)
    assert not np.isclose(np.angle(vortex),np.angle(plain))

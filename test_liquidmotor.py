from rocketpy.motors.LiquidMotor import Tank, LiquidMotor, MassBasedTank, UllageBasedTank, MassFlowRateBasedTank
from rocketpy.motors.Fluid import Fluid
from math import isclose
import numpy as np


# @PBales1
def test_mass_based_motor():
    lox = Fluid(name = "LOx", density = 1141.7, quality = 1.0) #Placeholder quality value
    propane = Fluid(name = "Propane", density = 493, quality = 1.0) #Placeholder quality value
    n2 = Fluid(name = "Nitrogen Gas", density = 51.75, quality = 1.0) #Placeholder quality value; density value may be estimate
    
    example_motor = MassBasedTank("Example Tank", 0.1540, 0.66, 0.7, "Placeholder", "Placeholder", lox, n2) 
    #Need docs to be pushed + tank dimension values


# @curtisjhu
def test_ullage_based_motor():
    lox = Fluid(name = "LOx", density = 1141, quality = 1.0) #Placeholder quality value
    propane = Fluid(name = "Propane", density = 493, quality = 1.0) #Placeholder quality value
    n2 = Fluid(name = "Nitrogen Gas", density = 51.75, quality = 1.0) #Placeholder quality value; density value may be estimate

    ullageData = []
    ullageTank = UllageBasedTank("Ullage Tank",  diameter=3, height=4, endcap="flat", gas=n2, liquid=lox, ullage=ullageData)
    
    ullageTank.centerOfMass()
    ullageTank.netMassFlowRate()
    ullageTank.mass()
    ullageTank.liquidVolume()

# @gautamsaiy
def test_mfr_tank_basic1():
    def test(t, a):
        for i in np.arange(0, 10, .2):
            assert isclose(t.getValue(i), a(i))
            # print(t.getValue(i), a(i))
    
    def test_nmfr():
        nmfr = lambda x: liquid_mass_flow_rate_in + gas_mass_flow_rate_in - liquid_mass_flow_rate_out - gas_mass_flow_rate_out
        test(t.netMassFlowRate(), nmfr)

    def test_mass():
        m = lambda x: (initial_liquid_mass + (liquid_mass_flow_rate_in - liquid_mass_flow_rate_out) * x) + \
            (initial_gas_mass + (gas_mass_flow_rate_in - gas_mass_flow_rate_out) * x)
        lm = t.mass()
        test(lm, m)

    def test_uh():
        uh = lambda x: ((initial_liquid_mass + (liquid_mass_flow_rate_in - liquid_mass_flow_rate_out) * x) / lox.density)
        test((t.initial_liquid_mass + t.liquid_mass_flow_rate_in - t.liquid_mass_flow_rate_out) / t.liquid.density, uh)


    tank_radius_function = {(0, 5): 1}
    lox = Fluid(name = "LOx", density = 1141, quality = 1.0) #Placeholder quality value
    n2 = Fluid(name = "Nitrogen Gas", density = 51.75, quality = 1.0) #Placeholder quality value; density value may be estimate
    initial_liquid_mass = 5
    initial_gas_mass = .05
    liquid_mass_flow_rate_in = .1
    gas_mass_flow_rate_in = .01
    liquid_mass_flow_rate_out = .2
    gas_mass_flow_rate_out = .02

    t = MassFlowRateBasedTank("Test Tank", tank_radius_function,
            initial_liquid_mass, initial_gas_mass, liquid_mass_flow_rate_in,
            gas_mass_flow_rate_in, liquid_mass_flow_rate_out, 
            gas_mass_flow_rate_out, lox, n2)
    # test_nmfr()
    # test_mass()
    test_uh()




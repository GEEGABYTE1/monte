import streamlit as st 
import numpy as np
from numpy.random import normal, choice
from time import process_time, time
import matplotlib.pyplot as plt
from rocketpy import Environment, SolidMotor, Rocket, Flight

st.set_page_config(page_title='Rocket Sim', layout='wide')
st.title("Rocket Flight Sim")

analysis_parameters = {
    "rocket_mass": (7.257, 0.001),
    "rocket_inertia_11": (3.675, 0.03675),
    "rocket_inertia_33": (0.007, 0.00007),
    "motor_dry_mass": (1.000, 0.001),
    "motor_inertia_11": (1.675, 0.01675),
    "motor_inertia_33": (0.003, 0.00003),
    "motor_dry_mass_position": (0.5, 0.001),
    "impulse": (1415.15, 35.3),
    "burn_time": (5.274, 1),
    "nozzle_radius": (21.642 / 1000, 0.5 / 1000),
    "throat_radius": (8 / 1000, 0.5 / 1000),
    "grain_separation": (6 / 1000, 1 / 1000),
    "grain_density": (1707, 50),
    "grain_outer_radius": (21.4 / 1000, 0.375 / 1000),
    "grain_initial_inner_radius": (9.65 / 1000, 0.375 / 1000),
    "grain_initial_height": (120 / 1000, 1 / 1000),
    "radius": (40.45 / 1000, 0.001),
    "nozzle_position": (-1.024, 0.001),
    "grains_center_of_mass_position": (-0.571, 0.001),
    "power_off_drag": (0.9081 / 1.05, 0.033),
    "power_on_drag": (0.9081 / 1.05, 0.033),
    "nose_length": (0.274, 0.001),
    "nose_distance_to_CM": (1.134, 0.001),
    "fin_span": (0.077, 0.0005),
    "fin_root_chord": (0.058, 0.0005),
    "fin_tip_chord": (0.018, 0.0005),
    "fin_distance_to_CM": (-0.906, 0.001),
    "inclination": (84.7, 1),
    "heading": (53, 2),
    "rail_length": (5.7, 0.0005),
    "ensemble_member": list(range(10)),
    "cd_s_drogue": (0.349 * 1.3, 0.07),
    "lag_rec": (1, 0.5),
    "lag_se": (0.73, 0.16),
}



def flight_settings(analysis_parameters, total_number):
    i = 0
    while i < total_number:
        flight_setting = {}
        for parameter_key, parameter_value in analysis_parameters.items():
            if isinstance(parameter_value, tuple):
                flight_setting[parameter_key] = normal(*parameter_value)
            else:
                flight_setting[parameter_key] = normal(*parameter_value)
        
        if flight_setting["lag_rec"] < 0 or flight_setting["lag_se"] < 0:
            continue 
        i += 1
        yield flight_setting


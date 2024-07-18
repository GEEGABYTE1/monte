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


def run_simulation(number_of_simulations):
    dispersion_general_results = []
    dispersion_results = {
        "out_of_rail_time": [],
        "out_of_rail_velocity": [],
        "apogee_time": [],
        "apogee_altitude": [],
        "apogee_x": [],
        "apogee_y": [],
        "impact_time": [],
        "impact_x": [],
        "impact_y": [],
        "impact_velocity": [],
        "initial_static_margin": [],
        "out_of_rail_static_margin": [],
        "final_static_margin": [],
        "number_of_events": [],
        "max_velocity": [],
        "drogue_triggerTime": [],
        "drogue_inflated_time": [],
        "drogue_inflated_velocity": [],
        "execution_time": []
    }
    initial_wall_time = time()
    initial_cpu_time = process_time() 

    Env = Environment(date=(2019, 8, 10, 21), latitude=-23.36311, longitude=-48.011389)
    Env.set_elevation(668)
    Env.max_expected_height = 1500 
    Env.set_atmospheric_model(
        type='Ensemble',
        file="monte_carlo_analysis_inputs/LASC2019_reanalysis.nc",
        dictionary="ECMWF"
    )

    def drogue_trigger(p, h, y):
        vertical_velocity = y[5]
        return True if vertical_velocity < 0 else False 
    
    for setting in flight_settings(analysis_parameters, number_of_simulations):
        start_time = process_time()
        try:
            Env.select_ensemble_member(setting["ensemble_member"])
            Keron = SolidMotor(
                thrust_source="monte_carlo_analysis_inputs/thrustCurve.csv",
                burn_time=5.274,
                reshape_thrust_curve=(setting["burn_time"], setting["impulse"]),
                nozzle_radius=setting["nozzle_radius"],
                throat_radius=setting["throat_radius"],
                grain_number=6,
                grain_separation=setting["grain_separation"],
                grain_density=setting["grain_density"],
                grain_outer_radius=setting["grain_outer_radius"],
                grain_initial_inner_radius=setting["grain_initial_inner_radius"],
                grain_initial_height=setting["grain_initial_height"],
                interpolation_method="linear",
                coordinate_system_orientation="nozzle_to_combustion_chamber",
                nozzle_position=setting["nozzle_position"],
                grains_center_of_mass_position=setting["grains_center_of_mass_position"],
                dry_mass=setting["motor_dry_mass"],
                dry_inertia=(
                    setting["motor_inertia_11"],
                    setting["motor_inertia_11"],
                    setting["motor_inertia_33"],
                ),
                center_of_dry_mass_position=setting["motor_dry_mass_position"],
            )
            Valetudo = Rocket(
                radius=setting["radius"],
                mass=setting["rocket_mass"],
                inertia=(
                    setting["rocket_inertia_11"],
                    setting["rocket_inertia_11"],
                    setting["rocket_inertia_33"],
                ),
                power_off_drag="monte_carlo_analysis_inputs/Cd_PowerOff.csv",
                power_on_drag="monte_carlo_analysis_inputs/Cd_PowerOn.csv",
                center_of_mass_without_motor=0,
                coordinate_system_orientation="tail_to_nose",
            )
            Valetudo.set_rail_buttons(0.224, -0.93, 30)
            Valetudo.add_motor(Keron, position=0)
            Valetudo.power_off_drag *= setting["power_off_drag"]
            Valetudo.power_on_drag *= setting["power_on_drag"]
            NoseCone = Valetudo.add_nose(
                length=setting["nose_length"],
                kind="vonKarman",
                position=setting["nose_distance_to_CM"] + setting["nose_length"],
            )
            FinSet = Valetudo.add_trapezoidal_fins(
                n=3,
                span=setting["fin_span"],
                root_chord=setting["fin_root_chord"],
                tip_chord=setting["fin_tip_chord"],
                position=setting["fin_distance_to_CM"],
                cant_angle=0,
                airfoil=None,
            )
            Drogue = Valetudo.add_parachute(
                "Drogue",
                cd_s=setting["cd_s_drogue"],
                trigger=drogue_trigger,
                sampling_rate=105,
                lag=setting["lag_rec"] + setting["lag_se"],
                noise=(0, 8.3, 0.5),
            )
            test_flight = Flight(
                rocket=Valetudo,
                environment=Env,
                rail_length=setting["rail_length"],
                inclination=setting["inclination"],
                heading=setting["heading"],
                max_time=600,
            )
            flight_data = test_flight.all_info()
            exec_time = process_time() - start_time
            flight_result = {
                "out_of_rail_time": flight_data.out_of_rail_time,
                "out_of_rail_velocity": flight_data.out_of_rail_velocity,
                "max_velocity": flight_data.speed.max,
                "apogee_time": flight_data.apogee_time,
                "apogee_altitude": flight_data.apogee,
                "apogee_x": flight_data.x_apogee,
                "apogee_y": flight_data.y_apogee,
                "impact_time": flight_data.t_final,
                "impact_x": flight_data.x_final,
                "impact_y": flight_data.y_final,
                "impact_velocity": flight_data.final_velocity,
                "initial_static_margin": test_flight.static_margin(0),
                "out_of_rail_static_margin": test_flight.static_margin(
                    flight_data.out_of_rail_time
                ),
                "final_static_margin": test_flight.static_margin(flight_data.apogee_time),
                "drogue_triggerTime": flight_data.drogue_deployment_time,
                "drogue_inflated_time": flight_data.drogue_inflation_time,
                "drogue_inflated_velocity": flight_data.drogue_inflation_velocity,
                "number_of_events": len(flight_data.time),
                "execution_time": exec_time,
            }
            for key in dispersion_results:
                dispersion_results[key].append(flight_result[key])
        except Exception as e:
            st.error(f"Simulation {setting} failed due to error: {str(e)}")
    
    return dispersion_results


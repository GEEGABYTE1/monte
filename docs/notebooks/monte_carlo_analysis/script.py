from time import process_time, time 
from rocketpy import Environment, SolidMotor, Rocket, Flight 

import numpy as np
from numpy.random import normal, choice 
from IPython.display import display

import matplotlib as mpl
import matplotlib.pyplot as plt
import streamlit as st

from imageio import imread 
from matplotlib.patches import Ellipse


mpl.rcParams["figure.figsize"] = [8, 5]
mpl.rcParams["figure.dpi"] = 120
mpl.rcParams["font.size"] = 14
mpl.rcParams["legend.fontsize"] = 14
mpl.rcParams["figure.titlesize"] = 14

st.title("Monte Carlo Rocket Sim GUI for RocketPy")

rocket_mass = 7.257
rocket_mass_std = 0.001

rocket_mass = st.sidebar.number_input("Rocket Mass Mean (kg)", value=7.257, step=0.01)
rocket_mass_std = st.sidebar.number_input("Rocket Mass Std Dev (kg)", value=0.001, step=0.01)
rocket_inertia_11 = st.sidebar.number_input("Rocket's Inertia Moment Perpendicular to its Axis (kg*m^2) ", value=3.675, step=0.01)
rocket_inertia_11_std = st.sidebar.number_input("Rocket's Inertia Moment Perpendicular to its Axis  Std (kg*m^2)", value=0.03675, step=0.01)
rocket_inertia_33 = st.sidebar.number_input("Rocket's Inertia Moment Relative to its Axis (kg*m^2) ", value=0.007, step=0.01)
rocket_inertia_33_std = st.sidebar.number_input("Rocket's Inertia Moment Relative to its Axis  Std (kg*m^2)", value=0.00007, step=0.01)
motor_dry_mass = st.sidebar.number_input("Motor Dry Mass without Propellant (kg)", value=1.000, step=0.01)
motor_dry_mass_std = st.sidebar.number_input("Motor Dry Mass without Propellant Std (kg)", value=0.001, step=0.01)
motor_inertia_11 = st.sidebar.number_input("Motor's Dry Inertia Moment Perpendicular to its Axis (kg*m^2)", value=1.675, step=0.01)
motor_intertia_11_std =  st.sidebar.number_input("Motor's Dry Inertia Moment Perpendicular to its Axis (kg*m^2) ", value=0.01675, step=0.01)
motor_inertia_33 = st.sidebar.number_input("Motor's Dry Inertia Moment Relative to its Axis (kg*m^2)", value=0.003, step=0.01)
motor_intertia_33_std =  st.sidebar.number_input("Motor's Dry Inertia Moment Relative to its Axis (kg*m^2) ", value=0.00003, step=0.01)
motor_dry_mass_position = st.sidebar.number_input("Distance between Rocket's Center of Dry Mass and motor's center of dry mass (m)".title(), value=0.5, step=0.01)
motor_dry_mass_position_std = st.sidebar.number_input("Distance between Rocket's Center of Dry Mass and motor's center of dry mass Std (m)".title(), value=0.001, step=0.01)
impulse = st.sidebar.number_input("Motor total impulse (N*s)".title(), value=1415.15, step=0.01)
impulse_std = st.sidebar.number_input("Motor total impulse Std (N*s)".title(), value=35.3, step=0.01)
burn_time = st.sidebar.number_input("Motor burn out time (s)".title(), value=5.274, step=0.01)
burn_time_std = st.sidebar.number_input("Motor burn out time (s) Std".title(), value=1.00, step=0.01)
nozzle_radius = st.sidebar.number_input("Motor nozzle radius (m)".title(), value=21.642, step=0.01)
nozzle_radius_std = st.sidebar.number_input("Motor nozzle radius Std (m)".title(), value=0.5, step=0.01)
throat_radius = st.sidebar.number_input("Motor's throat radius (m)".title(), value=8.00, step=0.01)
throat_radius_std = st.sidebar.number_input("Motor's throat radius Std(m)".title(), value=0.5, step=0.01)
grain_separation = st.sidebar.number_input("Motor's grain separation (axial distance between two grains) (m)".title(), value=6.00, step=0.01)
grain_separation_std = st.sidebar.number_input("Motor's grain separation (axial distance between two grains) Std (m)".title(), value=1.00, step=0.01)
grain_density = st.sidebar.number_input("Motor's grain density (kg/m^3)".title(), value=1707.00, step=0.01)
grain_density_std = st.sidebar.number_input("Motor's grain density Std (kg/m^3)".title(), value=50.00, step=0.01)
grain_outer_radius = st.sidebar.number_input("Motor's grain outer radius (m)".title(), value=21.4, step=0.01)
grain_outer_radius_std = st.sidebar.number_input("Motor's grain outer radius Std (m)".title(), value=0.375, step=0.01)
grain_initial_inner_radius = st.sidebar.number_input("Motor's grain inner radius (m)".title(), value=9.65, step=0.01)
grain_initial_inner_radius_std = st.sidebar.number_input("Motor's grain inner radius Std (m)".title(), value=0.375, step=0.01)
grain_initial_height = st.sidebar.number_input("Motor's grain height (m)".title(), value=120.00, step=0.01)
grain_initial_height_std = st.sidebar.number_input("Motor's grain height Std (m)".title(), value=1.00, step=0.01)
radius = st.sidebar.number_input("Rocket's radius (kg*m^2)".title(), value=40.45, step=0.01)
radius_std = st.sidebar.number_input("Rocket's radius Std (kg*m^2)".title(), value=0.001, step=0.01)
nozzle_position = st.sidebar.number_input("Distance between rocket's center of dry mass and nozzle exit plane (m) (negative)".title(), value=-1.024, step=0.01)
nozzle_position_std = st.sidebar.number_input("Distance between rocket's center of dry mass and nozzle exit plane Std (m) (negative)".title(), value=0.001, step=0.01)
power_off_drag = st.sidebar.number_input("Multiplier for rocket's drag curve. Usually has a mean value of 1 and a uncertainty of 5% to 10%".title(), value=0.9081, step=0.01)
power_off_drag_std = st.sidebar.number_input("Multiplier for rocket's drag curve Std".title(), value=0.033, step=0.01)
power_on_drag = st.sidebar.number_input("Multiplier for rocket's drag curve for ON drag".title(), value=0.90810, step=0.01)
power_on_drag_std = st.sidebar.number_input("Multiplier for rocket's drag curve for ON drag Std".title(), value=0.033, step=0.01)
nose_length = st.sidebar.number_input("Rocket's nose cone length (m)".title(), value=0.274, step=0.01)
nose_length_std = st.sidebar.number_input("Rocket's nose cone length Std (m)".title(), value=0.001, step=0.01)
nose_distance_to_CM = st.sidebar.number_input("Axial distance between rocket's center of dry mass and nearest point in its nose cone (m)".title(), value=1.134, step=0.01)
nose_distance_to_CM_std = st.sidebar.number_input("Axial distance between rocket's center of dry mass and nearest point in its nose cone Std (m)".title(), value=0.001, step=0.01)
fin_span  = st.sidebar.number_input("Fin span (m)".title(), value=0.077, step=0.01)
fin_span_std = st.sidebar.number_input("Fin span Std (m)".title(), value=0.0005, step=0.01)
fin_root_chord = st.sidebar.number_input("Fin Root Chord (m)".title(), value=0.058, step=0.01)
fin_root_chord_std = st.sidebar.number_input("Fin Root Chord Std ".title(), value=0.0005, step=0.01)
fin_tip_chord = st.sidebar.number_input("Fin Tip Chord (m)".title(), value=0.018, step=0.01)
fin_tip_chord_std = st.sidebar.number_input("Fin Tip Chord Std (m)".title(), value=0.0005, step=0.01)
fin_distance_to_CM = st.sidebar.number_input("Axial distance between rocket's center of dry mass and nearest point in its fin (m)".title(), value=-0.906, step=0.01)
fin_distance_to_CM_std = st.sidebar.number_input("Axial distance between rocket's center of dry mass and nearest point in its fin Std (m)".title(), value=0.001, step=0.01)
inclination = st.sidebar.number_input("Launch rail inclination angle relative to the horizontal plane (degrees)".title(), value=84.70, step=0.01)
inclination_std = st.sidebar.number_input("Launch rail inclination angle relative to the horizontal plane Std (degrees)".title(), value=1.00, step=0.01)
heading = st.sidebar.number_input("Launch rail heading relative to north (degrees)".title(), value=53.00, step=0.01)
heading_std = st.sidebar.number_input("Launch rail heading relative to north Std (degrees)".title(), value=2.00, step=0.01)
rail_length = st.sidebar.number_input("Launch rail length (m)".title(), value=5.7, step=0.01)
rail_length_std = st.sidebar.number_input("Launch rail length Std (m)".title(), value=0.0005, step=0.01)
cd_s_drogue = st.sidebar.number_input("Drag coefficient times reference area for the drogue chute (m^2)".title(), value=0.349, step=0.01)
cd_s_drogue_std = st.sidebar.number_input("Drag coefficient times reference area for the drogue chute Std (m^2)".title(), value=0.07, step=0.01)
lag_rec = st.sidebar.number_input("Time delay between parachute ejection signal is detected and parachute is inflated (s)".title(), value=1.00, step=0.01)
lag_rec_std = st.sidebar.number_input("Time delay between parachute ejection signal is detected and parachute is inflated Std (s)".title(), value=0.50, step=0.01)
lag_se = st.sidebar.number_input("Time delay between sensor signal is received and ejection signal is fired (s)".title(), value=0.73, step=0.01)
lag_se_std = st.sidebar.number_input("Time delay between sensor signal is received and ejection signal is fired Std (s)".title(), value=0.16, step=0.01)
grains_center_of_mass_position = st.sidebar.number_input("Distance between rocket's center of dry mass and and center of propellant mass (m) (negative)".title(), value=-0.73, step=0.01)
grains_center_of_mass_position_std = st.sidebar.number_input("Distance between rocket's center of dry mass and and center of propellant mass Std (m) (negative)".title(), value=0.001, step=0.01)
if st.button("Run Simulation"): 
        analysis_parameters = {
            # Mass Details
            # Rocket's dry mass without motor (kg) and its uncertainty (standard deviation)
            "rocket_mass": (rocket_mass, rocket_mass_std),
            # Rocket's inertia moment perpendicular to its axis (kg*m^2)
            "rocket_inertia_11": (rocket_inertia_11, rocket_inertia_11_std),
            # Rocket's inertia moment relative to its axis (kg*m^2)
            "rocket_inertia_33": (rocket_inertia_33, rocket_inertia_33_std),
            # Motors's dry mass without propellant (kg) and its uncertainty (standard deviation)
            "motor_dry_mass": (motor_dry_mass, motor_dry_mass_std),
            # Motor's dry inertia moment perpendicular to its axis (kg*m^2)
            "motor_inertia_11": (motor_inertia_11, motor_intertia_11_std),
            # Motors's dry inertia moment relative to its axis (kg*m^2)
            "motor_inertia_33": (motor_inertia_33, motor_intertia_33_std),
            # Distance between rocket's center of dry mass and motor's center of dry mass (m)
            "motor_dry_mass_position": (motor_dry_mass_position, motor_dry_mass_position_std),
            # Propulsion Details - run help(SolidMotor) for more information
            # Motor total impulse (N*s)
            "impulse": (impulse, impulse_std),
            # Motor burn out time (s)
            "burn_time": (burn_time, burn_time_std),
            # Motor's nozzle radius (m)
            "nozzle_radius": (nozzle_radius / 1000, nozzle_radius_std / 1000),
            # Motor's nozzle throat radius (m)
            "throat_radius": (throat_radius / 1000, throat_radius_std / 1000),
            # Motor's grain separation (axial distance between two grains) (m)
            "grain_separation": (grain_separation / 1000, grain_separation_std / 1000),
            # Motor's grain density (kg/m^3)
            "grain_density": (grain_density, grain_density_std),
            # Motor's grain outer radius (m)
            "grain_outer_radius": (grain_outer_radius / 1000, grain_outer_radius_std / 1000),
            # Motor's grain inner radius (m)
            "grain_initial_inner_radius": (grain_initial_inner_radius / 1000, grain_initial_inner_radius_std / 1000),
            # Motor's grain height (m)
            "grain_initial_height": (grain_initial_height / 1000, grain_initial_height_std / 1000),
            # Aerodynamic Details - run help(Rocket) for more information
            # Rocket's radius (kg*m^2)
            "radius": (radius / 1000, radius_std),
            # Distance between rocket's center of dry mass and nozzle exit plane (m) (negative)
            "nozzle_position": (nozzle_position, nozzle_position_std),
            # Distance between rocket's center of dry mass and and center of propellant mass (m) (negative)
            "grains_center_of_mass_position": (grains_center_of_mass_position, grains_center_of_mass_position_std),
            # Multiplier for rocket's drag curve. Usually has a mean value of 1 and a uncertainty of 5% to 10%
            "power_off_drag": (power_off_drag / 1.05, power_off_drag_std),
            # Multiplier for rocket's drag curve. Usually has a mean value of 1 and a uncertainty of 5% to 10%
            "power_on_drag": (power_on_drag / 1.05, power_on_drag_std),
            # Rocket's nose cone length (m)
            "nose_length": (nose_length, nose_length_std),
            # Axial distance between rocket's center of dry mass and nearest point in its nose cone (m)
            "nose_distance_to_CM": (nose_distance_to_CM, nose_distance_to_CM_std),
            # Fin span (m)
            "fin_span": (fin_span, fin_span_std),
            # Fin root chord (m)
            "fin_root_chord": (fin_root_chord, fin_root_chord_std),
            # Fin tip chord (m)
            "fin_tip_chord": (fin_tip_chord, fin_tip_chord_std),
            # Axial distance between rocket's center of dry mass and nearest point in its fin (m)
            "fin_distance_to_CM": (fin_distance_to_CM, fin_distance_to_CM_std),
            # Launch and Environment Details - run help(Environment) and help(Flight) for more information
            # Launch rail inclination angle relative to the horizontal plane (degrees)
            "inclination": (inclination, inclination_std),
            # Launch rail heading relative to north (degrees)
            "heading": (heading, heading_std),
            # Launch rail length (m)
            "rail_length": (rail_length, rail_length_std),
            # Members of the ensemble forecast to be used
            "ensemble_member": list(range(10)),
            # Parachute Details - run help(Rocket) for more information
            # Drag coefficient times reference area for the drogue chute (m^2)
            "cd_s_drogue": (cd_s_drogue * 1.3, cd_s_drogue_std),
            # Time delay between parachute ejection signal is detected and parachute is inflated (s)
            "lag_rec": (lag_rec, lag_rec_std),
            # Electronic Systems Details - run help(Rocket) for more information
            # Time delay between sensor signal is received and ejection signal is fired (s)
            "lag_se": (lag_se, lag_se_std),
        }

        def flight_settings(analysis_parameters, total_number):
            i = 0
            while i < total_number:
                # Generate a flight setting
                flight_setting = {}
                for parameter_key, parameter_value in analysis_parameters.items():
                    if type(parameter_value) is tuple:
                        flight_setting[parameter_key] = normal(*parameter_value)
                    else:
                        flight_setting[parameter_key] = choice(parameter_value)

                # Skip if certain values are negative, which happens due to the normal curve but isnt realistic
                if flight_setting["lag_rec"] < 0 or flight_setting["lag_se"] < 0:
                    continue

                # Update counter
                i += 1
                # Yield a flight setting
                yield flight_setting


        def export_flight_data(flight_setting, flight_data, exec_time):
            # Generate flight results
            flight_result = {
                "out_of_rail_time": flight_data.out_of_rail_time,
                "out_of_rail_velocity": flight_data.out_of_rail_velocity,
                "max_velocity": flight_data.speed.max,
                "apogee_time": flight_data.apogee_time,
                "apogee_altitude": flight_data.apogee - Env.elevation,
                "apogee_x": flight_data.apogee_x,
                "apogee_y": flight_data.apogee_y,
                "impact_time": flight_data.t_final,
                "impact_x": flight_data.x_impact,
                "impact_y": flight_data.y_impact,
                "impact_velocity": flight_data.impact_velocity,
                "initial_static_margin": flight_data.rocket.static_margin(0),
                "out_of_rail_static_margin": flight_data.rocket.static_margin(
                    flight_data.out_of_rail_time
                ),
                "final_static_margin": flight_data.rocket.static_margin(
                    flight_data.rocket.motor.burn_out_time
                ),
                "number_of_events": len(flight_data.parachute_events),
                "execution_time": exec_time,
            }

            # Take care of parachute results
            if len(flight_data.parachute_events) > 0:
                flight_result["drogue_triggerTime"] = flight_data.parachute_events[0][0]
                flight_result["drogue_inflated_time"] = (
                    flight_data.parachute_events[0][0] + flight_data.parachute_events[0][1].lag
                )
                flight_result["drogue_inflated_velocity"] = flight_data.speed(
                    flight_data.parachute_events[0][0] + flight_data.parachute_events[0][1].lag
                )
            else:
                flight_result["drogue_triggerTime"] = 0
                flight_result["drogue_inflated_time"] = 0
                flight_result["drogue_inflated_velocity"] = 0

            # Write flight setting and results to file
            dispersion_input_file.write(str(flight_setting) + "\n")
            dispersion_output_file.write(str(flight_result) + "\n")


        def export_flight_error(flight_setting):
            dispersion_error_file.write(str(flight_setting) + "\n")


        # Basic analysis info
        filename = "monte_carlo_analysis_outputs/valetudo_rocket_v0"
        number_of_simulations = 100

        # Create data files for inputs, outputs and error logging
        dispersion_error_file = open(str(filename) + ".disp_errors.txt", "w")
        dispersion_input_file = open(str(filename) + ".disp_inputs.txt", "w")
        dispersion_output_file = open(str(filename) + ".disp_outputs.txt", "w")

        # Initialize counter and timer
        i = 0

        initial_wall_time = time()
        initial_cpu_time = process_time()

        # Define basic Environment object
        Env = Environment(date=(2019, 8, 10, 21), latitude=-23.363611, longitude=-48.011389)
        Env.set_elevation(668)
        Env.max_expected_height = 1500
        Env.set_atmospheric_model(
            type="Ensemble",
            file="monte_carlo_analysis_inputs/LASC2019_reanalysis.nc",
            dictionary="ECMWF",
        )


        # Set up parachutes. This rocket, named Valetudo, only has a drogue chute.
        def drogue_trigger(p, h, y):
            # Check if rocket is going down, i.e. if it has passed the apogee
            vertical_velocity = y[5]
            # Return true to activate parachute once the vertical velocity is negative
            return True if vertical_velocity < 0 else False


        # Iterate over flight settings
        #out = display("Starting", display_id=True)

        for setting in flight_settings(analysis_parameters, number_of_simulations):
            start_time = process_time()
            i += 1

            # Update environment object
            Env.select_ensemble_member(setting["ensemble_member"])

            # Create motor
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
            # Create rocket
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

            # Edit rocket drag
            Valetudo.power_off_drag *= setting["power_off_drag"]
            Valetudo.power_on_drag *= setting["power_on_drag"]
            # Add rocket nose, fins and tail
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
            # Add parachute
            Drogue = Valetudo.add_parachute(
                "Drogue",
                cd_s=setting["cd_s_drogue"],
                trigger=drogue_trigger,
                sampling_rate=105,
                lag=setting["lag_rec"] + setting["lag_se"],
                noise=(0, 8.3, 0.5),
            )

            # Run trajectory simulation
            try:
                test_flight = Flight(
                    rocket=Valetudo,
                    environment=Env,
                    rail_length=setting["rail_length"],
                    inclination=setting["inclination"],
                    heading=setting["heading"],
                    max_time=600,
                )
                export_flight_data(setting, test_flight, process_time() - start_time)
            except Exception as E:
                print(E)
                export_flight_error(setting)

            # Register time
            #out.update(
            #  f"Curent iteration: {i:06d} | Average Time per Iteration: {(process_time() - initial_cpu_time)/i:2.6f} s"
            #)

        # Done

        ## Print and save total time
        final_string = f"Completed {i} iterations successfully. Total CPU time: {process_time() - initial_cpu_time} s. Total wall time: {time() - initial_wall_time} s"
        #out.update(final_string)
        dispersion_input_file.write(final_string + "\n")
        dispersion_output_file.write(final_string + "\n")
        dispersion_error_file.write(final_string + "\n")

        ## Close files
        dispersion_input_file.close()
        dispersion_output_file.close()
        dispersion_error_file.close()

        filename = "monte_carlo_analysis_outputs/valetudo_rocket_v0"

        # Initialize variable to store all results
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
            "execution_time": [],
        }

        # Get all dispersion results
        # Get file
        dispersion_output_file = open(str(filename) + ".disp_outputs.txt", "r+")

        # Read each line of the file and convert to dict
        for line in dispersion_output_file:
            # Skip comments lines
            if line[0] != "{":
                continue
            # Eval results and store them
            flight_result = eval(line)
            dispersion_general_results.append(flight_result)
            for parameter_key, parameter_value in flight_result.items():
                dispersion_results[parameter_key].append(parameter_value)

        # Close data file
        dispersion_output_file.close()



        #streamlit rendering




















        # Print number of flights simulated
        N = len(dispersion_general_results)
        st.write("Number of simulations: ", N)

        st.write(
            f'Out of Rail Time - Mean Value: {np.mean(dispersion_results["out_of_rail_time"]):0.3f} s'
        )
        st.write(
            f'Out of Rail Time - Standard Deviation: {np.std(dispersion_results["out_of_rail_time"]):0.3f} s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["out_of_rail_time"], bins=int(N**0.5))
        ax.set_title("Out of Rail Time")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Out of Rail Velocity - Mean Value: {np.mean(dispersion_results["out_of_rail_velocity"]):0.3f} m/s'
        )
        st.write(
            f'Out of Rail Velocity - Standard Deviation: {np.std(dispersion_results["out_of_rail_velocity"]):0.3f} m/s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["out_of_rail_velocity"], bins=int(N**0.5))
        ax.set_title("Out of Rail Velocity")
        ax.set_xlabel("Velocity (m/s)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Apogee Time - Mean Value: {np.mean(dispersion_results["apogee_time"]):0.3f} s'
        )
        st.write(
            f'Apogee Time - Standard Deviation: {np.std(dispersion_results["apogee_time"]):0.3f} s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["apogee_time"], bins=int(N**0.5))
        ax.set_title("Apogee Time")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Apogee Altitude - Mean Value: {np.mean(dispersion_results["apogee_altitude"]):0.3f} m'
        )
        st.write(
            f'Apogee Altitude - Standard Deviation: {np.std(dispersion_results["apogee_altitude"]):0.3f} m'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["apogee_altitude"], bins=int(N**0.5))
        ax.set_title("Apogee Altitude")
        ax.set_xlabel("Altitude (m)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Apogee X Position - Mean Value: {np.mean(dispersion_results["apogee_x"]):0.3f} m'
        )
        st.write(
            f'Apogee X Position - Standard Deviation: {np.std(dispersion_results["apogee_x"]):0.3f} m'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["apogee_x"], bins=int(N**0.5))
        ax.set_title("Apogee X Position")
        ax.set_xlabel("Apogee X Position (m)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Apogee Y Position - Mean Value: {np.mean(dispersion_results["apogee_y"]):0.3f} m'
        )
        st.write(
            f'Apogee Y Position - Standard Deviation: {np.std(dispersion_results["apogee_y"]):0.3f} m'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["apogee_y"], bins=int(N**0.5))
        ax.set_title("Apogee Y Position")
        ax.set_xlabel("Apogee Y Position (m)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Impact Time - Mean Value: {np.mean(dispersion_results["impact_time"]):0.3f} s'
        )
        st.write(
            f'Impact Time - Standard Deviation: {np.std(dispersion_results["impact_time"]):0.3f} s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["impact_time"], bins=int(N**0.5))
        ax.set_title("Impact Time")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Impact X Position - Mean Value: {np.mean(dispersion_results["impact_x"]):0.3f} m'
        )
        st.write(
            f'Impact X Position - Standard Deviation: {np.std(dispersion_results["impact_x"]):0.3f} m'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["impact_x"], bins=int(N**0.5))
        ax.set_title("Impact X Position")
        ax.set_xlabel("Impact X Position (m)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Impact Y Position - Mean Value: {np.mean(dispersion_results["impact_y"]):0.3f} m'
        )
        st.write(
            f'Impact Y Position - Standard Deviation: {np.std(dispersion_results["impact_y"]):0.3f} m'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["impact_y"], bins=int(N**0.5))
        ax.set_title("Impact Y Position")
        ax.set_xlabel("Impact Y Position (m)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Impact Velocity - Mean Value: {np.mean(dispersion_results["impact_velocity"]):0.3f} m/s'
        )
        st.write(
            f'Impact Velocity - Standard Deviation: {np.std(dispersion_results["impact_velocity"]):0.3f} m/s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["impact_velocity"], bins=int(N**0.5))
        ax.set_title("Impact Velocity")
        ax.set_xlim(-35, 0)
        ax.set_xlabel("Velocity (m/s)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Initial Static Margin - Mean Value: {np.mean(dispersion_results["initial_static_margin"]):0.3f} c'
        )
        st.write(
            f'Initial Static Margin - Standard Deviation: {np.std(dispersion_results["initial_static_margin"]):0.3f} c'
        )

        st.write(
            f'Out of Rail Static Margin - Mean Value: {np.mean(dispersion_results["out_of_rail_static_margin"]):0.3f} c'
        )
        st.write(
            f'Out of Rail Static Margin - Standard Deviation: {np.std(dispersion_results["out_of_rail_static_margin"]):0.3f} c'
        )

        st.write(
            f'Final Static Margin - Mean Value: {np.mean(dispersion_results["final_static_margin"]):0.3f} c'
        )
        st.write(
            f'Final Static Margin - Standard Deviation: {np.std(dispersion_results["final_static_margin"]):0.3f} c'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["initial_static_margin"], label="Initial", bins=int(N**0.5))
        ax.hist(dispersion_results["out_of_rail_static_margin"], label="Out of Rail", bins=int(N**0.5))
        ax.hist(dispersion_results["final_static_margin"], label="Final", bins=int(N**0.5))
        ax.legend()
        ax.set_title("Static Margin")
        ax.set_xlabel("Static Margin (c)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Maximum Velocity - Mean Value: {np.mean(dispersion_results["max_velocity"]):0.3f} m/s'
        )
        st.write(
            f'Maximum Velocity - Standard Deviation: {np.std(dispersion_results["max_velocity"]):0.3f} m/s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["max_velocity"], bins=int(N**0.5))
        ax.set_title("Maximum Velocity")
        ax.set_xlabel("Velocity (m/s)")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["number_of_events"])
        ax.set_title("Parachute Events")
        ax.set_xlabel("Number of Parachute Events")
        ax.set_ylabel("Number of Occurrences")
        st.pyplot(fig)

        st.write(
            f'Drogue Parachute Trigger Time - Mean Value: {np.mean(dispersion_results["drogue_triggerTime"]):0.3f} s'
        )
        st.write(
            f'Drogue Parachute Trigger Time - Standard Deviation: {np.std(dispersion_results["drogue_triggerTime"]):0.3f} s'
        )









        fig, ax = plt.subplots()

        plt.figure()
        ax.hist(dispersion_results["drogue_triggerTime"], bins=int(N**0.5))
        ax.set_title("Drogue Parachute Trigger Time")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Number of Occurences")
        st.pyplot(fig)

        st.write(
            f'Drogue Parachute Fully Inflated Time -         Mean Value: {np.mean(dispersion_results["drogue_inflated_time"]):0.3f} s'
        )
        st.write(
            f'Drogue Parachute Fully Inflated Time - Standard Deviation: {np.std(dispersion_results["drogue_inflated_time"]):0.3f} s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["drogue_inflated_time"], bins=int(N**0.5))
        ax.set_title("Drogue Parachute Fully Inflated Time")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Number of Occurences")
        st.pyplot(fig)

        st.write(
            f'Drogue Parachute Fully Inflated Velocity -         Mean Value: {np.mean(dispersion_results["drogue_inflated_velocity"]):0.3f} m/s'
        )
        st.write(
            f'Drogue Parachute Fully Inflated Velocity - Standard Deviation: {np.std(dispersion_results["drogue_inflated_velocity"]):0.3f} m/s'
        )

        fig, ax = plt.subplots()
        ax.hist(dispersion_results["drogue_inflated_velocity"], bins=int(N**0.5))
        ax.set_title("Drogue Parachute Fully Inflated Velocity")
        ax.set_xlabel("Velocity m/s)")
        ax.set_ylabel("Number of Occurences")
        st.pyplot(fig)



        img = imread("monte_carlo_analysis_inputs/Valetudo_basemap_final.jpg")
        apogee_x = np.array(dispersion_results["apogee_x"])
        apogee_y = np.array(dispersion_results["apogee_y"])
        impact_x = np.array(dispersion_results["impact_x"])
        impact_y = np.array(dispersion_results["impact_y"])

        def eigsorted(cov):
            vals, vecs = np.lingalg.eigh(cov)
            order = vals.argsort()[::-1]
            return vals[order], vecs[:, order]
        
        #error ellipses
        impactCov = np.cov(impact_x, impact_y)
        impactVals, impactVecs = eigsorted(impactCov)
        impactTheta = np.degrees(np.arctan2(*impactVecs[:, 0][::-1]))
        impactW, impactH = 2 * np.sqrt(impactVals)

        impact_ellipses = []
        for j in [1, 2, 3]:
            impactEll = Ellipse(
                xy=(np.mean(impact_x), np.mean(impact_y)),
                width=impactW * j,
                height=impactH * j,
                angle=impactTheta,
                color="black",
            )
            impactEll.set_facecolor((0, 0, 1, 0.2))
            impact_ellipses.append(impactEll)
            ax.add_artist(impactEll)
    
        apogeeCov = np.cov(apogee_x, apogee_y)
        apogeeVals, apogeeVecs = eigsorted(apogeeCov)
        apogeeTheta = np.degrees(np.arctan2(*apogeeVecs[:, 0][::-1]))
        apogeeW, apogeeH = 2 * np.sqrt(apogeeVals)


        #drawing 
        for j in [1, 2, 3]:
            apogeeEll = Ellipse(
                xy=(np.mean(apogee_x), np.mean(apogee_y)),
                width=apogeeW * j,
                height=apogeeH * j,
                angle=apogeeTheta,
                color="black",
            )
            apogeeEll.set_facecolor((0, 1, 0, 0.2))
            ax.add_artist(apogeeEll)


        
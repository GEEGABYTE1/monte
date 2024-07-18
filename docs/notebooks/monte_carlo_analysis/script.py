from time import process_time, time 
from rocketpy import Environment, SolidMotor, Rocket, Flight 

import numpy as np
from numpy.random import normal, choice 
from IPython.display import display

import matplotlib as mpl
import matplotlib.pyplot as plt
import streamlit as st


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


if st.button("Run Simulation"): 
        analysis_parameters = {
            # Mass Details
            # Rocket's dry mass without motor (kg) and its uncertainty (standard deviation)
            "rocket_mass": (7.257, 0.001),
            # Rocket's inertia moment perpendicular to its axis (kg*m^2)
            "rocket_inertia_11": (3.675, 0.03675),
            # Rocket's inertia moment relative to its axis (kg*m^2)
            "rocket_inertia_33": (0.007, 0.00007),
            # Motors's dry mass without propellant (kg) and its uncertainty (standard deviation)
            "motor_dry_mass": (1.000, 0.001),
            # Motor's dry inertia moment perpendicular to its axis (kg*m^2)
            "motor_inertia_11": (1.675, 0.01675),
            # Motors's dry inertia moment relative to its axis (kg*m^2)
            "motor_inertia_33": (0.003, 0.00003),
            # Distance between rocket's center of dry mass and motor's center of dry mass (m)
            "motor_dry_mass_position": (0.5, 0.001),
            # Propulsion Details - run help(SolidMotor) for more information
            # Motor total impulse (N*s)
            "impulse": (1415.15, 35.3),
            # Motor burn out time (s)
            "burn_time": (5.274, 1),
            # Motor's nozzle radius (m)
            "nozzle_radius": (21.642 / 1000, 0.5 / 1000),
            # Motor's nozzle throat radius (m)
            "throat_radius": (8 / 1000, 0.5 / 1000),
            # Motor's grain separation (axial distance between two grains) (m)
            "grain_separation": (6 / 1000, 1 / 1000),
            # Motor's grain density (kg/m^3)
            "grain_density": (1707, 50),
            # Motor's grain outer radius (m)
            "grain_outer_radius": (21.4 / 1000, 0.375 / 1000),
            # Motor's grain inner radius (m)
            "grain_initial_inner_radius": (9.65 / 1000, 0.375 / 1000),
            # Motor's grain height (m)
            "grain_initial_height": (120 / 1000, 1 / 1000),
            # Aerodynamic Details - run help(Rocket) for more information
            # Rocket's radius (kg*m^2)
            "radius": (40.45 / 1000, 0.001),
            # Distance between rocket's center of dry mass and nozzle exit plane (m) (negative)
            "nozzle_position": (-1.024, 0.001),
            # Distance between rocket's center of dry mass and and center of propellant mass (m) (negative)
            "grains_center_of_mass_position": (-0.571, 0.001),
            # Multiplier for rocket's drag curve. Usually has a mean value of 1 and a uncertainty of 5% to 10%
            "power_off_drag": (0.9081 / 1.05, 0.033),
            # Multiplier for rocket's drag curve. Usually has a mean value of 1 and a uncertainty of 5% to 10%
            "power_on_drag": (0.9081 / 1.05, 0.033),
            # Rocket's nose cone length (m)
            "nose_length": (0.274, 0.001),
            # Axial distance between rocket's center of dry mass and nearest point in its nose cone (m)
            "nose_distance_to_CM": (1.134, 0.001),
            # Fin span (m)
            "fin_span": (0.077, 0.0005),
            # Fin root chord (m)
            "fin_root_chord": (0.058, 0.0005),
            # Fin tip chord (m)
            "fin_tip_chord": (0.018, 0.0005),
            # Axial distance between rocket's center of dry mass and nearest point in its fin (m)
            "fin_distance_to_CM": (-0.906, 0.001),
            # Launch and Environment Details - run help(Environment) and help(Flight) for more information
            # Launch rail inclination angle relative to the horizontal plane (degrees)
            "inclination": (84.7, 1),
            # Launch rail heading relative to north (degrees)
            "heading": (53, 2),
            # Launch rail length (m)
            "rail_length": (5.7, 0.0005),
            # Members of the ensemble forecast to be used
            "ensemble_member": list(range(10)),
            # Parachute Details - run help(Rocket) for more information
            # Drag coefficient times reference area for the drogue chute (m^2)
            "cd_s_drogue": (0.349 * 1.3, 0.07),
            # Time delay between parachute ejection signal is detected and parachute is inflated (s)
            "lag_rec": (1, 0.5),
            # Electronic Systems Details - run help(Rocket) for more information
            # Time delay between sensor signal is received and ejection signal is fired (s)
            "lag_se": (0.73, 0.16),
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
        print("Number of simulations: ", N)


        print(
            f'Out of Rail Time -         Mean Value: {np.mean(dispersion_results["out_of_rail_time"]):0.3f} s'
        )
        print(
            f'Out of Rail Time - Standard Deviation: {np.std(dispersion_results["out_of_rail_time"]):0.3f} s'
        )

        plt.figure()
        plt.hist(dispersion_results["out_of_rail_time"], bins=int(N**0.5))
        plt.title("Out of Rail Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Out of Rail Velocity -         Mean Value: {np.mean(dispersion_results["out_of_rail_velocity"]):0.3f} m/s'
        )
        print(
            f'Out of Rail Velocity - Standard Deviation: {np.std(dispersion_results["out_of_rail_velocity"]):0.3f} m/s'
        )

        plt.figure()
        plt.hist(dispersion_results["out_of_rail_velocity"], bins=int(N**0.5))
        plt.title("Out of Rail Velocity")
        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Apogee Time -         Mean Value: {np.mean(dispersion_results["apogee_time"]):0.3f} s'
        )
        print(
            f'Apogee Time - Standard Deviation: {np.std(dispersion_results["apogee_time"]):0.3f} s'
        )

        plt.figure()
        plt.hist(dispersion_results["apogee_time"], bins=int(N**0.5))
        plt.title("Apogee Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Apogee Altitude -         Mean Value: {np.mean(dispersion_results["apogee_altitude"]):0.3f} m'
        )
        print(
            f'Apogee Altitude - Standard Deviation: {np.std(dispersion_results["apogee_altitude"]):0.3f} m'
        )

        plt.figure()
        plt.hist(dispersion_results["apogee_altitude"], bins=int(N**0.5))
        plt.title("Apogee Altitude")
        plt.xlabel("Altitude (m)")
        plt.ylabel("Number of Occurences")
        plt.show()

        # Real measured apogee for Valetudo = 860 m

        print(
            f'Apogee X Position -         Mean Value: {np.mean(dispersion_results["apogee_x"]):0.3f} m'
        )
        print(
            f'Apogee X Position - Standard Deviation: {np.std(dispersion_results["apogee_x"]):0.3f} m'
        )

        plt.figure()
        plt.hist(dispersion_results["apogee_x"], bins=int(N**0.5))
        plt.title("Apogee X Position")
        plt.xlabel("Apogee X Position (m)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Apogee Y Position -         Mean Value: {np.mean(dispersion_results["apogee_y"]):0.3f} m'
        )
        print(
            f'Apogee Y Position - Standard Deviation: {np.std(dispersion_results["apogee_y"]):0.3f} m'
        )

        plt.figure()
        plt.hist(dispersion_results["apogee_y"], bins=int(N**0.5))
        plt.title("Apogee Y Position")
        plt.xlabel("Apogee Y Position (m)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Impact Time -         Mean Value: {np.mean(dispersion_results["impact_time"]):0.3f} s'
        )
        print(
            f'Impact Time - Standard Deviation: {np.std(dispersion_results["impact_time"]):0.3f} s'
        )

        plt.figure()
        plt.hist(dispersion_results["impact_time"], bins=int(N**0.5))
        plt.title("Impact Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Impact X Position -         Mean Value: {np.mean(dispersion_results["impact_x"]):0.3f} m'
        )
        print(
            f'Impact X Position - Standard Deviation: {np.std(dispersion_results["impact_x"]):0.3f} m'
        )

        plt.figure()
        plt.hist(dispersion_results["impact_x"], bins=int(N**0.5))
        plt.title("Impact X Position")
        plt.xlabel("Impact X Position (m)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Impact Y Position -         Mean Value: {np.mean(dispersion_results["impact_y"]):0.3f} m'
        )
        print(
            f'Impact Y Position - Standard Deviation: {np.std(dispersion_results["impact_y"]):0.3f} m'
        )

        plt.figure()
        plt.hist(dispersion_results["impact_y"], bins=int(N**0.5))
        plt.title("Impact Y Position")
        plt.xlabel("Impact Y Position (m)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Impact Velocity -         Mean Value: {np.mean(dispersion_results["impact_velocity"]):0.3f} m/s'
        )
        print(
            f'Impact Velocity - Standard Deviation: {np.std(dispersion_results["impact_velocity"]):0.3f} m/s'
        )

        plt.figure()
        plt.hist(dispersion_results["impact_velocity"], bins=int(N**0.5))
        plt.title("Impact Velocity")
        # plt.grid()
        plt.xlim(-35, 0)
        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Initial Static Margin -             Mean Value: {np.mean(dispersion_results["initial_static_margin"]):0.3f} c'
        )
        print(
            f'Initial Static Margin -     Standard Deviation: {np.std(dispersion_results["initial_static_margin"]):0.3f} c'
        )

        print(
            f'Out of Rail Static Margin -         Mean Value: {np.mean(dispersion_results["out_of_rail_static_margin"]):0.3f} c'
        )
        print(
            f'Out of Rail Static Margin - Standard Deviation: {np.std(dispersion_results["out_of_rail_static_margin"]):0.3f} c'
        )

        print(
            f'Final Static Margin -               Mean Value: {np.mean(dispersion_results["final_static_margin"]):0.3f} c'
        )
        print(
            f'Final Static Margin -       Standard Deviation: {np.std(dispersion_results["final_static_margin"]):0.3f} c'
        )

        plt.figure()
        plt.hist(dispersion_results["initial_static_margin"], label="Initial", bins=int(N**0.5))
        plt.hist(
            dispersion_results["out_of_rail_static_margin"],
            label="Out of Rail",
            bins=int(N**0.5),
        )
        plt.hist(dispersion_results["final_static_margin"], label="Final", bins=int(N**0.5))
        plt.legend()
        plt.title("Static Margin")
        plt.xlabel("Static Margin (c)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Maximum Velocity -         Mean Value: {np.mean(dispersion_results["max_velocity"]):0.3f} m/s'
        )
        print(
            f'Maximum Velocity - Standard Deviation: {np.std(dispersion_results["max_velocity"]):0.3f} m/s'
        )

        plt.figure()
        plt.hist(dispersion_results["max_velocity"], bins=int(N**0.5))
        plt.title("Maximum Velocity")
        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        plt.figure()
        plt.hist(dispersion_results["number_of_events"])
        plt.title("Parachute Events")
        plt.xlabel("Number of Parachute Events")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Drogue Parachute Trigger Time -         Mean Value: {np.mean(dispersion_results["drogue_triggerTime"]):0.3f} s'
        )
        print(
            f'Drogue Parachute Trigger Time - Standard Deviation: {np.std(dispersion_results["drogue_triggerTime"]):0.3f} s'
        )

        plt.figure()
        plt.hist(dispersion_results["drogue_triggerTime"], bins=int(N**0.5))
        plt.title("Drogue Parachute Trigger Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Drogue Parachute Fully Inflated Time -         Mean Value: {np.mean(dispersion_results["drogue_inflated_time"]):0.3f} s'
        )
        print(
            f'Drogue Parachute Fully Inflated Time - Standard Deviation: {np.std(dispersion_results["drogue_inflated_time"]):0.3f} s'
        )

        plt.figure()
        plt.hist(dispersion_results["drogue_inflated_time"], bins=int(N**0.5))
        plt.title("Drogue Parachute Fully Inflated Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Number of Occurences")
        plt.show()

        print(
            f'Drogue Parachute Fully Inflated Velocity -         Mean Value: {np.mean(dispersion_results["drogue_inflated_velocity"]):0.3f} m/s'
        )
        print(
            f'Drogue Parachute Fully Inflated Velocity - Standard Deviation: {np.std(dispersion_results["drogue_inflated_velocity"]):0.3f} m/s'
        )

        plt.figure()
        plt.hist(dispersion_results["drogue_inflated_velocity"], bins=int(N**0.5))
        plt.title("Drogue Parachute Fully Inflated Velocity")
        plt.xlabel("Velocity m/s)")
        plt.ylabel("Number of Occurences")
        plt.show()
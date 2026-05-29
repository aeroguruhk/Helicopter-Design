# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 17:10:09 2026
Modifed on 24 May 2026with following fixes
1 fix unit consistency
2 add tail rotor power
3 improve induced power model
4 add tip Mach limits
@author: hknar
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import io
import base64
class HelicopterSweep:
    def __init__(self, m=None, passengers=None, R=None, main_blades=4, tail_blades=4):
        """
        Accepts single values or lists/arrays for m, passengers, or R.
        """
        self.m_values = np.atleast_1d(m) if m is not None else None
        self.passenger_values = np.atleast_1d(passengers) if passengers is not None else None
        self.R_values = np.atleast_1d(R) if R is not None else None
        self.main_blades = main_blades
        self.tail_blades = tail_blades

        # Storage for results
        self.results = []

        # Run calculations for each input
        self._compute_sweep()

    def _compute_single(self, m=0, passengers=0, R=0):
        """Compute helicopter parameters for a single input set."""
        heli = {}
        heli["MRotor"] = {"NBlades": self.main_blades}
        heli["TRotor"] = {"NBlades": self.tail_blades}
        heli["dimensions"] = {}

        # Case 1: mass given
        if m != 0:
            passengers = max(1, int(np.floor(np.log(m/1525)/0.0809)))
            R = 0.4885 * (m**0.308)
            MRarea= np.pi*R**2
            MDL = m*9.81/MRarea
            TR = 0.0886*0.5*m**0.393 # Tail rotor radius
            TRarea= np.pi*TR**2
            TDL= 1200.0
            heli["mass"], heli["passengers"] =m,passengers
            heli["MRotor"]["radius"],heli["MRotor"]["area"],heli["MRotor"]["DL"] = R, MRarea,MDL
            heli["TRotor"]["radius"],heli["TRotor"]["area"], heli["TRotor"]["TDL"] = TR, TRarea,TDL
        # Case 2: passengers given
        elif passengers != 0:
            m = 1525 * np.exp(0.0809 * passengers)
            R = 0.4885 * (m**0.308)
            MRarea= np.pi*R**2
            TR = 0.0886*0.5*m**0.393 #Tail rotor radius
            TRarea= np.pi*TR**2
            heli["mass"], heli["passengers"] =m,passengers
            heli["MRotor"]["radius"], heli["MRotor"]["area"] = R, MRarea
            heli["TRotor"]["radius"],heli["TRotor"]["area"] = TR, TRarea
        # Case 3: rotor radius given
        elif R != 0:
            m = (R / 0.4885)**(1/0.308)
            MRarea= np.pi*R**2
            TR = 0.0886*0.5*m**0.393 # Tail rotor radius
            TRarea= np.pi*TR**2
            passengers = max(1, int(np.floor(np.log(m/1525)/0.0809)))
            heli["mass"], heli["passengers"] =m,passengers
            heli["MRotor"]["radius"], heli["MRotor"]["area"] = R, MRarea
            heli["TRotor"]["radius"],heli["TRotor"]["area"] = TR, TRarea
        # General dimensions
         # General dimensions
        heli["dimensions"]["height"] = 0.642 * (2*R)**0.677
        heli["dimensions"]["length"] = 0.824 * (2*R)**1.056
        heli["dimensions"]["tiptotip"] = 1.09 * (2*R)**1.03
        heli["dimensions"]["width"] = 0.436 * (2*R)**0.697
        heli["dimensions"]["equivalentFlatPlateArea"] = (
            891.45 + np.sqrt(794683.1025 - 4412*(331.23 - m))) / 2206
            
        # Empirical geometric dimension estimates
        heli["dimensions"]["length"] = 2.0 * R
        heli["dimensions"]["width"] = 0.18 * (2.0 * R)
        heli["dimensions"]["height"] = 0.25 * (2.0 * R)
        heli["dimensions"]["tip_to_tip"] = 2.2 * R + TR

        # Main rotor
        heli["MRotor"]["angular_velocity"] = 280 / (2*R)**0.829 # radians/sec
        chord = 0.0108 * (m**0.539) / (heli["MRotor"]["NBlades"]**0.714)
        heli["MRotor"]["chord"] =chord
        heli["MRotor"]["blade"] = {"chord": chord, "AspectRatio": R/chord}
        heli["MRotor"]["solidity"] = chord * heli["MRotor"]["NBlades"] / (R*np.pi)
        
        heli["MRotor"]["mass"] = 0.45 * (m * MRarea)**0.5
        ######################################################################
        """
        Rad/sec = RPM*0.10472
        
        """
        # Tail rotor
        heli["TRotor"]["angular_velocity"] = 364 / (2*TR)**0.828 #radians/sec
        chord = 0.0058 * (m**0.506) / (heli["TRotor"]["NBlades"]**0.714)
        heli["TRotor"]["chord"] =chord
        heli["TRotor"]["blade"] = {"chord": chord, "AspectRatio": TR/chord}
        heli["TRotor"]["solidity"] = chord * heli["TRotor"]["NBlades"] / (TR*np.pi)
        heli["TRotor"]["mass"] = 0.025*m
        return heli
        #############################################################################
    def compute_main_rotor_mass(m, R, Nb):
        """
        SI-consistent conceptual main rotor mass model
        m  : helicopter mass (kg)
        R  : rotor radius (m)
        Nb : number of blades
        """

        A = np.pi * R**2

        # Base scaling (physics-based)
        Wmr = 0.45 * (m * A)**0.5

        # Blade count correction (weak effect)
        blade_factor = (Nb / 4.0)**0.3

        return Wmr * blade_factor
        ########################################################################
    
     ##############################################################################
    def _compute_sweep(self):
        """Compute results for all input values."""
        if self.m_values is not None:
            for m in self.m_values:
                self.results.append(self._compute_single(m=m))
        elif self.passenger_values is not None:
            for p in self.passenger_values:
                self.results.append(self._compute_single(passengers=p))
        elif self.R_values is not None:
            for r in self.R_values:
                self.results.append(self._compute_single(R=r))
    def compute_range(self, V_range, SFC=1.8e-7, fuel_fraction=0.3, rho=1.225, g=9.81, Cd0=0.011):
        for res in self.results:
            m = res["mass"]
            R = res["MRotor"]["radius"]
            A = np.pi * R**2
            Omega = res["MRotor"]["angular_velocity"]
            sigma = res["MRotor"]["solidity"]
            f = res["dimensions"]["equivalentFlatPlateArea"]
    
            W = m * g
            fuel_mass = fuel_fraction * m
    
            ranges = []   # ✅ moved inside loop
            endurances = []

            for V in V_range:
               # Pi = (W**1.5) / np.sqrt(2 * rho * A)
                Pi_hover = (W**1.5) / np.sqrt(2 * rho * A)
                mu = V / (Omega * R)
                Pi = Pi_hover / np.sqrt(1 + 4.65*mu**2)
                P0 = (sigma * Cd0 / 8) * rho * A * (Omega * R)**3
                Pp = 0.5 * rho * V**3 * f
                Ptr = 0.1 * (Pi + P0)
                Ptot = Pi + P0 + Pp +Ptr
                eta_trans =0.88
                Pengine = Ptot / eta_trans
                mdot_f = SFC * Pengine
                endurance = fuel_mass / mdot_f
                endurances.append(endurance / 3600)
    
                range_m = V * endurance
                ranges.append(range_m / 1000)

            res["Range"] = {   # ✅ now correctly assigned per helicopter
                "V": V_range,
                "Range_km": ranges,
                "Endurance_hours": endurances
            }

    def plot_range_vs_velocity(self, return_base64=False):
        plt.figure(figsize=(9,6))

        for res in self.results:
            V = res["Range"]["V"]
            R = res["Range"]["Range_km"]

            plt.plot(V, R, label=f"m={res['mass']:.0f} kg")

            # Mark max range point
            idx = np.argmax(R)
            plt.scatter(V[idx], R[idx])
            plt.text(V[idx], R[idx], f"Max @ {V[idx]:.1f} m/s", fontsize=8)

        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Range (km)")
        plt.title("Range vs Velocity")
        plt.legend()
        plt.grid(True)
        
        if return_base64:
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.getvalue()).decode('utf-8')
        else:
            plt.show()

    def compute_power_required(self, rho=1.225, g=9.81, Cd0=0.011, V_cruise=0, V_max=0):
        """Compute power required for all results."""
        for res in self.results:
            m = res["mass"]
            R = res["MRotor"]["radius"]
            A = np.pi * R**2
            Omega = res["MRotor"]["angular_velocity"]
            sigma = res["MRotor"]["solidity"]
            f = res["dimensions"]["equivalentFlatPlateArea"]

            W = m * g

            def calc_power(V):
                Pi_hover = (W**1.5) / np.sqrt(2 * rho * A)
                mu = V / (Omega * R)
                Pi = Pi_hover / np.sqrt(1 + 4.65*mu**2)
                P0 = (sigma * Cd0 / 8) * rho * A * (Omega * R)**3
                Pp = 0.5 * rho * V**3 * f
                Ptr = 0.1*(Pi+P0)
                return Pi, P0, Pp, Pi + P0 + Pp +Ptr

            # Evaluate power at hover, cruise, and max velocity
            powers = [calc_power(0), calc_power(V_cruise), calc_power(V_max)]
            # Find the condition with the maximum total power
            Pi, P0, Pp, Ptot = max(powers, key=lambda x: x[3])

            res["Power"] = {"Induced": Pi, "Profile": P0, "Parasite": Pp, "Total": Ptot}

    
    def plot_power(self, x_param="mass", components=("Induced","Profile","Parasite","Total")):
        """
        Plot power required vs input parameter for multiple components.
        components: tuple of power components to plot (default: all)
        """
        x_vals = [res[x_param] for res in self.results]

        plt.figure(figsize=(8,6))
        for comp in components:
            y_vals = [res["Power"][comp]/1000 for res in self.results]  # convert W to kW
            plt.plot(x_vals, y_vals, marker="o", label=comp)

        plt.xlabel(x_param)
        plt.ylabel("Power (kW)")
        plt.title(f"Power Components vs {x_param}")
        plt.legend()
        plt.grid(True)
        plt.show()
        ##################################################
    
    def plot_max_range_vs_mass(self):
        masses = []
        max_ranges = []
        opt_speeds = []

        for res in self.results:
            if "Range" not in res:
                raise ValueError("Run compute_range() before plotting.")

            V = np.array(res["Range"]["V"])
            R = np.array(res["Range"]["Range_km"])

            idx = np.argmax(R)

            masses.append(res["mass"])
            max_ranges.append(R[idx])
            opt_speeds.append(V[idx])

        plt.figure(figsize=(8,6))
        plt.plot(masses, max_ranges, marker='o')

        # Annotate optimum speeds
        for m, r, v in zip(masses, max_ranges, opt_speeds):
            plt.text(m, r, f"{v:.1f} m/s", fontsize=8)

        plt.xlabel("Mass (kg)")
        plt.ylabel("Maximum Range (km)")
        plt.title("Maximum Range vs Mass")
        plt.grid(True)
        plt.show()
    ####################################################
    def compute_optimum_speed(self, rho=1.225, g=9.81, Cd0=0.011):
        results_opt = []
        """
        Using our model:
        Induced power ≈ constant Pi
        Profile power ≈ constant
        Parasite power is proportional to V³
        so P total = A+ BV3, for max range minimize P/V
        V opt = (A/2B)raised to (1/3), here A = Pi+PO, B = 1/2*rho*f
        """
        for res in self.results:
            m = res["mass"]
            R = res["MRotor"]["radius"]
            A_disk = np.pi * R**2
            Omega = res["MRotor"]["angular_velocity"]
            sigma = res["MRotor"]["solidity"]
            f = res["dimensions"]["equivalentFlatPlateArea"]

            W = m * g

            # Constant terms
            #Pi = (W**1.5) / np.sqrt(2 * rho * A_disk)
            Pi_hover = (W**1.5) / np.sqrt(2 * rho * A_disk)
          # mu = V / (Omega * R)
           # Pi = Pi_hover / np.sqrt(1 + 4.65*mu**2)
            P0 = (sigma * Cd0 / 8) * rho * A_disk * (Omega * R)**3

            A_const = Pi_hover + P0
            B_const = 0.5 * rho * f

            V_opt = (A_const / (2 * B_const))**(1/3)

            results_opt.append({
                "mass": m,
                "V_opt": V_opt
            })

            # Store in object too
            res["V_opt"] = V_opt

        return results_opt
    #############################################################
    
    def export_to_csv(self, filename="helicopter_output.csv"):
        """Export computed results to CSV with extended rotor parameters."""

        fieldnames = [
            "mass",
            "passengers",

            # Main Rotor
            "MR_radius",
            "MR_area",
            "MR_disk_loading",
            "MR_solidity",
            "MR_blades",
            "MR_chord",
            "MR_mass",
            "MR_aspect_ratio",

            # Tail Rotor
            "TR_radius",
            "TR_area",
            "TR_disk_loading",
            "TR_solidity",
            "TR_blades",
            "TR_chord",
            "TR_aspect_ratio",
            "Emp_Weight_Total",
            "Emp_Weight_Total_kg"
        ]

        with open(filename, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()

            for heli in self.results:

                MR = heli.get("MRotor", {})
                TR = heli.get("TRotor", {})
                EW = heli.get("Weights_empirical", {})
                row = {
                    "mass": heli.get("mass"),
                    "passengers": heli.get("passengers"),

                    # Main Rotor
                    "MR_radius": MR.get("radius"),
                    "MR_area": MR.get("area"),
                    "MR_disk_loading": MR.get("MDL"),
                    "MR_solidity": MR.get("solidity"),
                    "MR_blades": MR.get("NBlades"),
                    "MR_chord": MR.get("chord"),
                    "MR_mass": MR.get("mass"),
                    "MR_aspect_ratio": MR.get("blade", {}).get("AspectRatio"),

                    # Tail Rotor
                    "TR_radius": TR.get("radius"),
                    "TR_area": TR.get("area"),
                    "TR_disk_loading": TR.get("TDL"),
                    "TR_solidity": TR.get("solidity"),
                    "TR_blades": TR.get("NBlades"),
                    "TR_chord": TR.get("chord"),
                    "TR_aspect_ratio": TR.get("blade", {}).get("AspectRatio"),
                    "Emp_Weight_Total": EW.get("Total"),
                    "Emp_Weight_Total_kg": EW.get("Total_kg"),
                }

                writer.writerow(row)

        print(f"CSV file '{filename}' written successfully.")
    ###########################################################
    def compute_empirical_weights(self, g=9.81):
        """
        Compute subsystem masses (kg) using 1983 empirical correlations (SI form).
        Adds a new dictionary: res["Weights_empirical"]
        """

        for res in self.results:

            m = res["mass"]                 # kg

            # --- Areas ---
            S_main = res["MRotor"]["area"]     # m²
            S_tail = res["TRotor"]["area"]     # m²
            S_body = res["dimensions"]["equivalentFlatPlateArea"]  # proxy

            # -------------------------------
            # 1. Main Rotor mass, M1 = Mblades + Mhub + Mcontrols
            M1 = 0.45 * (m * S_main)**0.5

            # -------------------------------
            # 2. Tail Rotor mass
            M2A = 0.025 * m
            """
            # Tail structure (approx form, constant can be tuned)
            M2B = 85.1 * np.log(S_tail) - 30

            # -------------------------------
            # 3. Body mass
            M3 = 3.79 * (S_body ** 1.917)

            # -------------------------------
            """
            # -------------------------------
            # 6. Propulsion

            # --- Engine (ASSUMPTION: proportional to power) ---
            # --- Engine (corrected SI form) ---
            if "Power" in res:
                P = res["Power"]["Total"]   # Watts
            else:
                P = 0

            """
            P_total in Watts
            returns kg
            """
            P_kW = P / 1000
            M6A = 0.30 * P_kW  

            # --- Drive system ---
            M6B = 76.5 * np.exp(7.43e-5 * S_main)

            # --- Fuel system ---
            fuel_mass = 0.3 * m   # consistent with fuel fraction model
            M6C = 0.10 * fuel_mass

            # -------------------------------
            # Total Empty Mass (kg)
            M_total = M1 + M2A + M6A + M6B + M6C

            res["Weights_empirical"] = {
                "MainRotor": M1,
                "TailRotor": M2A,
                "Power": P_kW,
                "Engine": M6A,
                "Drive_system": M6B,
                "FuelSystem": M6C,
                "Total": M_total,
                "Total_kg": M_total
            }
            print(res)
            #########################################
          
            
    def plot_optimum_speed_vs_mass(self):
        masses = []
        speeds = []

        for res in self.results:
            if "V_opt" not in res:
                raise ValueError("Run compute_optimum_speed() first.")

            masses.append(res["mass"])
            speeds.append(res["V_opt"])

        plt.figure(figsize=(8,6))
        plt.plot(masses, speeds, marker='o')

        plt.xlabel("Mass (kg)")
        plt.ylabel("Optimum Cruise Speed (m/s)")
        plt.title("Optimum Speed for Maximum Range")
        plt.grid(True)
        plt.show()
    ##################################################
    
    def plot_endurance_vs_velocity(self, V_range, SFC=1.8e-7, fuel_fraction=0.3,
                               rho=1.225, g=9.81, Cd0=0.011):

        plt.figure(figsize=(9,6))

        for res in self.results:
            m = res["mass"]
            R = res["MRotor"]["radius"]
            A = np.pi * R**2
            Omega = res["MRotor"]["angular_velocity"]
            sigma = res["MRotor"]["solidity"]
            f = res["dimensions"]["equivalentFlatPlateArea"]

            W = m * g
            fuel_mass = fuel_fraction * m

            endurance_vals = []

            for V in V_range:
                # Hover induced power
                Pi_hover = (W**1.5) / np.sqrt(2 * rho * A)

                # Advance ratio
                mu = V / (Omega * R)

                # Forward-flight induced power correction
                Pi = Pi_hover / np.sqrt(1 + 4.65*mu**2)
               # Pi = (W**1.5) / np.sqrt(2 * rho * A)
                P0 = (sigma * Cd0 / 8) * rho * A * (Omega * R)**3
                Pp = 0.5 * rho * V**3 * f
                Ptr = 0.1 * (Pi + P0)
                Ptot = Pi + P0 + Pp +Ptr
                eta_trans = 0.88
                Pengine = Ptot / eta_trans
                mdot_f = SFC * Pengine
                endurance = fuel_mass / mdot_f  # seconds

                endurance_vals.append(endurance / 3600)  # hours

            plt.plot(V_range, endurance_vals, label=f"m={m:.0f} kg")

            # Max endurance point
            idx = np.argmax(endurance_vals)
            plt.scatter(V_range[idx], endurance_vals[idx])
            plt.text(V_range[idx], endurance_vals[idx],
                     f"Max @ {V_range[idx]:.1f} m/s", fontsize=8)

        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Endurance (hours)")
        plt.title("Endurance vs Velocity (Loiter Condition)")
        plt.legend()
        plt.grid(True)
        plt.show()
    #######################################################
    def compute_fuel_fraction_for_range(self, target_range_km, V_cruise,
                                     SFC=1.8e-7, rho=1.225, g=9.81, Cd0=0.011):

        results_fuel = []

        for res in self.results:
            m = res["mass"]
            R = res["MRotor"]["radius"]
            A = np.pi * R**2
            Omega = res["MRotor"]["angular_velocity"]
            sigma = res["MRotor"]["solidity"]
            f = res["dimensions"]["equivalentFlatPlateArea"]

            W = m * g

            # Power at cruise
            Pi = (W**1.5) / np.sqrt(2 * rho * A)
            #Pi_hover = (W**1.5) / np.sqrt(2 * rho * A)
            #mu = V / (Omega * R)
            #Pi = Pi_hover / np.sqrt(1 + 4.65*mu**2)
            P0 = (sigma * Cd0 / 8) * rho * A * (Omega * R)**3
            Pp = 0.5 * rho * V_cruise**3 * f
            Ptr = 0.1 * (Pi + P0)
            Ptot = Pi + P0 + Pp +Ptr
            eta_trans = 0.88
            Pengine = Ptot / eta_trans
            mdot_f = SFC * Pengine

            # Required endurance
            target_range_m = target_range_km * 1000
            endurance_req = target_range_m / V_cruise

            fuel_mass = mdot_f * endurance_req
            fuel_fraction = fuel_mass / m

            res["fuel_fraction_required"] = fuel_fraction

            results_fuel.append({
                "mass": m,
                "fuel_fraction": fuel_fraction
            })

        return results_fuel
    ##################################################
    def plot_fuel_fraction_vs_mass(self):
        masses = []
        fuel_fracs = []

        for res in self.results:
            if "fuel_fraction_required" not in res:
                raise ValueError("Run compute_fuel_fraction_for_range() first.")

            masses.append(res["mass"])
            fuel_fracs.append(res["fuel_fraction_required"])

        plt.figure(figsize=(8,6))
        plt.plot(masses, fuel_fracs, marker='o')

        plt.xlabel("Mass (kg)")
        plt.ylabel("Fuel Fraction")
        plt.title("Fuel Fraction Required vs Mass")
        plt.grid(True)
        plt.show()
    #############################################
    def plot_constraint_diagram(self, rho=1.225, g=9.81, Cd0=0.011):

        DL_vals = []
        PL_vals = []

        for res in self.results:
            m = res["mass"]
            R = res["MRotor"]["radius"]
            A = np.pi * R**2

            Omega = res["MRotor"]["angular_velocity"]
            sigma = res["MRotor"]["solidity"]
            f = res["dimensions"]["equivalentFlatPlateArea"]

            W = m * g

            # Hover power (critical constraint)
            Pi = (W**1.5) / np.sqrt(2 * rho * A)
           # Pi_hover = (W**1.5) / np.sqrt(2 * rho * A)
            #mu = V / (Omega * R)
            #Pi = Pi_hover / np.sqrt(1 + 4.65*mu**2)
            P0 = (sigma * Cd0 / 8) * rho * A * (Omega * R)**3

            P_hover = Pi + P0

            DL = W / A
            PL = W / P_hover

            DL_vals.append(DL)
            PL_vals.append(PL)

        plt.figure(figsize=(8,6))
        plt.plot(DL_vals, PL_vals, marker='o')

        for i, res in enumerate(self.results):
            plt.text(DL_vals[i], PL_vals[i],
                     f"{res['mass']:.0f} kg", fontsize=8)

        plt.xlabel("Disk Loading (N/m²)")
        plt.ylabel("Power Loading (N/W)")
        plt.title("Constraint Diagram (Hover)")
        plt.grid(True)
        plt.show()
    ######################################
    def plot_power_vs_velocity(self, V_range, rho=1.225, g=9.81, Cd0=0.011,
                               reference_data=None, return_base64=False):
        """
        Plot power required vs forward velocity for each helicopter in the sweep.
        Optionally overlay reference data (dict of name: (V_array, P_array)).
        """
        plt.figure(figsize=(9,6))

        for res in self.results:
            m = res["mass"]
            R = res["MRotor"]["radius"]
            A = np.pi * R**2
            Omega = res["MRotor"]["angular_velocity"]
            sigma = res["MRotor"]["solidity"]
            f = res["dimensions"]["equivalentFlatPlateArea"]

            Pi_vals, P0_vals, Pp_vals, Ptot_vals = [], [], [], []

            for V in V_range:
                W = m * g
                #Pi = (W**1.5) / np.sqrt(2 * rho * A)
                Pi_hover = (W**1.5) / np.sqrt(2 * rho * A)
                mu = V / (Omega * R)
                Pi = Pi_hover / np.sqrt(1 + 4.65*mu**2)
                P0 = (sigma * Cd0 / 8) * rho * A * (Omega * R)**3
                Pp = 0.5 * rho * V**3 * f
                Ptr = 0.1 * (Pi + P0)
                Ptot = Pi + P0 + Pp +Ptr
                eta_trans = 0.88
                Pengine = Ptot / eta_trans
                Pi_vals.append(Pi/1000)
                P0_vals.append(P0/1000)
                Pp_vals.append(Pp/1000)
                Ptot_vals.append(Pengine/1000)

            # Plot total curve for this helicopter
            plt.plot(V_range, Ptot_vals, label=f"Total (m={m:.0f} kg)")
            # Mark minimum total power speed
            min_idx = np.argmin(Ptot_vals)
            plt.scatter(V_range[min_idx], Ptot_vals[min_idx], color='red', marker='x')
            plt.text(V_range[min_idx], Ptot_vals[min_idx],
                     f"Min @ {V_range[min_idx]:.1f} m/s", fontsize=8)

        # Overlay reference data if provided
        if reference_data:
            for name, (V_ref, P_ref) in reference_data.items():
                plt.plot(V_ref, P_ref, 'o--', label=f"{name} Reference")

        plt.xlabel("Forward Speed V (m/s)")
        plt.ylabel("Power (kW)")
        plt.title("Power Required vs Forward Speed")
        plt.legend()
        plt.grid(True)
        
        if return_base64:
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return base64.b64encode(buf.getvalue()).decode('utf-8')
        else:
            plt.show()

if __name__ == "__main__":
    # Sweep over mass values
    heli_sweep = HelicopterSweep(m=[6000, 8329,14969])
    heli_sweep.compute_power_required()
    V_range = np.linspace(0,100, 10)
    
    heli_sweep.compute_power_required()
    heli_sweep.compute_empirical_weights()
    
    for res in heli_sweep.results:
        print({
            "mass": res["mass"],
            "R": res["MRotor"]["radius"],
            "P_hover_kW": res["Power"]["Total"]/1000,
            "MR_mass": res["MRotor"]["mass"],
            "TR_mass": res["TRotor"].get("mass", None),
            "Engine_mass": res["Weights_empirical"]["Engine"],
            "Empty_mass": res["Weights_empirical"]["Total_kg"]
        })
    # NEW: Export
    heli_sweep.export_to_csv("output.csv")
    heli_sweep.compute_range(V_range)
    
    # Reference data (simplified, illustrative)
    V_ref = [0, 40, 80]
    P_huey = [900, 1100, 1500]   # UH-1 Huey published values (kW)
    P_blackhawk = [1600, 1800, 2200]  # UH-60 Black Hawk published values (kW)
    
    reference_data = {
        "UH-1 Huey": (V_ref, P_huey),
        "UH-60 Black Hawk": (V_ref, P_blackhawk)
    }
    
    # Plot with overlay
    
    V_range = np.linspace(0, 100, 10)
    
    heli_sweep.compute_range(V_range)   # MUST come first
    heli_sweep.plot_range_vs_velocity()
    heli_sweep.plot_max_range_vs_mass()
    
    heli_sweep.compute_optimum_speed()
    heli_sweep.plot_optimum_speed_vs_mass()
    
    heli_sweep.plot_power_vs_velocity(V_range, reference_data=reference_data)
    heli_sweep.compute_range(V_range)
    heli_sweep.plot_endurance_vs_velocity(V_range)
    
    heli_sweep.compute_fuel_fraction_for_range(target_range_km=500, V_cruise=60)
    heli_sweep.plot_fuel_fraction_vs_mass()
    
    heli_sweep.plot_constraint_diagram()

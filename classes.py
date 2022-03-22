import os
import re
import numpy as np
import math


class Config:
    def __init__(self, locals_config):
        """
        Contains the configuration variables of the run analysis
        Args:
            locals_config (dict): locals() variable set where values are extracted from
        """
        # Energy angular distributions
        self.plot_EADS = self.get_var_from_locals(locals_config, "plot_EADS")
        self.plot_EDFs = self.get_var_from_locals(locals_config, "plot_EDFs")
        self.plot_ADFs = self.get_var_from_locals(locals_config, "plot_ADFs")
        self.plot_EADS_species = self.get_var_from_locals(locals_config, "plot_EADS_species")
        self.IEAD_max_energy = self.get_var_from_locals(locals_config, "IEAD_max_energy")
        self.EEAD_max_energy = self.get_var_from_locals(locals_config, "EEAD_max_energy")

        # electric field XT plots
        self.plot_XT_Efield = self.get_var_from_locals(locals_config, "plot_XT_Efield")
        self.XT_colorbar = self.get_var_from_locals(locals_config, "XT_colorbar")
        self.XT_lower_half = self.get_var_from_locals(locals_config, "XT_lower_half")

        # probe potential over time
        self.plot_local_potential_over_time = self.get_var_from_locals(locals_config, "plot_local_potential_over_time")
        self.potential_over_time_locations = self.get_var_from_locals(locals_config, "potential_over_time_locations")
        self.potential_over_time_labels = self.get_var_from_locals(locals_config, "potential_over_time_labels")

        # config path
        self.dir_root = self.get_var_from_locals(locals_config, "dir_root")

    @staticmethod
    def get_var_from_locals(locals_config, variable):
        """
        fetch and return variable from locals() dictionary
        Args:
            locals_config (dict): dictionary containing config variables ( locals() )
            variable (str): name of variable to be fetched
        """
        if variable in locals_config:
            return locals_config[variable]
        else:
            return None

    def print_config(self):
        """
        prints information regarding the config object
        """
        print("Configuration:  ")
        for attribute in self.__dict__:
            print(f"{attribute}  = {self.__dict__[attribute]}")


class Case:
    def __init__(self, path):
        # file information
        self.path = path
        self.name = os.path.basename(path)
        self.path_nam = self.find_file_by_name(".nam", [])
        self.path_log = self.find_file_by_name(".log", ["runscript"])
        self.path_out = self.find_file_by_name(".out", ["gmon", "mcs.", "mesh."])
        self.path_movie1 = self.find_file_by_name("movie1.pdt")
        self.path_mesh = self.find_file_by_name("mesh.dat")
        self.path_pcmc = self.find_file_by_name("pcmc.prof")
        self.mesh = self.read_mesh_file()
        # pcmc data
        self.pcmc_species = []                                              # names of pcmc species
        self.pcmc_eads = []                                                 # 2D energy-angle distributions
        self.pcmc_edfs = []                                                 # 1D energy distributions
        self.pcmc_adfs = []                                                 # 1D angle distributions
        self.pcmc_max_angle = (0, 0)                                        # max pcmc angle tuple: (ions, neutrals)
        self.pcmc_max_energy = (0, 0)                                       # max pcmc energy tuple: (ions, neutrals)
        self.pcmc_mean_energy = []                                          # mean particle energies
        self.pcmc_mode_energy = []                                          # mode value particle energies
        # power and voltage setup
        self.irfpow = self.find_nam_parameter("IRFPOW", expected_type="int")         # adjust voltages for power target
        self.powerICP = 0
        self.powerCCP1 = 0
        self.powerCCP2 = 0
        self.rfpnorma = self.find_nam_parameter("RFPNORMA", expected_type='list_float')        # target powers
        self.icustom = self.find_nam_parameter("ICUSTOM", expected_type='list_int')            # use of custom waveforms
        self.contains_custom = self.icustom != [0] * len(self.icustom)                         # custom waveforms flag
        self.custom_phase = self.find_nam_parameter("CUSTOM_PHASE", expected_type='list_float')     # phase of harmonics
        self.custom_relharm = self.find_nam_parameter("CUSTOM_RELHARM", expected_type='list_int')  # orders of harmonics
        self.custom_relamp = self.find_nam_parameter("CUSTOM_RELAMP", expected_type='list_float')  # amplitude of harm.
        if self.custom_phase:
            self.cwaveform_phase = self.custom_phase[1]                         # main relative phase of harmonics
        else:
            self.cwaveform_phase = None
        # material data
        self.cwafer = self.find_nam_parameter("CWAFER", expected_type='list_str')         # wafer materials
        self.metal_labels = self.find_nam_parameter("CMETAL", expected_type='list_str')   # metal material labels
        # quantities after convergence
        self.final_voltages = self.get_final_voltages()                          # voltage amplitudes on last iteration
        self.dc_bias = eval(self.find_out_parameter("DC BIAS"))                  # DC self-bias on last iteration
        self.ne_ave = eval(self.find_out_parameter("AVERAGE ELECTRON DENSITY"))  # average n_e on last iteration
        self.restart = self.find_nam_parameter("IRESTART", expected_type='int')           # IRESTART .nam-value
        self.rffac = self.find_nam_parameter("RFFAC", expected_type='int')                # RFFAC .nam-value
        self.freq = self.find_nam_parameter("FREQ", expected_type='float')                # FREQ .nam-value
        # properties of movie file
        self.movie_I, self.movie_J, self.movie_num_zones = self.movie_find_dimensions()
        self.potential_xt = []                  # XT data of potential
        self.potential_over_time = []           # nested list of time varying potential and different locations
        self.potential_over_time_location = []  # sample locations of time varying potentials
        self.potential_over_time_labels = []    # sample label of time varying potentials

    def find_file_by_name(self, str_match, str_exclude=None):
        """
        searches case folder for file, matching base on filename
        check parent directory if not found

        Args:
            str_match: string to match
            str_exclude: string to exclude

        Returns:
            Path to file or "None" if not found
        """
        # type checking
        if str_exclude is None:
            str_exclude = []
        if type(str_exclude) == str:
            str_exclude = [str_exclude]
        # loop over files in case directory
        for name_file in os.listdir(self.path):
            if str_match in name_file.lower() and not any(x in name_file.lower() for x in str_exclude):
                return os.path.join(self.path, name_file)
        abspath = os.path.abspath(self.path)
        uppath = os.path.dirname(os.path.abspath(self.path))

        for name_file in os.listdir(uppath):
            if str_match in name_file.lower() and not any(x in name_file.lower() for x in str_exclude):
                return os.path.join(uppath, name_file)
        print(abspath)
        print(uppath)
        print(f"\033[91mcould not find file matching '{str_match}'\033[0m")

        return None

    def find_nam_parameter(self, str_match, expected_type="str"):
        """
        finds value of .nam parameter in .nam list file

        Args:
            str_match (str): string used for name matching
            expected_type(str): expected type of value,
                       supported values: int, float, str, int_list, float_list, str_list
        """
        parameter = ""
        with open(self.path_nam, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if re.search(rf' *{str_match} *=', line):
                    if "!" in line:
                        parameter = line[line.find("=")+1:line.find("!")]
                    else:
                        parameter = re.findall(r'= *(.*)(?=,*\s*$)', line)[0]
                    # remove leading and trailing whitespaces
                    parameter = parameter.strip()
                    # remove last ',' symbol if present
                    if parameter[-1] == ",":
                        parameter = parameter[:-1]
                    break
        # format results according to expected_type
        if parameter:
            if expected_type == 'int':
                return int(float(parameter))
            if expected_type == 'float':
                return float(parameter)
            if expected_type == 'str':
                return parameter
            if expected_type == 'list_int':
                parameter = [int(i) for i in list(eval(parameter))]
                return parameter
            if expected_type == 'list_float':
                return list(eval(parameter))
            if expected_type == 'list_str':
                return list(eval(parameter))
        return None

    def find_out_parameter(self, str_match):
        """
        finds and returns value of parameter in out file
        Args:
            str_match (str): match string
        """
        with open(self.path_out, 'r') as f:
            lines = f.readlines()
            for line in lines[::-1]:
                if re.search(rf' *{str_match} *', line):

                    parameter = re.findall(rf' *{str_match} *(.*) ', line)[0]
                    return parameter
        return None

    def get_final_voltages(self):
        """
        extracts final voltages from .out file
        Returns:
            list[float]: final voltage amplitudes
        """
        voltages = []
        for num_metal, metal in enumerate(self.metal_labels):
            # check if voltages are adjusted to match power target
            if self.find_nam_parameter("IRFPOW", expected_type='int'):
                parameter = 0
                with open(self.path_out, 'r') as f:
                    lines = f.readlines()
                    for line in lines[::-1]:
                        regex = rf' +{metal} +\d\.\d\d\dE\+\d\d  +(\d\.\d\d\dE\+\d\d) +'
                        if re.search(regex, line):
                            parameter = eval(re.findall(regex, line)[0])
                            break
                voltages.append(parameter)
            # if voltages are not adjusted use initial voltages
            else:
                voltages = self.find_nam_parameter("VRFM", expected_type='list_float')
                # remove commas and interpret as float
                #for i, s in enumerate(voltages):
                 #   voltages[i] = float(s.replace(",", ""))

        return voltages

    def read_mesh_file(self):
        """

        Returns:
            np.ndarray: 2D array containing the mesh geometry
        """
        if self.path_mesh is None:
            print("path to mesh file not found")
            print("mesh file not read")
            return None

        with open(self.path_mesh) as f:
            mesh = []
            lines = f.readlines()
            for line in lines:
                reg_phrase_end = r"^[\t\s]*\*+[\t\s]*"
                if re.search(reg_phrase_end, line):
                    break
                reg_phrase_start = r"[\t\s]*#*\*+"
                if re.search(reg_phrase_start, line):
                    continue
                mesh.append(list(line[1:].rstrip()))
        mesh = np.array(mesh)
        return mesh

    def read_pcmc_file(self, normalize_theta=True):
        """
        reads energy angular distribution from pcmc file
        Args:
            normalize_theta (bool): flag for the normalization based on solid angle

        Returns:
            array containing EADS and list of found species

        """
        # check if path to pcmc file is set. if not, return None
        if self.path_pcmc is None:
            print("pcmc file not found! IEAD ad EEAD plotting aborted")
            return None, None

        print(f"   Case: {self.name}")
        print(f"   File path: {self.path_pcmc}")
        print(f"   Species:\n     ", end='')

        # initialize lists
        pcmc_species = []
        pcmc_eads = []
        # find formatting
        with open(self.path_pcmc) as f:
            reg_phrase = r"[\t\s]*(\d+)[\t\s]+(\d+\.\d+E[+-]\d+)[\t\s]+(\d+\.\d+E[+-]\d+)"
            lines = f.readlines()
            num_lines_matched = 0
            for num_line, line in enumerate(lines):
                if re.search(reg_phrase, line):
                    num_lines_matched += 1
                    results = re.search(reg_phrase, line).groups()
                    if num_lines_matched == 1:
                        pcmc_angle_num_bins = eval(results[0])
                        pcmc_angle_max_neutrals = eval(results[1])
                        pcmc_angle_max_ions = eval(results[2])
                        self.pcmc_max_angle = (pcmc_angle_max_ions, pcmc_angle_max_neutrals)

                    if num_lines_matched == 2:
                        pcmc_energy_num_bins = eval(results[0])
                        pcmc_energy_max_neutrals = eval(results[1])
                        pcmc_energy_max_ions = eval(results[2])
                        self.pcmc_max_energy = (pcmc_energy_max_ions, pcmc_energy_max_neutrals)

                    if num_lines_matched >= 2:
                        pcmc_line_header_end = num_line+2
                        break

        # read data:
        # ----------------------------------------------------------------------
        # initialize lists
        with open(self.path_pcmc) as f:
            reg_phrase_exp = r"(\d+\.\d+E[+-]\d+)"
            num_species_found = 0
            # skip header lines
            for i in range(pcmc_line_header_end):
                next(f)
            # loop over lines
            for line in f:
                # find species line ( if line does not start with space )
                if line[0] != " ":
                    if num_species_found > 0:
                        pcmc_eads.append(np.array(rows))
                    name_species = line.strip()
                    pcmc_species.append(name_species)
                    num_species_found += 1
                    rows = []
                    row = []
                    print(f"{name_species},  ", end='')
                    continue

                # find exponential numbers in line
                numbers_found = re.findall(reg_phrase_exp, line)
                if numbers_found:
                    for n in numbers_found:
                        if len(row) < pcmc_energy_num_bins:
                            row.append(eval(n))
                        if len(row) == pcmc_energy_num_bins:
                            rows.append(row)
                            row = []
            print("\n")
            # add eads to list after last iteration
            pcmc_eads.append(np.array(rows))

            # add to Case object (modified)
            self.pcmc_species = pcmc_species

            # normalize magnitude to total flux
            normalize_magnitude = False
            if normalize_magnitude:
                for i, ead in enumerate(pcmc_eads):
                    pcmc_eads[i] = ead / sum(ead)

            # reshape
            for i, array in enumerate(pcmc_eads):
                height, width = array.shape
                array[0:int(height / 2), :] = np.flipud(array[0:int(height / 2), :])
                array = array.transpose()
                pcmc_eads[i] = array

        # normalize theta by solid angle
        if normalize_theta:
            dthi = pcmc_angle_max_ions/pcmc_angle_num_bins
            thim = np.zeros(pcmc_angle_num_bins)
            factors = []
            for i, j in enumerate(thim):
                thim[i] = (dthi*(i+1))
                factors.append(dthi*np.sin(thim[i]*3.141592/180.))
            factors = factors[::-1] + factors

            for i, ead in enumerate(pcmc_eads):
                for angle_bin in range(90):
                    ead[:, angle_bin] = ead[:, angle_bin] / factors[angle_bin]
                pcmc_eads[i] = ead

        self.pcmc_eads = pcmc_eads

        return pcmc_species, pcmc_eads

    def generate_edfs(self):
        """
        integrates Energy Angular Distributions to generate 1D energy distribution
        """
        for ead in self.pcmc_eads:
            edf = np.sum(ead, axis=1)

            if sum(edf) != 0:
                edf = edf/sum(edf)
            self.pcmc_edfs.append(edf)

    def get_mean_and_mode_energies(self):
        """
        determines mean particle energy of  pcmc species
        """

        # check if pcmc_edfs are found in case object
        if not self.pcmc_edfs:
            return

        # loop over edfs
        for i, edf in enumerate(self.pcmc_edfs):
            # generate energy scale
            scale_energy = np.linspace(0, float(self.pcmc_max_energy[0]), num=edf.shape[0])
            # compute mean energy
            energy_mean = 0
            for j, e in enumerate(scale_energy):
                energy_mean += e * edf[j]
            # normalize
            if np.sum(edf) != 0:
                energy_mean = energy_mean / np.sum(edf)
            # append to mean energies list
            self.pcmc_mean_energy.append(energy_mean)

            # find mode energy value
            bin_mode = np.argmax(edf)
            mode_energy = scale_energy[bin_mode]
            # append to mode energy list
            self.pcmc_mode_energy.append(mode_energy)

    def generate_adfs(self):
        """
        integrates Energy Angular Distributions to generate 1D angular distribution
        """
        for ead in self.pcmc_eads:
            adf = np.sum(ead, axis=1)
            adf = adf/sum(adf)
            self.pcmc_edfs.append(adf)

    def movie_find_dimensions(self):
        """
        finds dimensions of the movie file
        """
        # skip if no movie1 file was found
        if not self.path_movie1:
            return None, None, None

        num_zones = self.find_nam_parameter("IMOVIE_FRAMES", expected_type='int') + 1
        # in tecplot notation, array dimensions are generally I,J,K
        # here denoted as "idim", "jdim" and "kdim" in order to conform to code style conventions
        # find dimensions of Tecplot file
        # --------------------------------------------------------------------
        # find 'I' and 'J' variable list
        with open(self.path_movie1) as f:
            for line, row in enumerate(f):
                # find I
                if re.search(r'.*I= *([0-9]*),', row):
                    idim = int(re.findall(r'.*I= *([0-9]*),', row)[0])
                # find J
                if re.search(r'.*J= *([0-9]*),', row):
                    jdim = int(re.findall(r'.*J= *([0-9]*),', row)[0])
                    break

        # Find max T if not provided
        if num_zones:
            return idim, jdim, num_zones
        else:
            # loop over lines beginning from last
            for line in reversed(list(open(self.path_movie1))):
                # fast coarse match
                if re.search(r'.*?T.*?', line):
                    # precise match
                    regexpr = r'.*=\s*(\d+).*$'
                    if re.search(regexpr, line):
                        num_zones = int(re.findall(regexpr, line)[0])
                        return idim, jdim,  num_zones

    def movie_read_to_xt_array(self, quantity):
        """
        loads xt data of 'quantity into array' and returns it

        Args:
            quantity (str): Name of quantity to read
        """
        # skip if no movie1 file was found
        if not self.path_movie1:
            return None
        # find 'quantity' position in file
        pos_in_zone = None
        with open(self.path_movie1) as f:
            for line, row in enumerate(f):
                if quantity in row:
                    pos_in_zone = line - 2  # no of variable, R and Z omitted starting with 0
                    break
        # check if quantity was found
        if not pos_in_zone:
            print(f"could not find quantity '{quantity}' in movie file")
            return None

        # in tecplot notation, array dimensions are generally I,J,K
        # here denoted as "idim", "jdim" and "kdim" in order to conform to code style conventions

        idim = self.movie_I
        jdim = self.movie_J
        zones = self.movie_num_zones

        # read in data
        # --------------------------------------------------------------------
        lines_per_i = math.ceil(idim / 7)

        # find lines of zones
        row = []
        line_i_begin = 0
        line_i_end = 0
        array = np.zeros((jdim, idim, zones))
        zone = -1
        j = 0
        line_var_begin = 0
        line_var_end = 0
        with open(self.path_movie1, "r") as fp:
            for line, rowtext in enumerate(fp):

                if 'ZONE' in rowtext:
                    zone += 1
                    line_zone = line + 2

                    # because constants X and R are not repeated
                    line_var_begin = line_zone + (jdim * lines_per_i * pos_in_zone) + (j * lines_per_i)
                    if zone == 0:
                        line_var_begin = line_zone + (jdim * lines_per_i * (pos_in_zone + 2)) + (j * lines_per_i)
                    line_var_end = line_var_begin + lines_per_i * jdim
                    line_i_begin = line_var_begin
                    line_i_end = line_i_begin + lines_per_i

                if line_var_begin <= line < line_var_end:

                    if line_i_begin <= line < line_i_end:
                        row.extend([float(x) for x in rowtext.split()])

                    if line == line_i_end - 1:
                        array[j, :, zone] = row
                        j += 1
                        line_i_begin = line
                        line_i_end = line_i_begin + lines_per_i + 1
                        row = []
                    if j == jdim:
                        j = 0

        return array

    def get_local_potential_over_time(self, locs, labels):
        """
        fetches the local potential at point xy from movie file
        Args:
            locs (list of tuple): List of tuples of coordinates for probing
            labels (list of str): Labels of probe points
        """
        # call movie_read_to_xt_array if self.potential_xt is None r 0
        if not self.potential_xt:
            self.potential_xt = self.movie_read_to_xt_array("PPOT")

        # check if XT-potential was found
        if self.potential_xt is None:
            return None

        if type(locs) is tuple and len(locs) == 2:
            locs = [locs]

        for loc in locs:
            if loc[1] < self.potential_xt.shape[0] and loc[0] < self.potential_xt.shape[1]:
                self.potential_over_time.append(self.potential_xt[loc[1], loc[0], :])
        self.potential_over_time_labels = labels

        return


class Cases:
    def __init__(self):
        self.cases = []
        self.dir_source = None
        self.constant_power = None
        self.constant_phase = None

        self.mean_energies = []
        self.phases = []

    def __getitem__(self, i) -> Case:
        """

        Args:
            i (int): case index

        Returns:
            Case: 
        """
        return self.cases[i]

    def add_case(self, case):
        self.cases.append(case)

    def __len__(self):
        return len(self.cases)

    def scan_directory_for_cases(self, dir_root):
        """
        scans directory for cases directories.
        identified based on the presence of a .log file
        Args:
            dir_root (str): root directory for scanning
        """
        dirs_in_root = os.listdir(dir_root)

        print("scanning directory...")
        for dir_in_root in dirs_in_root:

            dir_sub = os.path.join(dir_root, dir_in_root)
            # check if file:
            if os.path.isfile(dir_sub):
                continue

            # check if directory
            if os.path.isdir(dir_sub):

                # if directory, find out if source or case directory or neither
                # -------------------------------------------------------------------
                paths_in_sub = os.listdir(dir_sub)
                # loop over files and folders
                for path_in_sub in paths_in_sub:

                    # check if its a file

                    if os.path.isfile(os.path.join(dir_sub, path_in_sub)):
                        # check for specific files to determine if it is a source file
                        # ---------------------------------------------------------------------
                        # if it's a source file label the directory as source directory

                        if "icp.f" in path_in_sub.lower() or "icp.for" in path_in_sub.lower():
                            self.dir_source = os.path.join(path_in_sub)
                            break

                        if ".log" in path_in_sub.lower():
                            path_case = os.path.join(dir_sub)
                            self.add_case(Case(path_case))
                            break

                    # if it is a directory, search directory
                    else:
                        paths_in_sub_sub = os.listdir(os.path.join(dir_sub, path_in_sub))
                        for path_in_sub_sub in paths_in_sub_sub:

                            if ".log" in path_in_sub_sub.lower():
                                path_case = os.path.join(dir_sub, path_in_sub)
                                self.add_case(Case(path_case))
                                break

    def get_value_pair(self, name_var1, name_var2, custom_waveform_only=False):
        """
        extracts value pair from cases and sorts by first value
        Args:
            name_var1 (str): name of first variable (x)
            name_var2 (str): name of first variable (y)
            custom_waveform_only (bool): only consider cases that contain custom waveforms

        Returns: nested list with sorted variables (x and y)

        """
        var1 = []
        var2 = []
        for case in self.cases:
            # check if custom waveform case:
            if custom_waveform_only and not case.contains_custom:
                continue
            var1.append(eval("case."+name_var1))
            var2.append(eval("case."+name_var2))
            # sort lists based on var1
            var1, var2 = (list(t) for t in zip(*sorted(zip(var1, var2))))
        return [var1, var2]

    def determine_sweep(self):

        powers = []
        metals = []
        phases = []
        const_power = True
        const_phase = True
        const_voltage = True

        for case in self.cases:
            powers.append(case.rfpnorma)
            metals.append(case.metal_labels)
            phases.append(case.cwaveform_phase)

        for num_case, case in enumerate(self.cases):
            if case.irfpow == 0:
                const_power = False
            else:
                const_voltage = False

            if powers[num_case] != powers[0]:
                const_power = False
            if phases[num_case] != phases[0]:
                const_phase = False

        # set object values
        self.constant_power = const_power
        self.constant_phase = const_phase

        # print sweep metrics
        print("\nanalyzing sweep metrics:\n-----------------------------------")
        print(powers)
        if const_power:
            print("constant power:")
            for i, p in enumerate(powers[0]):
                if p != 0:
                    print(f"   power on '{metals[0][i]}': {p} W")
        else:
            print("varied power")

        if const_phase:
            print(f"constant   phases : {phases[0]}] °")
        else:
            print("varied phase:")
            for i, p in enumerate(phases):
                print(f"   phase shift: {p} °")

    def print_info(self):
        """
        prints basic information about cases
        """
        print(f"{len(self.cases)} cases found:")
        for case in self.cases:
            print("\n   " + case.name)
            print("      DC-bias:", case.dc_bias)
            print("      custom waveforms:", case.contains_custom)
            print("      phase:", case.cwaveform_phase)

    def read_pcmc_file_all(self, config):
        """
        reads all pcmc data into case objects
        optional generation of 1D energy distribution and/or 1D angular distribution

        Args:
            config (Config) : config object
        """
        print("\nreading pcmc files:")
        print("-" * 80)
        for case in self.cases:
            if case.path_pcmc is not None:
                case.read_pcmc_file()
            if config.plot_EDFs:
                case.generate_edfs()
                case.get_mean_and_mode_energies()

    def load_sweep_data(self):
        for case in self.cases:
            # skip if case contains no tailored waveform
            if not case.contains_custom:
                continue
            # find mean energies corresponding to species to be plotted
            for i, species in enumerate(case.pcmc_species):
                if species == "ION-TOT":
                    self.mean_energies.append(case.pcmc_mean_energy[i])
                    break
            self.phases.append(case.custom_phase[1])

            # sort lists based on phase
            self.phases, self.mean_energies = (list(t) for t in zip(*sorted(zip(self.phases, self.mean_energies))))

    def save_data_to_csv(self, dir_root):

        # generate data folder if non existent
        path_data = os.path.join(dir_root, "Data")
        if not os.path.isdir(path_data):
            os.mkdir(path_data)
        print("saving data to file" )
        f = open(os.path.join(path_data, "data.csv"), "w")
        f.write("phase, mean energy\n")
        for i, phase in enumerate(self.phases):
            f.write(f"{phase}, {self.mean_energies[i]}\n")
        f.close()



import os
import re
import numpy as np
import pylab as pl


class Case:
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)
        self.path_nam = self.find_file_by_name(".nam", [])
        self.path_log = self.find_file_by_name(".log", ["runscript"])
        self.path_out = self.find_file_by_name(".out", ["gmon", "mcs.", "mesh."])
        self.path_movie1 = self.find_file_by_name("movie1.pdt")
        self.path_mesh = self.find_file_by_name("mesh.dat")
        self.path_pcmc = self.find_file_by_name("pcmc.prof")
        self.mesh = self.read_mesh_file(self.path_mesh)
        # pcmc data
        self.pcmc_species = []
        self.pcmc_eads = []
        self.pcmc_max_angle = (None, None)                                    # max pcmc angle tuple: (ions, neutrals)
        self.pcmc_max_energy = (None, None)                                  # max pcmc energy tuple: (ions, neutrals)
        self.irfpow = int(eval(self.find_nam_parameter("IRFPOW")))           # adjust voltages for power target
        self.powerICP = 0
        self.powerCCP1 = 0
        self.powerCCP2 = 0
        self.cwafer = list(eval(self.find_nam_parameter("CWAFER")))              # wafer materials
        self.metal_labels = list(eval(self.find_nam_parameter("CMETAL")))        # metal material labels
        self.final_voltages = self.get_final_voltages()                          # voltage amplitudes on last iteration
        self.rfpnorma = list(eval(self.find_nam_parameter("RFPNORMA")))
        self.dc_bias = eval(self.find_out_parameter("DC BIAS"))
        self.ne_ave = eval(self.find_out_parameter("AVERAGE ELECTRON DENSITY"))
        self.restart = int(self.find_nam_parameter("IRESTART"))
        self.icustom = list(eval(self.find_nam_parameter("ICUSTOM")))
        self.contains_custom = self.icustom != [0] * len(self.icustom)
        self.custom_phase = eval(self.find_nam_parameter("CUSTOM_PHASE"))
        self.custom_relharm = eval(self.find_nam_parameter("CUSTOM_RELHARM"))
        self.custom_relamp = eval(self.find_nam_parameter("CUSTOM_RELAMP"))
        self.cwaveform_phase = eval(self.find_nam_parameter("CUSTOM_PHASE"))[1]
        self.rffac = int(float(self.find_nam_parameter("RFFAC")))
        self.freq = float(self.find_nam_parameter("FREQ"))

    def find_file_by_name(self, str_match, str_exclude=[]):
        """
        searches case folder for file, matching base on filename
        check parent directory if not found

        Args:
            str_match: string to match
            str_exclude: string to exclude

        Returns:
            Path to file or "None" if not found

        """
        # loop over files in case directory
        for name_file in os.listdir(self.path):
            if str_match in name_file.lower() and not any(x in name_file.lower() for x in str_exclude):
                return os.path.join(self.path, name_file)
        abspath = os.path.abspath(self.path)
        uppath = os.path.dirname(os.path.abspath(self.path))
        print(abspath)
        print(uppath)
        for name_file in os.listdir(uppath):
            if str_match in name_file.lower() and not any(x in name_file.lower() for x in str_exclude):
                return os.path.join(uppath, name_file)

        print(f"\033[91mcould not find file matching '{str_match}'\033[0m")

        return None

    def find_nam_parameter(self, str_match, multiple_values=False):

        with open(self.path_nam, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if re.search(rf' *{str_match} *=', line):
                    if "!" in line:
                        parameter = re.findall(r'= *(.*)(?=,*\s*\!)', line)[0]
                    else:
                        parameter = re.findall(r'= *(.*)(?=,*\s*$)', line)[0]
                    # remove leading and trailing whitespaces
                    parameter = parameter.strip()
                    # remove last ',' symbol if present
                    if parameter[-1] == ",":
                        parameter = parameter[:-1]
                    return parameter
        return None

    def find_out_parameter(self, str_match):
        with open(self.path_out, 'r') as f:
            lines = f.readlines()
            for line in lines[::-1]:
                if re.search(rf' *{str_match} *', line):

                    parameter = re.findall(rf' *{str_match} *(.*) ', line)[0]
                    return parameter
        return None

    def get_final_voltages(self):
        voltages = []
        for num_metal, metal in enumerate(self.metal_labels):
            parameter = 0
            with open(self.path_out, 'r') as f:
                lines = f.readlines()
                for line in lines[::-1]:
                    regex = rf' +{metal} +\d\.\d\d\dE\+\d\d  +(\d\.\d\d\dE\+\d\d) +'
                    if re.search(regex, line):
                        parameter = eval(re.findall(regex, line)[0])
                        break
            voltages.append(parameter)
        return voltages

    def read_mesh_file(self, path_mesh):
        with open(path_mesh) as f:
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

        # check if path to pcmc file is set. if not, return None
        if self.path_pcmc is None:
            print("pcmc file not found! IEAD ad EEAD plotting aborted")
            return None, None

        print("reading pcmc file:")
        print("species:")

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
            #lines = f.readlines()[pcmc_line_header_end:]
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
                    print(f"   {name_species}")
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
            dthn = pcmc_angle_max_ions/pcmc_angle_num_bins
            THIM = np.zeros(pcmc_angle_num_bins)
            factors = []
            for i, j in enumerate(THIM):
                THIM[i] = (dthi*(i+1))
                factors.append(dthi*np.sin(THIM[i]*3.141592/180.))
            factors = factors[::-1] + factors

            for i, ead in enumerate(pcmc_eads):
                for angle_bin in range(90):
                    ead[:, angle_bin] = ead[:, angle_bin] / factors[angle_bin]
                pcmc_eads[i] = ead#/sum(ead)

        self.pcmc_eads = pcmc_eads

        return pcmc_species, pcmc_eads


class Cases:
    def __init__(self):
        self.cases = []
        self.dir_source = None
        self.constant_power = None
        self.constant_phase = None

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

        dirs_in_root = os.listdir(dir_root)

        print("scanning directory...")
        for dir_in_root in dirs_in_root:

            dir_sub = os.path.join(dir_root, dir_in_root)
            # check if file:
            string_dir_type = "{0:<17}".format("not recognised")
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


